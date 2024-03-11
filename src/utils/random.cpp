#include "random.h"
#include "geometry/primitive.h"

#include <cmath>
// todo: remove
#include <iostream>

const float step = 1e-4;

namespace rng
{

Engine get_generator() {
    std::random_device a;
    return {std::minstd_rand(a())};
}

}

RandomDistribution::RandomDistribution() {
    rng = rng::get_generator();
}

UniformDistribution::UniformDistribution(float min, float max) : dist(min, max), norm_dist(0, 1) {}

float UniformDistribution::sample() {
    return dist(rng);
}

float UniformDistribution::norm_sample() {
    return norm_dist(rng);
}

vector3f UniformDistribution::sphere_sample(vector3f /*point*/, vector3f normal) {
    vector3f a{norm_dist(rng), norm_dist(rng), norm_dist(rng)};
    normalize(a);
    if (dot(a, normal) < 0.f)
        a = -a;
    return a;
}

float UniformDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    return M_1_PIf32 / 2;
}

CosineWeightedDistribution::CosineWeightedDistribution() : uni_dist(0, 1) {}

vector3f CosineWeightedDistribution::sphere_sample(vector3f point, vector3f normal) {
    vector3f dir{};
    do {
        dir = uni_dist.sphere_sample(point, normal) + normal;
    } while (length(dir) < 1e-12);
    normalize(dir);
    return dir;
}

float CosineWeightedDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    return std::max(0.f, dot(direction, normal)) * M_1_PIf32;
}

LightDistribution::LightDistribution(const Primitive &object) : uni_dist(-1, 1) {
    primitive = &object;
}

vector3f LightDistribution::sphere_sample(vector3f point, vector3f normal) {
    vector3f res{};
    switch (primitive->type_) {
        case Primitive::Box: {
            const vector3f &box_size = primitive->param_;
            vector3f area{ 4.f * box_size.y * box_size.z, 4.f * box_size.x * box_size.z, 4.f * box_size.x * box_size.y};
            float full_area =area[0] + area[1] + area[2];

            const float sample = (uni_dist.sample() + 1.f) * full_area / 2;
            // choose face based on sample
            int face = -1;
            if (sample < area[0]) {
                face = 0;
            } else if (sample < area[0] + area[1]) {
                face = 1;
            } else if (sample < full_area) {
                face = 2;
            } else {
                throw std::runtime_error("Unexpected face value.");
            }

            const float side = (uni_dist.sample() < 0.f) ? -1.f : 1.f;

            const int face1 = (face + 1) % 3;
            const int face2 = (face + 2) % 3;
            res[face] = side * box_size[face];
            res[face1] = box_size[face1] * uni_dist.sample();
            res[face2] = box_size[face2] * uni_dist.sample();
            break;
        }
        case Primitive::Ellipsoid: {
            const vector3f &ellipsoid_radius = primitive->param_;
            // get point on sphere
            res = ::normal(vector3f{uni_dist.norm_sample(), uni_dist.norm_sample(), uni_dist.norm_sample()});
            res *= ellipsoid_radius;
            break;
        }
        case Primitive::Plane:
        default:
            throw std::runtime_error("rassert failed: light sample from Plane primitive");
    }
    res = primitive->to_global(res);
    return ::normal(res - point);
}

static float a = 0.f;
static int a_times = 0;

float LightDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    float probability[2] = {0, 0};
    Intersection intersection[2];

    intersection[0] = primitive->intersect({point, direction});
    if (!intersection[0]) {
        return 0.f;
    }
    intersection[1] = primitive->intersect({point + direction * intersection[0].distance + direction * step, direction});
    if (intersection[1]) {
        intersection[1].distance += intersection[0].distance;
    }


    switch (primitive->type_) {
        case Primitive::Box: {
            const vector3f &box_size = primitive->param_;
            float full_area_size = 8 * (box_size[0] * box_size[1] + box_size[1] * box_size[2] + box_size[2] * box_size[0]);
            probability[0] = 1 / full_area_size;
            probability[1] = probability[0];
            break;
        }
        case Primitive::Ellipsoid: {
            const vector3f &ellipsoid_radius = primitive->param_;

            for (int i = 0; i <= 1; ++i) {
                if (!intersection[i])
                    break;

                // get normal of point of sphere corresponding to the point of ellipsoid
                vector3f ell_point = primitive->to_local(point + direction * intersection[i].distance);
                vector3f sphere_normal = ell_point / ellipsoid_radius;

                float area_squared = 0.f;
                for (int f = 0; f < 3; ++f) {
                    const int f0 = (f) % 3;
                    const int f1 = (f + 1) % 3;
                    const int f2 = (f + 2) % 3;
                    area_squared += sphere_normal[f0] * sphere_normal[f0]
                                  * ellipsoid_radius[f1] * ellipsoid_radius[f1]
                                  * ellipsoid_radius[f2] * ellipsoid_radius[f2];
                }

                probability[i] = 0.25f * M_1_PIf32 / std::sqrt(area_squared);
            }
            break;
        }
        case Primitive::Plane:
        default:
            return 0;
    }

    float res = 0.f;
    for (int i = 0; i <= 1; ++i) {
        if (!intersection[i])
            break;
        res += std::abs(probability[i] * (intersection[i].distance * intersection[i].distance) / std::abs(dot(intersection[i].normal, direction)));
    }
    return res;
}


MixedDistribution::MixedDistribution() : uni(0, 1) {}

void MixedDistribution::add_dist(const random_distribution_sh_ptr& dist) {
    distributions.push_back(dist);
}

vector3f MixedDistribution::sphere_sample(vector3f point, vector3f normal) {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    int sample = std::floor(uni.sample() * (float)distributions.size());
    if (sample == distributions.size())
        sample -= 1;
    return distributions[sample]->sphere_sample(point, normal);
}

float MixedDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    float prob = 0.f;
    for (auto &d : distributions) {
        prob += d->pdf(point, normal, direction);
    }
    return prob / (float)distributions.size();
}

size_t MixedDistribution::get_size() const {
    return distributions.size();
}
