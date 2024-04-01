#include "random.h"
#include <geometry/primitive.h>

#include <cmath>
// todo: remove
#include <iostream>

const float step = 1e-4;

namespace rng
{

Engine get_generator(size_t seed) {
    std::random_device a;
    Engine rng(a());
    if (seed != 0)
        rng.seed(seed);
    return rng;
}

}

static uniform_float_d uniDist(-1.f, 1.f);
static normal_d normDist(0, 1.f);

float UniformDistribution::sample(Engine &rng) {
    return uniDist(rng);
}

float UniformDistribution::norm_sample(Engine &rng) {
    return normDist(rng);
}

vector3f UniformDistribution::uni_sphere_sample(vector3f /*point*/, vector3f normal, Engine &rng) {
    vector3f a{normDist(rng), normDist(rng), normDist(rng)};
    normalize(a);
    if (dot(a, normal) < 0.f)
        a = -a;
    return a;
}

vector3f UniformDistribution::sphere_sample(vector3f point, vector3f normal, Engine &rng) {
    return uni_sphere_sample(point, normal, rng);
}

float UniformDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    return M_1_PIf32 / 2;
}

CosineWeightedDistribution::CosineWeightedDistribution() {}

vector3f CosineWeightedDistribution::sphere_sample(vector3f point, vector3f normal, Engine &rng) {
    vector3f dir{};
    do {
        dir = UniformDistribution::uni_sphere_sample(point, normal, rng) + normal;
    } while (length(dir) < 1e-12);
    normalize(dir);
    return dir;
}

float CosineWeightedDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    return std::max(0.f, dot(direction, normal)) * M_1_PIf32;
}

LightDistribution::LightDistribution(const Primitive &object) {
    primitive = &object;
}

vector3f LightDistribution::sphere_sample(vector3f point, vector3f normal, Engine &rng) {
    vector3f res{};
    switch (primitive->type) {
        case Primitive::Box: {
            const auto &box_size = std::get<vector3f>(primitive->param_);
            vector3f area{4.f * box_size.y * box_size.z, 4.f * box_size.x * box_size.z, 4.f * box_size.x * box_size.y};
            float full_area = area[0] + area[1] + area[2];

            const float sample = (UniformDistribution::sample(rng) + 1.f) * full_area / 2;
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

            const float side = (UniformDistribution::sample(rng) < 0.f) ? -1.f : 1.f;

            const int face1 = (face + 1) % 3;
            const int face2 = (face + 2) % 3;
            res[face] = side * box_size[face];
            res[face1] = box_size[face1] * UniformDistribution::sample(rng);
            res[face2] = box_size[face2] * UniformDistribution::sample(rng);
            break;
        }
        case Primitive::Ellipsoid: {
            const auto &ellipsoid_radius = std::get<vector3f>(primitive->param_);
            // get point on sphere
            res = ::normal(vector3f{UniformDistribution::norm_sample(rng), UniformDistribution::norm_sample(rng), UniformDistribution::norm_sample(rng)});
            res *= ellipsoid_radius;
            break;
        }
        case Primitive::Triangle: {
            auto params = std::get<std::array<vector3f, 3>>(primitive->param_);
            vector3f &triangle_origin = params[0];
            vector3f &U = params[1];
            vector3f &V = params[2];
            //todo
            float u = (UniformDistribution::sample(rng) + 1.f) / 2;
            float v = (UniformDistribution::sample(rng) + 1.f) / 2;
            if (u + v > 1) {
                u = 1 - u;
                v = 1 - v;
            }
            res = triangle_origin + u * U + v * V;
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
    if (primitive->type != Primitive::Triangle) {
        intersection[1] = primitive->intersect(
                {point + direction * intersection[0].distance + direction * step, direction});
        if (intersection[1]) {
            intersection[1].distance += intersection[0].distance;
        }
    }


    switch (primitive->type) {
        case Primitive::Box: {
            const auto &box_size = std::get<vector3f>(primitive->param_);
            float full_area_size = 8 * (box_size[0] * box_size[1] + box_size[1] * box_size[2] + box_size[2] * box_size[0]);
            probability[0] = 1 / full_area_size;
            probability[1] = probability[0];
            break;
        }
        case Primitive::Ellipsoid: {
            const auto &ellipsoid_radius = std::get<vector3f>(primitive->param_);

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
        case Primitive::Triangle: {
            probability[0] = 1.f / primitive->cache.triangle_area;
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

void MixedDistribution::add_dist(const random_distribution_sh_ptr& dist) {
    distributions.push_back(dist);
    count.push_back(0);
}

vector3f MixedDistribution::sphere_sample(vector3f point, vector3f normal, Engine &rng) {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    int sample = std::floor((UniformDistribution::sample(rng) + 1.f) * 0.5f * (float)distributions.size());
    if (sample == distributions.size())
        sample -= 1;
    count[sample]++;
    return distributions[sample]->sphere_sample(point, normal, rng);
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

MixedDistribution::~MixedDistribution() {
//    std::cout << "Calls:\n";
//    for (int i = 0; i < distributions.size(); ++i) {
//        std::cout << "\t";
//        if (dynamic_cast<LightDistribution*>(distributions[i].get())) {
//            std::cout << "Light: ";
//        } else if (dynamic_cast<MixedDistribution*>(distributions[i].get())) {
//            std::cout << "Mix: ";
//        } else {
//            std::cout << "Cos: ";
//        }
//        std::cout << count[i] << std::endl;
//    }
}
