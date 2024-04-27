#include "random.h"
#include <geometry/primitive.h>

#include <cmath>

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

vector3f UniformDistribution::uni_sphere_sample(Engine &rng) {
    vector3f a{normDist(rng), normDist(rng), normDist(rng)};
    return ::normal(a);
}

vector3f UniformDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    vector3f sphere_sample = uni_sphere_sample(rng);
    if (dot(sphere_sample, normal) < 0.f)
        sphere_sample = -sphere_sample;
    return sphere_sample;
}

float UniformDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    return M_1_PIf32 / 2;
}

vector3f CosineWeightedDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    vector3f dir{};
    do {
        dir = UniformDistribution::uni_sphere_sample(rng) + normal;
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

vector3f LightDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    vector3f res{};

    auto params = primitive->position;
    float u = (UniformDistribution::sample(rng) + 1.f) / 2;
    float v = (UniformDistribution::sample(rng) + 1.f) / 2;
    if (u + v > 1) {
        u = 1 - u;
        v = 1 - v;
    }
    const vector3f &triangle_origin = params[0];
    const vector3f &U = params[1];
    const vector3f &V = params[2];

    res = triangle_origin + u * U + v * V;
    return ::normal(res - point);
}

float LightDistribution::pdfWithHint(vector3f point, vector3f direction, Intersection hint) {
    float probability = 0.f;
    Intersection intersection = hint;

    // it is guaranteed that cache was pre-calculated
    probability = 1.f / primitive->cache.triangle_area;

    return std::abs(probability * (intersection.distance * intersection.distance) / dot(intersection.normal, direction));
}

float LightDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    Intersection intersection;

    intersection= primitive->intersect({point, direction});
    if (!intersection) {
        return 0.f;
    }
    return pdfWithHint(point, direction, intersection);
}

vector3f VisibleNormalDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    // find rotation quaternion from local Z to absolute normal;
    vector3f Z{0.f, 0.f, 1.f};
    vector4f rotation = quat_from_two_vectors(normal, Z);
    vector3f Ve = rotate(eye_direction, rotation);

    vector3f Vh = ::normal(vector3f{Ve.x * alpha.x, Ve.y * alpha.y, Ve.z});
    float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
    vector3f T1 = lensq > 0 ? vector3f{-Vh.y, Vh.x, 0} * (1.f/ std::sqrt(lensq)) : vector3f{1,0,0};
    vector3f T2 = cross(Vh, T1);

    // get uniform sample on a disc
    vector2f t{};
    {
        int a = std::floor((UniformDistribution::sample(rng) + 2.f) * 2.f);
        vector2f u{1.f, 1.f};
        while (u.x * u.x + u.y * u.y > 1.f) {
            u.x = UniformDistribution::sample(rng);
            for (int i = -1; i < a; ++i) {
                u.y = UniformDistribution::sample(rng);
            }
        }
        t = u;
    }

    float s = 0.5f + 0.5f * Vh.z;
    t.y = (1.f - s) * std::sqrt(1.f - t.x * t.x) + s * t.y;

    vector3f Nh = t.x * T1 + t.y * T2 + std::sqrt(std::max(0.f, 1.f - t.x * t.x - t.y * t.y)) * Vh;
    vector3f Ne = ::normal(vector3f{alpha.x * Nh.x, alpha.y * Nh.y, std::max(0.f, Nh.z)});

//    rotation = quat_from_two_vectors(Z, normal);
    Ne = rotate(Ne, *rotation);

    return ::normal(2 * Ne * dot(Ne, eye_direction) - eye_direction);
}

float VisibleNormalDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    vector3f Z{0.f, 0.f, 1.f};
    vector4f rotation = quat_from_two_vectors(normal, Z);

    vector3f sample_normal = ::normal(direction + eye_direction);

    vector3f V = rotate(eye_direction, rotation);
    vector3f Ni = rotate(sample_normal, rotation);

    if (dot(V, Ni) < 0.f)
        return 0.f;

    float ax2 = alpha.x * alpha.x;
    float ay2 = alpha.y * alpha.y;

    float invD = M_PIf32 * alpha.x * alpha.y * std::pow(Ni.x * Ni.x / ax2 + Ni.y * Ni.y / ay2 + Ni.z * Ni.z, 2);
    float invG1 = 0.5f + 0.5f * std::sqrt(1.f + (ax2 * V.x * V.x + ay2 * V.y * V.y) / (V.z * V.z));

    return 1.f / (4.f * invD * invG1 * dot(V, Z));
}

void VisibleNormalDistribution::update(vector2f alpha_, vector3f eye_direction_) {
    alpha.x = std::max(0.05f, alpha_.x);
    alpha.y = std::max(0.05f, alpha_.y);
    eye_direction = normal(eye_direction_);
}

void MixedDistribution::add_distr(const random_distribution_sh_ptr& dist) {
    distributions.push_back(dist);
}

vector3f MixedDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    int sample = std::floor((UniformDistribution::sample(rng) + 1.f) * 0.5f * (float)distributions.size());
    if (sample == distributions.size())
        sample -= 1;
    return distributions[sample]->sample(point, normal, rng);
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

size_t MixedDistribution::size() const {
    return distributions.size();
}

ManyLightsDistribution::ManyLightsDistribution(const std::vector<primitive_sh_ptr>& primitives) {
    for (const auto& p : primitives) {
        if (!p->emissive())
            continue;
        objects.push_back(p);
    }
    if (objects.empty())
        return;
    bvh.buildBVH(objects);
    for (auto &o: objects) {
        distributions.emplace_back(*o);
    }
}

vector3f ManyLightsDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    int sample = std::floor((UniformDistribution::sample(rng) + 1.f) * 0.5f * (float)distributions.size());
    if (sample == distributions.size())
        sample -= 1;
    return distributions[sample].sample(point, normal, rng);
}

float ManyLightsDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    float prob = 0.f;
    std::vector<Intersection> intersections = bvh.intersectAll(objects, {point, direction});
    for (auto &i : intersections) {
        LightDistribution &d = distributions[i.object_id];
        prob += d.pdfWithHint(point, direction, i);
    }
    return prob / (float)distributions.size();
}

size_t ManyLightsDistribution::size() const {
    return distributions.size();
}

vector3f SceneDistribution::sample(vector3f point, vector3f normal, Engine &rng) {
    float sample = (UniformDistribution::sample(rng) + 1.f) * 3 * 0.5f;
    if (sample <= 1.f || !light.size())
        return cosine.sample(point, normal, rng);
    else if (sample <= 2.f)
        return light.sample(point, normal, rng);
    else
        return vndf.sample(point, normal, rng);
}

float SceneDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    if (!light.size())
        return cosine.pdf(point, normal, direction);
    return (cosine.pdf(point, normal, direction) + light.pdf(point, normal, direction) + vndf.pdf(point, normal, direction)) / 3;
//    return (cosine.pdf(point, normal, direction) + light.pdf(point, normal, direction)) / 2;
//    return cosine.pdf(point, normal, direction);
//    return vndf.pdf(point, normal, direction);
}

void SceneDistribution::update_vndf(float roughness2, vector3f direction) {
    vndf.update({roughness2, roughness2}, -direction);
}
