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
    float sample = UniformDistribution::sample(rng);
    if (sample < 0.f || !light.size())
        return cosine.sample(point, normal, rng);
    else
        return light.sample(point, normal, rng);
}

float SceneDistribution::pdf(vector3f point, vector3f normal, vector3f direction) {
    if (!light.size())
        return cosine.pdf(point, normal, direction);
    return (cosine.pdf(point, normal, direction) + light.pdf(point, normal, direction)) * 0.5f;
}
