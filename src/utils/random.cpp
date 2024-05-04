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


static uniform_float_d uniDist(-1.f, 1.f);
static normal_d normDist(0, 1.f);

float uniform::sample(Engine &rng) {
    return uniDist(rng);
}

vector3f uniform::sphere(Engine &rng) {
    vector3f a{normDist(rng), normDist(rng), normDist(rng)};
    return ::normal(a);
}

vector3f uniform::hemisphere(vector3f normal, Engine &rng) {
    vector3f sphere_sample = sphere(rng);
    if (dot(sphere_sample, normal) < 0.f)
        sphere_sample = -sphere_sample;
    return sphere_sample;
}

float uniform::pdf() {
    return M_1_PIf32 / 2;
}

float normal::sample(Engine &rng) {
    return normDist(rng);
}

vector3f cosine_weighted::sample(vector3f normal, Engine &rng) {
    vector3f dir{};
    do {
        dir = uniform::sphere(rng) + normal;
    } while (length(dir) < 1e-12);
    normalize(dir);
    return dir;
}

float cosine_weighted::pdf(vector3f normal, vector3f direction) {
    return std::max(0.f, dot(direction, normal)) * M_1_PIf32;
}

LightDistribution::LightDistribution(const Primitive &object) {
    primitive = &object;
}

vector3f LightDistribution::sample(vector3f point, Engine &rng) const {
    vector3f res{};

    auto params = primitive->position;
    float u = (uniform::sample(rng) + 1.f) / 2;
    float v = (uniform::sample(rng) + 1.f) / 2;
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

float LightDistribution::pdfWithHint(vector3f direction, Intersection hint) const {
    float probability = 0.f;
    Intersection intersection = hint;

    // it is guaranteed that cache was pre-calculated
    probability = 1.f / primitive->cache.triangle_area;

    return std::abs(probability * (intersection.distance * intersection.distance) / dot(intersection.normal, direction));
}

float LightDistribution::pdf(vector3f point, vector3f direction) const {
    Intersection intersection;

    intersection= primitive->intersect({point, direction});
    if (!intersection) {
        return 0.f;
    }
    return pdfWithHint(direction, intersection);
}

vector3f visible_normal::sample(vector3f normal, vector3f eye_direction, float alpha, Engine &rng) {
    // find rotation quaternion from local Z to absolute normal;
    vector3f Z{0.f, 0.f, 1.f};
    vector4f rotation = quat_from_two_vectors(normal, Z);
    vector3f Ve = rotate(eye_direction, rotation);

    vector3f Vh = ::normal(vector3f{Ve.x * alpha, Ve.y * alpha, Ve.z});
    float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
    vector3f T1 = lensq > 0 ? vector3f{-Vh.y, Vh.x, 0} * (1.f/ std::sqrt(lensq)) : vector3f{1,0,0};
    vector3f T2 = cross(Vh, T1);

    // get uniform sample on a disc
    vector2f t;
    {
        vector2f u;
        do {
            u.x = uniform::sample(rng);
            u.y = uniform::sample(rng);
        } while (u.x * u.x + u.y * u.y > 1.f);
        t = u;
    }

    float s = 0.5f + 0.5f * Vh.z;
    t.y = (1.f - s) * std::sqrt(1.f - t.x * t.x) + s * t.y;

    vector3f Nh = t.x * T1 + t.y * T2 + std::sqrt(std::max(0.f, 1.f - t.x * t.x - t.y * t.y)) * Vh;
    vector3f Ne = ::normal(vector3f{alpha * Nh.x, alpha * Nh.y, std::max(0.f, Nh.z)});

    Ne = rotate(Ne, *rotation);

    return ::normal(2 * Ne * dot(Ne, eye_direction) - eye_direction);
}

float visible_normal::pdf(vector3f normal, vector3f eye_direction, float alpha, vector3f direction) {
    vector3f Z{0.f, 0.f, 1.f};
    vector4f rotation = quat_from_two_vectors(normal, Z);

    vector3f sample_normal = ::normal(direction + eye_direction);

    vector3f V = rotate(eye_direction, rotation);
    vector3f Ni = rotate(sample_normal, rotation);

    if (dot(V, Ni) < 0.f)
        return 0.f;

    float alpha2 = alpha * alpha;

    float invD = M_PIf32 * alpha * alpha * std::pow(Ni.x * Ni.x / alpha2 + Ni.y * Ni.y / alpha2 + Ni.z * Ni.z, 2);
    float invG1 = 0.5f + 0.5f * std::sqrt(1.f + (alpha2 * V.x * V.x + alpha2 * V.y * V.y) / (V.z * V.z));

    return 1.f / (4.f * invD * invG1 * dot(V, Z));
}

ManyLightsDistribution::ManyLightsDistribution(const std::vector<Primitive>& primitives) {
    for (const auto& p : primitives) {
        if (!p.emissive())
            continue;
        objects.push_back(p);
    }
    if (objects.empty())
        return;
    bvh.buildBVH(objects);
    for (auto &o: objects) {
        distributions.emplace_back(o);
    }
}

vector3f ManyLightsDistribution::sample(vector3f point, Engine &rng) const {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    int sample = std::floor((uniform::sample(rng) + 1.f) * 0.5f * (float)distributions.size());
    if (sample == distributions.size())
        sample -= 1;
    return distributions[sample].sample(point, rng);
}

float ManyLightsDistribution::pdf(vector3f point, vector3f direction) const {
    if (distributions.empty())
        throw std::runtime_error("No distributions to sample from");
    float prob = 0.f;
    std::vector<Intersection> intersections = bvh.intersectAll(objects, {point, direction});
    for (auto &i : intersections) {
        prob += distributions[i.object_id].pdfWithHint(direction, i);
    }
    return prob / (float)distributions.size();
}

size_t ManyLightsDistribution::size() const {
    return distributions.size();
}

vector3f
SceneDistribution::sample(const vector3f &point, const vector3f &normal, const vector3f &eye_direction, float roughness2, Engine &rng) const {
    float sample = (uniform::sample(rng) + 1.f) * 3 * 0.5f;
//    float sample = 3.f;
    if (sample <= 1.f || !light.size())
        return cosine_weighted::sample(normal, rng);
    else if (sample <= 2.f)
        return light.sample(point, rng);
    else
        return visible_normal::sample(normal, eye_direction, roughness2, rng);
}

float
SceneDistribution::pdf(const vector3f &point, const vector3f &normal, const vector3f &eye_direction, float roughness2, const vector3f &direction) const {
    if (!light.size())
        return cosine_weighted::pdf(normal, direction);
    return (cosine_weighted::pdf(normal, direction)
                + light.pdf(point, direction)
                + visible_normal::pdf(normal, eye_direction, roughness2, direction)
            ) / 3;
}

}