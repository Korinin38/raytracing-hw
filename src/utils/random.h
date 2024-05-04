#pragma once

#include <geometry/primitive.h>
#include <core/bvh.h>
#include <utils/vector.h>

#include <random>
#include <utility>
#include <vector>
#include <memory>

typedef std::uniform_int_distribution<int> uniform_int_d;
typedef std::uniform_real_distribution<float> uniform_float_d;
typedef std::normal_distribution<float> normal_d;

//typedef pcg32 RandomGenerator;
typedef std::minstd_rand Engine;

namespace rng {
class MixedDistribution;
class SceneDistribution;
}
typedef std::shared_ptr<rng::MixedDistribution> mixed_distribution_sh_ptr;
typedef std::shared_ptr<rng::SceneDistribution> scene_distribution_sh_ptr;

namespace rng {
    Engine get_generator(size_t seed = 0);

namespace uniform {
    float sample(Engine &rng);
    vector3f sphere(Engine &rng);
    vector3f hemisphere(vector3f normal, Engine &rng);
    float pdf();
}

namespace normal {
    float sample(Engine &rng);
}

namespace cosine_weighted {
    static vector3f sample(vector3f normal, Engine &rng);
    static float pdf(vector3f normal, vector3f direction);
}

namespace visible_normal {
    static vector3f sample(vector3f normal, vector3f eye_direction, float alpha, Engine &rng);
    static float pdf(vector3f normal, vector3f eye_direction, float alpha, vector3f direction);
}

    class LightDistribution {
    public:
        explicit LightDistribution(const Primitive &object);
        vector3f sample(vector3f point, Engine &rng) const;
        float pdf(vector3f point, vector3f direction) const;
        inline float pdfWithHint(vector3f direction, Intersection hint) const;
    private:
        const Primitive *primitive;
    };

    class ManyLightsDistribution {
    public:
        explicit ManyLightsDistribution(const std::vector<Primitive>& primitives);
        vector3f sample(vector3f point, Engine &rng) const;
        float pdf(vector3f point, vector3f direction) const;
    //    void add_distr(const random_distribution_sh_ptr& dist);
    //    [[nodiscard]] size_t size() const;
        size_t size() const;
    private:
        std::vector<Primitive> objects;
        std::vector<LightDistribution> distributions;
        BVH bvh;
    };

    class SceneDistribution {
    public:
        explicit SceneDistribution(const std::vector<Primitive>& light_objects) : light(light_objects) {};
        vector3f sample(const vector3f &point, const vector3f &normal, const vector3f &eye_direction, float roughness2, Engine &rng) const;
        float pdf(const vector3f &point, const vector3f &normal, const vector3f &eye_direction, float roughness2, const vector3f &direction) const;
    private:
        ManyLightsDistribution light;
    };
}
