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

class RandomDistribution;
class MixedDistribution;
class SceneDistribution;
typedef std::shared_ptr<RandomDistribution> random_distribution_sh_ptr;
typedef std::shared_ptr<MixedDistribution> mixed_distribution_sh_ptr;
typedef std::shared_ptr<SceneDistribution> scene_distribution_sh_ptr;

class RandomDistribution {
public:
    // get direction in a hemisphere defined by normal
    virtual vector3f sample(vector3f point, vector3f normal, Engine &rng) = 0;
    // get pdf of direction that we got by sample
    virtual float pdf(vector3f point, vector3f normal, vector3f direction) = 0;
};

class UniformDistribution : public RandomDistribution {
public:
    static float sample(Engine &rng);
    static float norm_sample(Engine &rng);
    static vector3f uni_sphere_sample(Engine &rng);
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
};

class CosineWeightedDistribution : public RandomDistribution {
public:
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
};

class Primitive;
class Intersection;

class LightDistribution : public RandomDistribution {
public:
    explicit LightDistribution(const Primitive &object);
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
    inline float pdfWithHint(vector3f point, vector3f direction, Intersection hint);
private:
    const Primitive *primitive;
};

class MixedDistribution : public RandomDistribution {
public:
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
    void add_distr(const random_distribution_sh_ptr& dist);
    [[nodiscard]] size_t get_size() const;
private:
    std::vector<random_distribution_sh_ptr> distributions;
};

class ManyLightsDistribution : public RandomDistribution {
public:
    explicit ManyLightsDistribution(const std::vector<primitive_sh_ptr>& primitives);
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
//    void add_distr(const random_distribution_sh_ptr& dist);
//    [[nodiscard]] size_t get_size() const;
private:
    std::vector<primitive_sh_ptr> objects;
    std::vector<LightDistribution> distributions;
    BVH bvh;
};

class SceneDistribution : public RandomDistribution {
public:
    explicit SceneDistribution(const std::vector<primitive_sh_ptr>& light_objects) : light(light_objects) {};
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
private:
    ManyLightsDistribution light;
    CosineWeightedDistribution cosine;
};

namespace rng {
Engine get_generator(size_t seed = 0);
//vector3f get_sphere(Engine rng);
}
