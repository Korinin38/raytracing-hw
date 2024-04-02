#pragma once

#include <random>
#include "vector.h"
#include <vector>
#include <memory>

typedef std::uniform_int_distribution<int> uniform_int_d;
typedef std::uniform_real_distribution<float> uniform_float_d;
typedef std::normal_distribution<float> normal_d;

//typedef pcg32 RandomGenerator;
typedef std::minstd_rand Engine;

class RandomDistribution;
class MixedDistribution;
typedef std::shared_ptr<RandomDistribution> random_distribution_sh_ptr;
typedef std::shared_ptr<MixedDistribution> mixed_distribution_sh_ptr;

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

class LightDistribution : public RandomDistribution {
public:
    explicit LightDistribution(const Primitive &object);
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
private:
    const Primitive *primitive;
};

class MixedDistribution : public RandomDistribution {
public:
    vector3f sample(vector3f point, vector3f normal, Engine &rng) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
    void add_distr(const random_distribution_sh_ptr& dist);
    [[nodiscard]] size_t get_size() const;
    ~MixedDistribution();
private:
    std::vector<random_distribution_sh_ptr> distributions;
    std::vector<int> count;
};

namespace rng {
Engine get_generator(size_t seed = 0);
//vector3f get_sphere(Engine rng);
}
