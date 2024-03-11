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
    virtual vector3f sphere_sample(vector3f point, vector3f normal) = 0;
    // get pdf of direction that we got by sphere_sample
    virtual float pdf(vector3f point, vector3f normal, vector3f direction) = 0;
protected:
    RandomDistribution();
    Engine rng;
};

class UniformDistribution : public RandomDistribution {
public:
    UniformDistribution(float min, float max);
    float sample();
    float norm_sample();
    vector3f sphere_sample(vector3f point, vector3f normal) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
    uniform_float_d dist;
    normal_d norm_dist;
};

class CosineWeightedDistribution : public RandomDistribution {
public:
    CosineWeightedDistribution();
    vector3f sphere_sample(vector3f point, vector3f normal) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
private:
    UniformDistribution uni_dist;
};

class Primitive;

class LightDistribution : public RandomDistribution {
public:
    LightDistribution(const Primitive &object);
    vector3f sphere_sample(vector3f point, vector3f normal) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
private:
    const Primitive *primitive;
    UniformDistribution uni_dist;
};

class MixedDistribution : public RandomDistribution {
public:
    explicit MixedDistribution();
    vector3f sphere_sample(vector3f point, vector3f normal) override;
    float pdf(vector3f point, vector3f normal, vector3f direction) override;
    void add_dist(const random_distribution_sh_ptr& dist);
private:
    std::vector<random_distribution_sh_ptr> distributions;
    UniformDistribution uni;
};

namespace rng {
Engine get_generator();
//vector3f get_sphere(Engine rng);
}
