#pragma once

#include <random>
#include "vector.h"
#include "pcg/include/pcg_random.hpp"

typedef std::uniform_int_distribution<int> uniform_int_d;
typedef std::uniform_real_distribution<float> uniform_float_d;
typedef std::normal_distribution<float> normal_d;

typedef pcg32 RandomGenerator;

namespace rng {
RandomGenerator get_generator();
vector3f get_sphere(RandomGenerator rng);
}
