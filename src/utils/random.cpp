#include "random.h"

#include <cmath>

namespace rng
{

RandomGenerator get_generator() {
    return {pcg_extras::seed_seq_from<std::random_device>()};
}

static uniform_float_d angle_d(0, 2 * M_PI);
static uniform_float_d height_d(-1.f, 1.f);

vector3f get_sphere(RandomGenerator rng) {
    float theta = angle_d(rng);
    float height = height_d(rng);

    float proj = std::sqrt(1 - height * height);

    return {proj * std::cos(theta), proj * std::sin(theta), height};
}

}