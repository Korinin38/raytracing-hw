#include "random.h"

#include <cmath>

namespace rng
{

RandomGenerator get_generator() {
    return {std::minstd_rand()};
}

static uniform_float_d angle_d(0, 2 * M_PI);
static uniform_float_d height_d(-1.f, 1.f);
//static normal_d n_d(-1.f, 1.f);
//static normal_d height_d(-1.f, 1.f);

vector3f get_sphere(RandomGenerator rng) {

    float theta = angle_d(rng);
    float height = height_d(rng);

//    float x = height_d(rng);
//    float y = height_d(rng);
//    float z = height_d(rng);
//
//    while (x * x + y * y + z * z > 1) {
//        x = height_d(rng);
//        y = height_d(rng);
//        z = height_d(rng);
//    }
//
//    return {x, y, z};

    float proj = std::sqrt(1 - height * height);

    return {proj * std::cos(theta), proj * std::sin(theta), height};
}

}