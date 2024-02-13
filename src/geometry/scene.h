#pragma once
#include "primitive.h"
#include "utils/base.h"

#include <string>
#include <vector>

class Scene {
public:
    Scene(const std::string &filename);

    // pixel-wise dimensions
    int width_;
    int height_;

    vector3f bg_color_;
    vector3f camera_position_;

    // 0 - right
    // 1 - up
    // 2 - forward
    vector3f camera_axes_[3];

    // in radians
    float camera_fov_x_;

//    std::vector<Primitive> objects_;
    std::vector<primitive_uniq_ptr> objects_;
private:
    class SceneParser;
};