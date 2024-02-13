#pragma once

#include "utils/base.h"

#include <memory>

class Camera;

typedef std::unique_ptr<Camera> camera_uniq_ptr;

class Camera {
public:
    Camera(vector2i size, vector3f position, vector3f axes[3], float fov_x);

    // pixel-wise dimensions
    vector2i canvas_size_{};
    vector3f position_{};

    // 0 - right
    // 1 - up
    // 2 - forward
    vector3f axes_[3]{};

    // in radians
    vector2f fov_{};
};
