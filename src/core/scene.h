#pragma once

#include "geometry/primitive.h"
#include "core/camera.h"
#include "utils/base.h"

#include <string>
#include <vector>

class Scene {
public:
    Scene(const std::string &filename);

    camera_uniq_ptr camera_;
    vector3f bg_color_;
    std::vector<primitive_uniq_ptr> objects_;
private:
    class SceneParser;
};