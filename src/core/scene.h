#pragma once

#include "geometry/primitive.h"
#include "geometry/light_source.h"
#include "core/camera.h"
#include "utils/base.h"

#include <string>
#include <vector>

typedef bool (*ProgressFunc)(int progress, void *userData);

class Scene {
public:
    Scene(const std::string &filename);

    void render(ProgressFunc = nullptr) const;
    void draw_into(const std::string &filename) const;

    camera_uniq_ptr camera_;
    vector3f bg_color_;
    int ray_depth_ = 1;
    int samples_ = 16;
    vector3f ambient_{0.f, 0.f, 0.f};
    std::vector<primitive_sh_ptr> objects_;
    std::vector<light_source_sh_ptr> light_;

    Intersection intersect(Ray r, float max_distance = 1e9, bool no_light = false) const;
private:
    const float gamma = 1.f / 2.2f;
    const float step = 4e-4;
    class SceneParser;
};