#pragma once

#include <core/camera.h>
#include <core/bvh.h>
#include <geometry/primitive.h>
#include <utils/base.h>
#include <utils/timer.h>

#include <string>
#include <vector>

typedef std::function<bool(int, timer&)> ProgressFunc;

struct Scene {
public:
    Scene(camera_uniq_ptr &camera_, std::vector<primitive_sh_ptr> objects, vector3f bg_color, int ray_depth, int samples, float max_distance);

    void render(ProgressFunc = nullptr) const;
    void draw_into(const std::string &filename) const;

    camera_uniq_ptr camera;
    vector3f bg_color;
    int ray_depth;
    int samples;
    std::vector<primitive_sh_ptr> objects;
    BVH bvh;

    Intersection intersect(Ray r, Engine &rng, bool no_light = false) const;
private:
    const float max_distance;
    const float gamma_ = 1.f / 2.2f;
    mutable random_distribution_sh_ptr random_distributions_;
};