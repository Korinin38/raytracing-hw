#pragma once

#include <core/camera.h>
#include <core/bvh.h>
#include <geometry/primitive.h>
#include <geometry/light_source.h>
#include <utils/base.h>
#include <utils/timer.h>

#include <string>
#include <vector>

typedef std::function<bool(int, timer&)> ProgressFunc;

struct Scene {
public:
    Scene(camera_uniq_ptr &camera_, std::vector<primitive_sh_ptr> objects, vector3f bg_color = {}, int ray_depth = 6, int samples = 256, vector3f ambient = {});

    void render(ProgressFunc = nullptr) const;
    void draw_into(const std::string &filename) const;

    camera_uniq_ptr camera;
    vector3f bg_color;
    int ray_depth;
    int samples;
    vector3f ambient;
    std::vector<primitive_sh_ptr> objects;
    BVH bvh;
    std::vector<light_source_sh_ptr> light;

    Intersection intersect(Ray r, Engine &rng, float max_distance = 1e9, bool no_light = false) const;
private:
    const float gamma_ = 1.f / 2.2f;
    mutable random_distribution_sh_ptr random_distributions_;
};