#pragma once

#include <core/camera.h>
#include <core/bvh.h>
#include <geometry/primitive.h>
#include <geometry/light_source.h>
#include <utils/base.h>

#include <string>
#include <vector>

typedef bool (*ProgressFunc)(int progress, void *userData);

struct Scene {
public:
    Scene(const std::string &filename);

    void render(ProgressFunc = nullptr) const;
    void draw_into(const std::string &filename) const;

    camera_uniq_ptr camera;
    vector3f bg_color{};
    int ray_depth = 1;
    int samples = 16;
    vector3f ambient{0.f, 0.f, 0.f};
    std::vector<primitive_sh_ptr> objects;
    BVH bvh;
    std::vector<light_source_sh_ptr> light;

    Intersection intersect(Ray r, Engine &rng, float max_distance = 1e9, bool no_light = false) const;
private:
    const float gamma_ = 1.f / 2.2f;
    class SceneParser;
    mutable MixedDistribution random_distributions_;
};