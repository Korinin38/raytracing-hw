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
    Scene(camera_uniq_ptr &camera_, std::vector<Mesh> meshes_, std::vector<Primitive> objects_,
          std::vector<Texture> textures_, int ray_depth_, int samples_, float max_distance_);

    void render(ProgressFunc = nullptr) const;
    void draw_into(const std::string &filename) const;

    camera_uniq_ptr camera;
    vector3f bg_color;
    int ray_depth;
    int samples;
    std::vector<Mesh> meshes;
    std::vector<Primitive> objects;
    std::vector<Texture> textures;
    BVH bvh;

    Texture hdr;

    Intersection intersect(Ray r, Engine &rng, bool no_light = false) const;
private:
    const float max_distance;
    const float gamma_ = 1.f / 2.2f;
    scene_distribution_sh_ptr scene_distribution_;
};