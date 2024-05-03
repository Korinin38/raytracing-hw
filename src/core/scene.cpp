#include "scene.h"
#include <utils/base.h>
#include <utils/timer.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cmath>
#include <utility>

const float step = 1e-4;

void Scene::render(ProgressFunc callback) const {
    uniform_float_d offset(-0.5f, 0.5f);
    const unsigned int canvas_size = camera->canvas.height() * camera->canvas.width();
    std::vector<vector3f> sample_canvas(canvas_size, {0.f, 0.f, 0.f});

    int progress = 0;
    const double progress_coefficient = 100.f / (double) canvas_size;

    timer t;
    if (callback)
        callback(progress, t);
    else
        std::cout << "Render launched." << std::endl;

    #pragma omp parallel for shared(t, offset, sample_canvas, progress_coefficient) schedule(guided, 16) collapse(2)
    for (int j = 0; j < camera->canvas.height(); ++j) {
        for (int i = 0; i < camera->canvas.width(); ++i) {
            Engine rng = rng::get_generator(j * camera->canvas.width() + i);
            for (int s = 0; s < samples; ++s) {
                vector2f pix_offset{offset(rng), offset(rng)};
                vector2i pix_pos{i, j};
                Ray r = camera->cast_in_pixel(pix_pos, pix_offset);
                r.power = ray_depth;

                auto intersection = intersect(r, rng);
                sample_canvas[j * camera->canvas.width() + i] += intersection.color;
            }

            if (callback)
            #pragma omp critical
            {
                progress += 1;
                callback(std::floor(progress * progress_coefficient), t);
            }
        }
    }

    const float normalizer = 1.f / (float)samples;
    #pragma omp parallel for default(none) shared(normalizer, sample_canvas) collapse(2)
    for (int j = 0; j < camera->canvas.height(); ++j) {
        for (int i = 0; i < camera->canvas.width(); ++i) {
            vector3f color = sample_canvas[j * camera->canvas.width() + i];
            color *= normalizer;
            color = aces_tonemap(color);
            color = pow(color, gamma_);
            camera->canvas.set({i, j}, normal_to_ch8bit(color));
        }
    }
}

void Scene::draw_into(const std::string &filename) const {
    camera->canvas.write_to(filename);
}

Intersection Scene::intersect(Ray r, Engine &rng, bool no_color) const {
    if (r.power <= 0) {
        return {};
    }
    r.power -= 1;
    Intersection intersection;
    intersection.successful = false;
    intersection.distance = max_distance;
    intersection.color = bg_color;
    intersection.inside = false;
    intersection.object_id = -1;

    auto intersect_bvh = bvh.intersect(objects, r);

    if (intersect_bvh && intersect_bvh.distance < intersection.distance) {
        intersection = intersect_bvh;
        intersection.color = objects[intersection.object_id]->material.emission;
    }

    // simple check with no light or color
    if (no_color || intersection.object_id >= objects.size())
        return intersection;

    Primitive &obj = *objects[intersection.object_id].get();

    vector3f pos = r.origin + r.direction * intersection.distance;
    vector3f shading_normal = obj.get_shading_normal(intersection.local_coords);
    if (intersection.inside)
        shading_normal = -shading_normal;
//    vector3f shading_normal = intersection.normal;

    vector3f dir{};
    float pdf = 0.f;
    dir = scene_distribution_->sample(pos, shading_normal, -r.direction, obj.material.roughness2, rng);
    if (dot(dir, shading_normal) <= 0.f) {
        if (dot(dir, intersection.normal) <= 0.f) {
            return intersection;
        }
        shading_normal = intersection.normal;
    }
    pdf = scene_distribution_->pdf(pos, shading_normal, -r.direction, obj.material.roughness2, dir);
    if (pdf <= 0.f || isnanf(pdf))
        return intersection;
    Ray reflect_ray(pos + dir * step, dir);
    reflect_ray.power = r.power;
    auto reflect_inter = intersect(reflect_ray, rng);
    float coeff = 1 / pdf;

    vector3f half = normal(dir - r.direction);

    const float f0_base = 0.04f;

    float specular_visibility = smith_joint_masking_shadowing_function(obj.material.roughness2, shading_normal, -r.direction, dir)
                                   * (1.f / (4 * std::abs(dot(shading_normal, r.direction)) * std::abs(dot(shading_normal, dir))));
    float specular_brdf = ggx_microfacet_distribution(obj.material.roughness2, shading_normal, half) * specular_visibility;

    float VdotH = std::abs(dot(-r.direction, half));
    vector3f metal_brdf = specular_brdf * conductor_fresnel(obj.material.color, VdotH);
    vector3f diffuse_brdf = obj.material.color * M_1_PIf32;
    float dielectric_specular_coeff = conductor_fresnel(f0_base, VdotH);
    vector3f dielectric_brdf = diffuse_brdf * (1 - dielectric_specular_coeff) + vector3f{1.f, 1.f, 1.f} * specular_brdf * dielectric_specular_coeff;

    vector3f material = dielectric_brdf * (1 - obj.material.metallic) + metal_brdf * obj.material.metallic;

    intersection.color += reflect_inter.color * coeff * (material) * dot(dir, shading_normal) * obj.material.alpha;

    return intersection;

//    // refraction
//    if (!obj.transparent())
//        return intersection;
//    float cos_in = -dot(shading_normal, r.direction);
//    float eta_1 = 1.f;
//    float eta_2 = obj.material.ior;
//    if (intersection.inside)
//        std::swap(eta_1, eta_2);
//
//    float refractive_index = eta_1 / eta_2;
//    float sin_out = refractive_index * std::sqrt(1 - cos_in * cos_in);
//
//    if (sin_out >= 1.f) {
//        return intersection;
//    }
//    float cos_out = std::sqrt(1 - sin_out * sin_out);
//    float refr_coeff = eta_1 / eta_2;
//
//    vector3f refr_dir = refr_coeff * r.direction + (refr_coeff * cos_in - cos_out) * shading_normal;
//    normalize(refr_dir);
//    Ray refract_ray(pos + refr_dir * step, refr_dir);
//    refract_ray.power = r.power;
//    auto refract_inter = intersect(refract_ray, rng);
//
//    vector3f refract_color{};
//    if (refract_inter) {
//        refract_color = refract_inter.color;
//        if (!intersection.inside)
//            refract_color *= obj.material.color;
//    } else {
//        refract_color = bg_color;
//    }
//    intersection.color += refract_color * (1.f - obj.material.alpha);
//
//    return intersection;
}

Scene::Scene(camera_uniq_ptr &camera_, std::vector<primitive_sh_ptr> objects_, int ray_depth_, int samples_, float max_distance)
    : camera(std::move(camera_)),
      objects(std::move(objects_)),
      bg_color({0.f, 0.f, 0.f}),
      ray_depth(ray_depth_),
      samples(samples_),
      max_distance(max_distance)
{
    bvh.buildBVH(objects);
    std::cout << bvh.nodes.size() << " nodes in BVH." << std::endl;
    size_t max_node_count = 0;
    std::vector<int> node_obj_count(5, 0);
    for (auto n : bvh.nodes) {
        if (n.primitive_count < node_obj_count.size())
            ++node_obj_count[n.primitive_count];
        max_node_count = std::max(n.primitive_count, max_node_count);
    }
    std::cout << "Node object count:" << std::endl;
    for (int i = 0; i < node_obj_count.size(); ++i) {
        std::cout << "\t" << i<< ": " << node_obj_count[i] << std::endl;

    }
    std::cout << "At most " << max_node_count << " objects in single node." << std::endl;

    // NB: as objects with roughness < eps are breaking the view, we clamp it to 0.03f;
    const float ROUGHNESS_LIMIT = 0.03f;
    for (auto &o : objects)
        o->material.roughness2 = std::max(ROUGHNESS_LIMIT, o->material.roughness2);

    scene_distribution_ = std::make_shared<rng::SceneDistribution>(objects);
}
