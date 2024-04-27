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
    std::vector<vector3f> sample_canvas(camera->canvas.height() * camera->canvas.width(), {0.f, 0.f, 0.f});

    const unsigned int canvas_size = camera->canvas.height() * camera->canvas.width();
    int progress = 0;

    timer t;
    if (callback)
        callback(progress, t);
    else
        std::cout << "Render launched." << std::endl;

//    #pragma omp parallel for shared(t, offset, sample_canvas) schedule(guided, 16) collapse(2)
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
                callback(std::floor(100 * progress / (double)canvas_size), t);
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

//    size_t intersected_idx = -1;

    auto intersect_bvh = bvh.intersect(objects, r);

    if (intersect_bvh && intersect_bvh.distance < intersection.distance) {
        intersection = intersect_bvh;
        intersection.color = objects[intersection.object_id]->material.emission;
    }

    // simple check with no light or color
    if (no_color || intersection.object_id >= objects.size())
        return intersection;

    {
        Primitive &obj = *objects[intersection.object_id].get();

        vector3f pos = r.origin + r.direction * intersection.distance;
        vector3f shading_normal = obj.get_shading_normal(intersection.local_coords);
        if (intersection.inside)
            shading_normal = -shading_normal;

        vector3f dir{};
        float pdf = 0.f;
        float cos;
        random_distributions_->update_vndf(obj.material.roughness2, r.direction);
        dir = random_distributions_->sample(pos, shading_normal, rng);
        cos = dot(dir, intersection.normal);
        if (cos <= 0.f)
            return intersection;
        pdf = random_distributions_->pdf(pos, shading_normal, dir);
        if (pdf <= 0.f || isnanf(pdf))
            return intersection;
        Ray reflect_ray(pos + dir * step, dir);
        reflect_ray.power = r.power;
        auto reflect_inter = intersect(reflect_ray, rng);
        float coeff = 1 / pdf;

        vector3f half = normal(dir - r.direction);

        vector3f diffuse_color = lerp(reflect_inter.color, vector3f{0.f, 0.f, 0.f}, obj.material.metallic);
        const float f0_base = 0.04f;
        vector3f f0 = lerp(vector3f{f0_base, f0_base, f0_base}, reflect_inter.color, obj.material.metallic);
        float alpha = std::max(0.03f, obj.material.roughness2);
        vector3f F = f0 + (vector3f{1.f, 1.f, 1.f} - f0) * std::pow(1 - std::abs(dot(-r.direction, half)), 5.f);
        vector3f f_diffuse = (vector3f{1.f, 1.f, 1.f} - F) * diffuse_color * M_1_PI;

        float specular_visibility = smith_joint_masking_shadowing_function(alpha, shading_normal, half, -r.direction, dir)
                                       * (1.f / (4 * std::abs(dot(shading_normal, r.direction)) * std::abs(dot(shading_normal, dir))));
        vector3f f_specular = F * ggx_microfacet_distribution(alpha, shading_normal, half) * specular_visibility;

        intersection.color += obj.material.color * coeff * (f_diffuse + f_specular) * cos * obj.material.alpha;

        return intersection;

        // refraction
        if (!obj.transparent())
            return intersection;
        float cos_in = -dot(shading_normal, r.direction);
        float eta_1 = 1.f;
        float eta_2 = obj.material.ior;
        if (intersection.inside)
            std::swap(eta_1, eta_2);

        float refractive_index = eta_1 / eta_2;
        float sin_out = refractive_index * std::sqrt(1 - cos_in * cos_in);

        if (sin_out >= 1.f) {
            return intersection;
        }
        float cos_out = std::sqrt(1 - sin_out * sin_out);
        float refr_coeff = eta_1 / eta_2;

        vector3f refr_dir = refr_coeff * r.direction + (refr_coeff * cos_in - cos_out) * shading_normal;
        normalize(refr_dir);
        Ray refract_ray(pos + refr_dir * step, refr_dir);
        refract_ray.power = r.power;
        auto refract_inter = intersect(refract_ray, rng);

        vector3f refract_color{};
        if (refract_inter) {
            refract_color = refract_inter.color;
            if (!intersection.inside)
                refract_color *= obj.material.color;
        } else {
            refract_color = bg_color;
        }
        intersection.color += refract_color * (1.f - obj.material.alpha);
    }
    return intersection;
}

Scene::Scene(camera_uniq_ptr &camera_, std::vector<primitive_sh_ptr> objects_, vector3f bg_color_, int ray_depth_, int samples_, float max_distance)
    : camera(std::move(camera_)),
      objects(std::move(objects_)),
      bg_color(bg_color_),
      ray_depth(ray_depth_),
      samples(samples_),
      max_distance(max_distance)
{
    mixed_distribution_sh_ptr light_distr = std::make_shared<MixedDistribution>();
    for (const auto& o : objects) {
        if (o->emissive())
            light_distr->add_distr(std::make_shared<LightDistribution>(*o));
    }

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

    random_distributions_ = std::make_shared<SceneDistribution>(objects);
//    mixed_distribution_sh_ptr m = std::make_shared<MixedDistribution>();
//    m->add_distr(std::make_shared<UniformDistribution>());
//    m->add_distr(std::make_shared<CosineWeightedDistribution>());
//    if (light_distr->size() > 0)
//        m->add_distr(light_distr);
//    random_distributions_ = m;
}
