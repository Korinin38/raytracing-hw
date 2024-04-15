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

    #pragma omp parallel for shared(t, offset, sample_canvas) schedule(guided, 16) collapse(2)
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

Intersection Scene::intersect(Ray r, Engine &rng, float max_distance, bool no_color) const {
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

    for (int i = objects.size() - 1; i >= 0; --i) {
        if (objects[i]->type != Primitive::Plane)
            break;
        auto intersect_plane = objects[i]->intersect(r);
        if (intersect_plane && intersect_plane.distance < intersection.distance) {
            intersection = intersect_plane;
            intersection.object_id = i;
            intersection.color = objects[intersection.object_id]->material.emission;
        }
    }

    auto intersect_bvh = bvh.intersect(objects, r);

    if (intersect_bvh && intersect_bvh.distance < intersection.distance) {
        intersection = intersect_bvh;
        intersection.color = objects[intersection.object_id]->material.emission;
    }

//    intersection.color = objects[intersection.object_id]->emission;

//    for (int i = 0; i < objects.size(); ++i) {
//        auto intersect_obj = objects[i]->intersect(r);
//        if (!intersect_obj) continue;
//
//        if (intersect_obj.distance > intersection.distance) {
//            continue;
//        }
//        intersection.object_id = i;
//        intersection.successful = true;
//        intersection.inside = intersect_obj.inside;
//        intersection.distance = intersect_obj.distance;
//        intersection.normal = intersect_obj.normal;
//        intersection.color = objects[i]->emission;
//    }

    // simple check with no light or color
    if (no_color || intersection.object_id >= objects.size())
        return intersection;

    {
        Primitive &obj = *objects[intersection.object_id].get();

        vector3f pos = r.origin + r.direction * intersection.distance;

        switch(obj.material.type) {
            case Material::Type::Diffuse: {
                vector3f dir{};
                float pdf = 0.f;
                float cos;
                dir = random_distributions_->sample(pos, intersection.normal, rng);
                cos = dot(dir, intersection.normal);
                if (cos <= 0.f)
                    break;
                pdf = random_distributions_->pdf(pos, intersection.normal, dir);
                if (pdf <= 0.f || isnanf(pdf))
                    break;
                Ray reflect_ray(pos + dir * step, dir);
                reflect_ray.power = r.power;
                auto reflect_inter = intersect(reflect_ray, rng, max_distance);
                float coeff = M_1_PIf32 / pdf;
                intersection.color += obj.material.color * coeff * reflect_inter.color * cos;
                break;
            }
            case Material::Type::Dielectric: {
                float cos_in = -dot(intersection.normal, r.direction);
                float eta_1 = 1.f;
                float eta_2 = obj.material.ior;
                if (intersection.inside)
                    std::swap(eta_1, eta_2);
                float refractive_index = eta_1 / eta_2;
                float sin_out = refractive_index * std::sqrt(1 - cos_in * cos_in);

                float direction;
                float reflection_coefficient;
                if (sin_out >= 1.f) {
                    direction = 0;
                    reflection_coefficient = 1;
                } else {
                    float r0 = std::pow((eta_1 - eta_2) / (eta_1 + eta_2), 2.f);
                    reflection_coefficient = r0 + (1 - r0) * std::pow(1 - cos_in, 5.f);
                    uniform_float_d ray_chooser(0.f, 1.f);
                    direction = ray_chooser(rng);
                }

                if (direction < reflection_coefficient)
                {
                    // reflected
                    vector3f dir = r.direction - 2 * intersection.normal * dot(intersection.normal, r.direction);
                    normalize(dir);
                    Ray reflect_ray(pos + dir * step, dir);
                    reflect_ray.power = r.power;
                    auto reflect_inter = intersect(reflect_ray, rng, max_distance);
                    intersection.color += reflect_inter.color;
                } else {
                    // refracted
                    float cos_out = std::sqrt(1 - sin_out * sin_out);
                    float coeff = eta_1 / eta_2;

                    vector3f dir = coeff * r.direction + (coeff * cos_in - cos_out) * intersection.normal;
                    normalize(dir);
                    Ray reflect_ray(pos + dir * step, dir);
                    reflect_ray.power = r.power;
                    auto refract_inter = intersect(reflect_ray, rng, max_distance);

                    vector3f refract_color{};
                    if (refract_inter) {
                        refract_color = refract_inter.color;
                        if (!intersection.inside)
                            refract_color *= obj.material.color;
                    } else {
                        refract_color = bg_color;
                    }
                    intersection.color += refract_color;
                }
                break;
            }
            case Material::Type::Metallic: {
                vector3f dir = r.direction - 2 * intersection.normal * dot(intersection.normal, r.direction);
                normalize(dir);
                Ray reflect_ray(pos + dir * step, dir);
                reflect_ray.power = r.power;
                auto reflect_inter = intersect(reflect_ray, rng, max_distance);
                intersection.color += obj.material.color * reflect_inter.color;
                break;
            }
        }
    }
    return intersection;
}

Scene::Scene(camera_uniq_ptr &camera_, std::vector<primitive_sh_ptr> objects_, vector3f bg_color_, int ray_depth_, int samples_, vector3f ambient_, float max_distance)
    : camera(std::move(camera_)),
      objects(std::move(objects_)),
      bg_color(bg_color_),
      ray_depth(ray_depth_),
      samples(samples_),
      ambient(ambient_),
      max_distance(max_distance)
{
    mixed_distribution_sh_ptr light_distr = std::make_shared<MixedDistribution>();
    for (const auto& o : objects) {
        if (o->emissive() && o->type != Primitive::Plane)
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
////    m->add_distr(std::make_shared<UniformDistribution>());
//    m->add_distr(std::make_shared<CosineWeightedDistribution>());
//    if (light_distr->size() > 0)
//        m->add_distr(light_distr);
//    random_distributions_ = m;
}
