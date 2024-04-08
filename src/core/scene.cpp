#include "scene.h"
#include "utils/base.h"
#include "utils/timer.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cmath>

const float step = 1e-4;

class Scene::SceneParser {
public:
    SceneParser(Scene &scene, const std::string &filename);
};

enum ParseStage {
    UNKNOWN =           0,
    AMBIENT_LIGHT =     5,
    RAY_DEPTH =         6,
    NEW_PRIMITIVE =     7,
    NEW_LIGHT =         9,
    SAMPLES =           10,

    DIMENSIONS =        (1 << 0),
    BG_COLOR =          (1 << 1),
    CAMERA_POSITION =   (1 << 2),
    CAMERA_RIGHT =      (1 << 3),
    CAMERA_UP =         (1 << 4),
    CAMERA_FORWARD =    (1 << 5),
    CAMERA_FOV_X =      (1 << 6),
    READY =             (1 << 7) - 1
};

inline ParseStage get_parse_stage(const std::string& cmd) {
    if (cmd == "DIMENSIONS")        return DIMENSIONS; else
    if (cmd == "BG_COLOR")          return BG_COLOR; else
    if (cmd == "CAMERA_POSITION")   return CAMERA_POSITION; else
    if (cmd == "CAMERA_RIGHT")      return CAMERA_RIGHT; else
    if (cmd == "CAMERA_UP")         return CAMERA_UP; else
    if (cmd == "CAMERA_FORWARD")    return CAMERA_FORWARD; else
    if (cmd == "CAMERA_FOV_X")      return CAMERA_FOV_X; else
    if (cmd == "AMBIENT_LIGHT")     return AMBIENT_LIGHT; else
    if (cmd == "RAY_DEPTH")         return RAY_DEPTH; else
    if (cmd == "SAMPLES")           return SAMPLES; else
    if (cmd == "NEW_PRIMITIVE")     return NEW_PRIMITIVE; else
    if (cmd == "NEW_LIGHT")         return NEW_LIGHT; else
                                    return UNKNOWN;
}

Scene::Scene(const std::string &filename) {
    (SceneParser(*this, filename));

    mixed_distribution_sh_ptr light_distr = std::make_shared<MixedDistribution>();
    for (const auto& o : objects) {
        if (o->emissive() && o->type != Primitive::Plane)
            light_distr->add_distr(std::make_shared<LightDistribution>(*o));
    }

    bvh.buildBVH(objects);

    random_distributions_.add_distr(std::make_shared<CosineWeightedDistribution>());
//    random_distributions_.add_distr(std::make_shared<UniformDistribution>());
    if (light_distr->get_size() > 0)
        random_distributions_.add_distr(light_distr);
}

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

    #pragma omp parallel for shared(offset, sample_canvas) schedule(guided, 16) collapse(2)
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
            intersection.color = objects[intersection.object_id]->emission;
        }
    }

    auto intersect_bvh = bvh.intersect(objects, r);

    if (intersect_bvh && intersect_bvh.distance < intersection.distance) {
        intersection = intersect_bvh;
        intersection.color = objects[intersection.object_id]->emission;
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

        switch(obj.material) {
            case Primitive::Diffuse: {
                vector3f dir{};
                float pdf = 0.f;
                float cos;
                dir = random_distributions_.sample(pos, intersection.normal, rng);
                cos = dot(dir, intersection.normal);
                if (cos <= 0.f)
                    break;
                pdf = random_distributions_.pdf(pos, intersection.normal, dir);
                if (pdf <= 0.f || isnanf(pdf))
                    break;
                Ray reflect_ray(pos + dir * step, dir);
                reflect_ray.power = r.power;
                auto reflect_inter = intersect(reflect_ray, rng, max_distance);
                float coeff = M_1_PIf32 / pdf;
                intersection.color += obj.color * coeff * reflect_inter.color * cos;
                break;
            }
            case Primitive::Dielectric: {
                float cos_in = -dot(intersection.normal, r.direction);
                float eta_1 = 1.f;
                float eta_2 = obj.ior;
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
                            refract_color *= obj.color;
                    } else {
                        refract_color = bg_color;
                    }
                    intersection.color += refract_color;
                }
                break;
            }
            case Primitive::Metallic: {
                vector3f dir = r.direction - 2 * intersection.normal * dot(intersection.normal, r.direction);
                normalize(dir);
                Ray reflect_ray(pos + dir * step, dir);
                reflect_ray.power = r.power;
                auto reflect_inter = intersect(reflect_ray, rng, max_distance);
                intersection.color += obj.color * reflect_inter.color;
                break;
            }
        }
    }
    return intersection;
}

Scene::SceneParser::SceneParser(Scene &scene, const std::string &filename) {
    vector2i cam_canvas{};
    vector3f cam_position{};
    vector3f cam_axes[3];
    float    cam_fov_x;

    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("File open error");

    int parse_stages = 0;
    std::string line;

    // scene parameters
    while (std::getline(in, line)) {
        std::stringstream ss(line);

        std::string cmd;
        ss >> cmd;
        if (cmd.empty())
            continue;

        ParseStage stage = get_parse_stage(cmd);
        parse_stages |= stage;
        switch (stage) {
            case UNKNOWN:
                if (!scene.objects.empty() && scene.objects.back()->parse(line))
                    break;
                if (!scene.light.empty() && scene.light.back()->parse(line))
                    break;

                std::cout << "Warning: unknown command: " << cmd << std::endl;
                break;
            case DIMENSIONS:
                parse_stages |= DIMENSIONS;
                cam_canvas = vec2i_from_string(line, cmd.length() + 1);
                break;
            case BG_COLOR:
                parse_stages |= BG_COLOR;
                scene.bg_color = vec3f_from_string(line, cmd.length() + 1);
                break;
            case CAMERA_POSITION:
                parse_stages |= CAMERA_POSITION;
                cam_position = vec3f_from_string(line, cmd.length() + 1);
                break;
            case CAMERA_RIGHT:
                parse_stages |= CAMERA_RIGHT;
                cam_axes[0] = normal(vec3f_from_string(line, cmd.length() + 1));
                break;
            case CAMERA_UP:
                parse_stages |= CAMERA_UP;
                cam_axes[1] = normal(vec3f_from_string(line, cmd.length() + 1));
                break;
            case CAMERA_FORWARD:
                parse_stages |= CAMERA_FORWARD;
                cam_axes[2] = normal(vec3f_from_string(line, cmd.length() + 1));
                break;
            case CAMERA_FOV_X:
                parse_stages |= ParseStage::CAMERA_FOV_X;
                cam_fov_x = float_from_string(line, cmd.length() + 1);
                break;
            case AMBIENT_LIGHT:
                scene.ambient = vec3f_from_string(line, cmd.length() + 1);
                break;
            case RAY_DEPTH:
                scene.ray_depth = int_from_string(line, cmd.length() + 1);
                break;
            case SAMPLES:
                scene.samples = int_from_string(line, cmd.length() + 1);
                break;
            case NEW_PRIMITIVE:
                scene.objects.emplace_back(new Primitive());
                break;
            case NEW_LIGHT:
                scene.light.emplace_back(new LightSource());
                break;
            case READY:
                break;
        }
    }
    if (parse_stages != ParseStage::READY) {
        throw std::runtime_error("Wrong file format: " + filename);
    }
    scene.camera = std::make_unique<Camera>(cam_canvas, cam_position, cam_axes, cam_fov_x);
}
