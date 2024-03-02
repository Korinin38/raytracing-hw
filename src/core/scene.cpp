#include "scene.h"
#include "utils/base.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cmath>

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

    DIMENSIONS =        (1 << 0),
    BG_COLOR =          (1 << 1),
    CAMERA_POSITION =   (1 << 2),
    CAMERA_RIGHT =      (1 << 3),
    CAMERA_UP =         (1 << 4),
    CAMERA_FORWARD =    (1 << 5),
    CAMERA_FOV_X =      (1 << 6),
    READY =             (1 << 7) - 1
};

ParseStage get_parse_stage(std::string cmd) {
    if (cmd == "DIMENSIONS")        return DIMENSIONS; else
    if (cmd == "BG_COLOR")          return BG_COLOR; else
    if (cmd == "CAMERA_POSITION")   return CAMERA_POSITION; else
    if (cmd == "CAMERA_RIGHT")      return CAMERA_RIGHT; else
    if (cmd == "CAMERA_UP")         return CAMERA_UP; else
    if (cmd == "CAMERA_FORWARD")    return CAMERA_FORWARD; else
    if (cmd == "CAMERA_FOV_X")      return CAMERA_FOV_X; else
    if (cmd == "AMBIENT_LIGHT")     return AMBIENT_LIGHT; else
    if (cmd == "RAY_DEPTH")         return RAY_DEPTH; else
    if (cmd == "NEW_PRIMITIVE")     return NEW_PRIMITIVE; else
    if (cmd == "NEW_LIGHT")         return NEW_LIGHT; else
                                    return UNKNOWN;
}

Scene::Scene(const std::string &filename) {
    (SceneParser(*this, filename));
}

void Scene::render() const {
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < camera_->canvas_.height(); ++j) {
        for (int i = 0; i < camera_->canvas_.width(); ++i) {
            vector2i pix_pos{i, j};
            Ray r = camera_->cast_in_pixel(pix_pos);
            r.power = ray_depth_;
            vector3f color = bg_color_;

            if (i == 607 && j == 72) {
                int a = 1;
            }

            auto intersection = intersect(r);
            if (intersection) {
                color = intersection->color;
            }
            color = aces_tonemap(color);
            color = pow(color, gamma);
            camera_->canvas_.set(pix_pos, normal_to_ch8bit(color));
        }
    }
}

void Scene::draw_into(const std::string &filename) const {
    camera_->canvas_.write_to(filename);
}

std::optional<Intersection> Scene::intersect(Ray r, float max_distance, bool no_color) const {
    if (r.power <= 0) {
        Intersection blackness;
        blackness.color = {0.f};
        blackness.distance = max_distance;
        blackness.inside = false;
        blackness.normal = {1.f, 0.f, 0.f};
        return std::make_optional(blackness);
//        return {};
    }
    r.power -= 1;
    Intersection intersection;
    intersection.distance = max_distance;
    intersection.color = {0.f};
    intersection.inside = false;

    bool intersection_happened = false;
    for (const auto& o : objects_) {
        auto intersect_obj = o->intersect(r);
        if (!intersect_obj) continue;

        if (intersect_obj->distance > intersection.distance) {
            continue;
        }
        intersection_happened = true;
        intersection.inside = intersect_obj->inside;
        intersection.distance = intersect_obj->distance;
        intersection.normal = intersect_obj->normal;
        intersection.color = {0.f};

        // simple check with no light or color
        if (no_color)
            continue;

        vector3f pos = r.position + r.direction * intersect_obj->distance;

        switch(o->material_) {
            case Primitive::Diffuse:
                intersection.color = o->color_ * ambient_;
                break;
            case Primitive::Dielectric: {
                float cos_in = -dot(intersect_obj->normal, r.direction);
                float eta_1 = 1.f;
                float eta_2 = o->ior_;
                if (intersect_obj->inside)
                    std::swap(eta_1, eta_2);
                float refractive_index = eta_1 / eta_2;
                float sin_out = refractive_index * std::sqrt(1 - cos_in * cos_in);

                float reflection_coefficient;
                if (sin_out >= 1.f) {
                    reflection_coefficient = 0;
                } else {
                    float r0 = std::pow((eta_1 - eta_2) / (eta_1 + eta_2), 2.f);
                    reflection_coefficient = r0 + (1 - r0) * std::pow(1 - cos_in, 5.f);
                }

                // reflected
                {
                    vector3f dir = r.direction - 2 * intersect_obj->normal * dot(intersect_obj->normal, r.direction);
                    normalize(dir);
                    Ray reflect_ray(pos + dir * 1e-4, dir);
                    reflect_ray.power = r.power;
                    auto reflect_inter = intersect(reflect_ray, max_distance);
                    vector3f reflect_color{};
                    if (!reflect_inter) {
                         reflect_color = reflection_coefficient * bg_color_;
                    } else {
                        reflect_color = reflection_coefficient * reflect_inter->color;
                    }
                    intersection.color += reflect_color;
                }
                if (sin_out >= 1.f)
                    continue;
                // refracted
                {
                    // todo: still not ideal
                    float cos_out = std::sqrt(1 - sin_out * sin_out);
                    float coeff = eta_1 / eta_2;

                    vector3f dir = coeff * r.direction + (coeff * cos_in - cos_out) * intersect_obj->normal;
                    normalize(dir);
                    Ray reflect_ray(pos + dir * 1e-4, dir);
                    reflect_ray.power = r.power;
                    auto refract_inter = intersect(reflect_ray, max_distance);

                    vector3f refract_color{};
                    if (refract_inter) {
                        refract_color = (1 - reflection_coefficient) * refract_inter->color;
                        if (intersect_obj->inside)
                            refract_color *= o->color_;
                    } else {
                        if (refract_inter->inside)
                            refract_color = (1 - reflection_coefficient) * bg_color_;
                        else
                            refract_color = {0.f};
                    }
                    intersection.color += refract_color;
                }
                continue;
            }
            case Primitive::Metallic: {
                vector3f dir = r.direction - 2 * intersect_obj->normal * dot(intersect_obj->normal, r.direction);
                normalize(dir);
                Ray reflect_ray(pos + dir * 1e-4, dir);
                reflect_ray.power = r.power;
                auto reflect_inter = intersect(reflect_ray, max_distance);
                if (!reflect_inter) {
                    intersection.color += o->color_ * bg_color_;
                } else {
                    intersection.color += o->color_ * reflect_inter->color;
                }
                continue;
            }
        }

        for (const auto& l : light_) {
            vector3f dir{};
            float light_distance = max_distance;

            if (l->type_ == LightSource::Directional) {
                dir = l->direction_;
            } else {
                dir = l->position_ - pos;
                light_distance = length(dir);
                normalize(dir);
            }
            Ray light_ray(pos + dir * 1e-4, dir);

            auto light_inter = intersect(light_ray, light_distance, true);
            if (light_inter)
                continue;
            vector3f light_color = l->intensity_;
            if (l->type_ == LightSource::Positional) {
                light_color *= (1.f / (l->attenuation_[0]
                                                    + l->attenuation_[1] * light_distance
                                                    + l->attenuation_[2] * light_distance * light_distance));
            }

            light_color *= std::max(0.f, dot(dir, intersect_obj->normal));
            intersection.color += o->color_ * light_color;
        }
    }
    if (intersection_happened)
        return std::make_optional(intersection);
    else
        return {};
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
                if (!scene.objects_.empty() && scene.objects_.back()->parse(line))
                    break;
                if (!scene.light_.empty() && scene.light_.back()->parse(line))
                    break;

                std::cout << "Warning: unknown command: " << cmd << std::endl;
                break;
            case DIMENSIONS:
                parse_stages |= DIMENSIONS;
                cam_canvas = vec2i_from_string(line, cmd.length() + 1);
                break;
            case BG_COLOR:
                parse_stages |= BG_COLOR;
                scene.bg_color_ = vec3f_from_string(line, cmd.length() + 1);
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
                scene.ambient_ = vec3f_from_string(line, cmd.length() + 1);
                break;
            case RAY_DEPTH:
                scene.ray_depth_ = int_from_string(line, cmd.length() + 1);
                break;
            case NEW_PRIMITIVE:
                scene.objects_.emplace_back(new Primitive());
                break;
            case NEW_LIGHT:
                scene.light_.emplace_back(new LightSource());
                break;
            case READY:
                break;
        }
    }
    if (parse_stages != ParseStage::READY) {
        throw std::runtime_error("Wrong file format: " + filename);
    }
    scene.camera_ = std::make_unique<Camera>(cam_canvas, cam_position, cam_axes, cam_fov_x);
}
