#include "scene.h"
#include "utils/base.h"

#include <fstream>
#include <sstream>
#include <memory>

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

void Scene::render() {
    for (int j = 0; j < camera_->canvas_.height(); ++j) {
        for (int i = 0; i < camera_->canvas_.width(); ++i) {
            vector2i pix_pos{i, j};
            Ray r = camera_->cast_in_pixel(pix_pos);

            vector3f color = bg_color_;
            float distance = 1e9;

            for (const auto& o : objects_) {
                auto intersection = o->intersect(r);
                if (intersection) {
                    if (intersection.value().distance < distance) {
                        distance = intersection.value().distance;
                        color = o->color_;
                    }
                }
            }
            color = color + ambient_;
            color = aces_tonemap(color);
            color = pow(color, gamma);
            camera_->canvas_.set(pix_pos, normal_to_ch8bit(color));
        }
    }
}

void Scene::draw_into(const std::string &filename) const {
    camera_->canvas_.write_to(filename);
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

        ParseStage stage = get_parse_stage(cmd);
        parse_stages |= stage;
        switch (stage) {
            case UNKNOWN:
                if (!scene.objects_.empty() && !scene.objects_.back()->parse(line) ||
                    !scene.light_.empty() && !scene.light_.back()->parse(line))
                    break;
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
            case NEW_LIGHT:
                scene.light_.emplace_back(new Primitive());
        }
    }
    if (parse_stages != ParseStage::READY) {
        throw std::runtime_error("Wrong file format: " + filename);
    }
    scene.camera_ = std::make_unique<Camera>(cam_canvas, cam_position, cam_axes, cam_fov_x);
}

