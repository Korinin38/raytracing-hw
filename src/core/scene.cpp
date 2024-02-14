#include "scene.h"
#include "utils/base.h"
#include "geometry/objects.h"

#include <fstream>
#include <memory>

class Scene::SceneParser {
public:
    SceneParser(Scene &scene, const std::string &filename);
};

enum ParseStage {
    DIMENSIONS = (1 << 0),
    BG_COLOR = (1 << 1),
    CAMERA_POSITION = (1 << 2),
    CAMERA_RIGHT = (1 << 3),
    CAMERA_UP = (1 << 4),
    CAMERA_FORWARD = (1 << 5),
    CAMERA_FOV_X = (1 << 6),
    READY = (1 << 7) - 1
};

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
                auto intersection = o->intersects(r);
                if (intersection) {
                    if (intersection.value() < distance) {
                        distance = intersection.value();
                        color = o->color_;
                    }
                }
            }

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
        size_t space_pos = line.find(' ');

        std::string cmd;
        if (space_pos == std::string::npos)
            cmd = line;
        else
            cmd = line.substr(0, space_pos);

//         = line.substr(0, space_pos);

        if (cmd == "DIMENSIONS") {
            parse_stages |= ParseStage::DIMENSIONS;
            cam_canvas = vec2i_from_string(line, space_pos + 1);
        } else if (cmd == "BG_COLOR") {
            parse_stages |= ParseStage::BG_COLOR;
            scene.bg_color_ = vec3f_from_string(line, space_pos + 1);
        } else if (cmd == "CAMERA_POSITION") {
            parse_stages |= ParseStage::CAMERA_POSITION;
            cam_position = vec3f_from_string(line, space_pos + 1);
        } else if (cmd == "CAMERA_RIGHT") {
            parse_stages |= ParseStage::CAMERA_RIGHT;
            cam_axes[0] = normal(vec3f_from_string(line, space_pos + 1));
        } else if (cmd == "CAMERA_UP") {
            parse_stages |= ParseStage::CAMERA_UP;
            cam_axes[1] = normal(vec3f_from_string(line, space_pos + 1));
        } else if (cmd == "CAMERA_FORWARD") {
            parse_stages |= ParseStage::CAMERA_FORWARD;
            cam_axes[2] = normal(vec3f_from_string(line, space_pos + 1));
        } else if (cmd == "CAMERA_FOV_X") {
            parse_stages |= ParseStage::CAMERA_FOV_X;
            cam_fov_x = float_from_string(line, space_pos + 1);
        } else if (cmd == "NEW_PRIMITIVE") {
            break;
        }

    }
    if (parse_stages != ParseStage::READY) {
        throw std::runtime_error("Wrong file format: " + filename);
    }
    scene.camera_ = std::make_unique<Camera>(cam_canvas, cam_position, cam_axes, cam_fov_x);

    // primitive parameters
    while (std::getline(in, line)) {
        size_t space_pos = line.find(' ');
        if (space_pos == std::string::npos)
            continue;

        std::string cmd = line.substr(0, space_pos);

        if (cmd == "PLANE") {
            vector3f n = normal(vec3f_from_string(line, space_pos + 1));
            scene.objects_.emplace_back(new Plane(n));
            scene.objects_.back()->parse(in);
        } else if (cmd == "ELLIPSOID") {
            vector3f radius = vec3f_from_string(line, space_pos + 1);
            scene.objects_.emplace_back(new Ellipsoid(radius));
            scene.objects_.back()->parse(in);
        } else if (cmd == "BOX") {
            vector3f sizes = vec3f_from_string(line, space_pos + 1);
            scene.objects_.emplace_back(new Box(sizes));
            scene.objects_.back()->parse(in);
        }
    }
}

