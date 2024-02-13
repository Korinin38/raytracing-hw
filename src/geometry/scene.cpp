#include "scene.h"
#include "utils/base.h"
#include "objects.h"

#include <fstream>

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

Scene::SceneParser::SceneParser(Scene &scene, const std::string &filename) {
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
            vector2i dims = vec2i_from_string(line, space_pos + 1);
            scene.width_ = dims.x;
            scene.height_ = dims.y;
        } else if (cmd == "BG_COLOR") {
            parse_stages |= ParseStage::BG_COLOR;
            scene.bg_color_ = vec3f_from_string(line, space_pos + 1);
        } else if (cmd == "CAMERA_POSITION") {
            parse_stages |= ParseStage::CAMERA_POSITION;
            scene.camera_position_ = vec3f_from_string(line, space_pos + 1);
        } else if (cmd == "CAMERA_RIGHT") {
            parse_stages |= ParseStage::CAMERA_RIGHT;
            scene.camera_axes_[0] = normal(vec3f_from_string(line, space_pos + 1));
        } else if (cmd == "CAMERA_UP") {
            parse_stages |= ParseStage::CAMERA_UP;
            scene.camera_axes_[1] = normal(vec3f_from_string(line, space_pos + 1));
        } else if (cmd == "CAMERA_FORWARD") {
            parse_stages |= ParseStage::CAMERA_FORWARD;
            scene.camera_axes_[2] = normal(vec3f_from_string(line, space_pos + 1));
        } else if (cmd == "CAMERA_FOV_X") {
            parse_stages |= ParseStage::CAMERA_FOV_X;
            scene.camera_fov_x_ = float_from_string(line, space_pos + 1);
        } else if (cmd == "NEW_PRIMITIVE") {
            break;
        }

    }
    if (parse_stages != ParseStage::READY) {
        throw std::runtime_error("Wrong file format: " + filename);
    }

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
        }
        else if (cmd == "ELLIPSOID") {
            vector3f radius = vec3f_from_string(line, space_pos + 1);
            scene.objects_.emplace_back(new Ellipsoid(radius));
            scene.objects_.back()->parse(in);
        }
        if (cmd == "PLANE") {
            vector3f sizes = vec3f_from_string(line, space_pos + 1);
            scene.objects_.emplace_back(new Box(sizes));
            scene.objects_.back()->parse(in);
        }
    }
}

