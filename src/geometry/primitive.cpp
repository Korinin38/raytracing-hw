#include "primitive.h"

#include <string>

void Primitive::parse(std::ifstream &in) {
    std::string line;
    int parse_stages = 0;
    while (std::getline(in, line)) {
        size_t space_pos = line.find(' ');
        std::string cmd;

        if (space_pos == std::string::npos)
            cmd = line;
        else
            cmd = line.substr(0, space_pos);

        if (cmd == "POSITION") {
            position_ = vec3f_from_string(line, space_pos + 1);
        } else if (cmd == "ROTATION") {
            rotation_ = normal(vec4f_from_string(line, space_pos + 1));
        } else if (cmd == "COLOR") {
            parse_stages |= ParseStage::COLOR;
            color_ = vec3f_from_string(line, space_pos + 1);
        } else if (cmd == "NEW_PRIMITIVE") {
            break;
        }
    }
    if (parse_stages != ParseStage::READY) {
        throw std::runtime_error("Invalid primitive structure: required values missing");
    }
}

void Primitive::translateRay(Ray &ray) const {
    for (int i = 0; i < 3; ++i)
        ray.position[i] -= position_[i];

//    ray.direction = rotate(ray.direction, rotation_);
}

Ray::Ray(vector3f p, vector3f d) : position(p), direction(d) {
    normalize(direction);
}
