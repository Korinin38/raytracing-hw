#include "light_source.h"
#include <sstream>

bool LightSource::parse(const std::string &line) {
    std::stringstream ss(line);
    std::string cmd;

    ss >> cmd;

    if (cmd == "LIGHT_INTENSITY") {
        intensity_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "LIGHT_DIRECTION") {
        type_ = Directional;
        direction_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "LIGHT_POSITION") {
        type_ = Positional;
        position_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "LIGHT_ATTENUATION") {
        type_ = Positional;
        attenuation_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    }
    return false;
}
