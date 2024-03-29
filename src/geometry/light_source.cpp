#include "light_source.h"
#include <sstream>

bool LightSource::parse(const std::string &line) {
    std::stringstream ss(line);
    std::string cmd;

    ss >> cmd;

    if (cmd == "LIGHT_INTENSITY") {
		intensity = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "LIGHT_DIRECTION") {
		type = Directional;
		direction = normal(vec3f_from_string(line, cmd.length() + 1));
        return true;
    } else if (cmd == "LIGHT_POSITION") {
		type = Positional;
		position = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "LIGHT_ATTENUATION") {
		type = Positional;
		attenuation = vec3f_from_string(line, cmd.length() + 1);
        return true;
    }
    return false;
}
