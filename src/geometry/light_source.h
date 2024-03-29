#pragma once

#include "utils/base.h"

#include <memory>

struct LightSource;

typedef std::shared_ptr<LightSource> light_source_sh_ptr;

struct LightSource {
public:
    enum Type {
        Directional,
        Positional
    };
    Type type;
    vector3f intensity;
    vector3f direction;
    vector3f position;
    vector3f attenuation;

    bool parse(const std::string& line);
};
