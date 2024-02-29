#pragma once

#include "utils/base.h"

#include <memory>

class LightSource;

typedef std::shared_ptr<LightSource> light_source_sh_ptr;

class LightSource {
public:
    enum Type {
        Directional,
        Positional
    };
    Type type_;
    vector3f intensity_;
    vector3f direction_;
    vector3f position_;
    vector3f attenuation_;

    bool parse(const std::string& line);
};
