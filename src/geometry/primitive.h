#pragma once

#include "utils/base.h"

#include <fstream>
#include <memory>

class Primitive;

typedef std::unique_ptr<Primitive> primitive_uniq_ptr;

class Primitive {
public:
    vector3f position_ = {0, 0, 0};
    vector4f rotation_ = {0, 0, 0, 1};

    vector3f color_ = {};

    virtual void parse(std::ifstream &in);

private:
    enum ParseStage {
        COLOR = 0b1,
        READY = 0b1,
    };
};

class Ray {
public:
    // start position
    vector3f p;
    // direction
    vector3f d;
};