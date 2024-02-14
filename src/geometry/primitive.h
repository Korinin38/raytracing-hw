#pragma once

#include "utils/base.h"

#include <fstream>
#include <optional>
#include <memory>

class Primitive;
class Ray;

typedef std::shared_ptr<Primitive> primitive_sh_ptr;

class Primitive {
public:
    vector3f position_ = {0, 0, 0};
    vector4f rotation_ = {0, 0, 0, 1};

    vector3f color_ = {};

    virtual void parse(std::ifstream &in);
    [[nodiscard]] virtual std::optional<float> intersects(Ray ray) const = 0;

protected:
    // shift and rotate Ray to make itself behave like axis-aligned, in origin
    void translateRay(Ray &ray) const;

private:
    enum ParseStage {
        COLOR = 0b1,
        READY = 0b1,
    };
};

class Ray {
public:
    Ray(vector3f p, vector3f d);
    // start position
    vector3f position;
    // direction
    vector3f direction;
};