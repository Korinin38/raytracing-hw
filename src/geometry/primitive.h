#pragma once

#include "utils/base.h"

#include <fstream>
#include <optional>
#include <memory>

class Primitive;
class Ray;

typedef std::shared_ptr<Primitive> primitive_sh_ptr;

typedef struct {
    float distance;
    vector3f normal;
    bool inside;
} intersection;

class Primitive {
public:
    enum Type {
        Box,
        Ellipsoid,
        Plane
    };

    Type     type_;
    vector3f param_;
    vector3f position_ = {0, 0, 0};
    vector4f rotation_ = {0, 0, 0, 1};

    vector3f color_ = {0, 0, 0};

    Primitive() = default;
    bool parse(const std::string& line);
    [[nodiscard]] std::optional<intersection> intersect(Ray ray) const;

protected:
    // shift and rotate Ray to make itself behave like axis-aligned, in origin
    void translateRay(Ray &ray) const;
};

class Ray {
public:
    Ray(vector3f p, vector3f d);
    vector3f position;
    vector3f direction;
};