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
} Intersection;

class Primitive {
public:
    enum GeomType {
        Box,
        Ellipsoid,
        Plane
    };

    enum Material {
        Diffuse,
        Dielectric,
        Metallic
    };

    GeomType type_ = Plane;
    Material material_ = Diffuse;
    float ior_ = 1.0;
    vector3f param_ = {};
    vector3f position_ = {};
    vector4f rotation_ = {0, 0, 0,1};

    vector3f color_ = {0, 0, 0};

    Primitive() = default;
    bool parse(const std::string& line);
    [[nodiscard]] std::optional<Intersection> intersect(Ray ray) const;

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