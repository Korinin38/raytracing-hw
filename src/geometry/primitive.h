#pragma once

#include "utils/base.h"

#include <fstream>
#include <optional>
#include <memory>
#include <variant>

class Primitive;
class Ray;

typedef std::shared_ptr<Primitive> primitive_sh_ptr;

typedef struct {
    bool successful = false;
    float distance;
    vector3f normal;
    vector3f color;
    bool inside = false;
    operator bool() const;
} Intersection;

class Primitive {
public:
    enum GeomType {
        Box,
        Ellipsoid,
        Plane,
        Triangle
    };

    enum Material {
        Diffuse,
        Dielectric,
        Metallic
    };

    GeomType type_ = Plane;
    Material material_ = Diffuse;
    float ior_ = 1.0;
    std::variant<vector3f, std::tuple<vector3f, vector3f, vector3f>> param_ = {};
    vector3f position_ = {};
    vector4f rotation_ = {0, 0, 0,1};

    vector3f color_ = {0, 0, 0};
    vector3f emission_ = {0, 0, 0};

    Primitive() = default;
    bool parse(const std::string& line);
    bool emissive() const;
    [[nodiscard]] Intersection intersect(Ray ray) const;

    vector3f to_global(vector3f local) const;
    vector3f to_local(vector3f global) const;
protected:
    // shift and rotate Ray to make itself behave like axis-aligned, in origin
    void translateRay(Ray &ray) const;
};

class Ray {
public:
    Ray(vector3f p, vector3f d);
    vector3f position;
    vector3f direction;

    // amount of jumps possible
    int power = 1;
};