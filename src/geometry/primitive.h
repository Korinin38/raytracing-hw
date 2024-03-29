#pragma once

#include "utils/base.h"

#include <fstream>
#include <optional>
#include <memory>
#include <variant>

struct Primitive;
struct Ray;

typedef std::shared_ptr<Primitive> primitive_sh_ptr;

typedef struct {
    bool successful = false;
    float distance;
    vector3f normal;
    vector3f color;
    bool inside = false;
    operator bool() const;
} Intersection;

struct Primitive {
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

    GeomType type = Plane;
    Material material = Diffuse;
    float ior = 1.0;
    std::variant<vector3f, std::tuple<vector3f, vector3f, vector3f>> param_ = {};
    vector3f position = {};
    vector4f rotation = {0, 0, 0, 1};

    vector3f color = {0, 0, 0};
    vector3f emission = {0, 0, 0};

    Primitive() = default;
    bool parse(const std::string& line);
    bool emissive() const;
    [[nodiscard]] Intersection intersect(Ray ray) const;

    vector3f to_global(vector3f local) const;
    vector3f to_local(vector3f global) const;
protected:
    // shift and rotate Ray to make itself behave like axis-aligned, in origin
    void transformRay(Ray &ray) const;
};

struct Ray {
public:
    Ray(vector3f p, vector3f d);
    vector3f position;
    vector3f direction;

    // amount of jumps possible
    int power = 1;
};