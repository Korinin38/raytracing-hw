#pragma once

#include <utils/vector.h>

#include <fstream>
#include <optional>
#include <memory>
#include <variant>

struct AABB;
struct Primitive;
struct Ray;

typedef std::shared_ptr<Primitive> primitive_sh_ptr;

typedef struct AABB {
    vector3f min = get_max_vec3f();
    vector3f max = get_min_vec3f();

    void grow(vector3f p) {
        min = ::min(min, p);
        max = ::max(max, p);
    }

    void grow(AABB aabb) {
        min = ::min(min, aabb.min);
        max = ::max(max, aabb.max);
    }

    [[nodiscard]]
    bool empty() const {
        return (min.x > max.x);
    }

    [[nodiscard]]
    vector3f size() const {
        return max - min;
    }

    [[nodiscard]]
    vector3f center() const {
        return min + size() * 0.5f;
    }

    [[nodiscard]]
    float surface_area() const {
        vector3f size = this->size();
        return 2.f * (size.x * size.y + size.x * size.z + size.y * size.z);
    }

    bool in(vector3f point) const {
        bool a = true;
        for (int i = 0; i < 3; ++i) {
            a &= (point[i] >= min[i]);
            a &= (point[i] <= max[i]);
        }
        return a;
    }

    std::optional<float> intersect(Ray r) const;
} AABB;

typedef struct Intersection {
    bool successful = false;
    float distance;
    vector3f normal;
    vector3f color;
    bool inside = false;
    size_t object_id = -1;
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
    std::variant<vector3f, std::array<vector3f, 3>> param_ = {};
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

    AABB aabb() const;

    struct CalculationCache {
        vector3f triangle_normal{};
        float triangle_area = 0.f;
        vector3f box_inv_size{};
        AABB aabb;
    };
    mutable CalculationCache cache;
protected:
    // shift and rotate Ray to make itself behave like axis-aligned, in origin
    void transformRay(Ray &ray) const;
};

struct Ray {
public:
    Ray(vector3f p, vector3f d);
    vector3f origin;
    vector3f direction;

    // amount of jumps possible
    int power = 1;
};