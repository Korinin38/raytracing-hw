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
    vector2f local_coords;
    vector3f color;
    bool inside = false;
    size_t object_id = -1;
    operator bool() const;
} Intersection;

struct Material {
    enum Type {
        Diffuse,
        Dielectric,
        Metallic
    };

    Type type = Diffuse;
    float ior = 1.0;
    vector3f color = {0, 0, 0};
    vector3f emission = {0, 0, 0};
    float metallic = 1.f;
    float roughness = 1.f;
};

struct Primitive {
public:
    Material material;
    // 0: origin point
    // 1: U
    // 2: V
    vector3f position[3] = {};
    vector3f normal[3] = {};
    vector2f texcoord[3] = {};

    Primitive() = default;
    bool emissive() const;
    [[nodiscard]] Intersection intersect(Ray ray) const;

    vector3f get_geometric_normal() const;
    vector3f get_shading_normal(vector2f local_coords) const;

    vector3f to_global(vector3f local) const;
    vector3f to_local(vector3f global) const;

    AABB aabb() const;

    struct CalculationCache {
        vector3f triangle_normal{};
        float triangle_area = 0.f;
        AABB aabb;
    };
    mutable CalculationCache cache;
};

struct Ray {
public:
    Ray(vector3f p, vector3f d);
    vector3f origin;
    vector3f direction;
    // 1 / direction used for speedup
    vector3f inv_direction;

    // amount of jumps possible
    int power = 1;
};