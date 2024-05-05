#pragma once

#include <utils/vector.h>
#include <utils/matrix.h>

#include <fstream>
#include <optional>
#include <memory>
#include <variant>
#include <vector>

struct AABB;
struct Material;
struct Primitive;
struct Texture;
struct Ray;

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
    float ior = 1.f;
    float alpha = 1.f;
    vector3f base_color = {1.f, 1.f, 1.f};
    vector3f emission = {0, 0, 0};
    float metallic = 1.f;
    float roughness2 = 1.f; // alpha

    static const int NO_TEXTURE = -1;
    int base_color_i = NO_TEXTURE;
    int normal_i = NO_TEXTURE;
    int metallic_roughness_i = NO_TEXTURE;
    int emission_i = NO_TEXTURE;
};

struct Mesh {
    Material material;
    matrix4d normal_transform;
};

struct Primitive {
public:
    // 0: origin point
    // 1: U
    // 2: V
    vector3f position[3] = {};
    vector3f normal[3] = {};
    vector2f texcoord[3] = {};

    Primitive() = default;
    bool emissive() const;
    bool transparent() const;
    [[nodiscard]] Intersection intersect(Ray ray) const;

    vector3f get_geometric_normal() const;
    vector3f get_shading_normal(vector2f local_coords) const;

    vector2f get_texcoord(vector2f local_coords) const;

    vector3f get_color(vector2f local_coords) const;
    vector3f get_emission(vector2f local_coords) const;
    std::tuple<float, float> get_metallic_roughness(vector2f local_coords) const;

    AABB aabb() const;

    struct CalculationCache {
        vector3f triangle_normal{};
        float triangle_area = 0.f;
        AABB aabb;
    };

    mutable CalculationCache cache;

    // these are stored outside for shared usage
    Texture *textures;
    Mesh *mesh;
    int mesh_id;
};

struct Texture {
    std::vector<uint8_t> data;
    int width;
    int height;
    int channels;
    int bytes_per_channel;

    const unsigned char* get(int x, int y) const {
        x %= width;
        y %= height;
        return data.data() + (y * width + x) * channels * bytes_per_channel;
    }

    static float interpolate_bilinear(float a00, float a01, float a10, float a11, float dx, float dy) {
        return a00 * (1 - dx) * (1 - dy) + a01 * (1 - dx) * dy
             + a10 * dx * (1 - dy) + a11 * dx * dy;
    };
    // p in [0,1]x[0,1] (or, otherwise, transformed by modulo 1)
    // dst is valid buffer that contains at least (channels * bytes_per_channel) bytes of memory
    const vector4f sample(vector2f p) const {
        p.x -= std::floor(p.x);
        p.y -= std::floor(p.y);

        // expand to full canvas
        p.x *= (float)width;
        p.y *= (float)height;

        vector4f res;
        const uint8_t *ptrs = get(std::floor(p.x), std::floor(p.y));
        for (int i = 0; i < channels; ++i) {
            res[i] = ch8bit_to_normal(ptrs[i]);
        }
        return res;
    }

    // interpolates sample
    const vector4f sRGBA_sample(vector2f p) const {
        return interpolate_sample(p, true);
    }

    // interpolates sample
    const vector4f RGBA_sample(vector2f p) const {
        return interpolate_sample(p);
    }

private:
    const vector4f interpolate_sample(vector2f p, bool srgb = false) const {
        // repeat
        p.x -= std::floor(p.x);
        p.y -= std::floor(p.y);

        // expand to full canvas
        p.x *= (float)width;
        p.y *= (float)height;

        vector4f res;

        vector2i point = {(int)std::floor(p.x), (int)std::floor(p.y)};
        vector2f inter_delta = {p.x - std::floor(p.x), p.y - std::floor(p.y)};

        const uint8_t *ptrs[4];
        for (int dy = 0; dy < 2; ++dy)
            for (int dx = 0; dx < 2; ++dx) {
                ptrs[dx * 2 + dy] = get(point.x + dx, point.y + dy);
            }

        for (int i = 0; i < channels; ++i) {
            float vals[4];
            for (int dy = 0; dy < 2; ++dy)
                for (int dx = 0; dx < 2; ++dx) {
                    vals[dx * 2 + dy] = ch8bit_to_normal(ptrs[dx * 2 + dy][i]);
                    if (srgb)
                        vals[dx * 2 + dy] = std::pow(vals[dx * 2 + dy], 2.2f);
                }
            res[i] = interpolate_bilinear(vals[0], vals[1], vals[2], vals[3],
                                          inter_delta.x, inter_delta.y);
        }

        return res;
    }
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

void set_mesh_texture_primitive_correspondence(std::vector<Mesh> &meshes, std::vector<Primitive> &primitives, std::vector<Texture> &textures);