#include "primitive.h"

#include <vector>
#include <string>
#include <cmath>
#include <sstream>

Intersection::operator bool() const {
    return successful;
}

Ray::Ray(vector3f p, vector3f d) : origin(p), direction(d) {
    normalize(direction);
    inv_direction = vector3f{1.f, 1.f, 1.f} / direction;
}

Intersection Primitive::intersect(Ray ray) const {
    const vector3f &triangle_origin = position[0];
    const vector3f &U = position[1];
    const vector3f &V = position[2];

    vector3f cross_dir_V = (cross(ray.direction, V));
    float det = dot(U, cross_dir_V);
    if (-1e-6 < det && det < 1e-6) {
        return {}; // ray is parallel to the triangle
    }

    float inv_det = 1.f / det;
    vector3f s = ray.origin - triangle_origin;
    float u = inv_det * dot(s, cross_dir_V);
    if (u < 0 || u > 1)
        return {};


    vector3f cross_s_U = cross(s, U);
    float v = inv_det * dot(ray.direction, cross_s_U);
    if (v < 0 || u + v > 1)
        return {};

    float t = inv_det * dot(V, cross_s_U);
    if (t < 0.f)
        return {};

    Intersection intersection;

    intersection.successful = true;
    intersection.color = material.color;
    intersection.local_coords = {u, v};
    intersection.normal = get_geometric_normal();
    intersection.distance = t;
    if (dot(ray.direction, intersection.normal) > 0) {
        intersection.inside = true;
        intersection.normal = -intersection.normal;
    }

    return intersection;
}

bool Primitive::emissive() const {
    return (!material.emission.is_zero());
}

AABB Primitive::aabb() const {
    if (!cache.aabb.empty()) return cache.aabb;

    vector3f points[3];
    points[0] = position[0];
    points[1] = position[0] + position[1];
    points[2] = position[0] + position[2];
    for (auto &p : points) {
        cache.aabb.grow(p);
    }

    return cache.aabb;
}

vector3f Primitive::get_geometric_normal() const {
    if (cache.triangle_normal.is_zero()) {
        cache.triangle_normal = cross(position[1], position[2]);
        cache.triangle_area = 0.5f * length(cache.triangle_normal);
        normalize(cache.triangle_normal);
    }
    return cache.triangle_normal;
}

vector3f Primitive::get_shading_normal(vector2f local_coords) const {
    if (local_coords.x < 0 || local_coords.x > 1 || local_coords.y < 0 || local_coords.x + local_coords.y > 1)
        throw std::runtime_error("Invalid local coordinates");
    return ::normal((1 - local_coords.x - local_coords.y) * normal[0] + local_coords.x * normal[1] + local_coords.y * normal[2]);
}

bool Primitive::transparent() const {
    return (material.alpha < 1.f);
}

std::optional<float> AABB::intersect(Ray r) const {
    // got from https://github.com/erich666/GraphicsGems/blob/master/gems/RayBox.c
    // see: https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
    enum {
        DIMS = 3,
        Left = 0,
        Mid = 1,
        Right = 2
    };
    bool inside = true;
    int quadrant[DIMS];
    float candidatePlane[DIMS];

    for (int i = 0; i < DIMS; ++i) {
        if (r.origin[i] < min[i]) {
            quadrant[i] = Left;
            candidatePlane[i] = min[i];
            inside = false;
        } else if (r.origin[i] > max[i]) {
            quadrant[i] = Right;
            candidatePlane[i] = max[i];
            inside = false;
        } else {
            quadrant[i] = Mid;
        }
    }
    if (inside)
        return std::make_optional(0.f);

    float maxT[DIMS];

    // Calculate T distances to candidate planes
    for (int i = 0; i < DIMS; i++)
        if (quadrant[i] != Mid && r.direction[i] !=0.f)
            maxT[i] = (candidatePlane[i] - r.origin[i]) * r.inv_direction[i];
        else
            maxT[i] = -1.;


    // Get largest of the maxT's for final choice of intersection
    int whichPlane = 0;
    for (int i = 1; i < DIMS; ++i) {
        if (maxT[whichPlane] < maxT[i]) {
            whichPlane = i;
        }
    }

    if (maxT[whichPlane] < 0.f) return {};

    vector3f coord{};
    for (int i = 0; i < DIMS; i++) {
        if (whichPlane != i) {
            coord[i] = r.origin[i] + maxT[whichPlane] * r.direction[i];
            if (coord[i] < min[i] || coord[i] > max[i])
                return {};

        } else {
            coord[i] = candidatePlane[i];
        }
    }

    return std::make_optional(length(coord - r.origin));
}
