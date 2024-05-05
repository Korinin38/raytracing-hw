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
    intersection.local_coords = {u, v};
    intersection.color = get_color(intersection.local_coords);
    intersection.normal = get_geometric_normal();
    intersection.distance = t;
    if (dot(ray.direction, intersection.normal) > 0) {
        intersection.inside = true;
        intersection.normal = -intersection.normal;
    }

    return intersection;
}

bool Primitive::emissive() const {
    return (!mesh->material.emission.is_zero());
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

    vector3f local_z = ::normal((1 - local_coords.x - local_coords.y) * normal[0] + local_coords.x * normal[1] + local_coords.y * normal[2]);

    if (mesh->material.normal_i == Material::NO_TEXTURE)
        return local_z;

    vector3f tangents[3] = {tangent[0].reduce(),tangent[1].reduce(),tangent[2].reduce()};
    vector3f local_x = ::normal(multiplyVector(mesh->normal_transform, ::normal((1 - local_coords.x - local_coords.y) * tangents[0] + local_coords.x * tangents[1] + local_coords.y * tangents[2])));
    vector3f local_y = cross(local_z, local_x) * tangent[0].w;

    vector2f tc = get_texcoord(local_coords);
    vector3f sample = textures[mesh->material.normal_i].RGBA_sample(tc).reduce();
    vector3f local_normal = sample.add(-0.5f) * 2.f;

    local_normal = ::normal(local_x * local_normal.x + local_y * local_normal.y + local_z * local_normal.z);
    return local_normal;
}

vector2f Primitive::get_texcoord(vector2f local_coords) const {
    return (1 - local_coords.x - local_coords.y) * texcoord[0] + local_coords.x * texcoord[1] + local_coords.y * texcoord[2];
}

vector3f Primitive::get_color(vector2f local_coords) const {
    if (mesh->material.base_color_i == Material::NO_TEXTURE)
        return mesh->material.base_color;
    vector2f tc = get_texcoord(local_coords);

    vector3f color = textures[mesh->material.base_color_i].sRGBA_sample(tc).reduce();
    vector3f &factor = mesh->material.base_color;
    return color * factor;
}

vector3f Primitive::get_emission(vector2f local_coords) const {
    if (mesh->material.emission_i == Material::NO_TEXTURE)
        return mesh->material.emission;
    vector2f tc = get_texcoord(local_coords);

    vector3f emission = textures[mesh->material.emission_i].sRGBA_sample(tc).reduce();
    vector3f &factor = mesh->material.emission;
    return emission * factor;
}

std::tuple<float, float> Primitive::get_metallic_roughness(vector2f local_coords) const {
    if (mesh->material.metallic_roughness_i == Material::NO_TEXTURE)
        return {mesh->material.roughness2, mesh->material.metallic};
    vector2f tc = get_texcoord(local_coords);

    vector4f mr_texture = textures[mesh->material.metallic_roughness_i].RGBA_sample(tc);
    float r2 = mr_texture.y * mr_texture.y;
    float m = mr_texture.z;
    return {r2 * mesh->material.roughness2, m * mesh->material.metallic};
}

bool Primitive::transparent() const {
    return (mesh->material.alpha < 1.f);
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

void set_mesh_texture_primitive_correspondence(std::vector<Mesh> &meshes, std::vector<Primitive> &primitives,
                                               std::vector<Texture> &textures) {
    for (auto &p : primitives) {
        p.mesh = &(meshes[p.mesh_id]);
        p.textures = textures.data();
    }
}
