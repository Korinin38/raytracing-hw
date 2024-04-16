#include "primitive.h"

#include <vector>
#include <string>
#include <cmath>
#include <sstream>

Intersection::operator bool() const {
    return successful;
}

bool Primitive::parse(const std::string& line) {
    std::stringstream ss(line);
    std::string cmd;

    ss >> cmd;

    if (cmd == "BOX") {
		type = Box;
        param_[0] = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "ELLIPSOID") {
		type = Ellipsoid;
        param_[0] = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "PLANE") {
		type = Plane;
        param_[0] = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "TRIANGLE") {
		type = Triangle;
        vector3f p1{}, p2{}, p3{};
        ss >> p1.x >> p1.y >> p1.z;
        ss >> p2.x >> p2.y >> p2.z;
        ss >> p3.x >> p3.y >> p3.z;
        param_[0] = p1;
        param_[1] = p2 - p1;
        param_[2] = p3 - p1;
        return true;
    } else if (cmd == "POSITION") {
		position = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "ROTATION") {
		rotation = normal(vec4f_from_string(line, cmd.length() + 1));
        return true;
    } else if (cmd == "COLOR") {
		material.color = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "EMISSION") {
        material.emission = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "METALLIC") {
		material.type = Material::Type::Metallic;
        return true;
    } else if (cmd == "DIELECTRIC") {
		material.type = Material::Type::Dielectric;
        return true;
    } else if (cmd == "IOR") {
        material.ior = float_from_string(line, cmd.length() + 1);
        return true;
    }
    return false;
}

void Primitive::transformRay(Ray &ray) const {
    ray.origin = ray.origin - position;

    ray.origin = rotate(ray.origin, *rotation);
    ray.direction = rotate(ray.direction, *rotation);

    ray.inv_direction = vector3f{1.f, 1.f, 1.f} / ray.direction;
}

Ray::Ray(vector3f p, vector3f d) : origin(p), direction(d) {
    normalize(direction);
    inv_direction = vector3f{1.f, 1.f, 1.f} / direction;
}

Intersection Primitive::intersect(Ray ray) const {
	transformRay(ray);

    switch (type) {
        case Box: {
            const auto &size_ = param_[0];
            vector3f t1{};
            vector3f t2{};
            for (int i = 0; i < 3; ++i) {
                t1[i] = (-size_[i] - ray.origin[i]) * ray.inv_direction[i];
                t2[i] = (size_[i] - ray.origin[i]) * ray.inv_direction[i];
                if (t1[i] > t2[i])
                    std::swap(t1[i], t2[i]);
            }
            float t1_max = std::max(std::max(t1.x, t1.y), t1.z);
            float t2_min = std::min(std::min(t2.x, t2.y), t2.z);

            if (t1_max > t2_min || t2_min < 0)
                return {};

            Intersection intersection;
            intersection.successful = true;
            intersection.distance = t1_max;
            intersection.color = material.color;
            if (t1_max < 0) {
                intersection.inside = true;
                intersection.distance = t2_min;
            }

            // normal calculation
            {
                vector3f intersection_point = ray.origin + intersection.distance * ray.direction;
                if (cache.box_inv_size.is_zero()) {
                    cache.box_inv_size = (vector3f{1.f, 1.f, 1.f} / size_);
                }
                intersection.normal = intersection_point * cache.box_inv_size;

                float max_dist = 0.f;
                int max_idx = 0;
                for (int i = 0; i < 3; ++i) {
                    if (std::abs(intersection.normal[i]) >= max_dist) {
                        max_dist = std::abs(intersection.normal[i]);
                        max_idx = i;
                    }
                }
                for (int i = 0; i < 3; ++i) {
                    intersection.normal[i] = (i == max_idx) ? (intersection.normal[i] > 0.f ? 1.f : -1.f) : 0.f;
                }
                normalize(intersection.normal);

                if (intersection.inside)
                    intersection.normal = -intersection.normal;

                intersection.normal = rotate(intersection.normal, rotation);
            }
            return intersection;
        }
        case Plane: {
            const auto &normal_ = param_[0];

            float t = -dot(ray.origin, normal_) / dot(ray.direction, normal_);

            if (t < 0.f)
                return {};

            Intersection intersection;

            intersection.successful = true;
            intersection.color = material.color;
            intersection.normal = rotate(normal_, rotation);
            intersection.distance = t;
            if (dot(ray.direction, normal_) > 0) {
                intersection.inside = true;
                intersection.normal = -intersection.normal;
            }

            return intersection;
        }
        case Ellipsoid: {
            const auto &radius_ = param_[0];

            float a, b, c;
            vector3f o_r = ray.origin / radius_;
            vector3f d_r = ray.direction / radius_;
            a = dot(d_r, d_r);
            b = dot(o_r, d_r);
            c = dot(o_r, o_r);

            float h = b * b - a * (c - 1.f);
            if (h < 0.f)
                return {};

            h = std::sqrt(h);
            float t1 = (-b - h) / a;
            float t2 = (-b + h) / a;

            if (t2 < 0)
                return {};

            Intersection intersection;
            intersection.successful = true;
            intersection.color = material.color;
            intersection.distance = t1;

            if (t1 < 0) {
                intersection.inside = true;
                intersection.distance = t2;
            }

            // normal calculation
            {
                vector3f intersection_point = ray.origin + intersection.distance * ray.direction;
                intersection.normal = normal(intersection_point / (radius_ * radius_));
                intersection.normal = rotate(intersection.normal, rotation);
                if (intersection.inside)
                    intersection.normal = -intersection.normal;
            }

            return intersection;
        }
        case Triangle: {
            const vector3f &triangle_origin = param_[0];
            const vector3f &U = param_[1];
            const vector3f &V = param_[2];

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

            if (cache.triangle_normal.is_zero()) {
                cache.triangle_normal = cross(U, V);
                cache.triangle_area = 0.5f * length(cache.triangle_normal);
                normalize(cache.triangle_normal);
            }

            intersection.successful = true;
            intersection.color = material.color;
            intersection.normal = rotate(cache.triangle_normal, rotation);
            intersection.distance = t;
            if (dot(ray.direction, intersection.normal) > 0) {
                intersection.inside = true;
                intersection.normal = -intersection.normal;
            }

            return intersection;

        }
    }
    throw std::runtime_error("Reached unreachable code");
}

bool Primitive::emissive() const {
    return (!material.emission.is_zero());
}

vector3f Primitive::to_global(vector3f local) const {
    return rotate(local, rotation) + position;
}

vector3f Primitive::to_local(vector3f global) const {
    return rotate(global - position, *rotation);
}

AABB Primitive::aabb() const {
    if (!cache.aabb.empty()) return cache.aabb;
    switch (type) {
        case Ellipsoid: {
            // todo: make more tight AABB
            //  for now, use Box case
        }
        case Box: {
            auto size = param_[0];
            std::vector<vector3f> points;
            points.reserve(8);
            for (float sgn1 : {-1.f, 1.f})
            for (float sgn2 : {-1.f, 1.f})
            for (float sgn3 : {-1.f, 1.f}) {
                points.push_back({size.x * sgn1, size.y * sgn2, size.z * sgn3});
            }

            for (auto &p : points) {
                p = rotate(p,rotation);
                cache.aabb.grow(p + position);
            }
            break;
        }
        case Triangle: {
            vector3f points[3];
            points[0] = param_[0];
            points[1] = param_[0] + param_[1];
            points[2] = param_[0] + param_[2];
            for (auto &p : points) {
                cache.aabb.grow(rotate(p, rotation) + position);
            }
            break;
        }
        case Plane:
        default:
            throw std::runtime_error("Incorrect type to calculate AABB");
    }
    return cache.aabb;
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
