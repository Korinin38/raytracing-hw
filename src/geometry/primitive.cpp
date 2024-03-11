#include "primitive.h"

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
        type_ = Box;
        param_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "ELLIPSOID") {
        type_ = Ellipsoid;
        param_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "PLANE") {
        type_ = Plane;
        param_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "POSITION") {
        position_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "ROTATION") {
        rotation_ = normal(vec4f_from_string(line, cmd.length() + 1));
        return true;
    } else if (cmd == "COLOR") {
        color_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "EMISSION") {
        emission_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "METALLIC") {
        material_ = Metallic;
        return true;
    } else if (cmd == "DIELECTRIC") {
        material_ = Dielectric;
        return true;
    } else if (cmd == "IOR") {
        ior_ = float_from_string(line, cmd.length() + 1);
        return true;
    }
    return false;
}

void Primitive::translateRay(Ray &ray) const {
    for (int i = 0; i < 3; ++i)
        ray.position[i] -= position_[i];

    ray.position = rotate(ray.position, *rotation_);
    ray.direction = rotate(ray.direction, *rotation_);
}

Ray::Ray(vector3f p, vector3f d) : position(p), direction(d) {
    normalize(direction);
}

Intersection Primitive::intersect(Ray ray) const {
    translateRay(ray);

    switch (type_) {
        case Box: {
            const vector3f &size_ = param_;
            vector3f t1{};
            vector3f t2{};
            for (int i = 0; i < 3; ++i) {
                t1[i] = (-size_[i] - ray.position[i]) / ray.direction[i];
                t2[i] = (size_[i] - ray.position[i]) / ray.direction[i];
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
            intersection.color = color_;
            if (t1_max < 0) {
                intersection.inside = true;
                intersection.distance = t2_min;
            }

            // normal calculation
            {
                vector3f intersection_point = ray.position + intersection.distance * ray.direction;
                intersection.normal = intersection_point / size_;

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

                intersection.normal = rotate(intersection.normal, rotation_);
            }
            return intersection;
        }
        case Plane: {
            const vector3f &normal_ = param_;

            float t = -dot(ray.position, normal_) / dot(ray.direction, normal_);

            if (t < 0.f)
                return {};

            Intersection intersection;

            intersection.successful = true;
            intersection.color = color_;
            intersection.normal = rotate(normal_, rotation_);
            intersection.distance = t;
            if (dot(ray.direction, normal_) > 0) {
                intersection.inside = true;
                intersection.normal = -intersection.normal;
            }

            return intersection;
        }
        case Ellipsoid: {
            const vector3f &radius_ = param_;

            float a, b, c;
            vector3f o_r = ray.position / radius_;
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
            intersection.color = color_;
            intersection.distance = t1;

            if (t1 < 0) {
                intersection.inside = true;
                intersection.distance = t2;
            }

            // normal calculation
            {
                vector3f intersection_point = ray.position + intersection.distance * ray.direction;
                intersection.normal = normal(intersection_point / (radius_ * radius_));
                intersection.normal = rotate(intersection.normal, rotation_);
                if (intersection.inside)
                    intersection.normal = -intersection.normal;
            }

            return intersection;
        }
    }
}

bool Primitive::emissive() {
    vector3f zero{};
    return (emission_ != zero);
}

vector3f Primitive::to_global(vector3f local) const {
    return rotate(local, rotation_) + position_;
}

vector3f Primitive::to_local(vector3f global) const {
    return rotate(global - position_, *rotation_);
}
