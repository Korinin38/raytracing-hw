#include "primitive.h"

#include <string>
#include <cmath>
#include <sstream>

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

std::optional<intersection> Primitive::intersect(Ray ray) const {
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

            intersection i;

            if (t1_max < 0) {
                i.inside = true;
                i.distance = t2_min;
                // todo: normal
            }

            return std::make_optional(i);
        }
        case Plane: {
            const vector3f &normal_ = param_;

            float t = -dot(ray.position, normal_) / dot(ray.direction, normal_);
            intersection i;

            if (t < 0.f)
                return {};

            i.normal = normal_;
            i.distance = std::abs(t);

            return std::make_optional(i);
        }
        case Ellipsoid: {
            const vector3f &radius_ = param_;

            float a, b, c;
            vector3f o_div_r = ray.position / radius_;
            vector3f d_div_r = ray.direction / radius_;
            a = dot(d_div_r, d_div_r);
            b = dot(o_div_r, d_div_r);
            c = dot(o_div_r, o_div_r);

            float h = b * b - a * (c - 1.f);
            if (h < 0.f)
                return {};

            h = std::sqrt(h);
            float t1 = (-b - h) / a;
            float t2 = (-b + h) / a;

            if (t2 < 0)
                return {};
            intersection i;

            if (t1 < 0) {
                i.inside = true;
                i.distance = t2;
                // todo: normal
            }
            return std::make_optional(i);
        }
    }
}
