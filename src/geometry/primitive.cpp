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
		type = Box;
        param_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "ELLIPSOID") {
		type = Ellipsoid;
        param_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "PLANE") {
		type = Plane;
        param_ = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "TRIANGLE") {
		type = Triangle;
        vector3f p1, p2, p3;
        ss >> p1.x >> p1.y >> p1.z;
        ss >> p2.x >> p2.y >> p2.z;
        ss >> p3.x >> p3.y >> p3.z;
        param_ = std::tuple(p1, p2, p3);
        return true;
    } else if (cmd == "POSITION") {
		position = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "ROTATION") {
		rotation = normal(vec4f_from_string(line, cmd.length() + 1));
        return true;
    } else if (cmd == "COLOR") {
		color = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "EMISSION") {
		emission = vec3f_from_string(line, cmd.length() + 1);
        return true;
    } else if (cmd == "METALLIC") {
		material = Metallic;
        return true;
    } else if (cmd == "DIELECTRIC") {
		material = Dielectric;
        return true;
    } else if (cmd == "IOR") {
		ior = float_from_string(line, cmd.length() + 1);
        return true;
    }
    return false;
}

void Primitive::transformRay(Ray &ray) const {
	ray.position = ray.position - position;

    ray.position = rotate(ray.position, *rotation);
    ray.direction = rotate(ray.direction, *rotation);
}

Ray::Ray(vector3f p, vector3f d) : position(p), direction(d) {
    normalize(direction);
}

Intersection Primitive::intersect(Ray ray) const {
	transformRay(ray);

    switch (type) {
        case Box: {
            const auto &size_ = std::get<vector3f>(param_);
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
            intersection.color = color;
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

                intersection.normal = rotate(intersection.normal, rotation);
            }
            return intersection;
        }
        case Plane: {
            const auto &normal_ = std::get<vector3f>(param_);

            float t = -dot(ray.position, normal_) / dot(ray.direction, normal_);

            if (t < 0.f)
                return {};

            Intersection intersection;

            intersection.successful = true;
            intersection.color = color;
            intersection.normal = rotate(normal_, rotation);
            intersection.distance = t;
            if (dot(ray.direction, normal_) > 0) {
                intersection.inside = true;
                intersection.normal = -intersection.normal;
            }

            return intersection;
        }
        case Ellipsoid: {
            const auto &radius_ = std::get<vector3f>(param_);

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
            intersection.color = color;
            intersection.distance = t1;

            if (t1 < 0) {
                intersection.inside = true;
                intersection.distance = t2;
            }

            // normal calculation
            {
                vector3f intersection_point = ray.position + intersection.distance * ray.direction;
                intersection.normal = normal(intersection_point / (radius_ * radius_));
                intersection.normal = rotate(intersection.normal, rotation);
                if (intersection.inside)
                    intersection.normal = -intersection.normal;
            }

            return intersection;
        }
        case Triangle: {
            //todo
        }
    }
    throw std::runtime_error("Reached unreachable code");
}

bool Primitive::emissive() const {
    vector3f zero{};
    return (emission != zero);
}

vector3f Primitive::to_global(vector3f local) const {
    return rotate(local, rotation) + position;
}

vector3f Primitive::to_local(vector3f global) const {
    return rotate(global - position, *rotation);
}
