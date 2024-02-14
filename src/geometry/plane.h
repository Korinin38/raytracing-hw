#pragma once

#include "primitive.h"

class Plane : public Primitive {
public:
    Plane(vector3f n);

    vector3f normal_;

    [[nodiscard]] std::optional<float> intersects(Ray ray) const override;
};

Plane::Plane(vector3f n) : normal_(n) {}

std::optional<float> Plane::intersects(Ray ray) const {
    translateRay(ray);

    float t = - dot(ray.position, normal_) / dot(ray.direction, normal_);

    if (t < 0.f)
        return {};
    return std::make_optional(t);
}
