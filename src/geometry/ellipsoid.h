#pragma once

#include "primitive.h"

#include <cmath>

class Ellipsoid : public Primitive {
public:
    Ellipsoid(vector3f r);
    vector3f radius_;

    [[nodiscard]] std::optional<float> intersects(Ray ray) const override;
};

Ellipsoid::Ellipsoid(vector3f r) : radius_(r) {}

std::optional<float> Ellipsoid::intersects(Ray ray) const {
    translateRay(ray);

    float a, b, c;
    vector3f o_div_r = ray.position / radius_;
    vector3f d_div_r = ray.direction / radius_;
    a = dot(d_div_r, d_div_r);
    b = dot(o_div_r, d_div_r);
    c = dot(o_div_r, o_div_r);

    float h = b * b - a * (c - 1.f);
    if( h < 0.f )
        return {};

    h = std::sqrt(h);
    float t1 = (-b - h) / a;
    float t2 = (-b + h) / a;

    if (t2 < 0)
        return {};

    if (t1 < 0)
        t1 = t2;
    return std::make_optional(t1);
}
