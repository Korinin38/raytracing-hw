#pragma once

#include "primitive.h"

class Box : public Primitive {
public:
    Box(vector3f s);
    vector3f size_;
    [[nodiscard]] std::optional<float> intersects(Ray ray) const override;
};

Box::Box(vector3f s) : size_(s) {}

std::optional<float> Box::intersects(Ray ray) const {
    translateRay(ray);

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

    if (t1_max < 0)
        t1_max = t2_min;
    return std::make_optional(t1_max);
}
