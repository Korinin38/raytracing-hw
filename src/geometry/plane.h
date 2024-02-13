#pragma once
#include "primitive.h"

class Plane : public Primitive {
public:
    Plane(vector3f n);
    vector3f normal_;
};

Plane::Plane(vector3f n) : normal_(n) {}
