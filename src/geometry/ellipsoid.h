#pragma once
#include "primitive.h"

class Ellipsoid : public Primitive {
public:
    Ellipsoid(vector3f r);
    vector3f radius_;
};

Ellipsoid::Ellipsoid(vector3f r) : radius_(r) {}
