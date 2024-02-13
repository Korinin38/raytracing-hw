#pragma once
#include "primitive.h"

class Box : public Primitive {
public:
    Box(vector3f s);
    vector3f size_;
};

Box::Box(vector3f s) : size_(s) {}
