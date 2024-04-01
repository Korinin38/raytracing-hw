#pragma once

#include <utils/vector.h>
#include <geometry/primitive.h>
#include <vector>

const size_t invalidKey = std::numeric_limits<size_t>::max();

struct Node {
    AABB aabb;
    size_t left = invalidKey;
    size_t right = invalidKey;
    size_t first_primitive_id = invalidKey;
    size_t primitive_count = 0;
};

class BVH {
    std::vector<Node> nodes;
    size_t buildNode(std::vector<primitive_sh_ptr> &primitives, size_t first, size_t count);
    void buildBVH(std::vector<primitive_sh_ptr> &primitives) {
        buildNode(primitives, 0, primitives.size());
    }
};
