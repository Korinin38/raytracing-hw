#pragma once

#include <utils/vector.h>
#include <geometry/primitive.h>
#include <algorithm>
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
public:
    std::vector<Node> nodes;
    void buildBVH(std::vector<primitive_sh_ptr> &primitives);
    [[nodiscard]]
    Intersection intersect(const std::vector<primitive_sh_ptr> &primitives, Ray r, size_t node_id = 0) const;
private:

    struct StackBuildNode {
        StackBuildNode(size_t p, size_t f, size_t c) : place(p), first(f), count(c) {}
        size_t place;
        size_t first;
        size_t count;
    };
    size_t buildNode(std::vector<StackBuildNode> &nodes, std::vector<primitive_sh_ptr> &primitives);
};
