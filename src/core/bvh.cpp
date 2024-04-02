#include "bvh.h"
#include <algorithm>

size_t BVH::buildNode(std::vector<primitive_sh_ptr> &primitives, size_t first, size_t count) { // NOLINT(*-no-recursion)
    size_t place = nodes.size();
    nodes.emplace_back();
    Node &cur = nodes.back();

    auto begin = primitives.begin() + (ptrdiff_t)first;
    auto end = primitives.begin() + (ptrdiff_t)first + (ptrdiff_t)count;

    for (auto it = begin; it != end; ++it) {
        cur.aabb.grow(it->get()->aabb());
    }

    if (count <= 4) {
        cur.first_primitive_id = first;
        cur.primitive_count = count;
        return place;
    }

    vector3f aabb_size = cur.aabb.size();
    int aabb_widest_axis = 0;
    {
        float max_size = aabb_size.x;
        for (int i = 0; i < 3; ++i) {
            if (max_size < aabb_size[i]) {
                max_size = aabb_size[i];
                aabb_widest_axis = i;
            }
        }
    }

    // partition
    float middle = cur.aabb.min[aabb_widest_axis] + aabb_size[aabb_widest_axis] / 2;
    std::function<bool(primitive_sh_ptr)> pred = [aabb_widest_axis, middle](const primitive_sh_ptr &a) {
        return ((a.get())->aabb().max[aabb_widest_axis] < middle);
    };

    auto res = std::partition(begin, end, pred);
    if (res == begin || res == end) {
        cur.first_primitive_id = first;
        cur.primitive_count = count;
        return place;
    }

    cur.left = buildNode(primitives, begin - primitives.begin(), res - begin);
    cur.right = buildNode(primitives, res - primitives.begin(), primitives.end() - res);

    return place;
}

Intersection BVH::intersect(Ray r, size_t node_id) {
    Node & node = nodes[node_id];

    if (!node.aabb.intersect(r))
        return {};

    Intersection intersection;

    // todo
    if (node.left != invalidKey) {
//        Intersection a = intersect(r, )

    }

    return Intersection();
}
