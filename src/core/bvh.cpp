#include "bvh.h"
#include <algorithm>

namespace {
    std::vector<primitive_sh_ptr>::iterator buildHelperNaive(std::vector<primitive_sh_ptr> &primitives, AABB node_aabb, size_t place, size_t first, size_t count) {
        auto begin = primitives.begin() + (ptrdiff_t)first;
        auto end = primitives.begin() + (ptrdiff_t)first + (ptrdiff_t)count;

        vector3f aabb_size = node_aabb.size();
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
        float middle = node_aabb.min[aabb_widest_axis] + aabb_size[aabb_widest_axis] / 2;
        std::function<bool(primitive_sh_ptr)> pred = [aabb_widest_axis, middle](const primitive_sh_ptr &a) {
            return ((a.get())->aabb().max[aabb_widest_axis] < middle);
        };

        return std::partition(begin, end, pred);
    }

    std::vector<primitive_sh_ptr>::iterator buildHelperSAH(std::vector<primitive_sh_ptr> &primitives, AABB node_aabb, size_t place, size_t first, size_t count) {
        throw std::runtime_error("Not implemented");
    }
}

size_t BVH::buildNode(std::vector<StackBuildNode> &nodes_q, std::vector<primitive_sh_ptr> &primitives) { // NOLINT(*-no-recursion)
    StackBuildNode stack = nodes_q.back();
    nodes_q.pop_back();

    size_t place = stack.place;
    Node &cur = nodes[place];

    size_t first = stack.first;
    size_t count = stack.count;

    auto begin = primitives.begin() + (ptrdiff_t)first;
    auto end = primitives.begin() + (ptrdiff_t)first + (ptrdiff_t)count;

    if ((ptrdiff_t)first + (ptrdiff_t)count > primitives.size())
        throw std::runtime_error("What?");
    for (auto it = begin; it != end; ++it) {
        cur.aabb.grow(it->get()->aabb());
    }

    if (count <= 4) {
        cur.first_primitive_id = first;
        cur.primitive_count = count;
        return place;
    }

    auto partition = buildHelperNaive(primitives, cur.aabb, place, first, count);
//    auto partition = buildHelperSAH(primitives, cur.aabb, place, first, count);

    if (partition == begin || partition == end) {
        cur.first_primitive_id = first;
        cur.primitive_count = count;
        return place;
    }

    size_t l_fc = begin - primitives.begin();
    size_t l_cnt = partition - begin;
    size_t r_fc = partition - primitives.begin();
    size_t r_cnt = end - partition;

    nodes[place].left = nodes.size();
    nodes_q.emplace_back(nodes.size(), l_fc, l_cnt);
    nodes.emplace_back();
    nodes[place].right = nodes.size();
    nodes_q.emplace_back(nodes.size(), r_fc, r_cnt);
    nodes.emplace_back();

    return place;
}

Intersection BVH::intersect(const std::vector<primitive_sh_ptr> &primitives, Ray r, size_t node_id) const {
    const Node & node = nodes[node_id];

    if (!node.aabb.intersect(r))
        return {};

    Intersection intersection;
    intersection.distance = 1e9;

    // todo
    if (node.left != invalidKey) {
        Intersection a = intersect(primitives, r, node.left);
        if (a && a.distance < intersection.distance)
            intersection = a;
    }
    if (node.right != invalidKey) {
        Intersection a = intersect(primitives, r, node.right);
        if (a && a.distance < intersection.distance)
            intersection = a;
    }
    if (!node.primitive_count)
        return intersection;

    for (size_t i = node.first_primitive_id; i < node.first_primitive_id + node.primitive_count; ++i) {
        Intersection a = primitives[i]->intersect(r);
        if (a && a.distance < intersection.distance) {
            intersection = a;
            intersection.object_id = i;
        }
    }

    return intersection;
}
