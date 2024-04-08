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
        const int BIN_SIZE = 8;
        auto begin = primitives.begin() + (ptrdiff_t)first;
        auto end = primitives.begin() + (ptrdiff_t)first + (ptrdiff_t)count;

        int dim = 0;
        auto predicate = [dim](primitive_sh_ptr &a, primitive_sh_ptr &b) {
            return ((a->aabb().max[dim] + a->aabb().min[dim]) / 2 < (b->aabb().max[dim] + b->aabb().min[dim]) / 2);
        };

        const int bins_count = ((end - begin - 1) / BIN_SIZE) + 1;
        std::vector<AABB> sah_bins(bins_count);
        std::vector<AABB> sah_left_scan(bins_count);
        std::vector<AABB> sah_right_scan(bins_count);
        int last_bin_count = 0;

        float best = node_aabb.surface_area() * (float)count;
        auto result = begin;
        int best_dim = 2;

        for (dim = 0; dim < 3; ++dim) {
            std::sort(begin, end, predicate);

            #pragma omp parallel for shared(primitives)
            for (int bin = 0; bin < bins_count; ++bin) {
                if (bin == sah_bins.size() - 1)
                    ++last_bin_count;
                for (auto obj = begin + bin * BIN_SIZE; obj != end; obj++) {
                    sah_bins[bin].grow(obj->get()->aabb());
                }
            }

            sah_left_scan[0].grow(sah_bins[0]);
            for (int i = 1; i < bins_count; ++i) {
                sah_left_scan[i].grow(sah_left_scan[i - 1]);
                sah_left_scan[i].grow(sah_bins[i]);
            }

            sah_right_scan[bins_count - 1].grow(sah_bins[bins_count - 1]);
            for (int i = bins_count - 2; i >= 0; --i) {
                sah_right_scan[i].grow(sah_right_scan[i + 1]);
                sah_right_scan[i].grow(sah_bins[i]);
            }

            for (int i = 1; i < bins_count; ++i) {
                float ls = sah_left_scan[i].surface_area();
                float rs = sah_right_scan[i + 1].surface_area();

                float metric = ls * (float)(i * BIN_SIZE) + rs * (float)(count - i * BIN_SIZE);
                if (metric < best) {
                    best = metric;
                    result = begin + (i * BIN_SIZE);
                    best_dim = dim;
                }
            }
        }

        if (best_dim != 2) {
            dim = best_dim;
            std::sort(begin, end, predicate);
        }
        return result;

//        throw std::runtime_error("Not implemented");
//        return result;
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

//    if (count <= 1) {
    if (count <= 4) {
        cur.first_primitive_id = first;
        cur.primitive_count = count;
        return place;
    }

//    auto partition = buildHelperNaive(primitives, cur.aabb, place, first, count);
    auto partition = buildHelperSAH(primitives, cur.aabb, place, first, count);

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

void BVH::buildBVH(std::vector<primitive_sh_ptr> &primitives) {
    size_t npo = 0;
    for (auto &p : primitives) {
        if (p->type != Primitive::Plane)
            ++npo;
    }
    auto plane_it = std::partition(primitives.begin(), primitives.end(), [](primitive_sh_ptr &p) {return (p->type != Primitive::Plane);} );
    size_t non_plane_objects = plane_it - primitives.begin();

    std::vector<StackBuildNode> nodes_q;
    nodes_q.emplace_back(nodes.size(), 0, non_plane_objects);
    nodes.emplace_back();
    int a = 0;
    while (!nodes_q.empty()) {
        buildNode(nodes_q, primitives);
        a++;
    }
}
