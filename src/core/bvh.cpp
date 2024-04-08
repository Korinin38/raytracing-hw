#include "bvh.h"
#include <algorithm>

namespace {
    size_t buildHelperNaive(std::vector<primitive_sh_ptr> &primitives, AABB node_aabb, size_t first, size_t count, size_t &dim) {
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

        dim = aabb_widest_axis;
        // partition
        float middle = node_aabb.min[aabb_widest_axis] + aabb_size[aabb_widest_axis] / 2;
        std::function<bool(primitive_sh_ptr)> pred = [aabb_widest_axis, middle](const primitive_sh_ptr &a) {
            return ((a.get())->aabb().max[aabb_widest_axis] < middle);
        };

        auto part = std::partition(begin, end, pred);
        return (part - begin);
    }

    size_t buildHelperSAH(std::vector<primitive_sh_ptr> &primitives, AABB node_aabb, size_t first, size_t count, size_t &best_dim) {
        const int BIN_SIZE = 1;
        auto begin = primitives.begin() + (ptrdiff_t)first;
        auto end = primitives.begin() + (ptrdiff_t)first + (ptrdiff_t)count;


//        const int bins_count = ((count - 1) / BIN_SIZE) + 1;
        const int bins_count = count;
        std::vector<AABB> sah_bins(bins_count);
        std::vector<float> sah_left_scan(bins_count);
        std::vector<float> sah_right_scan(bins_count);

        float best = node_aabb.surface_area() * (float)count;
        auto result = begin;
        best_dim = 2;
        for (int dim = 0; dim < 3; ++dim) {
            auto predicate = [dim](primitive_sh_ptr &a, primitive_sh_ptr &b) {
                return ((a->aabb().max[dim] + a->aabb().min[dim]) < (b->aabb().max[dim] + b->aabb().min[dim]));
            };

            std::sort(begin, end, predicate);

//            #pragma omp parallel for shared(primitives)
            for (int bin = 0; bin < bins_count; ++bin) {
                for (auto obj = begin + bin * BIN_SIZE; obj != end; obj++) {
                    sah_bins[bin].grow(obj->get()->aabb());
                }
            }

            AABB tmp = sah_bins[0];
            sah_left_scan[0] = tmp.surface_area();
            for (int i = 1; i < bins_count; ++i) {
                tmp.grow(sah_bins[i]);
                sah_left_scan[i] = tmp.surface_area();
            }

            tmp = sah_bins[bins_count - 1];
            sah_right_scan[bins_count - 1] = tmp.surface_area();
            for (int i = bins_count - 2; i >= 0; --i) {
                tmp.grow(sah_bins[i]);
                sah_right_scan[i] = tmp.surface_area();
            }

            for (int i = 0; i < bins_count - 1; ++i) {
                float ls = sah_left_scan[i];
                float rs = sah_right_scan[i + 1];

                int left_obj_count = (i + 1) * BIN_SIZE;

                float metric = ls * (float)(left_obj_count) + rs * (float)(count - left_obj_count);
                if (metric < best) {
                    best = metric;
                    result = begin + left_obj_count;
                    best_dim = dim;
                }
            }
        }

        if (best_dim != 2) {
            auto predicate = [best_dim](primitive_sh_ptr &a, primitive_sh_ptr &b) {
                return ((a->aabb().max[best_dim] + a->aabb().min[best_dim]) < (b->aabb().max[best_dim] + b->aabb().min[best_dim]));
            };
            std::sort(begin, end, predicate);
        }
        return (result - begin);
    }
}

size_t BVH::buildNode(std::vector<StackBuildNode> &nodes_q, std::vector<primitive_sh_ptr> &primitives) { // NOLINT(*-no-recursion)
    StackBuildNode stack = nodes_q.back();
    nodes_q.pop_back();

    size_t place = stack.place;

    size_t first = stack.first;
    size_t count = stack.count;

    if (count == 0)
        throw std::runtime_error("What?");

    auto begin = primitives.begin() + (ptrdiff_t)first;
    auto end = primitives.begin() + (ptrdiff_t)first + (ptrdiff_t)count;

    if ((ptrdiff_t)first + (ptrdiff_t)count > primitives.size())
        throw std::runtime_error("What?");
    for (auto it = begin; it != end; ++it) {
        nodes[place].aabb.grow(it->get()->aabb());
    }

    if (count <= 4) {
//    if (count <= 1) {
        nodes[place].first_primitive_id = first;
        nodes[place].primitive_count = count;
        return place;
    }

//    auto partition = buildHelperNaive(primitives, cur.aabb, first, count);
//    auto partition = buildHelperSAH(primitives, cur.aabb, first, count);
    size_t res = buildHelperNaive(primitives, nodes[place].aabb, first, count, nodes[place].split_dim);
    auto partition = begin + res;
//    size_t res = buildHelperSAH(primitives, nodes[place].aabb, first, count, nodes[place].split_dim);
//    auto partition = begin + res;

    if (partition == begin || partition == end) {
        nodes[place].first_primitive_id = first;
        nodes[place].primitive_count = count;
        return place;
    }

    size_t l_fc = begin - primitives.begin();
    size_t l_cnt = partition - begin;
    size_t r_fc = partition - primitives.begin();
    size_t r_cnt = end - partition;

    nodes[place].left = nodes.size();
    nodes.emplace_back();
    nodes_q.emplace_back(nodes[place].left, l_fc, l_cnt);
    nodes[place].right = nodes.size();
    nodes.emplace_back();
    nodes_q.emplace_back(nodes[place].right, r_fc, r_cnt);

    return place;
}

Intersection BVH::intersect(const std::vector<primitive_sh_ptr> &primitives, Ray r, size_t node_id, bool early_out) const {
    const Node & node = nodes[node_id];

    if (!node.aabb.intersect(r))
        return {};

//    // todo: remove?
//    early_out = false;

    Intersection intersection;
    intersection.distance = 1e9;

    if (r.direction[node.split_dim] > 0) {
        if (node.left != invalidKey) {
            Intersection a = intersect(primitives, r, node.left);
            if (a && a.distance < intersection.distance)
                intersection = a;

        }
//        if (intersection && early_out)
//            return intersection;
        if (node.right != invalidKey) {
            Intersection a = intersect(primitives, r, node.right);
            if (a && a.distance < intersection.distance)
                intersection = a;
        }
    }
    else {
        if (node.right != invalidKey) {
            Intersection a = intersect(primitives, r, node.right);
            if (a && a.distance < intersection.distance)
                intersection = a;
        }
//        if (intersection && early_out)
//            return intersection;
        if (node.left != invalidKey) {
            Intersection a = intersect(primitives, r, node.left);
            if (a && a.distance < intersection.distance)
                intersection = a;
        }
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
