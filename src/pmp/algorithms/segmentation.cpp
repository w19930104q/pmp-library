// Copyright 2011-2026 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#include "pmp/algorithms/segmentation.h"

#include "pmp/algorithms/curvature.h"
#include "pmp/algorithms/utilities.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>

namespace pmp {

namespace {

// vertex classification based on curvature signs
enum : int
{
    FLAT = 0,
    CONVEX = 1,
    CONCAVE = 2,
    SADDLE = 3
};

// union-find with path compression
int find_root(std::vector<int>& parent, int i)
{
    int root = i;
    while (parent[root] != root)
        root = parent[root];
    while (parent[i] != root)
    {
        int next = parent[i];
        parent[i] = root;
        i = next;
    }
    return root;
}

} // namespace

void curvature_segmentation(SurfaceMesh& mesh, unsigned int min_segment_size,
                            int smoothing_steps)
{
    auto fseg = mesh.face_property<int>("f:segment", -1);

    if (mesh.n_vertices() == 0)
        return;

    // per-vertex mean and Gaussian curvature, estimated from the curvature
    // tensor over two-ring neighborhoods and smoothed for robustness
    const auto n = mesh.vertices_size();
    std::vector<Scalar> mean_curv(n, 0.0), gauss_curv(n, 0.0);

    curvature(mesh, Curvature::Mean, smoothing_steps, true, true);
    auto vcurv = mesh.get_vertex_property<Scalar>("v:curv");
    for (auto v : mesh.vertices())
        mean_curv[v.idx()] = vcurv[v];

    curvature(mesh, Curvature::Gauss, smoothing_steps, true, true);
    vcurv = mesh.get_vertex_property<Scalar>("v:curv");
    for (auto v : mesh.vertices())
        gauss_curv[v.idx()] = vcurv[v];

    mesh.remove_vertex_property(vcurv);

    // Estimate the noise floor of the curvature values from their
    // high-frequency residual (deviation from the neighborhood mean) using
    // the median absolute deviation, a robust estimator that is insensitive
    // to the small fraction of vertices on actual features.
    auto noise_level = [&](const std::vector<Scalar>& value) {
        std::vector<Scalar> residual;
        residual.reserve(mesh.n_vertices());
        for (auto v : mesh.vertices())
        {
            Scalar sum(0.0);
            unsigned int count = 0;
            for (auto vv : mesh.vertices(v))
            {
                sum += value[vv.idx()];
                ++count;
            }
            const Scalar r =
                count ? value[v.idx()] - sum / count : Scalar(0.0);
            residual.push_back(std::fabs(r));
        }
        auto mid = residual.begin() + residual.size() / 2;
        std::nth_element(residual.begin(), mid, residual.end());
        // median absolute deviation -> standard deviation of white noise
        return Scalar(1.4826) * *mid;
    };

    // Classification thresholds: the maximum of
    //  - a geometric floor corresponding to features on the scale of the
    //    bounding box,
    //  - a discretization floor: curvature responses below an angular scale
    //    of a few degrees per edge are within the discretization noise of
    //    any mesh and must not drive the classification,
    //  - three times the estimated noise level.
    // Adapting the thresholds to the noise floor keeps the classification
    // (and hence the segment boundaries) stable under small perturbations
    // of the input geometry.
    const auto bb = bounds(mesh);
    const Scalar diag = norm(bb.max() - bb.min());
    const Scalar safe_diag = std::max(diag, Scalar(1e-12));
    Scalar mean_edge(0.0);
    for (auto e : mesh.edges())
        mean_edge += distance(mesh.position(mesh.vertex(e, 0)),
                              mesh.position(mesh.vertex(e, 1)));
    mean_edge /= std::max(mesh.n_edges(), size_t(1));
    const Scalar safe_edge = std::max(mean_edge, Scalar(1e-12));
    const Scalar angle_floor(0.05); // ~3 degrees
    const Scalar eps_h =
        std::max({Scalar(1.0) / safe_diag, angle_floor / safe_edge,
                  Scalar(3.0) * noise_level(mean_curv)});
    const Scalar eps_k = std::max(
        {Scalar(1.0) / (safe_diag * safe_diag),
         angle_floor * angle_floor / (safe_edge * safe_edge),
         Scalar(3.0) * noise_level(gauss_curv)});

    // classify vertices by curvature sign
    std::vector<int> vclass(n, FLAT);
    for (auto v : mesh.vertices())
    {
        const Scalar h = mean_curv[v.idx()];
        const Scalar k = gauss_curv[v.idx()];
        if (k < -eps_k)
            vclass[v.idx()] = SADDLE;
        else if (h > eps_h)
            vclass[v.idx()] = CONVEX;
        else if (h < -eps_h)
            vclass[v.idx()] = CONCAVE;
        else
            vclass[v.idx()] = FLAT;
    }

    // Regularize the classification by mode filtering: replace each vertex
    // class by the majority class in its one-ring neighborhood (its own
    // class counted twice, ties keep the current class). This removes
    // isolated classification flips caused by noise and straightens class
    // boundaries, which is essential for stable segment boundaries.
    for (int iter = 0; iter < 2; ++iter)
    {
        std::vector<int> filtered = vclass;
        for (auto v : mesh.vertices())
        {
            int count[4] = {0, 0, 0, 0};
            count[vclass[v.idx()]] = 2;
            for (auto vv : mesh.vertices(v))
                ++count[vclass[vv.idx()]];
            int best = vclass[v.idx()];
            for (int c = 0; c < 4; ++c)
                if (count[c] > count[best])
                    best = c;
            filtered[v.idx()] = best;
        }
        vclass = filtered;
    }

    // region growing: connected components of equally classified vertices
    std::vector<int> vseg(n, -1);
    int n_segments = 0;
    std::vector<Vertex> stack;
    for (auto seed : mesh.vertices())
    {
        if (vseg[seed.idx()] != -1)
            continue;

        const int cls = vclass[seed.idx()];
        vseg[seed.idx()] = n_segments;
        stack.push_back(seed);
        while (!stack.empty())
        {
            const Vertex v = stack.back();
            stack.pop_back();
            for (auto vv : mesh.vertices(v))
            {
                if (vseg[vv.idx()] == -1 && vclass[vv.idx()] == cls)
                {
                    vseg[vv.idx()] = n_segments;
                    stack.push_back(vv);
                }
            }
        }
        ++n_segments;
    }

    // default minimum segment size, proportional to mesh resolution
    if (min_segment_size == 0)
        min_segment_size =
            std::max<unsigned int>(4, mesh.n_vertices() / 200);

    // merge small segments into the neighbor with the most similar mean
    // curvature (ties broken by largest shared border, then lowest ID)
    std::vector<int> parent(n_segments);
    std::iota(parent.begin(), parent.end(), 0);

    while (true)
    {
        // resolve segment IDs to current roots
        for (auto v : mesh.vertices())
            vseg[v.idx()] = find_root(parent, vseg[v.idx()]);

        // per-segment statistics: size, mean curvature sum, shared borders
        std::vector<unsigned int> size(n_segments, 0);
        std::vector<Scalar> sum_h(n_segments, 0.0);
        for (auto v : mesh.vertices())
        {
            ++size[vseg[v.idx()]];
            sum_h[vseg[v.idx()]] += mean_curv[v.idx()];
        }
        std::vector<std::map<int, unsigned int>> border(n_segments);
        for (auto e : mesh.edges())
        {
            const int s0 = vseg[mesh.vertex(e, 0).idx()];
            const int s1 = vseg[mesh.vertex(e, 1).idx()];
            if (s0 != s1)
            {
                ++border[s0][s1];
                ++border[s1][s0];
            }
        }

        // find smallest segment below the size threshold that has a neighbor
        int candidate = -1;
        for (int s = 0; s < n_segments; ++s)
        {
            if (parent[s] != s || size[s] >= min_segment_size ||
                border[s].empty())
                continue;
            if (candidate == -1 || size[s] < size[candidate])
                candidate = s;
        }
        if (candidate == -1)
            break;

        // find most similar neighbor
        const Scalar mean_h = sum_h[candidate] / size[candidate];
        int best = -1;
        Scalar best_diff = 0.0;
        unsigned int best_border = 0;
        for (const auto& [t, count] : border[candidate])
        {
            const Scalar diff =
                std::fabs(mean_h - sum_h[t] / size[t]);
            if (best == -1 || diff < best_diff ||
                (diff == best_diff &&
                 (count > best_border || (count == best_border && t < best))))
            {
                best = t;
                best_diff = diff;
                best_border = count;
            }
        }

        // merge candidate into best
        parent[candidate] = best;
    }

    // relabel segments consecutively, ordered by first vertex occurrence
    std::vector<int> relabel(n_segments, -1);
    int n_final = 0;
    for (auto v : mesh.vertices())
    {
        const int root = find_root(parent, vseg[v.idx()]);
        if (relabel[root] == -1)
            relabel[root] = n_final++;
        vseg[v.idx()] = relabel[root];
    }

    // transfer vertex segments to faces by majority vote
    for (auto f : mesh.faces())
    {
        std::map<int, unsigned int> count;
        for (auto v : mesh.vertices(f))
            ++count[vseg[v.idx()]];

        int best_seg = -1;
        unsigned int best_count = 0;
        for (const auto& [s, c] : count)
        {
            if (c > best_count)
            {
                best_seg = s;
                best_count = c;
            }
        }
        fseg[f] = best_seg;
    }
}

} // namespace pmp
