// Copyright 2011-2026 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"

#include <pmp/algorithms/segmentation.h>
#include <pmp/algorithms/shapes.h>
#include <pmp/algorithms/utilities.h>
#include <pmp/io/io.h>

#include <map>
#include <set>
#include <random>
#include <vector>

using namespace pmp;

namespace {

// number of distinct segment IDs in the face property
size_t count_segments(const SurfaceMesh& mesh)
{
    auto fseg = mesh.get_face_property<int>("f:segment");
    std::map<int, int> ids;
    for (auto f : mesh.faces())
        ++ids[fseg[f]];
    return ids.size();
}

// fraction of faces whose segment changed after matching segments by
// maximum overlap between the two labelings; faces within two rings of a
// segment boundary in either labeling are excluded (a boundary on a
// discrete mesh can only shift by whole face rings)
double interior_mismatch_fraction(const SurfaceMesh& mesh,
                                  const std::vector<int>& labels_a,
                                  const std::vector<int>& labels_b)
{
    // overlap counts between segments of a and b
    std::map<std::pair<int, int>, size_t> overlap;
    for (auto f : mesh.faces())
        ++overlap[{labels_a[f.idx()], labels_b[f.idx()]}];

    // match each segment of a to the segment of b with maximum overlap
    std::map<int, int> match;
    for (const auto& [key, count] : overlap)
    {
        auto it = match.find(key.first);
        if (it == match.end() || count > overlap.at({key.first, it->second}))
            match[key.first] = key.second;
    }

    // mark faces within two rings of a segment boundary in either labeling
    // (a boundary on a discrete mesh can only shift by whole face rings, so
    // faces in that zone cannot be expected to keep their label)
    std::vector<char> near_boundary(mesh.faces_size(), 0);
    for (const auto* labels : {&labels_a, &labels_b})
    {
        for (auto f : mesh.faces())
        {
            for (auto h : mesh.halfedges(f))
            {
                auto of = mesh.face(mesh.opposite_halfedge(h));
                if (of.is_valid() && (*labels)[of.idx()] != (*labels)[f.idx()])
                {
                    near_boundary[f.idx()] = 1;
                    for (auto hh : mesh.halfedges(of))
                    {
                        auto oof = mesh.face(mesh.opposite_halfedge(hh));
                        if (oof.is_valid())
                            near_boundary[oof.idx()] = 1;
                    }
                }
            }
        }
    }

    size_t mismatches = 0, interior = 0;
    for (auto f : mesh.faces())
    {
        if (near_boundary[f.idx()])
            continue;
        ++interior;
        if (match[labels_a[f.idx()]] != labels_b[f.idx()])
            ++mismatches;
    }
    return interior ? static_cast<double>(mismatches) / interior : 0.0;
}

// edges on a segment boundary for a given labeling
std::set<Edge> boundary_edges(const SurfaceMesh& mesh,
                              const std::vector<int>& labels)
{
    std::set<Edge> edges;
    for (auto e : mesh.edges())
    {
        auto h0 = mesh.halfedge(e, 0);
        auto h1 = mesh.halfedge(e, 1);
        auto f0 = mesh.face(h0);
        auto f1 = mesh.face(h1);
        if (f0.is_valid() && f1.is_valid() &&
            labels[f0.idx()] != labels[f1.idx()])
            edges.insert(e);
    }
    return edges;
}

// fraction of boundary edges in `from` that have no boundary edge in `to`
// incident to one of their vertices (one-ring tolerance)
double displaced_boundary_fraction(const SurfaceMesh& mesh,
                                   const std::set<Edge>& from,
                                   const std::set<Edge>& to)
{
    if (from.empty())
        return 0.0;

    std::vector<char> near_vertex(mesh.vertices_size(), 0);
    for (auto e : to)
    {
        near_vertex[mesh.vertex(e, 0).idx()] = 1;
        near_vertex[mesh.vertex(e, 1).idx()] = 1;
    }

    size_t displaced = 0;
    for (auto e : from)
    {
        if (!near_vertex[mesh.vertex(e, 0).idx()] &&
            !near_vertex[mesh.vertex(e, 1).idx()])
            ++displaced;
    }
    return static_cast<double>(displaced) / from.size();
}

} // namespace

TEST(SegmentationTest, empty_mesh)
{
    SurfaceMesh mesh;
    curvature_segmentation(mesh);
    EXPECT_TRUE(mesh.has_face_property("f:segment"));
}

TEST(SegmentationTest, sphere_single_segment)
{
    // a sphere is uniformly convex: expect a single segment
    SurfaceMesh mesh = icosphere(5);
    curvature_segmentation(mesh);
    auto fseg = mesh.get_face_property<int>("f:segment");
    ASSERT_TRUE(fseg);
    for (auto f : mesh.faces())
        EXPECT_EQ(fseg[f], 0);
    EXPECT_EQ(count_segments(mesh), 1u);
}

TEST(SegmentationTest, fandisk)
{
    SurfaceMesh mesh;
    read(mesh, "data/off/fandisk.off");
    curvature_segmentation(mesh);

    auto fseg = mesh.get_face_property<int>("f:segment");
    ASSERT_TRUE(fseg);

    // all faces assigned, and more than one segment found
    for (auto f : mesh.faces())
        EXPECT_GE(fseg[f], 0);
    EXPECT_GT(count_segments(mesh), 3u);
}

TEST(SegmentationTest, bunny)
{
    SurfaceMesh mesh;
    read(mesh, "data/off/bunny.off");
    curvature_segmentation(mesh);

    auto fseg = mesh.get_face_property<int>("f:segment");
    ASSERT_TRUE(fseg);
    for (auto f : mesh.faces())
        EXPECT_GE(fseg[f], 0);
    EXPECT_GT(count_segments(mesh), 3u);
}

TEST(SegmentationTest, stability)
{
    SurfaceMesh mesh;
    read(mesh, "data/off/fandisk.off");

    curvature_segmentation(mesh);
    auto fseg = mesh.get_face_property<int>("f:segment");
    std::vector<int> labels_a(mesh.faces_size());
    for (auto f : mesh.faces())
        labels_a[f.idx()] = fseg[f];

    // perturb vertices with noise below 0.1% of the bbox diagonal
    const Scalar diag = norm(bounds(mesh).max() - bounds(mesh).min());
    std::mt19937 rng(42);
    std::uniform_real_distribution<Scalar> uniform(-1.0, 1.0);
    for (auto v : mesh.vertices())
    {
        mesh.position(v) +=
            Scalar(0.0005) * diag *
            Point(uniform(rng), uniform(rng), uniform(rng));
    }

    curvature_segmentation(mesh);
    fseg = mesh.get_face_property<int>("f:segment");
    std::vector<int> labels_b(mesh.faces_size());
    for (auto f : mesh.faces())
        labels_b[f.idx()] = fseg[f];

    // segment interiors must be unaffected by the perturbation
    const double interior_mismatch =
        interior_mismatch_fraction(mesh, labels_a, labels_b);
    EXPECT_LT(interior_mismatch, 0.05);

    // segment boundaries must stay in place: on a discrete mesh a boundary
    // can only shift by whole face rings, so each boundary edge must have a
    // counterpart within one vertex ring in the other labeling (thresholds
    // validated to hold with margin for different noise seeds)
    const auto boundary_a = boundary_edges(mesh, labels_a);
    const auto boundary_b = boundary_edges(mesh, labels_b);
    EXPECT_LT(displaced_boundary_fraction(mesh, boundary_a, boundary_b),
              0.10);
    EXPECT_LT(displaced_boundary_fraction(mesh, boundary_b, boundary_a),
              0.10);
}
