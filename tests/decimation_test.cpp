// Copyright 2017-2021 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"

#include <cmath>

#include "pmp/algorithms/decimation.h"
#include "pmp/algorithms/features.h"
#include "pmp/algorithms/shapes.h"
#include "pmp/algorithms/subdivision.h"
#include "pmp/algorithms/triangulation.h"
#include "helpers.h"

using namespace pmp;

// plain simplification test
TEST(DecimationTest, simplification)
{
    auto mesh = vertex_onering();
    decimate(mesh, 3);
    EXPECT_EQ(mesh.n_vertices(), size_t(3));
}

// simplify mesh with boundary: boundary has to be preserved
TEST(DecimationTest, simplification_boundary)
{
    auto mesh = plane(10);
    triangulate(mesh);
    ASSERT_TRUE(mesh.is_triangle_mesh());

    decimate(mesh, 20);

    // all boundary vertices have to stay on the boundary of the unit square
    auto vpoint = mesh.vertex_property<Point>("v:point");
    size_t n_boundary = 0;
    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
        {
            ++n_boundary;
            const Point& p = vpoint[v];
            const bool on_boundary = p[0] < 1e-6 || p[0] > 1.0 - 1e-6 ||
                                     p[1] < 1e-6 || p[1] > 1.0 - 1e-6;
            EXPECT_TRUE(on_boundary);
            EXPECT_NEAR(p[2], 0.0, 1e-6);
        }
    }
    EXPECT_GT(n_boundary, size_t(0));
}

// simplify curved mesh with boundary: boundary has to be preserved
TEST(DecimationTest, simplification_boundary_curved)
{
    // paraboloid patch over the unit square, boundary on the square border
    auto mesh = plane(10);
    triangulate(mesh);
    auto vpoint = mesh.vertex_property<Point>("v:point");
    for (auto v : mesh.vertices())
    {
        const Point& p = vpoint[v];
        vpoint[v] = Point(p[0], p[1], p[0] * p[0] + p[1] * p[1]);
    }

    size_t n_boundary_before = 0;
    for (auto v : mesh.vertices())
        if (mesh.is_boundary(v))
            ++n_boundary_before;

    decimate(mesh, 40);
    ASSERT_EQ(mesh.n_vertices(), size_t(40));

    // all boundary vertices have to stay on the boundary of the unit square
    size_t n_boundary = 0;
    size_t n_constrained = 0;
    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
        {
            ++n_boundary;
            const Point& p = vpoint[v];
            const bool on_boundary = p[0] < 1e-6 || p[0] > 1.0 - 1e-6 ||
                                     p[1] < 1e-6 || p[1] > 1.0 - 1e-6;
            EXPECT_TRUE(on_boundary);

            // vertex constrained to a boundary chord does not lie on the
            // paraboloid surface anymore
            if (std::abs(p[2] - (p[0] * p[0] + p[1] * p[1])) > 1e-6)
                ++n_constrained;
        }
    }

    // the boundary is simplified but not removed
    EXPECT_GT(n_boundary, size_t(0));
    EXPECT_LT(n_boundary, n_boundary_before);

    // boundary collapses constrained vertices to the boundary
    EXPECT_GT(n_constrained, size_t(0));
}

// simplify with feature edge preservation enabled
TEST(DecimationTest, simplification_with_features)
{
    auto mesh = hexahedron();
    triangulate(mesh);
    detect_features(mesh, 45);
    loop_subdivision(mesh);
    decimate(mesh, 8, 0, 0, 0, 45);
    EXPECT_EQ(mesh.n_vertices(), size_t(8));
}

// simplify with respect to texture coordinates and seams
TEST(DecimationTest, simplification_texture_mesh)
{
    SurfaceMesh mesh = texture_seams_mesh();

    // if the test mesh does not have texcoords,
    // this test won't work
    ASSERT_TRUE(mesh.has_halfedge_property("h:tex"));

    auto seam_threshold = 1.0e-2;
    decimate(mesh, mesh.n_vertices() - 4,
             10.0,           // aspect ratio
             0.0,            // edge length
             0.0,            // max valence
             135.0,          // normal deviation
             0.0,            // Hausdorff
             seam_threshold, // seam threshold
             1);             // seam angle deviation

    size_t seam_edges = 0;
    auto texcoords = mesh.get_halfedge_property<TexCoord>("h:tex");
    auto seams = mesh.add_edge_property<bool>("e:seam");
    if (texcoords)
    {
        for (auto e : mesh.edges())
        {
            // texcoords are stored in halfedge pointing towards a vertex
            Halfedge h0 = mesh.halfedge(e, 0);
            Halfedge h1 = mesh.halfedge(e, 1);     //opposite halfedge
            Halfedge h0p = mesh.prev_halfedge(h0); // start point edge 0
            Halfedge h1p = mesh.prev_halfedge(h1); // start point edge 1

            // if start or end points differ more than seam_threshold
            // the corresponding edge is a texture seam
            if (norm(texcoords[h1] - texcoords[h0p]) > seam_threshold ||
                norm(texcoords[h0] - texcoords[h1p]) > seam_threshold)
            {
                seam_edges++;
                seams[e] = true;
            }
        }
    }

    // test loop case 2
    auto se = mesh.find_edge(Vertex(4), Vertex(11));
    // test loop case 1
    auto se2 = mesh.find_edge(Vertex(4), Vertex(0));

    EXPECT_EQ(mesh.n_vertices(), size_t(12));
    EXPECT_EQ(seam_edges, size_t(13));
    EXPECT_EQ(seams[se], 1);
    EXPECT_EQ(seams[se2], 1);
}
