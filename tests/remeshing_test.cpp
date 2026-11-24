// Copyright 2017-2021 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"

#include "pmp/algorithms/remeshing.h"
#include "pmp/algorithms/features.h"
#include "pmp/algorithms/shapes.h"
#include "pmp/algorithms/triangulation.h"
#include "pmp/algorithms/utilities.h"

#include "helpers.h"

using namespace pmp;

// adaptive remeshing
TEST(RemeshingTest, adaptive_remeshing_with_features)
{
    auto mesh = cylinder();
    triangulate(mesh);
    detect_features(mesh, 25);
    auto bb = bounds(mesh).size();
    adaptive_remeshing(mesh, 0.001 * bb, // min length
                       1.0 * bb,         // max length
                       0.001 * bb);      // approx. error
    EXPECT_EQ(mesh.n_vertices(), 6u);
}

TEST(RemeshingTest, adaptive_remeshing_with_boundary)
{
    // mesh with boundary
    auto mesh = open_cone();
    auto bb = bounds(mesh).size();
    adaptive_remeshing(mesh, 0.01 * bb, // min length
                       1.0 * bb,        // max length
                       0.01 * bb);      // approx. error
    EXPECT_EQ(mesh.n_vertices(), size_t(107));
}

TEST(RemeshingTest, adaptive_remeshing_with_selection)
{
    auto mesh = icosphere(1);

    // select half the vertices
    auto selected = mesh.add_vertex_property<bool>("v:selected");
    for (auto v : mesh.vertices())
        if (mesh.position(v)[1] > 0)
            selected[v] = true;

    auto bb = bounds(mesh).size();
    adaptive_remeshing(mesh, 0.01 * bb, // min length
                       1.0 * bb,        // max length
                       0.01 * bb);      // approx. error
    EXPECT_EQ(mesh.n_vertices(), size_t(58));
}

TEST(RemeshingTest, uniform_remeshing)
{
    auto mesh = open_cone();
    uniform_remeshing(mesh, 0.5);
    EXPECT_EQ(mesh.n_vertices(), size_t(40));
}

TEST(RemeshingTest, uniform_remeshing_skinny_triangles)
{
    // anisotropically sampled heightfield: triangle aspect ratio > 10
    SurfaceMesh mesh;
    const int nx = 40, ny = 4;
    const double lx = 4.0, ly = 4.0;
    std::vector<Vertex> verts;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
        {
            const double x = lx * i / nx;
            const double y = ly * j / ny;
            verts.push_back(mesh.add_vertex(
                Point(x, y, 0.3 * std::sin(2.0 * x) * std::cos(2.0 * y))));
        }
    auto idx = [&](int i, int j) { return verts[j * (nx + 1) + i]; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            mesh.add_triangle(idx(i, j), idx(i + 1, j), idx(i, j + 1));
            mesh.add_triangle(idx(i + 1, j), idx(i + 1, j + 1),
                              idx(i, j + 1));
        }

    uniform_remeshing(mesh, 0.6, 10, /*use_projection=*/true);

    // all triangles should have reasonable angles after remeshing
    auto points = mesh.get_vertex_property<Point>("v:point");
    for (auto f : mesh.faces())
    {
        Point p[3];
        int k = 0;
        for (auto v : mesh.vertices(f))
            p[k++] = points[v];
        for (int i = 0; i < 3; ++i)
        {
            const auto d0 = normalize(p[(i + 1) % 3] - p[i]);
            const auto d1 = normalize(p[(i + 2) % 3] - p[i]);
            // min. angle of 20 degrees
            EXPECT_LT(dot(d0, d1), 0.94);
        }
    }
}
