// Copyright 2024 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"

#include "pmp/algorithms/hole_filling.h"
#include "pmp/algorithms/hole_filling_intersection.h"
#include "pmp/surface_mesh.h"

#include "helpers.h"

#include <algorithm>
#include <array>
#include <numbers>
#include <vector>

using namespace pmp;

namespace {

// Collect every triangle of the mesh as (points, vertices).
struct TriangleSoup
{
    std::vector<std::array<Point, 3>> points;
    std::vector<std::array<Vertex, 3>> vertices;
};

TriangleSoup triangle_soup(const SurfaceMesh& mesh)
{
    TriangleSoup soup;
    for (auto f : mesh.faces())
    {
        std::vector<Vertex> fv;
        for (auto v : mesh.vertices(f))
            fv.push_back(v);
        for (size_t k = 1; k + 1 < fv.size(); ++k)
        {
            soup.points.push_back({mesh.position(fv[0]), mesh.position(fv[k]),
                                   mesh.position(fv[k + 1])});
            soup.vertices.push_back({fv[0], fv[k], fv[k + 1]});
        }
    }
    return soup;
}

// Count intersections between triangles that share no vertex. Connectivity
// adjacent triangles always meet along an edge/vertex and are not counted.
int count_self_intersections(const SurfaceMesh& mesh)
{
    const auto soup = triangle_soup(mesh);
    int hits = 0;
    for (size_t i = 0; i < soup.points.size(); ++i)
        for (size_t j = i + 1; j < soup.points.size(); ++j)
        {
            bool shared = false;
            for (auto a : soup.vertices[i])
                for (auto b : soup.vertices[j])
                    if (a == b)
                        shared = true;
            if (!shared &&
                interdetail::triangles_intersect(soup.points[i],
                                                 soup.points[j]))
                ++hits;
        }
    return hits;
}

bool is_manifold_mesh(const SurfaceMesh& mesh)
{
    for (auto e : mesh.edges())
        if (!mesh.face(mesh.halfedge(e, 0)).is_valid() ||
            !mesh.face(mesh.halfedge(e, 1)).is_valid())
            return false;
    for (auto v : mesh.vertices())
        if (!mesh.is_manifold(v))
            return false;
    return true;
}

int count_boundary_edges(const SurfaceMesh& mesh)
{
    int n = 0;
    for (auto e : mesh.edges())
        if (mesh.is_boundary(e))
            ++n;
    return n;
}

Halfedge find_boundary(const SurfaceMesh& mesh)
{
    for (auto h : mesh.halfedges())
        if (mesh.is_boundary(h))
            return h;
    return {};
}

std::vector<Halfedge> boundary_loops(SurfaceMesh& mesh)
{
    auto visited =
        mesh.add_halfedge_property<bool>("HoleFillingTest:visited", false);
    std::vector<Halfedge> loops;
    for (auto h : mesh.halfedges())
    {
        if (!mesh.is_boundary(h) || visited[h])
            continue;
        loops.push_back(h);
        Halfedge x = h;
        do
        {
            visited[x] = true;
            x = mesh.next_halfedge(x);
        } while (x != h);
    }
    mesh.remove_halfedge_property(visited);
    return loops;
}

// Disk whose inner rim is a folded 3D loop. The planar projection of the rim
// self-intersects (angle phi=t+a sin 2t with a>0.5) while the loop stays
// embedded in 3D thanks to the antisymmetric z=b sin 2t lift. A transition
// band relaxes the fold before joining perfect concentric circles, and a
// planar cap closes the disk. The whole input is embedded (verified).
SurfaceMesh folded_hole_disk(int n = 96, double a = 1.0, double b = 1.2)
{
    SurfaceMesh m;
    const int L = 8, outer = 6, nplanar = 4;
    const double dr = 0.12, band_drop = 0.8, rout = 3.0, base_z = -2.0;
    const double r_band = 1.0 + L * dr;
    const int nrings = (L + 1) + outer + nplanar;
    std::vector<std::vector<Vertex>> v(nrings,
                                       std::vector<Vertex>(n));

    for (int i = 0; i < n; ++i)
    {
        const double t = 2.0 * std::numbers::pi * i / n;
        const double phi = t + a * std::sin(2.0 * t);
        v[0][i] = m.add_vertex(
            Point(std::cos(phi), std::sin(phi), b * std::sin(2.0 * t)));
    }
    for (int k = 1; k <= L; ++k)
    {
        const double s = static_cast<double>(k) / L;
        const double amp = 1.0 - s;
        const double r = 1.0 + k * dr;
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            const double phi = t + amp * a * std::sin(2.0 * t);
            const double z = amp * b * std::sin(2.0 * t) - s * band_drop;
            v[k][i] = m.add_vertex(
                Point(r * std::cos(phi), r * std::sin(phi), z));
        }
    }
    const int base = L + 1;
    for (int q = 0; q < outer; ++q)
    {
        const double s = static_cast<double>(q + 1) / outer;
        const double r = r_band + (rout - r_band) * s;
        const double z = -band_drop + s * (base_z + band_drop);
        const int k = base + q;
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            v[k][i] =
                m.add_vertex(Point(r * std::cos(t), r * std::sin(t), z));
        }
    }
    for (int q = 0; q < nplanar; ++q)
    {
        const double rr =
            rout * (1.0 - static_cast<double>(q + 1) / (nplanar + 1));
        const int k = base + outer + q;
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            v[k][i] = m.add_vertex(
                Point(rr * std::cos(t), rr * std::sin(t), base_z));
        }
    }
    for (int k = 0; k + 1 < nrings; ++k)
        for (int i = 0; i < n; ++i)
        {
            const int j = (i + 1) % n;
            m.add_triangle(v[k][i], v[k][j], v[k + 1][j]);
            m.add_triangle(v[k][i], v[k + 1][j], v[k + 1][i]);
        }
    const Vertex c = m.add_vertex(Point(0, 0, base_z));
    for (int i = 0; i < n; ++i)
    {
        const int j = (i + 1) % n;
        m.add_triangle(v[nrings - 1][i], v[nrings - 1][j], c);
    }
    return m;
}

// Perforated disk in the yz-plane at x=px, closed on its outward side by a
// hemispherical cap. Its only boundary is the facing circular hole.
SurfaceMesh facing_perforated_disk(int n, double px, int sign, double r = 1.0,
                                   double R = 1.7)
{
    SurfaceMesh m;
    const int nann = 5, ncap = 5;
    const int nrings = nann + ncap;
    std::vector<std::vector<Vertex>> v(nrings,
                                       std::vector<Vertex>(n));
    for (int k = 0; k < nann; ++k)
    {
        const double s = static_cast<double>(k) / (nann - 1);
        const double rr = r + (R - r) * s;
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            v[k][i] = m.add_vertex(
                Point(px, rr * std::cos(t), rr * std::sin(t)));
        }
    }
    for (int q = 0; q < ncap; ++q)
    {
        const double theta =
            std::numbers::pi / 2.0 *
            (1.0 - static_cast<double>(q + 1) / (ncap + 1));
        const double x = px + sign * R * std::cos(theta);
        const double rad = R * std::sin(theta);
        const int k = nann + q;
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            v[k][i] = m.add_vertex(
                Point(x, rad * std::cos(t), rad * std::sin(t)));
        }
    }
    for (int k = 0; k + 1 < nrings; ++k)
        for (int i = 0; i < n; ++i)
        {
            const int j = (i + 1) % n;
            m.add_triangle(v[k][i], v[k][j], v[k + 1][j]);
            m.add_triangle(v[k][i], v[k + 1][j], v[k + 1][i]);
        }
    const Vertex c = m.add_vertex(Point(px + sign * R, 0, 0));
    for (int i = 0; i < n; ++i)
    {
        const int j = (i + 1) % n;
        m.add_triangle(v[nrings - 1][i], v[nrings - 1][j], c);
    }
    return m;
}

void append_mesh(SurfaceMesh& dst, const SurfaceMesh& src)
{
    std::vector<Vertex> remap(src.n_vertices());
    for (auto v : src.vertices())
        remap[v.idx()] = dst.add_vertex(src.position(v));
    for (auto f : src.faces())
    {
        std::vector<Vertex> fv;
        for (auto v : src.vertices(f))
            fv.push_back(remap[v.idx()]);
        dst.add_face(fv);
    }
}

} // namespace

TEST(HoleFillingTest, open_cone)
{
    auto mesh = open_cone();
    Halfedge h = find_boundary(mesh);
    EXPECT_TRUE(h.is_valid());
    fill_hole(mesh, h);
    h = find_boundary(mesh);
    EXPECT_FALSE(h.is_valid());
}

// A folded/taco-like boundary whose best-plane projection self-intersects.
// The filling must still close the hole and stay embedded and manifold.
TEST(HoleFillingTest, folded_self_intersecting_boundary)
{
    auto mesh = folded_hole_disk();
    EXPECT_EQ(count_self_intersections(mesh), 0); // input is embedded

    fill_hole(mesh, find_boundary(mesh));

    EXPECT_EQ(count_boundary_edges(mesh), 0);
    EXPECT_TRUE(is_manifold_mesh(mesh));
    EXPECT_EQ(count_self_intersections(mesh), 0);
}

// Two holes facing each other across a narrow gap. Each filling patch must
// avoid the opposite surface and the other patch.
TEST(HoleFillingTest, adjacent_facing_holes)
{
    const double gap = 0.2;
    const double px = gap / 2.0;
    SurfaceMesh mesh = facing_perforated_disk(48, -px, -1);
    append_mesh(mesh, facing_perforated_disk(48, px, 1));
    EXPECT_EQ(count_self_intersections(mesh), 0);

    const auto loops = boundary_loops(mesh);
    ASSERT_EQ(loops.size(), 2u);
    for (auto h : loops)
        fill_hole(mesh, h);

    EXPECT_EQ(count_boundary_edges(mesh), 0);
    EXPECT_TRUE(is_manifold_mesh(mesh));
    EXPECT_EQ(count_self_intersections(mesh), 0);
}

// The geometric predicates backing the post-fill intersection check, and the
// boundary unfolding.
TEST(HoleFillingTest, intersection_predicates)
{
    using interdetail::triangles_intersect;
    {
        std::array<Point, 3> A = {Point(-1, 0, 0), Point(1, 0, 0),
                                  Point(0, 1, 0)};
        std::array<Point, 3> B = {Point(0, -1, 0.5), Point(0, 1, 0.5),
                                  Point(0, 0, -0.5)};
        EXPECT_TRUE(triangles_intersect(A, B));
    }
    {
        std::array<Point, 3> A = {Point(0, 0, 0), Point(1, 0, 0),
                                  Point(0, 1, 0)};
        std::array<Point, 3> B = {Point(5, 0, 0), Point(6, 0, 0),
                                  Point(5, 1, 0)};
        EXPECT_FALSE(triangles_intersect(A, B));
    }
    {
        std::array<Point, 3> A = {Point(0, 0, 0), Point(2, 0, 0),
                                  Point(0, 2, 0)};
        std::array<Point, 3> B = {Point(1, 1, 0), Point(3, 1, 0),
                                  Point(1, 3, 0)};
        EXPECT_TRUE(triangles_intersect(A, B));
    }
    {
        // degenerate triangle never reports an intersection
        std::array<Point, 3> A = {Point(0, 0, 0), Point(1, 0, 0),
                                  Point(2, 0, 0)};
        std::array<Point, 3> B = {Point(0, -1, 0.5), Point(0, 1, 0.5),
                                  Point(0, 0, -0.5)};
        EXPECT_FALSE(triangles_intersect(A, B));
    }
}

// Count proper crossings between non-adjacent edges of a closed 2D polyline.
static int polygon_crossings(const std::vector<Vector<Scalar, 2>>& q)
{
    const int n = static_cast<int>(q.size());
    int hits = 0;
    for (int i = 0; i < n; ++i)
    {
        const int i2 = (i + 1) % n;
        for (int j = i + 2; j < n; ++j)
        {
            if (i == 0 && j == n - 1)
                continue;
            if (interdetail::segments_cross2(q[i], q[i2], q[j],
                                             q[(j + 1) % n], Scalar(1e-10)))
                ++hits;
        }
    }
    return hits;
}

TEST(HoleFillingTest, unfold_simple_polygon)
{
    using interdetail::unfold_simple_polygon;
    using Vec2 = Vector<Scalar, 2>;

    // 1. A self-intersecting planar projection (the taco fold:
    //    phi = t + a sin 2t with a > 0.5) must unfold to a simple polygon.
    {
        const int n = 64;
        const double a = 1.0;
        std::vector<Vec2> p(n);
        std::vector<Scalar> lengths(n);
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            const double phi = t + a * std::sin(2.0 * t);
            p[i] = Vec2(std::cos(phi), std::sin(phi));
        }
        for (int i = 0; i < n; ++i)
            lengths[i] = norm(p[(i + 1) % n] - p[i]);
        ASSERT_GT(polygon_crossings(p), 0); // the input really is folded

        const auto q = unfold_simple_polygon(p, lengths);

        ASSERT_EQ(q.size(), static_cast<size_t>(n));
        EXPECT_EQ(polygon_crossings(q), 0);
        // A folded loop is not developable with every edge length preserved.
        // The linear closure ramp perturbs individual edges (the tiny ones
        // near self-touching points most strongly), so the robust guarantees
        // are: exact closure, no collapsed edge, and a preserved total length.
        Scalar mean_len = 0;
        for (int i = 0; i < n; ++i)
            mean_len += lengths[i];
        mean_len /= n;
        Scalar in_perimeter = 0, out_perimeter = 0;
        for (int i = 0; i < n; ++i)
        {
            in_perimeter += lengths[i];
            const Scalar out_len = norm(q[(i + 1) % n] - q[i]);
            out_perimeter += out_len;
            EXPECT_GT(out_len, Scalar(1e-3) * mean_len); // no collapsed edge
        }
        // the closure ramp conserves the perimeter to within a few percent
        EXPECT_NEAR(out_perimeter, in_perimeter, Scalar(0.05) * in_perimeter);
        // the closure ramp brings the two ends within one edge length
        const Scalar closure_gap = norm(q[n - 1] - q[0]);
        EXPECT_LT(closure_gap, mean_len);
        // total turning of the result is one full counterclockwise turn
        Scalar total = 0;
        for (int i = 0; i < n; ++i)
        {
            const Vec2 e_in = q[i] - q[(i - 1 + n) % n];
            const Vec2 e_out = q[(i + 1) % n] - q[i];
            total += std::atan2(e_in[0] * e_out[1] - e_in[1] * e_out[0],
                                dot(e_in, e_out));
        }
        EXPECT_NEAR(total, Scalar(2) * std::numbers::pi, Scalar(1e-4));
    }

    // 2. A simple circle is already convex: it must remain simple and keep
    //    both its edge lengths and (nearly uniform) turning angles.
    {
        const int n = 48;
        std::vector<Vec2> p(n);
        std::vector<Scalar> lengths(n);
        for (int i = 0; i < n; ++i)
        {
            const double t = 2.0 * std::numbers::pi * i / n;
            p[i] = Vec2(std::cos(t), std::sin(t));
        }
        for (int i = 0; i < n; ++i)
            lengths[i] = norm(p[(i + 1) % n] - p[i]);

        const auto q = unfold_simple_polygon(p, lengths);

        EXPECT_EQ(polygon_crossings(q), 0);
        Scalar circle_mean = 0;
        for (int i = 0; i < n; ++i)
            circle_mean += lengths[i];
        circle_mean /= n;
        // the already-convex circle keeps every walked edge within 5% of its
        // length (the closing edge q[n-1]->q[0] absorbs the closure residual)
        for (int i = 0; i < n - 1; ++i)
            EXPECT_NEAR(norm(q[i + 1] - q[i]), lengths[i],
                        Scalar(0.05) * circle_mean);
        // and its ends are brought together to within an edge length
        EXPECT_LT(norm(q[n - 1] - q[0]), circle_mean);
    }
}
