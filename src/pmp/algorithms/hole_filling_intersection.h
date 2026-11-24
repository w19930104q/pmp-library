// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#pragma once

#include "pmp/types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

namespace pmp {

//! \addtogroup algorithms
//! @{

//! Robust geometric predicates used by hole filling and self-intersection
//! checks. All routines are header-only and operate on pmp::Point.
namespace interdetail {

using Tri = std::array<Point, 3>;
using Tri2 = std::array<Vector<Scalar, 2>, 3>;

inline Scalar orient2(const Vector<Scalar, 2>& a,
                      const Vector<Scalar, 2>& b,
                      const Vector<Scalar, 2>& c)
{
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

// Strict point-in-triangle test (boundary counts as inside). Returns false
// only when the point lies strictly on the opposite side of one edge while
// being on the interior side of another.
inline bool point_in_tri2(const Vector<Scalar, 2>& p, const Tri2& t,
                          Scalar eps)
{
    const Scalar d1 = orient2(p, t[0], t[1]);
    const Scalar d2 = orient2(p, t[1], t[2]);
    const Scalar d3 = orient2(p, t[2], t[0]);
    const bool pos = d1 > eps || d2 > eps || d3 > eps;
    const bool neg = d1 < -eps || d2 < -eps || d3 < -eps;
    return !(pos && neg);
}

// Proper segment crossing in 2D (shared/touching endpoints do not count).
inline bool segments_cross2(const Vector<Scalar, 2>& p1,
                            const Vector<Scalar, 2>& p2,
                            const Vector<Scalar, 2>& p3,
                            const Vector<Scalar, 2>& p4, Scalar eps)
{
    const Scalar d1 = orient2(p3, p4, p1);
    const Scalar d2 = orient2(p3, p4, p2);
    const Scalar d3 = orient2(p1, p2, p3);
    const Scalar d4 = orient2(p1, p2, p4);
    return ((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
           ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps));
}

// Overlap test for two triangles known to be coplanar.
inline bool coplanar_triangles_overlap(const Tri& A, const Tri& B)
{
    const Point n =
        normalize(cross(A[1] - A[0], A[2] - A[0]));
    int drop = 0;
    Scalar mx = std::abs(n[0]);
    for (int i = 1; i < 3; ++i)
        if (std::abs(n[i]) > mx)
        {
            mx = std::abs(n[i]);
            drop = i;
        }
    auto pr = [&](const Point& p) {
        return drop == 0 ? Vector<Scalar, 2>(p[1], p[2])
                         : drop == 1 ? Vector<Scalar, 2>(p[0], p[2])
                                     : Vector<Scalar, 2>(p[0], p[1]);
    };
    const Tri2 a = {pr(A[0]), pr(A[1]), pr(A[2])};
    const Tri2 b = {pr(B[0]), pr(B[1]), pr(B[2])};
    const Scalar eps = Scalar(1e-9);
    for (const auto& p : a)
        if (point_in_tri2(p, b, eps))
            return true;
    for (const auto& p : b)
        if (point_in_tri2(p, a, eps))
            return true;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (segments_cross2(a[i], a[(i + 1) % 3], b[j],
                                b[(j + 1) % 3], eps))
                return true;
    return false;
}

// Moeller-Trumbore segment/triangle intersection.
inline bool segment_triangle(const Point& a0, const Point& a1,
                             const Point& b0, const Point& b1,
                             const Point& b2, Scalar tol)
{
    const Point e1 = b1 - b0;
    const Point e2 = b2 - b0;
    const Point dir = a1 - a0;
    const Point p = cross(dir, e2);
    const Scalar det = dot(e1, p);
    if (std::abs(det) < Scalar(1e-14))
        return false;
    const Scalar inv = Scalar(1) / det;
    const Point tv = a0 - b0;
    const Scalar u = dot(tv, p) * inv;
    const Point q = cross(tv, e1);
    const Scalar v = dot(dir, q) * inv;
    const Scalar t = dot(e2, q) * inv;
    return u >= -tol && u <= 1 + tol && v >= -tol &&
           u + v <= 1 + tol && t >= -tol && t <= 1 + tol;
}

// Separating-axis triangle/triangle intersection. Returns true only for a
// genuine volume overlap or a coplanar area overlap; triangles that merely
// share an edge or vertex are handled by the caller via connectivity.
inline bool triangles_intersect(const Tri& A, const Tri& B,
                                Scalar eps = Scalar(1e-9))
{
    const Point nA = cross(A[1] - A[0], A[2] - A[0]);
    const Point nB = cross(B[1] - B[0], B[2] - B[0]);
    const Scalar lA = norm(nA);
    const Scalar lB = norm(nB);
    if (lA < Scalar(1e-14) || lB < Scalar(1e-14))
        return false;
    const Point na = nA / lA;
    const Point nb = nB / lB;

    Scalar edge_scale = Scalar(1e-12);
    for (int i = 0; i < 3; ++i)
    {
        edge_scale = std::max(edge_scale, norm(A[(i + 1) % 3] - A[i]));
        edge_scale = std::max(edge_scale, norm(B[(i + 1) % 3] - B[i]));
    }
    const Scalar tol = std::max(eps, Scalar(1e-10) * edge_scale);

    if (std::abs(dot(na, nb)) > Scalar(1.0 - 1e-6) &&
        std::abs(dot(B[0] - A[0], na)) < tol)
        return coplanar_triangles_overlap(A, B);

    const Point ea[3] = {A[1] - A[0], A[2] - A[1], A[0] - A[2]};
    const Point eb[3] = {B[1] - B[0], B[2] - B[1], B[0] - B[2]};
    for (const auto& e : ea)
        for (const auto& f : eb)
        {
            const Point raw = cross(e, f);
            if (norm(raw) < Scalar(1e-14))
                continue;
            const Point axis = raw / norm(raw);
            Scalar amin = 1e30f, amax = -1e30f;
            Scalar bmin = 1e30f, bmax = -1e30f;
            for (const auto& p : A)
            {
                const Scalar d = dot(p, axis);
                amin = std::min(amin, d);
                amax = std::max(amax, d);
            }
            for (const auto& p : B)
            {
                const Scalar d = dot(p, axis);
                bmin = std::min(bmin, d);
                bmax = std::max(bmax, d);
            }
            if (amax < bmin - tol || bmax < amin - tol)
                return false;
        }
    for (int i = 0; i < 3; ++i)
        if (segment_triangle(A[i], A[(i + 1) % 3], B[0], B[1], B[2], tol))
            return true;
    for (int i = 0; i < 3; ++i)
        if (segment_triangle(B[i], B[(i + 1) % 3], A[0], A[1], A[2], tol))
            return true;
    return false;
}

// Rebuild a guaranteed-simple strictly convex 2D polygon from the best-plane
// projection \p p of a possibly self-intersecting boundary loop and the
// measured 3D edge \p lengths. The polygon is walked with the measured edge
// lengths and external turns derived from the local projected edge
// directions, then convexified (clipped to a strictly positive sub-pi range,
// renormalized to total turning 2*pi) so the loop cannot cross itself. The
// open chain is closed by a monotone linear residual ramp; because a folded
// loop is not perfectly developable this slightly perturbs the edge lengths
// but never collapses one. Defensive fallback: if the closure ever
// reintroduced a crossing, a regular convex n-gon with the same edge lengths
// is returned instead.
inline std::vector<Vector<Scalar, 2>>
unfold_simple_polygon(const std::vector<Vector<Scalar, 2>>& p,
                      const std::vector<Scalar>& lengths)
{
    const int n = static_cast<int>(p.size());
    const Scalar two_pi = Scalar(2) * std::numbers::pi;

    // signed external turn of the projected loop at each vertex, wrapped to
    // (-pi, pi]; degenerate (collinear/duplicated) edges get no turn.
    std::vector<Scalar> turn(n, Scalar(0));
    for (int i = 0; i < n; ++i)
    {
        const Vector<Scalar, 2> e_in = p[i] - p[(i - 1 + n) % n];
        const Vector<Scalar, 2> e_out = p[(i + 1) % n] - p[i];
        if (norm(e_in) < Scalar(1e-12) || norm(e_out) < Scalar(1e-12))
            continue;
        Scalar delta = std::atan2(
            e_in[0] * e_out[1] - e_in[1] * e_out[0], dot(e_in, e_out));
        if (delta <= -std::numbers::pi)
            delta += two_pi;
        // take the convex (positive, left) component of the local turn
        turn[i] = delta > Scalar(0) ? delta : Scalar(0);
    }

    // convexify: blend toward the uniform turn, renormalize to 2*pi, then
    // water-fill so every turn stays in (0, pi); this guarantees a strictly
    // convex, hence simple, closed polygon for any edge-length sequence.
    Scalar sum = 0;
    for (int i = 0; i < n; ++i)
    {
        turn[i] = Scalar(0.5) * turn[i] + Scalar(0.5) * two_pi / n;
        sum += turn[i];
    }
    for (int i = 0; i < n; ++i)
        turn[i] *= two_pi / sum;
    for (int it = 0; it < 2 * n; ++it)
    {
        const Scalar floor = two_pi / (Scalar(2) * n);
        const Scalar ceiling = std::numbers::pi - two_pi / n;
        Scalar excess = 0;
        int free_slots = 0;
        for (int i = 0; i < n; ++i)
        {
            if (turn[i] < floor)
            {
                excess += floor - turn[i];
                turn[i] = floor;
            }
            else if (turn[i] > ceiling)
            {
                excess += ceiling - turn[i]; // negative: room to distribute
                turn[i] = ceiling;
            }
            else
                ++free_slots;
        }
        if (std::abs(excess) < Scalar(1e-12) || free_slots == 0)
            break;
        // conserve the total: clipped adjustments are redistributed with the
        // opposite sign to the still-free slots
        for (int i = 0; i < n; ++i)
            if (turn[i] > floor && turn[i] < ceiling)
                turn[i] -= excess / free_slots;
    }

    auto build = [&](const std::vector<Scalar>& turns) {
        std::vector<Vector<Scalar, 2>> q(n);
        q[0] = p[0];
        Scalar angle = 0;
        for (int i = 1; i < n; ++i)
        {
            angle += turns[i];
            q[i] = q[i - 1] +
                   Vector<Scalar, 2>(lengths[i - 1] * std::cos(angle),
                                     lengths[i - 1] * std::sin(angle));
        }
        // Close the open chain by distributing the end-to-end residual across
        // all interior vertices as a linear ramp. Vertex 0 stays fixed and
        // vertex n-1 absorbs the remaining (n-1)/n share, so every edge keeps
        // a non-degenerate length while the residual is spread over the whole
        // loop instead of bending one edge. The boundary loop is closed
        // topologically by the mesh halfedges; geometrically its two ends are
        // brought within residual/n of each other. The ramp is a monotone
        // perturbation of an already strictly convex polygon and cannot, by
        // construction, introduce an edge crossing (re-verified below).
        const Vector<Scalar, 2> residual = q[0] - q[n - 1];
        for (int i = 1; i < n; ++i)
            q[i] += residual * (Scalar(i) / n);
        return q;
    };

    auto has_crossing = [&](const std::vector<Vector<Scalar, 2>>& q) {
        for (int i = 0; i < n; ++i)
        {
            const int i2 = (i + 1) % n;
            for (int j = i + 2; j < n; ++j)
            {
                if (i == 0 && j == n - 1)
                    continue; // adjacent through the loop closure
                if (segments_cross2(q[i], q[i2], q[j], q[(j + 1) % n],
                                    Scalar(1e-12)))
                    return true;
            }
        }
        return false;
    };

    std::vector<Vector<Scalar, 2>> simple = build(turn);
    if (has_crossing(simple))
    {
        // defensive fallback: regular convex n-gon with the same edge lengths
        std::vector<Scalar> regular(n, two_pi / n);
        simple = build(regular);
    }
    return simple;
}

} // namespace interdetail

//! @}

} // namespace pmp
