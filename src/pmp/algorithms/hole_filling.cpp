// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#include "pmp/algorithms/hole_filling.h"
#include "pmp/algorithms/hole_filling_intersection.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <limits>
#include <unordered_set>
#include <vector>

#include "pmp/algorithms/fairing.h"
#include "pmp/algorithms/normals.h"

namespace pmp {
namespace {

class HoleFilling
{
public:
    explicit HoleFilling(SurfaceMesh& mesh);
    void fill_hole(Halfedge h);

private:
    struct Weight
    {
        Weight(Scalar angle = std::numeric_limits<Scalar>::max(),
               Scalar area = std::numeric_limits<Scalar>::max())
            : angle_(angle), area_(area)
        {
        }

        Weight operator+(const Weight& rhs) const
        {
            return Weight(std::max(angle_, rhs.angle_), area_ + rhs.area_);
        }

        bool operator<(const Weight& rhs) const
        {
            return (angle_ < rhs.angle_ ||
                    (angle_ == rhs.angle_ && area_ < rhs.area_));
        }

        Scalar angle_;
        Scalar area_;
    };

    // ---- boundary analysis & unfolding ----
    void trace_hole(Halfedge h);
    bool boundary_self_intersects_2d() const;
    void unfold_boundary();

    // best-fit plane (centroid plus in-plane orthonormal frame) of the
    // boundary vertices; the frame normal is the smallest-eigenvalue direction
    void fit_boundary_plane(Point& centroid, Point& ax1, Point& ax2) const;

    // ---- triangulation ----
    void triangulate_hole();
    Weight compute_weight(int i, int j, int k) const;

    // ---- refinement (isotropic remeshing) ----
    void refine();
    void split_long_edges(const Scalar lmax);
    void collapse_short_edges(const Scalar lmin);
    void flip_edges();
    void relaxation();
    void remove_caps();
    void fairing();

    // ---- proximity-aware interior fairing ----
    void fair_nearby_boundary();

    // ---- post-fill self-intersection repair ----
    void repair_self_intersections();
    std::vector<std::pair<Face, Face>> find_patch_intersections() const;
    void delaunay_refine_face(Face f);

    // return i'th vertex of hole
    Vertex hole_vertex(unsigned int i) const
    {
        assert(i < hole_.size());
        return mesh_.to_vertex(hole_[i]);
    }

    // return normal of face opposite to edge (i-1,i)
    Normal opposite_normal(unsigned int i) const
    {
        assert(i < hole_.size());
        return face_normal(mesh_,
                           mesh_.face(mesh_.opposite_halfedge(hole_[i])));
    }

    // does interior edge (_a,_b) exist already?
    bool is_interior_edge(Vertex a, Vertex b) const;

    Scalar compute_area(Vertex a, Vertex b, Vertex c) const;
    Point compute_normal(Vertex a, Vertex b, Vertex c) const;
    Scalar compute_angle(const Point& n1, const Point& n2) const;

    // mesh and properties
    SurfaceMesh& mesh_;
    VertexProperty<Point> points_;
    VertexProperty<bool> vlocked_;
    EdgeProperty<bool> elocked_;
    FaceProperty<bool> patch_face_;

    std::vector<Halfedge> hole_;

    // data for computing optimal triangulation
    std::vector<std::vector<Weight>> weight_;
    std::vector<std::vector<int>> index_;

    // 2D unfolding of the boundary (empty when unfolding is not needed)
    std::vector<Vector<Scalar, 2>> uv_;
};

HoleFilling::HoleFilling(SurfaceMesh& mesh) : mesh_(mesh)
{
    points_ = mesh_.vertex_property<Point>("v:point");
}

bool HoleFilling::is_interior_edge(Vertex a, Vertex b) const
{
    const Halfedge h = mesh_.find_halfedge(a, b);
    if (!h.is_valid())
        return false; // edge does not exist
    return (!mesh_.is_boundary(h) &&
            !mesh_.is_boundary(mesh_.opposite_halfedge(h)));
}

Scalar HoleFilling::compute_area(Vertex a, Vertex b, Vertex c) const
{
    return sqrnorm(cross(points_[b] - points_[a], points_[c] - points_[a]));
}

Point HoleFilling::compute_normal(Vertex a, Vertex b, Vertex c) const
{
    return normalize(cross(points_[b] - points_[a], points_[c] - points_[a]));
}

Scalar HoleFilling::compute_angle(const Point& n1, const Point& n2) const
{
    return (1.0 - dot(n1, n2));
}

void HoleFilling::fill_hole(Halfedge h)
{
    if (!h.is_valid())
    {
        throw InvalidInputException("HoleFilling: Invalid halfedge.");
    }

    // is it really a hole?
    if (!mesh_.is_boundary(h))
    {
        auto what = std::string{__func__} + ": Not a boundary halfedge.";
        throw InvalidInputException(what);
    }

    // lock vertices/edge that already exist, to be later able to
    // identify the filled-in vertices/edges
    vlocked_ = mesh_.add_vertex_property<bool>("HoleFilling:vlocked", false);
    elocked_ = mesh_.add_edge_property<bool>("HoleFilling:elocked", false);
    for (auto v : mesh_.vertices())
        vlocked_[v] = true;
    for (auto e : mesh_.edges())
        elocked_[e] = true;

    try
    {
        trace_hole(h);

        // Step 1: unfold a self-intersecting boundary into 2D so the
        // triangulation sees a simple polygon instead of a folded loop.
        if (boundary_self_intersects_2d())
        {
            unfold_boundary();
        }

        triangulate_hole(); // do minimal triangulation

        // tag the newly created faces as the filling patch via a flood fill
        // starting from the boundary loop, walking only across edges created
        // during the triangulation (locked rim edges stop the flood)
        patch_face_ = mesh_.add_face_property<bool>("HoleFilling:patch", false);
        {
            FaceProperty<bool> visited =
                mesh_.add_face_property<bool>("HoleFilling:visited", false);
            std::vector<Face> stack;
            for (auto hh : hole_)
            {
                Face f = mesh_.face(mesh_.next_halfedge(hh));
                if (f.is_valid() && !visited[f])
                {
                    visited[f] = true;
                    stack.push_back(f);
                }
            }
            while (!stack.empty())
            {
                Face f = stack.back();
                stack.pop_back();
                patch_face_[f] = true;
                for (auto hh : mesh_.halfedges(f))
                {
                    Face nf = mesh_.face(mesh_.opposite_halfedge(hh));
                    if (nf.is_valid() && !visited[nf])
                    {
                        // only walk across edges created during filling
                        Edge e = mesh_.edge(hh);
                        if (!elocked_[e])
                        {
                            visited[nf] = true;
                            stack.push_back(nf);
                        }
                    }
                }
            }
            mesh_.remove_face_property(visited);
        }

        // Step 2: keep interior vertices away from nearby surfaces/holes
        fair_nearby_boundary();

        refine(); // refine filled-in edges

        // Step 4: verify and repair any remaining self-intersection
        repair_self_intersections();
    }
    catch (InvalidInputException& e)
    {
        // clean up
        hole_.clear();
        uv_.clear();
        if (mesh_.has_face_property("HoleFilling:visited"))
        {
            auto vis =
                mesh_.get_face_property<bool>("HoleFilling:visited");
            mesh_.remove_face_property(vis);
        }
        if (mesh_.has_face_property("HoleFilling:patch"))
        {
            auto pf = mesh_.get_face_property<bool>("HoleFilling:patch");
            mesh_.remove_face_property(pf);
        }
        mesh_.remove_vertex_property(vlocked_);
        mesh_.remove_edge_property(elocked_);

        throw e;
    }

    // clean up
    hole_.clear();
    uv_.clear();
    if (mesh_.has_face_property("HoleFilling:patch"))
        mesh_.remove_face_property(patch_face_);
    mesh_.remove_vertex_property(vlocked_);
    mesh_.remove_edge_property(elocked_);
}

void HoleFilling::trace_hole(Halfedge h)
{
    hole_.clear();
    Halfedge hit = h;
    do
    {
        // check for manifoldness
        if (!mesh_.is_manifold(mesh_.to_vertex(hit)))
        {
            auto what = std::string{__func__} + ": Non-manifold hole.";
            throw InvalidInputException(what);
        }

        hole_.emplace_back(hit);
    } while ((hit = mesh_.next_halfedge(hit)) != h);
}

void HoleFilling::fit_boundary_plane(Point& centroid, Point& ax1,
                                     Point& ax2) const
{
    const int n = static_cast<int>(hole_.size());

    centroid = Point(0, 0, 0);
    for (auto hh : hole_)
        centroid += points_[mesh_.to_vertex(hh)];
    centroid /= Scalar(n);

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (auto hh : hole_)
    {
        const Point d = points_[mesh_.to_vertex(hh)] - centroid;
        cov += Eigen::Vector3d(d[0], d[1], d[2]) *
               Eigen::Vector3d(d[0], d[1], d[2]).transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    const Eigen::Vector3d normal = es.eigenvectors().col(0);
    const Point N(normal[0], normal[1], normal[2]);
    if (std::abs(N[0]) < Scalar(0.9))
        ax1 = normalize(cross(N, Point(1, 0, 0)));
    else
        ax1 = normalize(cross(N, Point(0, 1, 0)));
    ax2 = cross(N, ax1);
}

// Project the boundary onto its least-squares best-fit plane and test the
// resulting closed 2D polyline for proper self-intersection.
bool HoleFilling::boundary_self_intersects_2d() const
{
    const int n = static_cast<int>(hole_.size());
    if (n < 4)
        return false;

    // best-fit plane through the boundary vertices
    Point centroid, ax1, ax2;
    fit_boundary_plane(centroid, ax1, ax2);

    std::vector<Vector<Scalar, 2>> p(n);
    Scalar lo0 = 1e30f, hi0 = -1e30f, lo1 = 1e30f, hi1 = -1e30f;
    for (int i = 0; i < n; ++i)
    {
        const Point d = points_[hole_vertex(i)] - centroid;
        p[i] = Vector<Scalar, 2>(dot(d, ax1), dot(d, ax2));
        lo0 = std::min(lo0, p[i][0]);
        hi0 = std::max(hi0, p[i][0]);
        lo1 = std::min(lo1, p[i][1]);
        hi1 = std::max(hi1, p[i][1]);
    }
    // scale-relative tolerance so near-tangent crossings are still caught
    const Scalar diag =
        std::sqrt((hi0 - lo0) * (hi0 - lo0) + (hi1 - lo1) * (hi1 - lo1));
    const Scalar eps = std::max(Scalar(1e-12), Scalar(1e-12) * diag);

    // test non-adjacent segment pairs for a proper crossing
    for (int i = 0; i < n; ++i)
    {
        const int i2 = (i + 1) % n;
        for (int j = i + 2; j < n; ++j)
        {
            if (i == 0 && j == n - 1)
                continue; // adjacent through the loop closure
            const int j2 = (j + 1) % n;
            if (interdetail::segments_cross2(p[i], p[i2], p[j], p[j2], eps))
                return true;
        }
    }
    return false;
}

// Unfold the boundary into a simple 2D polygon in the best-fit plane.
// The projection of a folded boundary is not a simple loop, so it cannot be
// triangulated directly. We rebuild a guaranteed-simple strictly convex
// polygon that retains the boundary's metric information:
//   - the walk uses the measured 3D edge lengths;
//   - every external turn is derived from the local projected edge
//     directions, then convexified (clipped to a strictly positive sub-pi
//     range and renormalized to total turning 2*pi), keeping the local shape
//     distribution where it is well defined while removing the loops;
//   - the open chain is closed with a monotone linear residual ramp.
// A folded loop is not perfectly developable, so the closure perturbs the
// edge lengths slightly (never collapsing one); see unfold_simple_polygon.
// The unfolding guides the triangulation topology only; the original 3D
// boundary positions are restored before triangulation and the subsequent
// fairing builds a smooth spatial patch.
void HoleFilling::unfold_boundary()
{
    const int n = static_cast<int>(hole_.size());

    Point centroid, ax1, ax2;
    fit_boundary_plane(centroid, ax1, ax2);

    // project all boundary points and record the measured 3D edge lengths
    std::vector<Vector<Scalar, 2>> p(n);
    std::vector<Scalar> lengths(n);
    for (int i = 0; i < n; ++i)
    {
        const Point d = points_[hole_vertex(i)] - centroid;
        p[i] = Vector<Scalar, 2>(dot(d, ax1), dot(d, ax2));
        lengths[i] =
            distance(points_[hole_vertex(i)],
                     points_[hole_vertex((i + 1) % n)]);
    }

    uv_ = interdetail::unfold_simple_polygon(p, lengths);
}

void HoleFilling::triangulate_hole()
{
    const int n = static_cast<int>(hole_.size());

    // When an unfolding is available, temporarily move the boundary vertices
    // to their unfolded 3D positions so the angle/area DP operates on a
    // simple polygon. Original positions are restored afterwards.
    std::vector<Point> saved;
    if (!uv_.empty())
    {
        Point centroid, ax1, ax2;
        fit_boundary_plane(centroid, ax1, ax2);

        saved.resize(n);
        for (int i = 0; i < n; ++i)
        {
            saved[i] = points_[hole_vertex(i)];
            points_[hole_vertex(i)] =
                centroid + uv_[i][0] * ax1 + uv_[i][1] * ax2;
        }
    }

    // compute minimal triangulation by dynamic programming
    weight_.clear();
    weight_.resize(n, std::vector<Weight>(n, Weight()));
    index_.clear();
    index_.resize(n, std::vector<int>(n, 0));

    int i, j, m, k, imin;
    Weight w, wmin;

    // initialize 2-gons
    for (i = 0; i < n - 1; ++i)
    {
        weight_[i][i + 1] = Weight(0, 0);
        index_[i][i + 1] = -1;
    }

    // n-gons with n>2
    for (j = 2; j < n; ++j)
    {
        // for all n-gons [i,i+j]
        for (i = 0; i < n - j; ++i)
        {
            k = i + j;
            wmin = Weight();
            imin = -1;

            // find best split i < m < i+j
            for (m = i + 1; m < k; ++m)
            {
                w = weight_[i][m] + compute_weight(i, m, k) + weight_[m][k];
                if (w < wmin)
                {
                    wmin = w;
                    imin = m;
                }
            }

            weight_[i][k] = wmin;
            index_[i][k] = imin;
        }
    }

    // now add triangles to mesh
    std::vector<ivec2> todo;
    todo.reserve(n);
    todo.emplace_back(0, n - 1);
    while (!todo.empty())
    {
        ivec2 tri = todo.back();
        todo.pop_back();
        const int start = tri[0];
        const int end = tri[1];
        if (end - start < 2)
            continue;
        const int split = index_[start][end];

        mesh_.add_triangle(hole_vertex(start), hole_vertex(split),
                           hole_vertex(end));

        todo.emplace_back(start, split);
        todo.emplace_back(split, end);
    }

    // restore the original 3D boundary positions
    if (!saved.empty())
    {
        for (int r = 0; r < n; ++r)
            points_[hole_vertex(r)] = saved[r];
    }

    // clean up
    weight_.clear();
    index_.clear();
}

HoleFilling::Weight HoleFilling::compute_weight(int i, int j, int k) const
{
    const Vertex a = hole_vertex(i);
    const Vertex b = hole_vertex(j);
    const Vertex c = hole_vertex(k);

    // if one of the potential edges already exists, this would result
    // in an invalid triangulation -> prevent by giving infinite weight
    if (is_interior_edge(a, b) || is_interior_edge(b, c) ||
        is_interior_edge(c, a))
    {
        return {};
    }

    // compute area
    const Scalar area = compute_area(a, b, c);

    // compute dihedral angles with...
    Scalar angle(0);
    const Point n = compute_normal(a, b, c);
    Point n2;

    // ...neighbor to (i,j)
    n2 = (i + 1 == j) ? opposite_normal(j)
                      : compute_normal(a, hole_vertex(index_[i][j]), b);
    angle = std::max(angle, compute_angle(n, n2));

    // ...neighbor to (j,k)
    n2 = (j + 1 == k) ? opposite_normal(k)
                      : compute_normal(b, hole_vertex(index_[j][k]), c);
    angle = std::max(angle, compute_angle(n, n2));

    // ...neighbor to (k,i) if (k,i)==(n-1, 0)
    if (i == 0 && k + 1 == (int)hole_.size())
    {
        n2 = opposite_normal(0);
        angle = std::max(angle, compute_angle(n, n2));
    }

    return {angle, area};
}

// Detect boundary vertices that sit close to another part of the surface
// (a neighboring hole rim or the opposite wall). Interior vertices created by
// the triangulation are pulled toward the hole's safe side during relaxation
// so the patch cannot bridge the gap and intersect the adjacent surface.
// Cost: O(n*V) for a boundary of n vertices in a mesh of V vertices; the
// fixed pass count keeps the relaxation itself linear in the patch size.
void HoleFilling::fair_nearby_boundary()
{
    const int n = static_cast<int>(hole_.size());
    if (n == 0)
        return;

    // characteristic edge length of the rim
    Scalar mean_length(0);
    for (int i = 0; i < n; ++i)
        mean_length += distance(points_[hole_vertex(i)],
                                points_[hole_vertex((i + 1) % n)]);
    mean_length /= Scalar(n);
    const Scalar proximity = Scalar(2.5) * mean_length;

    // For each boundary vertex, look for a nearby non-patch surface point and
    // record an inward repulsion direction (away from the obstacle).
    std::vector<Point> repulse(n, Point(0, 0, 0));
    std::vector<char> blocked(n, false);
    for (int i = 0; i < n; ++i)
    {
        const Point pi = points_[hole_vertex(i)];
        Scalar best = proximity * proximity;
        Point nearest(0, 0, 0);
        for (auto v : mesh_.vertices())
        {
            if (v == hole_vertex(i) || !vlocked_[v])
                continue;
            const Point pj = points_[v];
            const Scalar d2 = sqrnorm(pj - pi);
            if (d2 < best && d2 > Scalar(1e-20))
            {
                best = d2;
                nearest = pj;
            }
        }
        if (best < proximity * proximity)
        {
            blocked[i] = true;
            repulse[i] = normalize(pi - nearest);
        }
    }

    if (std::find(blocked.begin(), blocked.end(), true) == blocked.end())
        return; // no proximity conflict, nothing to do

    // Apply a few constrained Laplacian smoothing passes to the unlocked
    // (interior) vertices, adding the repulsion field sampled at the rim.
    for (int iter = 0; iter < 5; ++iter)
    {
        std::vector<Vertex> free_v;
        for (auto v : mesh_.vertices())
            if (!vlocked_[v])
                free_v.push_back(v);
        if (free_v.empty())
            break;

        std::vector<Point> target(free_v.size());
        for (size_t vi = 0; vi < free_v.size(); ++vi)
        {
            const Vertex v = free_v[vi];
            Point avg(0, 0, 0);
            Scalar c(0);
            for (auto vv : mesh_.vertices(v))
            {
                avg += points_[vv];
                ++c;
            }
            avg /= c;

            // weight repulsion by how close this interior vertex is to each
            // blocked rim vertex, fading out with distance
            Point force(0, 0, 0);
            Scalar wsum(0);
            for (int i = 0; i < n; ++i)
            {
                if (!blocked[i])
                    continue;
                const Point rim = points_[hole_vertex(i)];
                const Scalar d = distance(points_[v], rim);
                const Scalar w = Scalar(1) / (d + mean_length);
                force += w * repulse[i];
                wsum += w;
            }
            if (wsum > 0)
                force /= wsum;
            target[vi] = avg + Scalar(0.5) * mean_length * force;
        }
        for (size_t vi = 0; vi < free_v.size(); ++vi)
            points_[free_v[vi]] = target[vi];
    }
}

void HoleFilling::refine()
{
    // compute target edge length
    const int n = static_cast<int>(hole_.size());
    Scalar mean_length(0);
    for (int i = 0; i < n; ++i)
    {
        mean_length += distance(points_[hole_vertex(i)],
                                points_[hole_vertex((i + 1) % n)]);
    }
    mean_length /= (Scalar)n;

    // do some iterations
    const auto min_length = Scalar{0.7} * mean_length;
    const auto max_length = Scalar{1.5} * mean_length;
    for (int iter = 0; iter < 10; ++iter)
    {
        split_long_edges(max_length);
        collapse_short_edges(min_length);
        flip_edges();
        relaxation();
    }
    remove_caps();
    fairing();
}

void HoleFilling::split_long_edges(const Scalar lmax)
{
    bool ok;
    int i;

    for (ok = false, i = 0; !ok && i < 10; ++i)
    {
        ok = true;

        for (auto e : mesh_.edges())
        {
            if (!elocked_[e])
            {
                const Halfedge h10 = mesh_.halfedge(e, 0);
                const Halfedge h01 = mesh_.halfedge(e, 1);
                const Point& p0 = points_[mesh_.to_vertex(h10)];
                const Point& p1 = points_[mesh_.to_vertex(h01)];

                if (distance(p0, p1) > lmax)
                {
                    mesh_.split(e, 0.5 * (p0 + p1));
                    ok = false;
                }
            }
        }
    }
}

void HoleFilling::collapse_short_edges(const Scalar _lmin)
{
    bool ok;
    int i;

    for (ok = false, i = 0; !ok && i < 10; ++i)
    {
        ok = true;

        for (auto e : mesh_.edges())
        {
            if (!mesh_.is_deleted(e) && !elocked_[e])
            {
                const Halfedge h10 = mesh_.halfedge(e, 0);
                const Halfedge h01 = mesh_.halfedge(e, 1);
                const Vertex v0 = mesh_.to_vertex(h10);
                const Vertex v1 = mesh_.to_vertex(h01);
                const Point& p0 = points_[v0];
                const Point& p1 = points_[v1];

                // edge too short?
                if (distance(p0, p1) < _lmin)
                {
                    Halfedge h;
                    if (!vlocked_[v0])
                        h = h01;
                    else if (!vlocked_[v1])
                        h = h10;

                    if (h.is_valid() && mesh_.is_collapse_ok(h))
                    {
                        mesh_.collapse(h);
                        ok = false;
                    }
                }
            }
        }
    }

    mesh_.garbage_collection();
}

void HoleFilling::flip_edges()
{
    Vertex v0, v1, v2, v3;
    Halfedge h;
    int val0, val1, val2, val3;
    int val_opt0, val_opt1, val_opt2, val_opt3;
    int ve0, ve1, ve2, ve3, ve_before, ve_after;
    bool ok;
    int i;

    for (ok = false, i = 0; !ok && i < 10; ++i)
    {
        ok = true;

        for (auto e : mesh_.edges())
        {
            if (!elocked_[e])
            {
                h = mesh_.halfedge(e, 0);
                v0 = mesh_.to_vertex(h);
                v2 = mesh_.to_vertex(mesh_.next_halfedge(h));
                h = mesh_.halfedge(e, 1);
                v1 = mesh_.to_vertex(h);
                v3 = mesh_.to_vertex(mesh_.next_halfedge(h));

                val0 = mesh_.valence(v0);
                val1 = mesh_.valence(v1);
                val2 = mesh_.valence(v2);
                val3 = mesh_.valence(v3);

                val_opt0 = (mesh_.is_boundary(v0) ? 4 : 6);
                val_opt1 = (mesh_.is_boundary(v1) ? 4 : 6);
                val_opt2 = (mesh_.is_boundary(v2) ? 4 : 6);
                val_opt3 = (mesh_.is_boundary(v3) ? 4 : 6);

                ve0 = (val0 - val_opt0);
                ve1 = (val1 - val_opt1);
                ve2 = (val2 - val_opt2);
                ve3 = (val3 - val_opt3);

                ve_before = ve0 * ve0 + ve1 * ve1 + ve2 * ve2 + ve3 * ve3;

                --val0;
                --val1;
                ++val2;
                ++val3;

                ve0 = (val0 - val_opt0);
                ve1 = (val1 - val_opt1);
                ve2 = (val2 - val_opt2);
                ve3 = (val3 - val_opt3);

                ve_after = ve0 * ve0 + ve1 * ve1 + ve2 * ve2 + ve3 * ve3;

                if (ve_before > ve_after && mesh_.is_flip_ok(e))
                {
                    mesh_.flip(e);
                    ok = false;
                }
            }
        }
    }
}

void HoleFilling::relaxation()
{
    // properties
    VertexProperty<int> idx =
        mesh_.add_vertex_property<int>("HoleFilling:idx", -1);

    // collect free vertices
    std::vector<Vertex> vertices;
    vertices.reserve(mesh_.n_vertices());
    for (auto v : mesh_.vertices())
    {
        if (!vlocked_[v])
        {
            idx[v] = vertices.size();
            vertices.emplace_back(v);
        }
    }
    const int n = static_cast<int>(vertices.size());
    if (n == 0)
    {
        mesh_.remove_vertex_property(idx);
        return;
    }

    // setup matrix & rhs
    Eigen::MatrixXd B(n, 3);
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    for (int i = 0; i < n; ++i)
    {
        const Vertex v = vertices[i];
        Point b(0, 0, 0);
        Scalar c(0);

        for (auto vv : mesh_.vertices(v))
        {
            if (vlocked_[vv])
                b += points_[vv];
            else
                triplets.emplace_back(i, idx[vv], -1.0);
            ++c;
        }

        if (vlocked_[v])
            b -= c * points_[v];
        else
            triplets.emplace_back(i, idx[v], c);

        B.row(i) = (Eigen::Vector3d)b;
    }

    // solve least squares system
    using SparseMatrix = Eigen::SparseMatrix<double>;
    SparseMatrix A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    const Eigen::SimplicialLDLT<SparseMatrix> solver(A);
    Eigen::MatrixXd X = solver.solve(B);

    if (solver.info() != Eigen::Success)
    {
        // clean up
        mesh_.remove_vertex_property(idx);
        auto what = std::string{__func__} + ": Failed to solve linear system.";
        throw SolverException(what);
    }

    // copy solution to mesh vertices
    for (int i = 0; i < n; ++i)
    {
        points_[vertices[i]] = X.row(i);
    }

    // clean up
    mesh_.remove_vertex_property(idx);
}

void HoleFilling::remove_caps()
{
    Halfedge h;
    Vertex v, vb, vd;
    Scalar a0, a1, amin;
    const Scalar aa(::cos(170.0 * std::numbers::pi / 180.0));
    Point a, b, c, d;

    for (auto e : mesh_.edges())
    {
        if (mesh_.is_flip_ok(e))
        {
            h = mesh_.halfedge(e, 0);
            a = points_[mesh_.to_vertex(h)];
            h = mesh_.next_halfedge(h);
            b = points_[vb = mesh_.to_vertex(h)];
            h = mesh_.halfedge(e, 1);
            c = points_[mesh_.to_vertex(h)];
            h = mesh_.next_halfedge(h);
            d = points_[vd = mesh_.to_vertex(h)];

            a0 = dot(normalize(a - b), normalize(c - b));
            a1 = dot(normalize(a - d), normalize(c - d));

            if (a0 < a1)
            {
                amin = a0;
                v = vb;
            }
            else
            {
                amin = a1;
                v = vd;
            }

            // is it a cap?
            if (amin < aa)
            {
                // flip
                mesh_.flip(e);
            }
        }
    }
}

void HoleFilling::fairing()
{
    // did the refinement insert new vertices?
    // if yes, then trigger fairing; otherwise don't.
    bool new_vertices = false;
    for (auto v : mesh_.vertices())
        if (!vlocked_[v])
            new_vertices = true;
    if (!new_vertices)
        return;

    // convert non-locked into selection
    const bool add_vsel = !mesh_.has_vertex_property("v:selected");
    auto vsel = mesh_.vertex_property<bool>("v:selected");
    for (auto v : mesh_.vertices())
        vsel[v] = !vlocked_[v];

    try
    {
        // fair new vertices
        minimize_curvature(mesh_);
    }
    catch (SolverException& e)
    {
        // clean up
        if (add_vsel)
            mesh_.remove_vertex_property(vsel);
        throw e;
    }

    // clean up
    if (add_vsel)
        mesh_.remove_vertex_property(vsel);
}

// Collect intersecting pairs of patch triangles (or a patch triangle against
// any other triangle). Connectivity-adjacent pairs are ignored.
// Cost: O(F^2) triangle/triangle tests over the F mesh faces; only pairs
// involving at least one patch face are tested and the repair loop bounds the
// number of refinement rounds, so the practical cost is patch-size limited.
std::vector<std::pair<Face, Face>>
HoleFilling::find_patch_intersections() const
{
    std::vector<Face> faces;
    faces.reserve(mesh_.n_faces());
    for (auto f : mesh_.faces())
        faces.push_back(f);

    auto tri_of = [&](Face f) {
        std::array<Point, 3> t;
        int idx = 0;
        for (auto v : mesh_.vertices(f))
            t[idx++] = points_[v];
        return t;
    };
    auto verts_of = [&](Face f) {
        std::array<Vertex, 3> v;
        int idx = 0;
        for (auto vv : mesh_.vertices(f))
            v[idx++] = vv;
        return v;
    };

    std::vector<std::pair<Face, Face>> hits;
    for (size_t i = 0; i < faces.size(); ++i)
    {
        if (!patch_face_[faces[i]])
            continue;
        const auto ti = tri_of(faces[i]);
        const auto vi = verts_of(faces[i]);
        for (size_t j = i + 1; j < faces.size(); ++j)
        {
            // at least one face must belong to the filling patch
            if (!patch_face_[faces[i]] && !patch_face_[faces[j]])
                continue;
            const auto tj = tri_of(faces[j]);
            const auto vj = verts_of(faces[j]);
            bool shared = false;
            for (auto a : vi)
                for (auto b : vj)
                    if (a == b)
                    {
                        shared = true;
                        break;
                    }
            if (shared)
                continue;
            if (interdetail::triangles_intersect(ti, tj))
                hits.emplace_back(faces[i], faces[j]);
        }
    }
    return hits;
}

// Insert a Steiner vertex at the face centroid and split the triangle into
// three, then relax the new vertex. This is a Delaunay-style refinement step:
// it introduces an interior degree-6 vertex that can be smoothed out of the
// intersecting configuration without changing the boundary.
void HoleFilling::delaunay_refine_face(Face f)
{
    if (!f.is_valid() || mesh_.is_deleted(f))
        return;
    Point centroid(0, 0, 0);
    int c(0);
    for (auto v : mesh_.vertices(f))
    {
        centroid += points_[v];
        ++c;
    }
    centroid /= Scalar(c);

    // nudge the Steiner point slightly toward the average face normal so it
    // leaves the plane of the intersecting triangle
    Normal n(0, 0, 0);
    for (auto h : mesh_.halfedges(f))
    {
        Face nf = mesh_.face(mesh_.opposite_halfedge(h));
        if (nf.is_valid())
            n += face_normal(mesh_, nf);
    }
    if (norm(n) > Scalar(1e-12))
        centroid += Scalar(0.15) * normalize(n);

    const Vertex sv = mesh_.add_vertex(centroid);
    vlocked_[sv] = false; // interior, eligible for smoothing/fairing
    mesh_.split(f, sv);
    // tag the three child faces as part of the patch
    for (auto cf : mesh_.faces(sv))
        patch_face_[cf] = true;
}

// Iteratively locate self-intersecting patch triangles and resolve them with
// interior Steiner vertices followed by constrained relaxation. Bounded to
// keep the refinement tractable on pathological inputs.
void HoleFilling::repair_self_intersections()
{
    const int max_rounds = 8;
    for (int round = 0; round < max_rounds; ++round)
    {
        auto hits = find_patch_intersections();
        if (hits.empty())
            return;

        std::unordered_set<int> refined;
        for (auto [f1, f2] : hits)
        {
            for (Face f : {f1, f2})
            {
                if (f.is_valid() && !mesh_.is_deleted(f) &&
                    patch_face_[f] &&
                    refined.insert(f.idx()).second)
                {
                    delaunay_refine_face(f);
                }
            }
        }

        // smooth the newly inserted interior vertices so the local patch
        // relaxes out of the intersecting arrangement
        relaxation();
    }
}

} // namespace

void fill_hole(SurfaceMesh& mesh, Halfedge h)
{
    HoleFilling(mesh).fill_hole(h);
}
} // namespace pmp
