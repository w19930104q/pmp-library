// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#pragma once

#include "pmp/surface_mesh.h"

namespace pmp {

//! \brief Fill the hole specified by halfedge \p h
//! \details Close holes (boundary loops of manifold vertices) robustly,
//! including complex, non-simple holes. The pipeline is:
//! 1. Analyze the hole boundary. If its 3D->2D projection self-intersects
//!    (a folded/taco-like boundary), unfold the boundary into a locally
//!    consistent 2D parameterization before triangulation.
//! 2. When the boundary runs close to another hole or the opposite surface,
//!    fair the interpolated interior vertices so the filling patch does not
//!    punch through the neighboring surface.
//! 3. Fill the hole with an angle/area-minimizing triangulation (Liepa's
//!    dynamic programming), followed by isotropic remeshing.
//! 4. Verify that the new triangles do not self-intersect. Any remaining
//!    intersection is removed by inserting interior Steiner vertices with
//!    Delaunay-style refinement.
//! The result is always manifold (every edge has two incident faces) and,
//! for embedded input meshes, free of self-intersections.
//! See \cite liepa_2003_filling for the underlying triangulation.
//! \pre The specified halfedge is valid.
//! \pre The specified halfedge is a boundary halfedge.
//! \pre The specified halfedge is not adjacent to a non-manifold hole.
//! \throw InvalidInputException in case on of the input preconditions is violated
//! \note This algorithm works on general polygon meshes.
//! \ingroup algorithms
void fill_hole(SurfaceMesh& mesh, Halfedge h);

} // namespace pmp
