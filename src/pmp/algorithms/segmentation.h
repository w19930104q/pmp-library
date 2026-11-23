// Copyright 2011-2026 the Polygon Mesh Processing Library developers.
// SPDX-License-Identifier: MIT

#pragma once

#include "pmp/surface_mesh.h"

namespace pmp {

//! \brief Segment the mesh into regions of similar curvature.
//! \details Computes per-vertex mean and Gaussian curvature (see curvature())
//! and classifies each vertex as convex, concave, flat, or saddle based on
//! the curvature signs. The classification is regularized by mode filtering
//! and connected regions of the same class are found by region growing,
//! resulting in an initial over-segmentation. Small segments are then merged
//! into the adjacent segment with the most similar mean curvature. The
//! resulting per-face segment IDs are stored in the face property
//! "f:segment" (IDs are consecutive, starting at 0).
//!
//! The classification thresholds adapt to the bounding box diagonal, the
//! mesh resolution, and a robust estimate of the curvature noise level, and
//! the curvature values are smoothed before classification. This makes the
//! segmentation robust against small perturbations of the input geometry.
//!
//! \param mesh The input mesh. Works on general polygon meshes.
//! \param min_segment_size Segments consisting of fewer vertices are merged
//!        into a neighboring segment. If 0, a value proportional to the mesh
//!        resolution is chosen automatically.
//! \param smoothing_steps Number of curvature smoothing steps. Increase for
//!        noisier meshes.
//! \ingroup algorithms
void curvature_segmentation(SurfaceMesh& mesh,
                            unsigned int min_segment_size = 0,
                            int smoothing_steps = 3);

} // namespace pmp
