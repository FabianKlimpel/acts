// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <cmath>
#include <functional>
#include <variant>

namespace Acts {
	
/// @brief These functions perform the transport of a covariance matrix using
/// given Jacobians. The required data is provided by the stepper object
/// with some additional data. Since this is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a seperate file.
namespace detail {

using Jacobian = std::variant<BoundMatrix, FreeToBoundMatrix, BoundToFreeMatrix, FreeMatrix>;

/// Create and return the bound state at the current position
///
/// @brief It does not check if the transported state is at the surface, this
/// needs to be guaranteed by the propagator
///
/// @param [in] geoContext The geometry context
/// @param [in, out] covarianceMatrix The covariance matrix of the state
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] transportJacobian Global jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] covTransport Decision whether the covariance transport should be
/// performed
/// @param [in] accumulatedPath Propagated distance
/// @param [in] surface Target surface on which the state is represented
///
/// @return A bound state:
///   - the parameters at the surface
///   - the stepwise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
std::tuple<BoundParameters, JacobianToBoundPars, double> boundState(
    std::reference_wrapper<const GeometryContext> geoContext,
    BoundSymMatrix& covarianceMatrix, BoundMatrix& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal, const FreeVector& parameters,
    bool covTransport, double accumulatedPath, const Surface& surface);

/// Create and return a curvilinear state at the current position
///
/// @brief This creates a curvilinear state.
///
/// @param [in, out] covarianceMatrix The covariance matrix of the state
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] transportJacobian Global jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] covTransport Decision whether the covariance transport should be
/// performed
/// @param [in] accumulatedPath Propagated distance
///
/// @return A curvilinear state:
///   - the curvilinear parameters at given position
///   - the stepweise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
std::tuple<CurvilinearParameters, JacobianToBoundPars, double> curvilinearState(
    BoundSymMatrix& covarianceMatrix, BoundMatrix& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal, const FreeVector& parameters,
    bool covTransport, double accumulatedPath);

  /// Create and return a free state at the current position
  ///
  /// @brief This creates a free state.
  ///
    /// @tparam result_t Defines the return type
  /// @param [in] state State that will be presented as @c FreeState
  /// @param [in] reinitialize Boolean flag whether reinitialization is needed,
  /// i.e. if this is an intermediate state of a larger propagation
  ///
  /// @return A free state:
  ///   - the free parameters at given position
  ///   - the stepweise jacobian towards it (from last location)
  ///   - and the path length (from start - for ordering)
  std::tuple<FreeParameters, Jacobian, double>
freeState(StepperState& state);

/// @brief Method for on-demand transport of the covariance to a new frame at current position in parameter space. It treats different scenarios:
/// - from local to local
/// - from local to global
/// - from global to local
/// - from global to global
///
/// @param [in] geoContext The geometry context
/// @param [in, out] covarianceMatrix The covariance matrix of the state
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] transportJacobian Global jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] surface is the surface to which the covariance is
///        forwarded to
/// @note No check is done if the position is actually on the surface
void covarianceTransport(
    std::reference_wrapper<const GeometryContext> geoContext,
    BoundSymMatrix& covarianceMatrix, BoundMatrix& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& jacobianLocalToGlobal, const FreeVector& parameters,
    const Surface& surface);

/// @brief Method for on-demand transport of the covariance to a new frame at
/// current position in parameter space
///
/// @param [in, out] covarianceMatrix The covariance matrix of the state
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] transportJacobian Global jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
/// @param [in] toLocal Specifies whether the target is bound or free
void covarianceTransport(BoundSymMatrix& covarianceMatrix,
                         BoundMatrix& jacobian, FreeMatrix& transportJacobian,
                         FreeVector& derivatives,
                         BoundToFreeMatrix& jacobianLocalToGlobal,
                         const Vector3D& direction, bool toLocal = true);
}  // namespace detail
}  // namespace Acts
