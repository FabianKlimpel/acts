// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <cmath>

// Acts includes
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterTypes.hpp"

namespace Acts {
enum ParDef : unsigned int {
  eLOC_0    = 0,  ///< first coordinate in local surface frame
  eLOC_1    = 1,  ///< second coordinate in local surface frame
  eLOC_R    = eLOC_0,
  eLOC_PHI  = eLOC_1,
  eLOC_RPHI = eLOC_0,
  eLOC_Z    = eLOC_1,
  eLOC_X    = eLOC_0,
  eLOC_Y    = eLOC_1,
  eLOC_D0   = eLOC_0,
  eLOC_Z0   = eLOC_1,
  ePHI      = 2,  ///< phi direction of momentum in global frame
  eTHETA    = 3,  ///< theta direction of momentum in global frame
  eQOP = 4,  ///< charge/momentum for charged tracks, for neutral tracks it is
             /// 1/momentum
  eT = 5,
  TrackParsDim
};

constexpr unsigned int GlobalParsDim = 8;

using ParID_t    = ParDef;
using ParValue_t = double;

using TrackVector         = ActsVector<ParValue_t, TrackParsDim>;
using TrackRowVector      = ActsRowVector<ParValue_t, TrackParsDim>;
using GlobalVector        = ActsVector<ParValue_t, GlobalParsDim>;
using TrackToGlobalMatrix = ActsMatrix<ParValue_t, GlobalParsDim, TrackParsDim>;
using GlobalToTrackMatrix = ActsMatrix<ParValue_t, TrackParsDim, GlobalParsDim>;
using TrackMatrix         = ActsMatrix<ParValue_t, TrackParsDim, TrackParsDim>;
using TrackSymMatrix      = ActsSymMatrix<ParValue_t, TrackParsDim>;
using GlobalMatrix = ActsMatrix<ParValue_t, GlobalParsDim, GlobalParsDim>;

template <ParID_t>
struct par_type;

template <ParID_t par>
using par_type_t = typename par_type<par>::type;

template <>
struct par_type<ParDef::eLOC_0>
{
  using type = local_parameter;
};

template <>
struct par_type<ParDef::eLOC_1>
{
  using type = local_parameter;
};

template <>
struct par_type<ParDef::ePHI>
{
  static constexpr double
  pMin()
  {
    return -M_PI;
  }
  static constexpr double
  pMax()
  {
    return M_PI;
  }
  using type = cyclic_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eTHETA>
{
  static constexpr double
  pMin()
  {
    return 0;
  }
  static constexpr double
  pMax()
  {
    return M_PI;
  }
  using type = bound_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eQOP>
{
  using type = unbound_parameter;
};

template <>
struct par_type<ParDef::eT>
{
  using type = unbound_parameter;
};
}  // namespace Acts