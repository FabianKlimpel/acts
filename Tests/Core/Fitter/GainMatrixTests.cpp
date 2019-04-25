// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE GainMatrix Tests
#include <boost/optional/optional_io.hpp>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <memory>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Fitter/GainMatrixUpdator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {

  using Identifier = unsigned long int;
  using Jacobian   = BoundParameters::CovMatrix_t;

  template <ParID_t... params>
  using MeasurementType = Measurement<Identifier, params...>;
  using TrackState      = TrackState<Identifier, BoundParameters>;

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  BOOST_AUTO_TEST_CASE(gain_matrix_updator)
  {
    // Make dummy measurement
    auto cylinder = Surface::makeShared<CylinderSurface>(nullptr, 3, 10);

    ActsSymMatrixD<2> cov;
    cov << 0.04, 0, 0, 0.1;
    TrackState mState(MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1>(
        cylinder, 0, std::move(cov), -0.1, 0.45));

    // Make dummy track parameter
    ActsSymMatrixD<5> covTrk_def;
    covTrk_def << 0.08, 0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 1;
    TrackSymMatrix covTrk;
    covTrk.block<5, 5>(0, 0) = covTrk_def;

    TrackVector parValues = TrackVector::Zero();
    parValues(0)          = 0.3;
    parValues(1)          = 0.5;
    parValues(2)          = 0.5 * M_PI;
    parValues(3)          = 0.3 * M_PI;
    parValues(4)          = 0.01 * M_PI;
    BoundParameters pars(
        tgContext,
        std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
        parValues,
        cylinder);

    // "update" track state with "prediction"
    mState.parameter.predicted  = std::move(pars);
    mState.parameter.jacobian   = Jacobian::Identity();
    mState.parameter.pathLength = 0.;

    // Gain matrix update and filtered state
    GainMatrixUpdator<BoundParameters> gmu;

    BOOST_CHECK(!mState.parameter.filtered);
    BOOST_CHECK(!mState.measurement.calibrated);
    BOOST_CHECK(gmu(tgContext, mState));
    // filtered is set now
    BOOST_CHECK(!!mState.parameter.filtered);
    // measurement was calibrated
    BOOST_CHECK(!!mState.measurement.calibrated);
    // ref surface is same on measurements and parameters
    BOOST_CHECK_EQUAL(
        MeasurementHelpers::getSurface(*mState.measurement.calibrated),
        cylinder.get());
    BOOST_CHECK_EQUAL(&(*mState.parameter.filtered).referenceSurface(),
                      cylinder.get());

    // assert contents of mState were updated
  }

}  // namespace Test
}  // namespace Acts
