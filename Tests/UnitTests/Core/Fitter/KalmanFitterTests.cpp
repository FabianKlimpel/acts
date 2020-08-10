// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Fitter/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <boost/math/distributions/chi_squared.hpp>

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

struct SourceLink
{
  using MeasurementType = std::variant<Measurement<SourceLink, BoundParametersIndices, eLOC_0>, 
						  Measurement<SourceLink, BoundParametersIndices,eLOC_1>, 
						  Measurement<SourceLink, BoundParametersIndices,eLOC_0, eLOC_1>, 
						  Measurement<SourceLink, FreeParametersIndices, eFreePos0, eFreePos1, eFreePos2>>;

  const MeasurementType* meas{nullptr};

  bool operator==(const SourceLink& rhs) const {return meas == rhs.meas;}

  const GeometryObject& referenceObject() const {return *MeasurementHelpers::getObject(*meas);}

  const MeasurementType& operator*() const {return *meas;}
};

// A few initialisations and definitionas
//~ using SourceLink = MinimalSourceLink;
using Covariance = BoundSymMatrix;

using Resolution = std::pair<ParID_t, double>;
using ElementResolution = std::vector<Resolution>;
using VolumeResolution = std::map<GeometryID::Value, ElementResolution>;
using DetectorResolution = std::map<GeometryID::Value, VolumeResolution>;
using DebugOutput = DebugOutputActor;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

ActsSymMatrixD<1> cov1D;
SymMatrix2D cov2D;

bool debugMode = false;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();
CalibrationContext calContext = CalibrationContext();

template <ParID_t... params>
using MeasurementType =
    Measurement<SourceLink, BoundParametersIndices, params...>;

/// @brief This struct creates FittableMeasurements on the
/// detector surfaces, according to the given smearing xxparameters
///
struct MeasurementCreator {
  /// @brief Constructor
  MeasurementCreator() = default;

  /// The detector resolution
  DetectorResolution detectorResolution;

  struct this_result {
    // The measurements
    std::vector<SourceLink::MeasurementType> measurements;

    // The outliers
    std::vector<SourceLink::MeasurementType> outliers;
    
    std::vector<SourceLink::MeasurementType> freeMeasurements;
  };

  using result_type = this_result;

  /// @brief Operater that is callable by an ActionList. The function collects
  /// the surfaces
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] result Vector of matching surfaces
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // monitor the current surface
    auto surface = state.navigation.currentSurface;
    if (surface and surface->associatedDetectorElement()) {
      auto geoID = surface->geoID();
      auto volumeID = geoID.volume();
      auto layerID = geoID.layer();
      // find volume and layer information for this
      auto vResolution = detectorResolution.find(volumeID);
      if (vResolution != detectorResolution.end()) {
        // find layer resolutions
        auto lResolution = vResolution->second.find(layerID);
        if (lResolution != vResolution->second.end()) {
          // Apply global to local
          Acts::Vector2D lPos;
          surface->globalToLocal(state.geoContext,
                                 stepper.position(state.stepping),
                                 stepper.direction(state.stepping), lPos);
          if (lResolution->second.size() == 1) {
            double sp = lResolution->second[0].second;
            cov1D << sp * sp;
            double dp = sp * gauss(generator);
            if (lResolution->second[0].first == eLOC_0) {
              // push back & move a LOC_0 measurement
              MeasurementType<eLOC_0> m0(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_0] + dp);
              result.measurements.push_back(std::move(m0));
              // push back & move a LOC_0 outlier
              MeasurementType<eLOC_0> o0(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_0] + sp * 10);
              result.outliers.push_back(std::move(o0));
            } else {
              // push back & move a LOC_1 measurement
              MeasurementType<eLOC_1> m1(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_1] + dp);
              result.measurements.push_back(std::move(m1));
              // push back & move a LOC_1 outlier
              MeasurementType<eLOC_1> o1(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_1] + sp * 10);
              result.outliers.push_back(std::move(o1));
            }
          } else if (lResolution->second.size() == 2) {
            // Create the measurment and move it
            double sx = lResolution->second[eLOC_0].second;
            double sy = lResolution->second[eLOC_1].second;
            cov2D << sx * sx, 0., 0., sy * sy;
            double dx = sx * gauss(generator);
            double dy = sy * gauss(generator);
            // push back & move a LOC_0, LOC_1 measurement
            MeasurementType<eLOC_0, eLOC_1> m01(surface->getSharedPtr(), {},
                                                cov2D, lPos[eLOC_0] + dx,
                                                lPos[eLOC_1] + dy);
            result.measurements.push_back(std::move(m01));
            // push back & move a LOC_0, LOC_1 outlier
            MeasurementType<eLOC_0, eLOC_1> o01(surface->getSharedPtr(), {},
                                                cov2D, lPos[eLOC_0] + sx * 10,
                                                lPos[eLOC_1] + sy * 10);
            result.outliers.push_back(std::move(o01));
          }
        }
      }
    }
    if(!surface && state.navigation.currentVolume)
    {
		const double resX = 200_um;
		const double resY = 200_um;
		const double resZ = 200_um;
		SymMatrix3D cov;
		cov << resX * resX, 0, 0,
			   0, resY * resY, 0,
			   0, 0, resZ * resZ;
		
		const double dx = resX * gauss(generator);
		const double dy = resY * gauss(generator);
		const double dz = resZ * gauss(generator);
		
		const Vector3D pos = stepper.position(state.stepping);
		
		Measurement<SourceLink, FreeParametersIndices, eFreePos0, eFreePos1, eFreePos2> meas(
				std::shared_ptr<const Volume>(state.navigation.currentVolume), {}, cov, pos.x() + dx, pos.y() + dy, pos.z() + dz);
		result.freeMeasurements.push_back(std::move(meas));
	}
  }
};

double dX, dY;
Vector3D pos;
const Surface* sur;

struct MeasurementCalibrator {
  template <typename parameters_t>
  SourceLink::MeasurementType operator()(
      const SourceLink& sourceLink,
      const parameters_t& /* parameters */) const {
    return *sourceLink;
  }
};

///
/// @brief Simplified material interaction effect by pure gaussian
/// deflection
///
struct MaterialScattering {
  /// @brief Constructor
  MaterialScattering() = default;

  /// @brief Main action list call operator for the scattering on material
  ///
  /// @todo deal momentum in a gaussian way properly
  ///
  /// @tparam propagator_state_t State of the propagator
  /// @param stepper_t Type of the stepper
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper) const {
    // Check if there is a surface with material and a covariance is existing
    if (state.navigation.currentSurface &&
        state.navigation.currentSurface->surfaceMaterial() &&
        state.stepping.cov != typename stepper_t::Covariance::Zero()) {
      // Sample angles
      std::normal_distribution<double> scatterAngle(
          0., 0.017);  //< \approx 1 degree
      double dPhi = scatterAngle(generator), dTheta = scatterAngle(generator);

      // Update the covariance
      state.stepping.cov(ePHI, ePHI) += dPhi * dPhi;
      state.stepping.cov(eTHETA, eTHETA) += dTheta * dTheta;

      // Update the angles
      auto direction = stepper.direction(state.stepping);
      double theta = std::acos(direction.z());
      double phi = std::atan2(direction.y(), direction.x());

      state.stepping.update(
          stepper.position(state.stepping),
          {std::sin(theta + dTheta) * std::cos(phi + dPhi),
           std::sin(theta + dTheta) * std::sin(phi + dPhi),
           std::cos(theta + dTheta)},
          std::max(stepper.momentum(state.stepping) -
                       std::abs(gauss(generator)) * UnitConstants::MeV,
                   0.));
    }
  }
};

struct MinimalOutlierFinder {
  /// The measurement significance criteria
  double measurementSignificanceCutoff = 0;
  /// The chi2 round-off error
  double chi2Tolerance = 10e-5;

  /// @brief Public call mimicking an outlier finder
  ///
  /// @tparam track_state_t Type of the track state
  ///
  /// @param state The track state to investigate
  ///
  /// @return Whether it's outlier or not
  template <typename track_state_t>
  bool operator()(const track_state_t& state) const {
    // Can't determine if it's an outlier if no calibrated measurement or no
    // predicted parameters
    if (not state.hasCalibrated() or not state.hasBoundPredicted()) { // TODO: free case
      return false;
    }

    // The predicted parameters coefficients
    const auto& predicted = state.boundPredicted();
    // The predicted parameters covariance
    const auto& predicted_covariance = state.boundPredictedCovariance();

    // Calculate the chi2 using predicted parameters and calibrated measurement
    double chi2 = std::numeric_limits<double>::max();
    visit_measurement(
        state.calibrated(), state.calibratedCovariance(),
        state.calibratedSize(),
        [&](const auto calibrated, const auto calibrated_covariance) {
          constexpr size_t measdim = decltype(calibrated)::RowsAtCompileTime;
          using par_t = ActsVectorD<measdim>;

          // Take the projector (measurement mapping function)
          const ActsMatrixD<measdim, eBoundParametersSize> H =
              state.projector()
                  .template topLeftCorner<measdim, eBoundParametersSize>();

          // Calculate the residual
          const par_t residual = calibrated - H * predicted;

          // Calculate the chi2
          chi2 = (residual.transpose() *
                  ((calibrated_covariance +
                    H * predicted_covariance * H.transpose()))
                      .inverse() *
                  residual)
                     .eval()(0, 0);
        });

    // In case the chi2 is too small
    if (std::abs(chi2) < chi2Tolerance) {
      return false;
    }
    // The chisq distribution
    boost::math::chi_squared chiDist(state.calibratedSize());
    // The p-Value
    double pValue = 1 - boost::math::cdf(chiDist, chi2);
    // If pValue is NOT significant enough => outlier
    return pValue > measurementSignificanceCutoff ? false : true;
  }
};

///
/// @brief Unit test for Kalman fitter with measurements along the x-axis
///
BOOST_AUTO_TEST_CASE(kalman_fitter_zero_field) {
  // Build detector
  CubicTrackingGeometry cGeometry(tgContext);
  auto detector = cGeometry();

  // Build navigator for the measurement creatoin
  Navigator mNavigator(detector);
  mNavigator.resolvePassive = false;
  mNavigator.resolveMaterial = true;
  mNavigator.resolveSensitive = true;

  // Use straingt line stepper to create the measurements
  StraightLineStepper mStepper;

  // Define the measurement propagator
  using MeasurementPropagator = Propagator<StraightLineStepper, Navigator>;

  // Build propagator for the measurement creation
  MeasurementPropagator mPropagator(mStepper, mNavigator);
  Vector3D mPos(-3_m, 0., 0.), mMom(1_GeV, 0., 0);
  SingleCurvilinearTrackParameters<NeutralPolicy> mStart(std::nullopt, mPos,
                                                         mMom, 42_ns);

  // Create action list for the measurement creation
  using MeasurementActions = ActionList<MeasurementCreator, DebugOutput>;
  using MeasurementAborters = AbortList<EndOfWorldReached>;

  auto pixelResX = Resolution(eLOC_0, 25_um);
  auto pixelResY = Resolution(eLOC_1, 50_um);
  auto stripResX = Resolution(eLOC_0, 100_um);
  auto stripResY = Resolution(eLOC_1, 150_um);

  ElementResolution pixelElementRes = {pixelResX, pixelResY};
  ElementResolution stripElementResI = {stripResX};
  ElementResolution stripElementResO = {stripResY};

  VolumeResolution pixelVolumeRes;
  pixelVolumeRes[2] = pixelElementRes;
  pixelVolumeRes[4] = pixelElementRes;

  VolumeResolution stripVolumeRes;
  stripVolumeRes[2] = stripElementResI;
  stripVolumeRes[4] = stripElementResO;
  stripVolumeRes[6] = stripElementResI;
  stripVolumeRes[8] = stripElementResO;

  DetectorResolution detRes;
  detRes[2] = pixelVolumeRes;
  detRes[3] = stripVolumeRes;

  // Set options for propagator
  PropagatorOptions<MeasurementActions, MeasurementAborters> mOptions(
      tgContext, mfContext);
  mOptions.debug = debugMode;
  mOptions.maxStepSize = 10_cm;
  auto& mCreator = mOptions.actionList.get<MeasurementCreator>();
  mCreator.detectorResolution = detRes;

  // Launch and collect - the measurements
  auto mResult = mPropagator.propagate(mStart, mOptions).value();
  if (debugMode) {
    const auto debugString =
        mResult.template get<DebugOutput::result_type>().debugString;
    std::cout << ">>>> Measurement creation: " << std::endl;
    std::cout << debugString;
  }

  // Extract measurements from result of propagation.
  // This vector owns the measurements
  std::vector<SourceLink::MeasurementType> measurements = std::move(
      mResult.template get<MeasurementCreator::result_type>().measurements);
  BOOST_CHECK_EQUAL(measurements.size(), 6u);

  // Make a vector of source links as input to the KF
  std::vector<SourceLink> sourcelinks;
  std::transform(measurements.begin(), measurements.end(),
                 std::back_inserter(sourcelinks),
                 [](const auto& m) { return SourceLink{&m}; });

std::vector<SourceLink::MeasurementType> freeMeasurements = std::move(
      mResult.template get<MeasurementCreator::result_type>().freeMeasurements);
      
std::vector<SourceLink> freeSourcelinks;
std::transform(freeMeasurements.begin(), freeMeasurements.end(),
                 std::back_inserter(freeSourcelinks),
                 [](const auto& m) { return SourceLink{&m}; }); 
                      
  // The KalmanFitter - we use the eigen stepper for covariance transport
  // Build navigator for the measurement creatoin
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3D(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  // Set initial parameters for the particle track
  Covariance cov;
  cov << 1000_um, 0., 0., 0., 0., 0., 0., 1000_um, 0., 0., 0., 0., 0., 0., 0.05,
      0., 0., 0., 0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0., 0.01, 0., 0., 0.,
      0., 0., 0., 1.;

  Vector3D rPos(-3_m, 10_um * gauss(generator), 100_um * gauss(generator));
  Vector3D rMom(1_GeV, 0.025_GeV * gauss(generator),
                0.025_GeV * gauss(generator));

  SingleCurvilinearTrackParameters<ChargedPolicy> rStart(cov, rPos, rMom, 1.,
                                                         42.);

  const Surface* rSurface = &rStart.referenceSurface();

  using Updater = GainMatrixUpdater;
  using Smoother = GainMatrixSmoother;
  using KalmanFitter =
      KalmanFitter<RecoPropagator>;

  MinimalOutlierFinder outlierFinder;
  outlierFinder.measurementSignificanceCutoff = 0.05;

  KalmanFitter kFitter(rPropagator,
                       getDefaultLogger("KalmanFilter", Logging::FATAL));

  KalmanFitterOptions<MinimalOutlierFinder> kfOptions(
      tgContext, mfContext, calContext, outlierFinder, rSurface);

  // Fit the track
  auto fitRes = kFitter.fit<Updater, Smoother, MeasurementCalibrator>(sourcelinks, rStart, kfOptions, freeSourcelinks);
  BOOST_CHECK(fitRes.ok());
  auto& fittedTrack = *fitRes;
  
  auto state = fittedTrack.fittedStates.getTrackState(fittedTrack.trackTip);
  std::cout << state.hasBoundFiltered() << " " << state.hasFreeFiltered() << std::endl;
  while(state.hasFreeFiltered())
	state = fittedTrack.fittedStates.getTrackState(state.previous());
if(state.hasBoundFiltered())
        std::cout << state.boundFiltered() << "\n" << state.boundFilteredCovariance() << std::endl;
  //~ auto fittedParameters = fittedTrack.fittedParameters.value();
std::cout << "start destruction" << std::endl;
fittedTrack.currentFreeMeasurements.clear();
std::cout << "clear " << &fittedTrack.currentFreeMeasurements << " | " << &fittedTrack.currentVolume << " | " << fittedTrack.currentVolume << std::endl;
fittedTrack.~KalmanFitterResult<SourceLink>(); // TODO: this line breaks
std::cout << "fitted track destroyed" << std::endl;
  //~ // Calculate global track parameters covariance matrix
  //~ const auto& [trackParamsCov, stateRowIndices] =
      //~ detail::globalTrackParametersCovariance(fittedTrack.fittedStates,
                                              //~ fittedTrack.trackTip);

  //~ // Check the size of the global track parameters size
  //~ BOOST_CHECK_EQUAL(stateRowIndices.size(), 6);
  //~ BOOST_CHECK_EQUAL(stateRowIndices.at(fittedTrack.trackTip), 30);
  //~ BOOST_CHECK_EQUAL(trackParamsCov.rows(), 6 * eBoundParametersSize);

  //~ // Make sure it is deterministic
  //~ fitRes = kFitter.fit<Updater, Smoother, MeasurementCalibrator>(sourcelinks, rStart, kfOptions);
  //~ BOOST_CHECK(fitRes.ok());
  //~ auto& fittedAgainTrack = *fitRes;
  //~ auto fittedAgainParameters = fittedAgainTrack.fittedParameters.value();

  //~ CHECK_CLOSE_REL(fittedParameters.parameters().template head<5>(),
                  //~ fittedAgainParameters.parameters().template head<5>(), 1e-5);
  //~ CHECK_CLOSE_ABS(fittedParameters.parameters().template tail<1>(),
                  //~ fittedAgainParameters.parameters().template tail<1>(), 1e-5);

  //~ // Change the order of the sourcelinks
  //~ std::vector<SourceLink> shuffledMeasurements = {
      //~ sourcelinks[3], sourcelinks[2], sourcelinks[1],
      //~ sourcelinks[4], sourcelinks[5], sourcelinks[0]};

  //~ // Make sure it works for shuffled measurements as well
  //~ fitRes = kFitter.fit<Updater, Smoother, MeasurementCalibrator>(shuffledMeasurements, rStart, kfOptions);
  //~ BOOST_CHECK(fitRes.ok());
  //~ auto& fittedShuffledTrack = *fitRes;
  //~ auto fittedShuffledParameters = fittedShuffledTrack.fittedParameters.value();

  //~ CHECK_CLOSE_REL(fittedParameters.parameters().template head<5>(),
                  //~ fittedShuffledParameters.parameters().template head<5>(),
                  //~ 1e-5);
  //~ CHECK_CLOSE_ABS(fittedParameters.parameters().template tail<1>(),
                  //~ fittedShuffledParameters.parameters().template tail<1>(),
                  //~ 1e-5);

  //~ // Remove one measurement and find a hole
  //~ std::vector<SourceLink> measurementsWithHole = {
      //~ sourcelinks[0], sourcelinks[1], sourcelinks[2], sourcelinks[4],
      //~ sourcelinks[5]};

  //~ // Make sure it works for shuffled measurements as well
  //~ fitRes = kFitter.fit<Updater, Smoother, MeasurementCalibrator>(measurementsWithHole, rStart, kfOptions);
  //~ BOOST_CHECK(fitRes.ok());
  //~ auto& fittedWithHoleTrack = *fitRes;
  //~ auto fittedWithHoleParameters = fittedWithHoleTrack.fittedParameters.value();

  //~ // Calculate global track parameters covariance matrix
  //~ const auto& [holeTrackTrackParamsCov, holeTrackStateRowIndices] =
      //~ detail::globalTrackParametersCovariance(fittedWithHoleTrack.fittedStates,
                                              //~ fittedWithHoleTrack.trackTip);

  //~ // Check the size of the global track parameters size
  //~ BOOST_CHECK_EQUAL(holeTrackStateRowIndices.size(), 6);
  //~ BOOST_CHECK_EQUAL(holeTrackStateRowIndices.at(fittedWithHoleTrack.trackTip),
                    //~ 30);
  //~ BOOST_CHECK_EQUAL(holeTrackTrackParamsCov.rows(), 6 * eBoundParametersSize);

  //~ // Count one hole
  //~ BOOST_CHECK_EQUAL(fittedWithHoleTrack.missedActiveSurfaces.size(), 1u);
  //~ // And the parameters should be different
  //~ //
  //~ // BOOST_CHECK(!Acts::Test::checkCloseRel(fittedParameters.parameters().template
  //~ // head<5>(), ~ fittedWithHoleParameters.parameters().template head<5>(), ~
  //~ // 1e-6));
  //~ std::cout << &(fittedParameters) << std::endl;
  //~ std::cout << &(fittedWithHoleParameters) << std::endl;
  //~ BOOST_CHECK(!Acts::Test::checkCloseRel(fittedParameters.parameters(),
                                         //~ fittedWithHoleParameters.parameters(),
                                         //~ 1e-6));

  //~ // Run KF fit in backward filtering mode
  //~ kfOptions.backwardFiltering = true;
  //~ // Fit the track
  //~ fitRes = kFitter.fit<Updater, Smoother, MeasurementCalibrator>(sourcelinks, rStart, kfOptions);
  //~ BOOST_CHECK(fitRes.ok());
  //~ auto fittedWithBwdFiltering = *fitRes;
  //~ // Check the filtering and smoothing status flag
  //~ BOOST_CHECK(fittedWithBwdFiltering.forwardFiltered);
  //~ BOOST_CHECK(not fittedWithBwdFiltering.smoothed);

  //~ // Count the number of 'smoothed' states
  //~ auto trackTip = fittedWithBwdFiltering.trackTip;
  //~ auto mj = fittedWithBwdFiltering.fittedStates;
  //~ size_t nSmoothed = 0;
  //~ mj.visitBackwards(trackTip, [&](const auto& state) {
    //~ if (state.hasBoundSmoothed())
      //~ nSmoothed++;
  //~ });
  //~ BOOST_CHECK_EQUAL(nSmoothed, 6u);

  //~ // Extract outliers from result of propagation.
  //~ // This vector owns the outliers
  //~ std::vector<SourceLink::MeasurementType> outliers = std::move(
      //~ mResult.template get<MeasurementCreator::result_type>().outliers);

  //~ // Replace one measurement with outlier
  //~ std::vector<SourceLink> measurementsWithOneOutlier = {
      //~ sourcelinks[0],           sourcelinks[1], sourcelinks[2],
      //~ SourceLink{&outliers[3]}, sourcelinks[4], sourcelinks[5]};

  //~ // Make sure it works with one outlier
  //~ fitRes = kFitter.fit<Updater, Smoother, MeasurementCalibrator>(measurementsWithOneOutlier, rStart, kfOptions);
  //~ BOOST_CHECK(fitRes.ok());
  //~ auto& fittedWithOneOutlier = *fitRes;

  //~ // Count the number of outliers
  //~ trackTip = fittedWithOneOutlier.trackTip;
  //~ mj = fittedWithOneOutlier.fittedStates;
  //~ size_t nOutliers = 0;
  //~ mj.visitBackwards(trackTip, [&](const auto& state) {
    //~ auto typeFlags = state.typeFlags();
    //~ if (typeFlags.test(TrackStateFlag::OutlierFlag)) {
      //~ nOutliers++;
    //~ }
  //~ });
  //~ BOOST_CHECK_EQUAL(nOutliers, 1u);
std::cout << "fit done" << std::endl;
kFitter.~KalmanFitter();
std::cout << "deconstructed" << std::endl;
}

}  // namespace Test
}  // namespace Acts
