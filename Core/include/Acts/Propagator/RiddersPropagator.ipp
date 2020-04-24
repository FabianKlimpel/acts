// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename propagator_t>
template <typename return_parameters_t, typename parameters_t,
          typename propagator_options_t>
auto Acts::RiddersPropagator<propagator_t>::propagate(
    const parameters_t& start, const propagator_options_t& options) const
    -> Result<action_list_t_result_t<
        return_parameters_t, typename propagator_options_t::action_list_type>> {
  using ThisResult = Result<action_list_t_result_t<
      CurvilinearParameters, typename propagator_options_t::action_list_type>>;

  // Propagate the nominal parameters
  auto nominalRet = m_propagator.propagate(start, options);
  if (not nominalRet.ok()) {
    return ThisResult::failure(nominalRet.error());
  }

  // Extract results from the nominal propagation
  auto nominalResult = std::move(nominalRet).value();
  const auto& nominalParameters& nominalParameters =
      nominalResult.endParameters->parameters();
  // Use the curvilinear surface of the propagated parameters as target
  const Surface& surface = nominalResult.endParameters->referenceSurface();

  // Steps for estimating derivatives
  std::vector<double> deviations = {-4e-3, -2e-3, 2e-3, 4e-3};

  // Allow larger distances for the oscillation
  propagator_options_t opts = options;

  // Test whether we want to end up bound
  if constexpr (return_parameters_t::is_local_representation) {
    // Allow larger distances for the oscillation
    propagator_options_t opts = options;
    opts.pathLimit *= 2.;

    // Pick the surface of the propagation as target
    const Surface& surface = nominalResult.endParameters->referenceSurface();

    // Case I: We start and end bound
    if constexpr (parameters_t::is_local_representation) {
      // Derivations of each parameter around the nominal parameters
      std::array<std::vector<BoundVector>, eBoundParametersSize> derivatives;

      // Wiggle each dimension individually
      for (unsigned int i = 0; i < eBoundParametersSize; i++) {
        derivatives[i] = wiggleDimension(opts, start, i, nominalParameters,
                                         deviations, surface);
      }

      // Exchange the result by Ridders Covariance
      const FullParameterSet& parSet =
          nominalResult.endParameters->getParameterSet();
      FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
      if (start.covariance()) {
        mParSet->setCovariance(std::get<BoundSymMatrix>(calculateCovariance(
            derivatives, *start.covariance(), deviations, Vector3D())));
      }
    }
    // Case II: We start free and end bound
    else {
      // Derivations of each parameter around the nominal parameters
      std::array<std::vector<BoundVector>, 7> derivatives;

      // Wiggle each dimension individually
      for (unsigned int i = 0; i < 7; i++) {
        derivatives[i] = wiggleDimension(opts, start, i, nominalParameters,
                                         deviations, surface);
      }

      // Exchange the result by Ridders Covariance
      const FullParameterSet& parSet =
          nominalResult.endParameters->getParameterSet();
      FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
      if (start.covariance()) {
        mParSet->setCovariance(std::get<BoundSymMatrix>(
            calculateCovariance(derivatives, *start.covariance(), deviations,
                                start.parameters().template segment<3>(4))));
      }
    }
  }
  // TODO: Case III and IV need some propagation distance settings, otherwise
  // the result should be uncomparable
  else {
    // Case III: We start bound and end free
    if constexpr (parameters_t::is_local_representation) {
      // Derivations of each parameter around the nominal parameters
      std::array<std::vector<FreeVector>, eBoundParametersSize> derivatives;

      // Wiggle each dimension individually
      for (unsigned int i = 0; i < eBoundParametersSize; i++) {
        derivatives[i] =
            wiggleDimension(options, start, i, nominalParameters, deviations);
      }

      // Exchange the result by Ridders Covariance
      FreeParameters* parameters =
          const_cast<FreeParameters*>(nominalResult.endParameters.get());
      if (start.covariance()) {
        parameters->covariance(std::get<FreeSymMatrix>(calculateCovariance(
            derivatives, *start.covariance(), deviations, Vector3D())));
      }
      nominalResult.endParameters =
          std::make_unique<const FreeParameters>(*parameters);
    }
    // Case IV: We start and end free
    else {
      // Derivations of each parameter around the nominal parameters
      std::array<std::vector<FreeVector>, 7> derivatives;

      // Wiggle each dimension individually
      for (unsigned int i = 0; i < 7; i++) {
        derivatives[i] =
            wiggleDimension(options, start, i, nominalParameters, deviations);
      }

      // Exchange the result by Ridders Covariance
      FreeParameters* parameters =
          const_cast<FreeParameters*>(nominalResult.endParameters.get());
      if (start.covariance()) {
        parameters->covariance(std::get<FreeSymMatrix>(
            calculateCovariance(derivatives, *start.covariance(), deviations,
                                start.parameters().template segment<3>(4))));
      }
      nominalResult.endParameters =
          std::make_unique<const FreeParameters>(*parameters);
    }
  }

  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
template <typename parameters_t, typename propagator_options_t>
auto Acts::RiddersPropagator<propagator_t>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const
    -> Result<action_list_t_result_t<
        BoundParameters, typename propagator_options_t::action_list_type>> {
  using ThisResult = Result<action_list_t_result_t<
      BoundParameters, typename propagator_options_t::action_list_type>>;

  // Propagate the nominal parameters
  auto nominalRet = m_propagator.propagate(start, target, options);
  if (not nominalRet.ok()) {
    return ThisResult::failure(nominalRet.error());
  }

  // Extract results from the nominal propagation
  auto nominalResult = std::move(nominalRet).value();
  const BoundVector& nominalParameters =
      nominalResult.endParameters->parameters();

  // Steps for estimating derivatives
  std::vector<double> deviations = {-4e-4, -2e-4, 2e-4, 4e-4};
  if (target.type() == Surface::Disc) {
    deviations = {{-3e-4, -1e-4, 1e-4, 3e-4}};
  }

  // - for planar surfaces the dest surface is a perfect destination
  // surface for the numerical propagation, as reference frame
  // aligns with the referenceSurface.transform().rotation() at
  // at any given time
  //
  // - for straw & cylinder, where the error is given
  // in the reference frame that re-aligns with a slightly different
  // intersection solution

  // Allow larger distances for the oscillation
  propagator_options_t opts = options;
  
  	  // Case I: We start bound
	  if constexpr (parameters_t::is_local_representation) {
		  opts.pathLimit *= 2.;

		  // Derivations of each parameter around the nominal parameters
		  std::array<std::vector<BoundVector>, eBoundParametersSize> derivatives;

		  // Wiggle each dimension individually
		  for (unsigned int i = 0; i < eBoundParametersSize; i++) {
			derivatives[i] =
				wiggleDimension(opts, start, i, nominalParameters, deviations, target);
		  }
		  // Exchange the result by Ridders Covariance
		  const FullParameterSet& parSet =
			  nominalResult.endParameters->getParameterSet();
		  FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
		  if (start.covariance()) {
			// Test if target is disc - this may lead to inconsistent results
			if (target.type() == Surface::Disc) {
			  for (const std::vector<BoundVector>& deriv : derivatives) {
				if (inconsistentDerivativesOnDisc(deriv)) {
				  // Set covariance to zero and return
				  // TODO: This should be changed to indicate that something went
				  // wrong
				  mParSet->setCovariance(BoundSymMatrix::Zero());
				  return std::move(nominalResult);
				}
			  }
			}
			mParSet->setCovariance(std::get<BoundSymMatrix>(
				calculateCovariance(derivatives, *start.covariance(), deviations, Vector3D())));
		  }
		}  
	// Case II: We start free
	  else
	  {
			// Derivations of each parameter around the nominal parameters
			std::array<std::vector<BoundVector>, 7>
				derivatives;
				
			// Wiggle each dimension individually
			for (unsigned int i = 0; i < 7; i++) {
			  derivatives[i] = wiggleDimension(
				  opts, start, i, nominalParameters,
				  deviations, target);
			}

			// Exchange the result by Ridders Covariance
			const FullParameterSet& parSet =
				nominalResult.endParameters->getParameterSet();
			FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
			if (start.covariance()) {			  
			// Test if target is disc - this may lead to inconsistent results
			if (target.type() == Surface::Disc) {
			  for (const std::vector<BoundVector>& deriv : derivatives) {
				if (inconsistentDerivativesOnDisc(deriv)) {
				  // Set covariance to zero and return
				  // TODO: This should be changed to indicate that something went
				  // wrong
				  mParSet->setCovariance(BoundSymMatrix::Zero());
				  return std::move(nominalResult);
				}
			  }
			}
			  mParSet->setCovariance(std::get<BoundSymMatrix>(
				  calculateCovariance(derivatives, *start.covariance(), deviations, start.parameters().template segment<3>(4))));
			} 
	  }
  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
bool Acts::RiddersPropagator<propagator_t>::inconsistentDerivativesOnDisc(
    const std::vector<Acts::BoundVector>& derivatives) const {
  // Test each component with each other
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    bool jumpedAngle = true;
    for (unsigned int j = 0; j < derivatives.size(); j++) {
      // If there is at least one with a similar angle then it seems to work
      // properly
      if (i != j &&
          std::abs(derivatives[i](1) - derivatives[j](1)) < 0.5 * M_PI) {
        jumpedAngle = false;
        break;
      }
    }
    // Break if a jump was detected
    if (jumpedAngle) {
      return true;
    }
  }
  return false;
}

template <typename propagator_t>
template <typename options_t, typename start_parameters_t>
std::vector<Acts::BoundVector>
Acts::RiddersPropagator<propagator_t>::wiggleDimension(
    const options_t& options, const start_parameters_t& startPars,
    const unsigned int param, const BoundVector& nominal,
    const std::vector<double>& deviations, const Surface& target) const {
  // Storage of the results
  std::vector<BoundVector> derivatives;
  derivatives.reserve(deviations.size());
  for (double h : deviations) {
    start_parameters_t tp =
        wiggleStartVector(options.geoContext, h, param, startPars);

    const auto& r = m_propagator.propagate(tp, target, options).value();
    // Collect the slope
    derivatives.push_back((r.endParameters->parameters() - nominal) / h);
std::cout << r.endParameters->parameters().transpose() << std::endl;
    // Correct angular results
    if constexpr (start_parameters_t::is_local_representation) {
      // Correct for a possible variation of phi around
      if (param == 2) {
        double phi0 = nominal(Acts::ePHI);
        double phi1 = r.endParameters->parameters()(Acts::ePHI);
        if (std::abs(phi1 + 2. * M_PI - phi0) < std::abs(phi1 - phi0))
          derivatives.back()[Acts::ePHI] = (phi1 + 2. * M_PI - phi0) / h;
        else if (std::abs(phi1 - 2. * M_PI - phi0) < std::abs(phi1 - phi0))
          derivatives.back()[Acts::ePHI] = (phi1 - 2. * M_PI - phi0) / h;
      }
    } else {
      if (param == 4 || param == 5) {
        double phi0 = nominal(Acts::ePHI);
        double phi1 = r.endParameters->parameters()(Acts::ePHI);
        if (std::abs(phi1 + 2. * M_PI - phi0) < std::abs(phi1 - phi0))
          derivatives.back()[Acts::ePHI] = (phi1 + 2. * M_PI - phi0) / h;
        else if (std::abs(phi1 - 2. * M_PI - phi0) < std::abs(phi1 - phi0))
          derivatives.back()[Acts::ePHI] = (phi1 - 2. * M_PI - phi0) / h;

        double theta0 = nominal(Acts::eTHETA);
        double theta1 = r.endParameters->parameters()(Acts::eTHETA);
        if (std::abs(theta1 + M_PI - theta0) < std::abs(theta1 - theta0))
          derivatives.back()[Acts::eTHETA] = (theta1 + M_PI - theta0) / h;
        else if (std::abs(theta1 - M_PI - theta0) < std::abs(theta1 - theta0))
          derivatives.back()[Acts::eTHETA] = (theta1 - M_PI - theta0) / h;
      }
    }
  }

  return derivatives;
}

template <typename propagator_t>
template <typename options_t, typename parameters_t>
std::vector<Acts::FreeVector>
Acts::RiddersPropagator<propagator_t>::wiggleDimension(
    const options_t& options, const parameters_t& startPars,
    const unsigned int param, const Acts::FreeVector& nominal,
    const std::vector<double>& deviations) const {
  // Storage of the results
  std::vector<FreeVector> derivatives;
  derivatives.reserve(deviations.size());
  for (double h : deviations) {
    parameters_t tp =
        wiggleStartVector(options.geoContext, h, param, startPars);

    const auto& r =
        m_propagator.template propagate<FreeParameters>(tp, options).value();

    // Collect the slope
    derivatives.push_back((r.endParameters->parameters() - nominal) / h);
  }

  return derivatives;
}

template <typename propagator_t>
template <typename parameters_t>
parameters_t Acts::RiddersPropagator<propagator_t>::wiggleStartVector(
    std::reference_wrapper<const GeometryContext> geoContext, double h,
    const unsigned int param, parameters_t tp) const {
  if constexpr (parameters_t::is_local_representation) {
    wiggleBoundStartVector(geoContext, h, param, tp);
  } else {
    wiggleFreeStartVector(geoContext, h, param, tp);
  }
  return tp;
}

template <typename propagator_t>
template <typename parameters_t>
void Acts::RiddersPropagator<propagator_t>::wiggleBoundStartVector(
    std::reference_wrapper<const GeometryContext> geoContext, double h,
    const unsigned int param, parameters_t& tp) const {
  // Treatment for theta
  if (param == eTHETA) {
    const double current_theta = tp.template get<eTHETA>();
    if (current_theta + h > M_PI) {
      h = M_PI - current_theta;
    }
    if (current_theta + h < 0) {
      h = -current_theta;
    }
  }

  // Modify start parameter
  switch (param) {
    case 0: {
      tp.template set<eLOC_0>(geoContext, tp.template get<eLOC_0>() + h);
      break;
    }
    case 1: {
      tp.template set<eLOC_1>(geoContext, tp.template get<eLOC_1>() + h);
      break;
    }
    case 2: {
      tp.template set<ePHI>(geoContext, tp.template get<ePHI>() + h);
      break;
    }
    case 3: {
      tp.template set<eTHETA>(geoContext, tp.template get<eTHETA>() + h);
      break;
    }
    case 4: {
      tp.template set<eQOP>(geoContext, tp.template get<eQOP>() + h);
      break;
    }
    case 5: {
      tp.template set<eT>(geoContext, tp.template get<eT>() + h);
      break;
    }
  }
}

template <typename propagator_t>
template <typename parameters_t>
void Acts::RiddersPropagator<propagator_t>::wiggleFreeStartVector(
    std::reference_wrapper<const GeometryContext> geoContext, double h,
    const unsigned int param, parameters_t& tp) const {
  // Modify start parameter
  switch (param) {
    case 0: {
      tp.template set<0>(geoContext, tp.template get<0>() + h);
      break;
    }
    case 1: {
      tp.template set<1>(geoContext, tp.template get<1>() + h);
      break;
    }
    case 2: {
      tp.template set<2>(geoContext, tp.template get<2>() + h);
      break;
    }
    case 3: {
      tp.template set<3>(geoContext, tp.template get<3>() + h);
      break;
    }
    case 4: {
      const double phi = std::atan2(tp.template get<5>(), tp.template get<4>());
      const double theta = std::acos(tp.template get<6>());
      tp.template set<4>(geoContext, std::cos(phi + h) * std::sin(theta));
      tp.template set<5>(geoContext, std::sin(phi + h) * std::sin(theta));
      break;
    }
    case 5: {
      const double phi = std::atan2(tp.template get<5>(), tp.template get<4>());
      const double theta = std::acos(tp.template get<6>());
      tp.template set<4>(geoContext, std::cos(phi) * std::sin(theta + h));
      tp.template set<5>(geoContext, std::sin(phi) * std::sin(theta + h));
      tp.template set<6>(geoContext, std::cos(theta + h));
      break;
    }
    case 6: {
      tp.template set<7>(geoContext, tp.template get<7>() + h);
      break;
    }
  }
}

template <typename propagator_t>
Acts::ActsMatrixD<7, 8>
Acts::RiddersPropagator<propagator_t>::anglesToDirectionsJacobian(
    const Vector3D dir) const {
  ActsMatrixD<7, 8> jacobian = ActsMatrixD<7, 8>::Zero();
  const double x = dir(0);  // == cos(phi) * sin(theta)
  const double y = dir(1);  // == sin(phi) * sin(theta)
  const double z = dir(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  jacobian(0, 0) = 1.;
  jacobian(1, 1) = 1.;
  jacobian(2, 2) = 1.;
  jacobian(3, 3) = 1.;
  jacobian(6, 7) = 1.;

  jacobian(4, 4) = -sinPhi * invSinTheta;
  jacobian(4, 5) = cosPhi * invSinTheta;
  jacobian(5, 4) = cosPhi * cosTheta;
  jacobian(5, 5) = sinPhi * cosTheta;
  jacobian(5, 6) = -invSinTheta * (1. - cosTheta * cosTheta);

  return jacobian;
}

template <typename propagator_t>
auto Acts::RiddersPropagator<propagator_t>::calculateCovariance(
    const std::array<std::vector<BoundVector>, Acts::eBoundParametersSize>&
        derivatives,
    const std::variant<Acts::BoundSymMatrix, Acts::FreeSymMatrix>& startCov,
    const std::vector<double>& deviations, const Vector3D /*unused*/) const
    -> const Covariance {
  BoundMatrix jacobian;
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    jacobian.col(i) = fitLinear(derivatives[i], deviations);
  }
  return BoundSymMatrix(jacobian * std::get<Acts::BoundSymMatrix>(startCov) *
                        jacobian.transpose());
}

template <typename propagator_t>
auto Acts::RiddersPropagator<propagator_t>::calculateCovariance(
    const std::array<std::vector<BoundVector>, 7>& derivatives,
    const std::variant<Acts::BoundSymMatrix, Acts::FreeSymMatrix>& startCov,
    const std::vector<double>& deviations, const Vector3D direction) const
    -> const Covariance {
  ActsMatrixD<eBoundParametersSize, 7> jac;
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    jac.col(i) = fitLinear(derivatives[i], deviations);
  }
  FreeToBoundMatrix jacobian = jac * anglesToDirectionsJacobian(direction);
  return BoundSymMatrix(jacobian * std::get<Acts::FreeSymMatrix>(startCov) *
                        jacobian.transpose());
}

template <typename propagator_t>
auto Acts::RiddersPropagator<propagator_t>::calculateCovariance(
    const std::array<std::vector<FreeVector>, Acts::eBoundParametersSize>&
        derivatives,
    const std::variant<Acts::BoundSymMatrix, Acts::FreeSymMatrix>& startCov,
    const std::vector<double>& deviations, const Vector3D /*unused*/) const
    -> const Covariance {
  BoundToFreeMatrix jacobian;
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    jacobian.col(i) = fitLinear(derivatives[i], deviations);
  }
  return FreeSymMatrix(jacobian * std::get<Acts::BoundSymMatrix>(startCov) *
                       jacobian.transpose());
}

template <typename propagator_t>
auto Acts::RiddersPropagator<propagator_t>::calculateCovariance(
    const std::array<std::vector<FreeVector>, 7>& derivatives,
    const std::variant<Acts::BoundSymMatrix, Acts::FreeSymMatrix>& startCov,
    const std::vector<double>& deviations, const Vector3D direction) const
    -> const Covariance {
  ActsMatrixD<8, 7> jac;
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    jac.col(i) = fitLinear(derivatives[i], deviations);
  }
  FreeMatrix jacobian = jac * anglesToDirectionsJacobian(direction);
  return FreeSymMatrix(jacobian * std::get<Acts::FreeSymMatrix>(startCov) *
                       jacobian.transpose());
}

template <typename propagator_t>
template <typename vector_t>
vector_t Acts::RiddersPropagator<propagator_t>::fitLinear(
    const std::vector<vector_t>& values,
    const std::vector<double>& deviations) const {
  vector_t A;
  vector_t C;
  A.setZero();
  C.setZero();
  double B = 0;
  double D = 0;
  const unsigned int N = deviations.size();

  for (unsigned int i = 0; i < N; ++i) {
    A += deviations.at(i) * values.at(i);
    B += deviations.at(i);
    C += values.at(i);
    D += deviations.at(i) * deviations.at(i);
  }

  vector_t b = (N * A - B * C) / (N * D - B * B);
  vector_t a = (C - B * b) / N;

  return a;
}
