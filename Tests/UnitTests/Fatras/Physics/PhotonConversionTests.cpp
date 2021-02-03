// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Common.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/PhotonConversion/PhotonConversion.hpp"

#include <limits>
#include <random>

#include "Dataset.hpp"

using Generator = std::ranlux48;

BOOST_AUTO_TEST_SUITE(FatrasPhotonConversion)

BOOST_DATA_TEST_CASE(PhotonConversion, Dataset::parametersPhotonConversion, phi,
                     lambda, seed) {
  using Scalar = ActsFatras::PhotonConversion::Scalar;
  using namespace Acts::UnitLiterals;

  Generator gen(seed);

  /// Produce not a photon
  ActsFatras::Particle particle =
      Dataset::makeParticle(Acts::PdgParticle::eElectron, phi, lambda, 1_GeV);
  ActsFatras::Particle particleInit = particle;

  ActsFatras::PhotonConversion pc;

  // No limits should be set
  std::pair<Scalar, Scalar> limits;
  limits = pc.generatePathLimits(gen, particle);
  BOOST_CHECK_EQUAL(limits.first, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(limits.second, std::numeric_limits<Scalar>::infinity());

  // No particles should be generated
  std::vector<ActsFatras::Particle> generated;
  BOOST_CHECK(!pc.run(gen, particle, generated));
  BOOST_CHECK(generated.empty());
  // Particle shouldn't be modified
  BOOST_CHECK_EQUAL(particle.fourPosition(), particleInit.fourPosition());
  BOOST_CHECK_EQUAL(particle.fourMomentum(), particleInit.fourMomentum());
  BOOST_CHECK_EQUAL(particle.process(), particleInit.process());
  BOOST_CHECK_EQUAL(particle.properTime(), particleInit.properTime());
  BOOST_CHECK_EQUAL(particle.pathInX0(), particleInit.pathInX0());
  BOOST_CHECK_EQUAL(particle.pathInL0(), particleInit.pathInL0());

  /// Produce a dead photon
  particle = Dataset::makeParticle(Acts::PdgParticle::eGamma, phi, lambda, 0);
  particleInit = particle;

  // No limits should be set - momentum too low
  limits = pc.generatePathLimits(gen, particle);
  BOOST_CHECK_EQUAL(limits.first, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(limits.second, std::numeric_limits<Scalar>::infinity());

  // No particles should be generated - momentum too low
  generated.clear();
  BOOST_CHECK(!pc.run(gen, particle, generated));
  BOOST_CHECK(generated.empty());
  // Particle shouldn't be modified
  BOOST_CHECK_EQUAL(particle.fourPosition(), particleInit.fourPosition());
  BOOST_CHECK_EQUAL(particle.fourMomentum(), particleInit.fourMomentum());
  BOOST_CHECK_EQUAL(particle.process(), particleInit.process());
  BOOST_CHECK_EQUAL(particle.properTime(), particleInit.properTime());
  BOOST_CHECK_EQUAL(particle.pathInX0(), particleInit.pathInX0());
  BOOST_CHECK_EQUAL(particle.pathInL0(), particleInit.pathInL0());

  /// Produce a low momentum photon
  particle =
      Dataset::makeParticle(Acts::PdgParticle::eGamma, phi, lambda, 1_keV);
  particleInit = particle;

  // No limits should be set - momentum too low
  limits = pc.generatePathLimits(gen, particle);
  BOOST_CHECK_EQUAL(limits.first, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(limits.second, std::numeric_limits<Scalar>::infinity());

  // No particles should be generated - momentum too low
  generated.clear();
  BOOST_CHECK(!pc.run(gen, particle, generated));
  BOOST_CHECK(generated.empty());
  // Particle shouldn't be modified
  BOOST_CHECK_EQUAL(particle.fourPosition(), particleInit.fourPosition());
  BOOST_CHECK_EQUAL(particle.fourMomentum(), particleInit.fourMomentum());
  BOOST_CHECK_EQUAL(particle.process(), particleInit.process());
  BOOST_CHECK_EQUAL(particle.properTime(), particleInit.properTime());
  BOOST_CHECK_EQUAL(particle.pathInX0(), particleInit.pathInX0());
  BOOST_CHECK_EQUAL(particle.pathInL0(), particleInit.pathInL0());

  /// Produce a high momentum photon
  particle =
      Dataset::makeParticle(Acts::PdgParticle::eGamma, phi, lambda, 1_GeV);
  particleInit = particle;

  // No limits should be set - momentum too low
  limits = pc.generatePathLimits(gen, particle);
  BOOST_CHECK_NE(limits.first, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(limits.second, std::numeric_limits<Scalar>::infinity());

  // No particles should be generated - momentum too low
  generated.clear();
  BOOST_CHECK(pc.run(gen, particle, generated));
  BOOST_CHECK_EQUAL(generated.size(), 2);

  // Test the children
  BOOST_CHECK((generated[0].pdg() == Acts::PdgParticle::eElectron) ||
              (generated[0].pdg() == Acts::PdgParticle::ePositron));
  BOOST_CHECK((generated[1].pdg() == Acts::PdgParticle::eElectron) ||
              (generated[1].pdg() == Acts::PdgParticle::ePositron));
  BOOST_CHECK_NE(generated[0].pdg(), generated[1].pdg());
  BOOST_CHECK_NE(generated[0].fourMomentum(), Acts::Vector4::Zero());
  BOOST_CHECK_NE(generated[1].fourMomentum(), Acts::Vector4::Zero());

  // Test for similar invariant masses
  auto momSum = generated[0].fourMomentum() + generated[1].fourMomentum();
  auto momVector = momSum.template segment<3>(Acts::eMom0);
  auto sSum = sqrt(momSum[Acts::eEnergy] * momSum[Acts::eEnergy] -
                   momVector.norm() * momVector.norm());

  auto sParticle =
      sqrt(particleInit.energy() * particleInit.energy() -
           particleInit.absoluteMomentum() * particleInit.absoluteMomentum());

  CHECK_CLOSE_OR_SMALL(sSum, sParticle, 1e-2, 1e-2);
}

BOOST_AUTO_TEST_SUITE_END()
