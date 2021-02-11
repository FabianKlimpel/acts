// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialProperties.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <array>
#include <random>

namespace ActsFatras {

/// Simulate electron energy loss using the Bethe-Heitler description.
///
/// Bethe-Heitler for electron bremsstrahlung description as described here:
/// "A Gaussian-mixture approximation of the Bethe–Heitler model of electron
/// energy loss by bremsstrahlung" R. Frühwirth
struct BetheHeitler {
  //~ using Scalar = Particle::Scalar;
  //~ using Vector3 = Particle::Vector3;
  //~ using phi = Acts::detail::VectorHelpers::phi;
  //~ using theta = Acts::detail::VectorHelpers::theta;

  /// A scaling factor to
  double scaleFactor = 1.;

  //~ // Simplified angle evaluation
  //~ bool uniformHertzDipoleAngle = false;
  
//~ // electron in original state, energy loss sampled
//~ Particle bremPhoton(const Particle &particle, Scalar gammaE, Scalar rndPsi, Scalar rndTheta1, Scalar rndTheta2, Scalar rndTheta3) const
//~ {
  //~ // ------------------------------------------------------
  //~ // simple approach
  //~ // (a) simulate theta uniform within the opening angle of the relativistic Hertz dipole
  //~ //      theta_max = 1/gamma
  //~ // (b)Following the Geant4 approximation from L. Urban -> encapsulate that later
  //~ //      the azimutal angle
    
  //~ Scalar psi    =  2. * M_PI * rndPsi;
    
  //~ // the start of the equation
  //~ Scalar theta = 0.;
  //~ Scalar E = particle.energy();

  //~ if (uniformHertzDipoleAngle) {
    //~ // the simplest simulation
    //~ theta = m/E * rndTheta1;  
  //~ } else {
    //~ // -----> 
    //~ theta = m/E;
    //~ // follow 
    //~ constexpr Scalar a = 0.625; // 5/8
    //~ Scalar u =  -log(rndTheta2 * rndTheta3) / a;
    //~ theta *= (rndTheta1 < 0.25 ) ? u : u / 3.; // 9./(9.+27) = 0.25
  //~ }
  
  //~ EX_MSG_VERBOSE("[ brem ]", "BremPhoton", "", "Simulated angle to electron    = " << theta << "." );
  
  //~ Vector3 particleDirection = particle.unitDirection();
  //~ Vector3 photonDirection = particleDirection;
  
  //~ // construct the combined rotation to the scattered direction
  //~ Acts::RotationMatrix3 rotation(
      //~ // rotation of the scattering deflector axis relative to the reference
      //~ Acts::AngleAxis3(psi, particleDirection) *
      //~ // rotation by the scattering angle around the deflector axis
      //~ Acts::AngleAxis3(theta, Acts::makeCurvilinearUnitU(particleDirection)));
  //~ photonDirection.applyOnTheLeft(rotation);
     
     //~ Particle photon(particle.barcode().makeDescendant(0), Acts::PdgParticle::eGamma);
     //~ photon.setProcess(eBremsstrahlung).setPosition4(particle.fourPosition()).setDirection(photonDirection).setAbsMomentum(gammaE);
  //~ return bremPhoton;
//~ } 



  /// Simulate energy loss and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Empty secondaries containers.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::array<Particle, 0> operator()(generator_t &generator,
                                     const Acts::MaterialProperties &slab,
                                     Particle &particle) const {
    // Take a random gamma-distributed value - depending on t/X0
    std::gamma_distribution<double> gDist(slab.thicknessInX0() / std::log(2.0),
                                          1.0);

    const auto u = gDist(generator);
    const auto z = std::exp(-u);
    const auto sampledEnergyLoss =
        std::abs(scaleFactor * particle.energy() * (z - 1.));

//~ std::uniform_distribution<Scalar> uDist(0., 1.);
//~ Particle photon = bremPhoton(particle, sampledEnergyLoss,uDist(generator), uDist(generator), uDist(generator), uDist(generator));
  //~ // recoil / save input momentum for validation
  //~ particle.setDirection(particle.unitDirection() * particle.absoluteMomentum() - gammaE * photonDirection);
  
    // apply the energy loss
    particle.correctEnergy(-sampledEnergyLoss);

    // TODO return the lost energy as a photon
    return {};
  }
};

}  // namespace ActsFatras
