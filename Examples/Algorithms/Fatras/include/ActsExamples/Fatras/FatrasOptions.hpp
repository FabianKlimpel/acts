// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsFatras/Kernel/Process.hpp"
#include "ActsFatras/Physics/StandardPhysicsLists.hpp"

#include <utility>

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Options {

/// Add Fatras options.
///
/// @param desc The options description to add options to
void addFatrasOptions(Description& desc);

/// Read Fatras options to create the algorithm config.
///
/// @tparam simulator_t type of the simulation kernel
/// @param vars         the variables to read from
/// @param simulator    the simulation kernel
template <typename simulator_t>
typename FatrasAlgorithm<simulator_t>::Config readFatrasConfig(
    const Variables& variables, simulator_t&& simulator) {
  using namespace Acts::UnitLiterals;
  using Config = typename FatrasAlgorithm<simulator_t>::Config;
  using PMin = ActsFatras::Min<ActsFatras::Casts::P>;

  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("FatrasOptions", Acts::Logging::INFO))

  Config cfg(std::forward<simulator_t>(simulator));

  // simulation particle cuts both for input and physics lists output
  const auto pmin = variables["fatras-pmin-gev"].as<double>() * 1_GeV;
  cfg.simulator.selectCharged.template get<PMin>().valMin = pmin;
  cfg.simulator.selectNeutral.template get<PMin>().valMin = pmin;
  cfg.simulator.charged.physics =
      ActsFatras::makeChargedElectroMagneticPhysicsList(pmin);
  cfg.simulator.neutral.physics =
      ActsFatras::makeNeutralPhysicsList(pmin);

  // all physics process are enabled by default
  if (not variables["fatras-em-scattering"].as<bool>()) {
    cfg.simulator.charged.physics
        .template disable<ActsFatras::detail::StandardScattering>();
    cfg.simulator.neutral.physics.template disable<ActsFatras::detail::StandardScattering>();
  }
  if (not variables["fatras-em-ionisation"].as<bool>()) {
    cfg.simulator.charged.physics
        .template disable<ActsFatras::detail::StandardBetheBloch>();
  }
  if (not variables["fatras-em-radiation"].as<bool>()) {
    cfg.simulator.charged.physics
        .template disable<ActsFatras::detail::StandardBetheHeitler>();
  }

  // select hit surfaces for charged particles
  const std::string hits = variables["fatras-hits"].as<std::string>();
  if (hits == "sensitive") {
    ACTS_DEBUG("Configure hits on sensitive surfaces");
    cfg.simulator.charged.selectHitSurface.sensitive = true;
    cfg.simulator.charged.selectHitSurface.material = false;
    cfg.simulator.charged.selectHitSurface.passive = false;
    
    cfg.simulator.neutral.selectHitSurface.sensitive = true;
    cfg.simulator.neutral.selectHitSurface.material = false;
    cfg.simulator.neutral.selectHitSurface.passive = false;
  } else if (hits == "material") {
    ACTS_DEBUG("Configure hits on material surfaces");
    cfg.simulator.charged.selectHitSurface.sensitive = false;
    cfg.simulator.charged.selectHitSurface.material = true;
    cfg.simulator.charged.selectHitSurface.passive = false;
    
    cfg.simulator.neutral.selectHitSurface.sensitive = false;
    cfg.simulator.neutral.selectHitSurface.material = true;
    cfg.simulator.neutral.selectHitSurface.passive = false;
  } else if (hits == "all") {
    ACTS_DEBUG("Configure hits on all surfaces");
    cfg.simulator.charged.selectHitSurface.sensitive = false;
    cfg.simulator.charged.selectHitSurface.material = false;
    // this includes sensitive and material surfaces
    cfg.simulator.charged.selectHitSurface.passive = true;
    
    cfg.simulator.neutral.selectHitSurface.sensitive = false;
    cfg.simulator.neutral.selectHitSurface.material = false;
    // this includes sensitive and material surfaces
    cfg.simulator.neutral.selectHitSurface.passive = true;
  } else {
    ACTS_WARNING("Invalid hits configuration '" << hits << "'");
    // none or unknown type -> record nothing
    cfg.simulator.charged.selectHitSurface.sensitive = false;
    cfg.simulator.charged.selectHitSurface.material = false;
    cfg.simulator.charged.selectHitSurface.passive = false;
    
    cfg.simulator.neutral.selectHitSurface.sensitive = false;
    cfg.simulator.neutral.selectHitSurface.material = false;
    cfg.simulator.neutral.selectHitSurface.passive = false;
  }

  return cfg;
}

}  // namespace Options
}  // namespace ActsExamples
