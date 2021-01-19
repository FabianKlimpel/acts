// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/StandardInteractions.hpp"

ActsFatras::StandardChargedElectroMagneticInteractions
ActsFatras::makeStandardChargedElectroMagneticInteractions(double pMin) {
  StandardChargedElectroMagneticInteractions pl;
  pl.get<detail::StandardBetheBloch>().selectOutputParticle.valMin = pMin;
  pl.get<detail::StandardBetheHeitler>().selectOutputParticle.valMin = pMin;
  pl.get<detail::StandardBetheHeitler>().selectChildParticle.valMin = pMin;
  return pl;
}

ActsFatras::NeutralPhysicsList
ActsFatras::makeNeutralPhysicsList(double minimumAbsMomentum) {
  NeutralPhysicsList pl;
  pl.get<detail::ParametrizedNuclearInteraction>().selectInput.valMin = minimumAbsMomentum;
  return pl;
}
