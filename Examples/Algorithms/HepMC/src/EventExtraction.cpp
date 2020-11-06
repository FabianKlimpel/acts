// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/HepMC/EventExtraction.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Particle.hpp"
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

namespace {
/// Stores the initial properties of a particle, the properties before the interaction and the particle properties after the interaction
using Particles = std::tuple<ActsExamples::SimParticle, ActsExamples::SimParticle,
                            std::vector<ActsExamples::SimParticle>>;

/// @brief This method searches for an outgoing particle from a vertex
///
/// @param [in] vertex The vertex
/// @param [in] id The track ID of the particle
///
/// @return The particle pointer if found, else nullptr
HepMC3::ConstGenParticlePtr processParticle(HepMC3::ConstGenVertexPtr vertex,
                                            const int id) {
  // Loop over all outgoing particles
  for (const auto& particle : vertex->particles_out()) {
    const int trackid =
        particle->attribute<HepMC3::IntAttribute>("TrackID")->value();
    // Compare ID
    if (trackid == id)
      return particle;
  }
  return nullptr;
}

/// @brief This method collects the material in X_0 and L_0 a particle has
/// passed from its creation up to a certain vertex.
///
/// @param [in, out] particle The particle that get the passed material attached
/// @param [in] vertex The end vertex of the collection
/// @param [in] id The track ID
void passedMaterial(ActsExamples::SimParticle& particle,
                    const HepMC3::ConstGenVertexPtr& vertex, const int id) {
  double x0 = 0.;
  double l0 = 0.;
  HepMC3::ConstGenParticlePtr currentParticle = nullptr;
  HepMC3::ConstGenVertexPtr currentVertex = vertex;
  // Loop backwards and test whether the track still exists
  while (currentVertex && !currentVertex->particles_in().empty() &&
         currentVertex->particles_in()[0]->attribute<HepMC3::IntAttribute>(
             "TrackID") &&
         currentVertex->particles_in()[0]
                 ->attribute<HepMC3::IntAttribute>("TrackID")
                 ->value() == id) {
    // Get the step length
    currentParticle = currentVertex->particles_in()[0];
    const double stepLength =
        currentParticle->attribute<HepMC3::DoubleAttribute>("StepLength")
            ->value();
    // Add the passed material
    x0 +=
        stepLength /
        currentParticle->attribute<HepMC3::DoubleAttribute>("NextX0")->value();
    l0 +=
        stepLength /
        currentParticle->attribute<HepMC3::DoubleAttribute>("NextL0")->value();
    currentVertex = currentParticle->production_vertex();
  }
  // Assign the passed material to the particle
  particle.setMaterialPassed(x0, l0);
}

/// @brief This function collects outgoing particles from a vertex while keeping
/// track of the future of the ingoing particle.
///
/// @param [in] vertex The vertex
/// @param [in] trackID The track ID of the ingoing particle
///
/// @return Vector containing the outgoing particles from a vertex
std::vector<ActsExamples::SimParticle> outgoingParticles(
    HepMC3::ConstGenVertexPtr vertex, const int trackID) {
  std::vector<ActsExamples::SimParticle> finalStateParticles;

  // Identify the ingoing particle in the outgoing particles
  HepMC3::ConstGenParticlePtr procPart = processParticle(vertex, trackID);

  // Test whether this particle survives or dies
  HepMC3::ConstGenVertexPtr endVertex = procPart->end_vertex();
  if (endVertex
          ->attribute<HepMC3::StringAttribute>("NextProcessOf-" +
                                               std::to_string(trackID))
          ->value() != "Death") {
    // Store the particle if it survives
    finalStateParticles.push_back(
        ActsExamples::HepMC3Particle::particle(procPart));
  } else {
    // Store the leftovers if it dies
    for (const HepMC3::ConstGenParticlePtr procPartOut :
         endVertex->particles_out())
      if (procPartOut->attribute<HepMC3::IntAttribute>("TrackID")->value() ==
              trackID &&
          procPartOut->end_vertex()) {
        for (const HepMC3::ConstGenParticlePtr dyingPartOut :
             procPartOut->end_vertex()->particles_out())
          finalStateParticles.push_back(
              ActsExamples::HepMC3Particle::particle(dyingPartOut));
      }
  }

  // Record the particles produced in this process
  const std::vector<std::string> attributes = endVertex->attribute_names();
  for (const auto& att : attributes) {
    // Search for initial parameters
    if (att.find("InitialParametersOf") != std::string::npos) {
      const std::vector<double> mom4 =
          endVertex->attribute<HepMC3::VectorDoubleAttribute>(att)->value();
      const HepMC3::FourVector& pos4 = endVertex->position();
      const int id = stoi(att.substr(att.find("-") + 1));
      HepMC3::ConstGenParticlePtr genParticle = processParticle(endVertex, id);
      ActsFatras::Barcode barcode = ActsFatras::Barcode().setParticle(id);
      auto pid = static_cast<Acts::PdgParticle>(genParticle->pid());

      // Build an Acts particle out of the data
      ActsExamples::SimParticle simParticle(barcode, pid);
      simParticle.setPosition4(pos4.x(), pos4.y(), pos4.z(), pos4.t());
      Acts::Vector3D mom3(mom4[0], mom4[1], mom4[2]);
      simParticle.setDirection(mom3.normalized());
      simParticle.setAbsMomentum(mom3.norm());

      // Store the particle
      finalStateParticles.push_back(simParticle);
    }
  }

  return finalStateParticles;
}

/// @brief This method filters and sorts the recorded interactions.
///
/// @param [in, out] interactions The recorded interactions
/// @param [in] cfg Configuration of the filtering
void filterAndSort(std::vector<Particles>& interactions,
                   const ActsExamples::EventExtraction::Config& cfg) {
  for (Particles& interaction : interactions) {
    for (auto cit = std::get<2>(interaction).cbegin();
         cit != std::get<2>(interaction).cend();) {
      // Test whether a particle fulfills the conditions
      if (cit->pdg() < cfg.minAbsPdg || cit->pdg() > cfg.maxAbsPdg ||
          cit->absMomentum() < cfg.pMin)
        std::get<2>(interaction).erase(cit);
      else
        cit++;
    }
  }

  // Sort the particles based on their momentum
  for (Particles& interaction : interactions) {
    std::sort(std::get<2>(interaction).begin(), std::get<2>(interaction).end(),
              [](ActsExamples::SimParticle& a, ActsExamples::SimParticle& b) {
                return a.absMomentum() > b.absMomentum();
              });
  }
}
}  // namespace

ActsExamples::EventExtraction::~EventExtraction() {}

ActsExamples::EventExtraction::EventExtraction(
    ActsExamples::EventExtraction::Config&& cnf, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("EventExtraction", level),
      m_cfg(std::move(cnf)) {
  if (m_cfg.inputEvents.empty()) {
    throw std::invalid_argument("Missing input event collection");
  }
  if (m_cfg.outputEventFraction.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
  if (m_cfg.extractionProcess.empty()) {
    throw std::invalid_argument("Missing extraction process");
  }
}

ActsExamples::ProcessCode ActsExamples::EventExtraction::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // Retrieve the initial particles
  const auto events =
      context.eventStore.get<std::vector<HepMC3::GenEvent>>(m_cfg.inputEvents);

  std::vector<Particles> fractions;
  for (const HepMC3::GenEvent& event : events) {
    // Fast exit
    if (event.particles().empty() || event.vertices().empty())
      break;

    // Get the initial particle
    HepMC3::ConstGenParticlePtr initialParticle = event.particles()[0];
    ActsExamples::SimParticle simParticle =
        HepMC3Particle::particle(initialParticle);

    // Get the final state particles
    ActsExamples::SimParticle particleToInteraction;
    std::vector<ActsExamples::SimParticle> finalStateParticles;
    // Search the process vertex
    for (const auto& vertex : event.vertices()) {
      const std::vector<std::string> attributes = vertex->attribute_names();
      for (const auto& attribute : attributes) {
        if (vertex->attribute_as_string(attribute).find(
                m_cfg.extractionProcess) != std::string::npos) {
          const int procID = stoi(attribute.substr(attribute.find("-") + 1));
          // Get the particle before the interaction
          particleToInteraction = HepMC3Particle::particle(vertex->particles_in()[0]);
          // Attach passed material to the particle
          passedMaterial(particleToInteraction, vertex, procID);
          // Record the final state particles
          finalStateParticles = outgoingParticles(vertex, procID);
        }
      }
    }
    fractions.push_back(std::make_tuple(simParticle, particleToInteraction, finalStateParticles));
  }

  // Filter and sort the record
  filterAndSort(fractions, m_cfg);

  ACTS_INFO(events.size() << " processed");

  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.outputEventFraction, std::move(fractions));

  return ActsExamples::ProcessCode::SUCCESS;
}
