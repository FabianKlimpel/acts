// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4HepMC/EventRecording.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

#include <iostream>
#include <stdexcept>

#include <FTFP_BERT.hh>
#include <HepMC3/GenParticle.h>
#include <HepMC3/WriterAscii.h>

#include "EventAction.hpp"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

#include "ActsExamples/Utilities/Paths.hpp"

namespace {

ActsExamples::SimParticleContainer
finalStateParticles(const std::vector<HepMC3::GenEvent>& events) {
	ActsExamples::SimParticleContainer finalStates;
	unsigned int eventCounter = 0;
	for(const auto& evt : events)
	{
		unsigned int counter = 0;
		int currentTrackID = 0;
		for(const auto& part : evt.particles())
		{
			//~ if(part->end_vertex() == nullptr && part->production_vertex()->id() != 0)
			if(part->attribute<HepMC3::IntAttribute>("TrackID") != nullptr)
			{
				HepMC3::FourVector mom = part->momentum();
				Acts::Vector3 mom3{mom.x(), mom.y(), mom.z()};
				mom3 *= Acts::UnitConstants::GeV;
				if(mom3.norm() < 50. * Acts::UnitConstants::MeV)
					continue;
				
				int pid = part->pid();
				if(currentTrackID != part->attribute<HepMC3::IntAttribute>("TrackID")->value())
				{
					currentTrackID = part->attribute<HepMC3::IntAttribute>("TrackID")->value();
					counter = 0;
				}
				ActsFatras::Barcode barcode;
				barcode.setVertexPrimary(eventCounter);
				barcode.setVertexSecondary(currentTrackID);
				barcode.setSubParticle(counter++);
				ActsFatras::Particle particle(barcode, static_cast<Acts::PdgParticle>(pid));
				
				auto vtx = part->production_vertex();
				HepMC3::FourVector pos = vtx->position();
				particle.setPosition4(pos.x(), pos.y(), pos.z(), pos.t());
				
				particle.setDirection(mom3.normalized());
				particle.setAbsoluteMomentum(mom3.norm());

				if(part->attribute<HepMC3::DoubleAttribute>("ProperTime") != nullptr)
					particle.setProperTime(part->attribute<HepMC3::DoubleAttribute>("ProperTime")->value() * Acts::UnitConstants::s);

				finalStates.insert(particle);
			}
		}
		eventCounter++;
	}
	
	return finalStates;
}
	
}

ActsExamples::EventRecording::~EventRecording() {
  m_runManager = nullptr;
}

ActsExamples::EventRecording::EventRecording(
    ActsExamples::EventRecording::Config&& cnf, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("EventRecording", level),
      m_cfg(std::move(cnf)),
      m_runManager(std::make_unique<G4RunManager>()) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (m_cfg.outputHepMcTracks.empty()) {
    throw std::invalid_argument("Missing output event collection");
  }
  if (!m_cfg.detectorConstruction) {
    throw std::invalid_argument("Missing detector construction object");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(m_cfg.detectorConstruction.release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new ActsExamples::RunAction());
  m_runManager->SetUserAction(
      new ActsExamples::EventAction(m_cfg.processesCombine));
  m_runManager->SetUserAction(
      new ActsExamples::PrimaryGeneratorAction(m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(
      new ActsExamples::SteppingAction(m_cfg.processesReject));
  m_runManager->Initialize();
}

ActsExamples::ProcessCode ActsExamples::EventRecording::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  auto path = perEventFilepath("", "2.hepmc3",
                               context.eventNumber);
  ACTS_DEBUG("Attempting to write event to " << path);
  HepMC3::WriterAscii writer(path);
                              
  // Retrieve the initial particles
  const auto initialParticles =
      context.eventStore.get<ActsExamples::SimParticleContainer>(
          m_cfg.inputParticles);

  // Storage of events that will be produced
  std::vector<HepMC3::GenEvent> events;
  events.reserve(initialParticles.size());

  for (const auto& part : initialParticles) {
    // Prepare the particle gun
    ActsExamples::PrimaryGeneratorAction::instance()->prepareParticleGun(part);

    // Begin with the simulation
    m_runManager->BeamOn(1);

    // Test if the event was aborted
    if (SteppingAction::instance()->eventAborted()) {
      continue;
    }

    // Set event start time
    HepMC3::GenEvent event = ActsExamples::EventAction::instance()->event();
    HepMC3::FourVector shift(0., 0., 0., part.time() / Acts::UnitConstants::mm);
    event.shift_position_by(shift);

    // Set beam particle properties
    const Acts::Vector4 momentum4 =
        part.fourMomentum() / Acts::UnitConstants::GeV;
    HepMC3::FourVector beamMom4(momentum4[0], momentum4[1], momentum4[2],
                                momentum4[3]);
    auto beamParticle = event.particles()[0];
    beamParticle->set_momentum(beamMom4);
    beamParticle->set_pid(part.pdg());

    if (m_cfg.processSelect.empty()) {
		// Remove vertices without outgoing particles
        for (auto it = event.vertices().crbegin();
             it != event.vertices().crend(); it++) {
          if ((*it)->particles_out().empty()) {
            event.remove_vertex(*it);
          }
        }
      // Store the result
      events.push_back(std::move(event));
std::cout << "Writing Event " << events.size() << std::endl;
      writer.write_event(events.back());
std::cout << "Writing Done " << std::endl;      
      if (writer.failed())
		  return ActsExamples::ProcessCode::ABORT;
    } else {
      bool storeEvent = false;
      // Test if the event has a process of interest in it
      for (const auto& vertex : event.vertices()) {
        if (vertex->id() == -1) {
          vertex->add_particle_in(beamParticle);
        }
        const std::vector<std::string> vertexAttributes =
            vertex->attribute_names();
        for (const auto& att : vertexAttributes) {
          if ((vertex->attribute_as_string(att).find(m_cfg.processSelect) !=
               std::string::npos) &&
              !vertex->particles_in().empty() &&
              vertex->particles_in()[0]->attribute<HepMC3::IntAttribute>(
                  "TrackID") &&
              vertex->particles_in()[0]
                      ->attribute<HepMC3::IntAttribute>("TrackID")
                      ->value() == 1) {
            storeEvent = true;
            break;
          }
        }
        if (storeEvent) {
          break;
        }
      }
      // Store the result
      if (storeEvent) {
        // Remove vertices without outgoing particles
        for (auto it = event.vertices().crbegin();
             it != event.vertices().crend(); it++) {
          if ((*it)->particles_out().empty()) {
            event.remove_vertex(*it);
          }
        }
        events.push_back(std::move(event));
std::cout << "Writing Event " << events.size() << std::endl;
      writer.write_event(events.back());
std::cout << "Writing Done " << std::endl;      
      if (writer.failed())
		  return ActsExamples::ProcessCode::ABORT;
      }
    }
  }

  ACTS_INFO(initialParticles.size() << " initial particles provided");
  ACTS_INFO(events.size() << " tracks generated");

  SimParticleContainer finalState = finalStateParticles(events);
  ACTS_INFO(finalState.size() << " final state particles found");
  
std::cout << "Closing the writer" << std::endl;
writer.close();
std::cout << "Writer closed" << std::endl;
  
  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.outputHepMcTracks, std::move(events));
  //~ context.eventStore.add(m_cfg.outputHepMcTracks, std::move(finalState));

  return ActsExamples::ProcessCode::SUCCESS;
}
