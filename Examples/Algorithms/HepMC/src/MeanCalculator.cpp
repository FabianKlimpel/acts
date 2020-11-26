// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/HepMC/MeanCalculator.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Particle.hpp"

#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include "Acts/Propagator/Navigator.hpp"

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

//~ #include "Acts/Utilities/Definitions.hpp"

namespace {

//~ /// @brief This method searches for an outgoing particle from a vertex
//~ ///
//~ /// @param [in] vertex The vertex
//~ /// @param [in] id The track ID of the particle
//~ ///
//~ /// @return The particle pointer if found, else nullptr
//~ HepMC3::ConstGenParticlePtr searchProcessParticleById(
    //~ HepMC3::ConstGenVertexPtr vertex, const int id) {
  //~ // Loop over all outgoing particles
  //~ for (const auto& particle : vertex->particles_out()) {
    //~ const int trackid =
        //~ particle->attribute<HepMC3::IntAttribute>("TrackID")->value();
    //~ // Compare ID
    //~ if (trackid == id) {
      //~ return particle;
    //~ }
  //~ }
  //~ return nullptr;
//~ }

//~ /// @brief This method collects the material in X_0 and L_0 a particle has
//~ /// passed from its creation up to a certain vertex.
//~ ///
//~ /// @param [in] vertex The end vertex of the collection
//~ /// @param [in] id The track ID
//~ /// @param [in, out] particle The particle that get the passed material attached
//~ void setPassedMaterial(const HepMC3::ConstGenVertexPtr& vertex, const int id,
                       //~ ActsExamples::SimParticle& particle) {
  //~ double x0 = 0.;
  //~ double l0 = 0.;
  //~ HepMC3::ConstGenParticlePtr currentParticle = nullptr;
  //~ HepMC3::ConstGenVertexPtr currentVertex = vertex;
  //~ // Loop backwards and test whether the track still exists
  //~ while (currentVertex && !currentVertex->particles_in().empty() &&
         //~ currentVertex->particles_in()[0]->attribute<HepMC3::IntAttribute>(
             //~ "TrackID") &&
         //~ currentVertex->particles_in()[0]
                 //~ ->attribute<HepMC3::IntAttribute>("TrackID")
                 //~ ->value() == id) {
    //~ // Get the step length
    //~ currentParticle = currentVertex->particles_in()[0];
    //~ const double stepLength =
        //~ currentParticle->attribute<HepMC3::DoubleAttribute>("StepLength")
            //~ ->value();
    //~ // Add the passed material
    //~ x0 +=
        //~ stepLength /
        //~ currentParticle->attribute<HepMC3::DoubleAttribute>("NextX0")->value();
    //~ l0 +=
        //~ stepLength /
        //~ currentParticle->attribute<HepMC3::DoubleAttribute>("NextL0")->value();
    //~ currentVertex = currentParticle->production_vertex();
  //~ }
  //~ // Assign the passed material to the particle
  //~ particle.setMaterialPassed(x0, l0);
//~ }

//~ /// @brief This function collects outgoing particles from a vertex while keeping
//~ /// track of the future of the ingoing particle.
//~ ///
//~ /// @param [in] vertex The vertex
//~ /// @param [in] trackID The track ID of the ingoing particle
//~ ///
//~ /// @return Vector containing the outgoing particles from a vertex
//~ std::vector<ActsExamples::SimParticle> selectOutgoingParticles(
    //~ HepMC3::ConstGenVertexPtr vertex, const int trackID) {
  //~ std::vector<ActsExamples::SimParticle> finalStateParticles;

  //~ // Identify the ingoing particle in the outgoing particles
  //~ HepMC3::ConstGenParticlePtr procPart =
      //~ searchProcessParticleById(vertex, trackID);

  //~ // Test whether this particle survives or dies
  //~ HepMC3::ConstGenVertexPtr endVertex = procPart->end_vertex();
  //~ if (endVertex
          //~ ->attribute<HepMC3::StringAttribute>("NextProcessOf-" +
                                               //~ std::to_string(trackID))
          //~ ->value() != "Death") {
    //~ // Store the particle if it survives
    //~ finalStateParticles.push_back(
        //~ ActsExamples::HepMC3Particle::particle(procPart));
  //~ } else {
    //~ // Store the leftovers if it dies
    //~ for (const HepMC3::ConstGenParticlePtr& procPartOut :
         //~ endVertex->particles_out())
      //~ if (procPartOut->attribute<HepMC3::IntAttribute>("TrackID")->value() ==
              //~ trackID &&
          //~ procPartOut->end_vertex()) {
        //~ for (const HepMC3::ConstGenParticlePtr dyingPartOut :
             //~ procPartOut->end_vertex()->particles_out()) {
          //~ finalStateParticles.push_back(
              //~ ActsExamples::HepMC3Particle::particle(dyingPartOut));
        //~ }
      //~ }
  //~ }

  //~ // Record the particles produced in this process
  //~ const std::vector<std::string> attributes = endVertex->attribute_names();
  //~ for (const auto& att : attributes) {
    //~ // Search for initial parameters
    //~ if (att.find("InitialParametersOf") != std::string::npos) {
      //~ const std::vector<double> mom4 =
          //~ endVertex->attribute<HepMC3::VectorDoubleAttribute>(att)->value();
      //~ const HepMC3::FourVector& pos4 = endVertex->position();
      //~ const int id = stoi(att.substr(att.find("-") + 1));
      //~ HepMC3::ConstGenParticlePtr genParticle =
          //~ searchProcessParticleById(endVertex, id);
      //~ ActsFatras::Barcode barcode = ActsFatras::Barcode().setParticle(id);
      //~ auto pid = static_cast<Acts::PdgParticle>(genParticle->pid());

      //~ // Build an Acts particle out of the data
      //~ ActsExamples::SimParticle simParticle(barcode, pid);
      //~ simParticle.setPosition4(pos4.x(), pos4.y(), pos4.z(), pos4.t());
      //~ Acts::Vector3D mom3(mom4[0], mom4[1], mom4[2]);
      //~ simParticle.setDirection(mom3.normalized());
      //~ simParticle.setAbsMomentum(mom3.norm());

      //~ // Store the particle
      //~ finalStateParticles.push_back(simParticle);
    //~ }
  //~ }

  //~ return finalStateParticles;
//~ }

//~ /// @brief This method filters and sorts the recorded interactions.
//~ ///
//~ /// @param [in] cfg Configuration of the filtering
//~ /// @param [in, out] interactions The recorded interactions
//~ void filterAndSort(
    //~ const ActsExamples::MeanCalculator::Config& cfg,
    //~ std::vector<ActsExamples::ExtractedSimulationProcess>& interactions) {
  //~ for (auto& interaction : interactions) {
    //~ for (auto cit = interaction.after.cbegin();
         //~ cit != interaction.after.cend();) {
      //~ // Test whether a particle fulfills the conditions
      //~ if (cit->pdg() < cfg.absPdgMin || cit->pdg() > cfg.absPdgMax ||
          //~ cit->absMomentum() < cfg.pMin) {
        //~ interaction.after.erase(cit);
      //~ } else {
        //~ cit++;
      //~ }
    //~ }
  //~ }

  //~ // Sort the particles based on their momentum
  //~ for (auto& interaction : interactions) {
    //~ std::sort(interaction.after.begin(), interaction.after.end(),
              //~ [](ActsExamples::SimParticle& a, ActsExamples::SimParticle& b) {
                //~ return a.absMomentum() > b.absMomentum();
              //~ });
  //~ }
//~ }

std::vector<ActsExamples::SimParticle>
collectG4Steps(const HepMC3::GenEvent& event, int trackID) {
	std::vector<ActsExamples::SimParticle> g4Steps;
	for (const auto& vertex : event.vertices()) {
		const auto material = vertex->attribute<HepMC3::StringAttribute>("Material")->value();
		if(material == "NoMaterial" || material == "Vacuum" || material == "Air")
			continue;
		for(const auto& particle :vertex->particles_out()) {
			if(particle->attribute<HepMC3::IntAttribute>("TrackID")->value() == trackID)
			{
				const auto& posVtx = vertex->position();
				const Acts::Vector4D position(posVtx.x(), posVtx.y(), posVtx.z(), posVtx.t());
				
				const auto& momPart = particle->momentum();
				const Acts::Vector3D momentum(momPart.x(), momPart.y(), momPart.z());
				
				ActsExamples::SimParticle g4Particle;
				g4Particle.setPosition4(position).setDirection(momentum.normalized()).setAbsMomentum(momentum.norm());
				g4Steps.push_back(g4Particle);
				break;
			}
		}
	}
	return g4Steps;
}

Acts::BoundVector
findClosestPoint(const std::vector<ActsExamples::SimParticle>& g4Steps, std::shared_ptr<const Acts::Surface> surface, const Acts::GeometryContext& gctx) {
	std::vector<std::pair<double, Acts::BoundVector>> pathLengthPosition;
	Acts::BoundaryCheck bCheck(false);
	for(const ActsExamples::SimParticle& g4Step : g4Steps)
	{
		const Acts::SurfaceIntersection intersection = surface->intersect(gctx, g4Step.position(), g4Step.unitDirection(), bCheck);
		if(intersection)
		{				
			const Acts::Vector3D pos3 = intersection.intersection.position;
			const Acts::Vector3D mom3 = g4Step.momentum4().template segment<3>(Acts::eMom0);
			//~ const Acts::Vector2D pos2 = surface->globalToLocal(gctx, pos3, mom3).value(); // Result<Vector2D>
			Acts::FreeVector freeParams;
			freeParams[eFreePos0] = pos3[eX];
			freeParams[eFreePos1] = pos3[eY];
			freeParams[eFreePos2] = pos3[eZ];
			freeParams[eFreeTime] = g4Step.time();
			freeParams[eFreeDir0] = g4Step.unitDirection()[eMom0];
			freeParams[eFreeDir1] = g4Step.unitDirection()[eMom0];
			freeParams[eFreeDir2] = g4Step.unitDirection()[eMom0];
			freeParams[eFreeQOverP] = (g4Step.charge() == 0. ? 1. : g4Step.charge()) / g4Step.absMomentum();
			
			Acts::BoundVector params = transformFreeToBoundParameters(freeParams, surface, gctx);
			
			const double pathLength = intersection.intersection.pathLength;
			pathLengthPosition.push_back(std::make_pair(pathLength, params));
			// TODO: test posG4 value
		}
	}
	const auto closest = std::min_element(pathLengthPosition.begin(), pathLengthPosition.end(), 
			[&](const std::pair<double, Acts::Vector2D>& pos1, const std::pair<double, Acts::Vector2D>& pos2) 
				{ return  pos1.first < pos2.first; });
	return closest->second;
}

Acts::BoundVector
mean(const std::vector<Acts::BoundVector>& positions) {
	
	Acts::BoundVector mean = Acts::BoundVector::Zero();
	for(const Acts::BoundVector& position : positions)
	{
		mean += position;
	}
	mean /= positions.size();
	return mean;
}

void
plotMean(const std::vector<std::pair<Acts::BoundVector, Acts::BoundVector>>& meanProbG4) {
	
	// TODO: should rather do this for all events once
	// TODO: split the plot in e.g. vs. r, phi, z, eta, ...
}

}  // namespace

ActsExamples::MeanCalculator::~MeanCalculator() {}

ActsExamples::MeanCalculator::MeanCalculator(
    ActsExamples::MeanCalculator::Config&& cnf,
    Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("MeanCalculator", level),
      m_cfg(std::move(cnf)) {
  if (m_cfg.inputEvents.empty()) {
    throw std::invalid_argument("Missing input event collection");
  }
  if (m_cfg.inputParticles.empty()) {
	throw std::invalid_argument("Missing input event collection");
  }
}

ActsExamples::ProcessCode ActsExamples::MeanCalculator::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // Retrieve the initial particles
  const auto events =
      context.eventStore.get<std::vector<HepMC3::GenEvent>>(m_cfg.inputEvents);

  // Retrieve the initial particles
  const auto initialParticles =
      context.eventStore.get<ActsExamples::SimParticleContainer>(
          m_cfg.inputParticles);
  
  // The stepper
  Acts::NullBField bfield;  
  Acts::EigenStepper stepper(bfield);
  
  // The Navigator
  Acts::Navigator navigator(m_cfg.trackingGeometry);
  
  // The Propagator
  Acts::Propagator propagator(stepper, navigator);
  
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::PropagatorOptions<Acts::ActionList<Acts::detail::SteppingLogger>> options(gctx, mctx, Acts::getDummyLogger());
  
  for(const ActsExamples::SimParticle& initialParticle : initialParticles)
  {
	  // Propagate the mean
	Acts::CurvilinearTrackParameters mean(initialParticle.position4(), initialParticle.unitDirection(), initialParticle.charge(), initialParticle.absMomentum());
	const auto& result = propagator.propagate(mean, options).value(); //result.ok()
	const auto stepperLog = result.get<typename Acts::detail::SteppingLogger::result_type>();
	
	// Walk over each step
	for(const auto& step : stepperLog.steps)
	{
		// Only care about surfaces
		if(!step.surface)
			continue;
		// Calculate the value of the mean on the surface
		// TODO: This transformation should be transformationFreeToBound(...)
		const BoundVector localPropagatedMean = step.surface->globalToLocal(gctx, step.position, step.momentum); // Result<Vector2D>
		
		// Now find the corresponding G4 steps
		std::vector<Acts::BoundVector> localG4Params;
		localG4Positions.reserve(events.size());
		
		for(const auto& event : events)
		{
			// Fast continue
			if(event.vertices().empty())
				continue;
			// Get the track ID that we follow
			const int trackID = event.vertices()[0]->particles_out()[0]->attribute<HepMC3::IntAttribute>("TrackID")->value();
			// The storage of each step
			std::vector<ActsExamples::SimParticle> g4Steps = collectG4Steps(event, trackID);
			
			localG4Params.push_back(findClosestPoint(g4Steps, step.surface, gctx));
		}
		
		const Acts::BoundVector meanG4 = mean(localG4Params);
	}
	// TODO: plot either here or below
	}
                   
  //~ std::vector<ActsExamples::ExtractedSimulationProcess> fractions;
  //~ for (const HepMC3::GenEvent& event : events) {
    //~ // Fast exit
    //~ if (event.particles().empty() || event.vertices().empty()) {
      //~ break;
    //~ }

    //~ // Get the initial particle
    //~ HepMC3::ConstGenParticlePtr initialParticle = event.particles()[0];
    //~ ActsExamples::SimParticle simParticle =
        //~ HepMC3Particle::particle(initialParticle);

    //~ // Get the final state particles
    //~ ActsExamples::SimParticle particleToInteraction;
    //~ std::vector<ActsExamples::SimParticle> finalStateParticles;
    //~ // Search the process vertex
    //~ bool vertexFound = false;
    //~ for (const auto& vertex : event.vertices()) {
      //~ const std::vector<std::string> attributes = vertex->attribute_names();
      //~ for (const auto& attribute : attributes) {
        //~ if (vertex->attribute_as_string(attribute).find(
                //~ m_cfg.extractionProcess) != std::string::npos) {
          //~ const int procID = stoi(attribute.substr(attribute.find("-") + 1));
          //~ // Get the particle before the interaction
          //~ particleToInteraction =
              //~ HepMC3Particle::particle(vertex->particles_in()[0]);
          //~ // Attach passed material to the particle
          //~ setPassedMaterial(vertex, procID, particleToInteraction);
          //~ // Record the final state particles
          //~ finalStateParticles = selectOutgoingParticles(vertex, procID);
          //~ vertexFound = true;
          //~ break;
        //~ }
      //~ }
      //~ if (vertexFound) {
        //~ break;
      //~ }
    //~ }
    //~ fractions.push_back(ActsExamples::ExtractedSimulationProcess{
        //~ simParticle, particleToInteraction, finalStateParticles});
  //~ }

  //~ // Filter and sort the record
  //~ filterAndSort(m_cfg, fractions);

  //~ ACTS_INFO(events.size() << " processed");

  //~ // Write the recorded material to the event store
  //~ context.eventStore.add(m_cfg.outputSimulationProcesses, std::move(fractions));

  return ActsExamples::ProcessCode::SUCCESS;
}