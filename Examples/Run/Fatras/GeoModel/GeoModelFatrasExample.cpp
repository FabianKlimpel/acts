// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Geant4HepMC/EventRecording.hpp"
#include "ActsExamples/Geant4/Geant4Options.hpp"

#include "ActsExamples/Geant4GeoModel/GdmlDetectorConstruction.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

#include "Fatras.hpp"
#include <boost/program_options.hpp>

namespace {

using namespace ActsExamples;

// collection names
static constexpr const char* kParticlesInput = "particles_input";
static constexpr const char* kParticlesSelection = "particles_selection";
static constexpr const char* kParticlesInitial = "particles_initial";
static constexpr const char* kParticlesFinal = "particles_final";
static constexpr const char* kSimHits = "simhits";

// input handling
void addInputOptions(ActsExamples::Options::Description& desc) {
  ActsExamples::Options::addParticleGunOptions(desc);
  /// TODO(add options for pythia)
  //~ ActsExamples::ParticleSelector::addOptions(desc);
}

void setupInput(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers) {
  auto logLevel = Options::readLogLevel(vars);

	/// TODO(use either pythia or pgun)
    // generate particle input from a particle gun
    auto gen = Options::readParticleGunOptions(vars);
    gen.outputParticles = kParticlesInput;
    gen.randomNumbers = randomNumbers;
    sequencer.addReader(std::make_shared<EventGenerator>(gen, logLevel));

  // add additional particle selection
  //~ auto select = ActsExamples::ParticleSelector::readConfig(vars);
  ActsExamples::ParticleSelector::Config select;
  select.inputParticles = kParticlesInput;
  select.outputParticles = kParticlesSelection;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::ParticleSelector>(select, logLevel));
}

// output handling
// output options are just the common output options
void setupOutput(const ActsExamples::Options::Variables& vars,
                 ActsExamples::Sequencer& sequencer) {
  auto logLevel = Options::readLogLevel(vars);
  auto outputDir =
      ensureWritableDirectory(vars["output-dir"].template as<std::string>());

    // write simulated particle initial states
    RootParticleWriter::Config writeInitial;
    writeInitial.inputParticles = kParticlesInitial;
    writeInitial.filePath = joinPaths(outputDir, "particles_initial.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(writeInitial, logLevel));

    // write simulated particle final states
    RootParticleWriter::Config writeFinal;
    writeFinal.inputParticles = kParticlesFinal;
    writeFinal.filePath = joinPaths(outputDir, "particles_final.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(writeFinal, logLevel));
        
    // write simulated hits
    RootSimHitWriter::Config writeHits;
    writeHits.inputSimHits = kSimHits;
    writeHits.filePath = joinPaths(outputDir, "hits.root");
    sequencer.addWriter(
        std::make_shared<RootSimHitWriter>(writeHits, logLevel));
}

// simulation handling

void setupSimulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  auto logLevel = Options::readLogLevel(vars);
  auto fatrasCfg = FatrasAlgorithm::readConfig(vars);
  fatrasCfg.inputParticles = kParticlesSelection;
  fatrasCfg.outputParticlesInitial = kParticlesInitial;
  fatrasCfg.outputParticlesFinal = kParticlesFinal;
  fatrasCfg.outputSimHits = kSimHits;
  fatrasCfg.randomNumbers = randomNumbers;
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.magneticField = ActsExamples::Options::readMagneticField(vars);

  sequencer.addAlgorithm(
      std::make_shared<FatrasAlgorithm>(std::move(fatrasCfg), logLevel));
}
}  // namespace

int main(int argc, char* argv[]) {
	auto detector = std::make_shared<TGeoDetector>();
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  // add general and detector-specific geometry options
  Options::addGeometryOptions(desc);
  detector->addOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addMagneticFieldOptions(desc);
  // algorithm-specific options
  FatrasAlgorithm::addOptions(desc);
  Options::addGeant4Options(desc);
  desc.add_options()(
      "gdml-to-gm-plugin",
      boost::program_options::value<std::string>()->default_value(""),
      "Path to libGDMLtoGM.so. Note: The gdml file has to be named "
      "gdmlfile.xml.");
      
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic services
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // setup sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));
  // setup detector geometry and material and the magnetic field
  auto [trackingGeometry, contextDecorators] = Geometry::build(vars, *detector);
std::cout << "TRACKING GEOMETRY: " << trackingGeometry << std::endl;       
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }
  // setup algorithm chain
  setupInput(vars, sequencer, randomNumbers);
  setupOutput(vars, sequencer);
  
  /// TODO: G4 BField
  // Prepare the detector
  auto gdmlFile = vars["gdml-to-gm-plugin"].as<std::string>();
  auto g4detector = std::make_unique<GdmlDetectorConstruction>(gdmlFile);
  
  EventRecording::Config erConfig;
  erConfig.inputParticles = kParticlesSelection;
  erConfig.outputHepMcTracks = "geant-outcome-tracks";
  erConfig.detectorConstruction = std::move(g4detector);
  erConfig.seed1 = vars["g4-rnd-seed1"].as<unsigned int>();
  erConfig.seed2 = vars["g4-rnd-seed2"].as<unsigned int>();
  //~ erConfig.processesReject = {"pi+Inelastic", "pi-Inelastic"}; 
  //~ erConfig.processesReject = {"kaon0SInelastic"}; 
  //~ erConfig.processSelect = "Decay"; 

    // write simulated particle final states
    RootParticleWriter::Config writeFinal;
    writeFinal.inputParticles = erConfig.outputHepMcTracks;
    writeFinal.filePath = joinPaths(ensureWritableDirectory(vars["output-dir"].template as<std::string>()), "particles_final_geant4.root");
    
  auto logLevel = Options::readLogLevel(vars);
  sequencer.addAlgorithm(
      std::make_shared<EventRecording>(std::move(erConfig), logLevel));
  sequencer.addWriter(
        std::make_shared<RootParticleWriter>(writeFinal, logLevel));
  setupSimulation(vars, sequencer, randomNumbers, trackingGeometry);
    
  // run the simulation
  return sequencer.run();
}