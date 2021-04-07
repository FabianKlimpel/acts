// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3ProcessExtractor.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Options.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Io/NuclearInteractions/RootNuclearInteractionParametersWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"

namespace {
void
addMyOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("multiplicity-Min",
      value<unsigned int>()->default_value({}),
      "Minimum multiplicity.");
  opt("multiplicity-Max",
      value<unsigned int>()->default_value({}),
      "Maximum multiplicity.");
  opt("record-Soft",
      value<bool>()->default_value({}),
      "Record soft.");
  opt("write-hist",
      value<bool>()->default_value({}),
      "Write histograms.");
  opt("num-simulated-events",
      value<unsigned int>()->default_value({}),
      "Number of simulated events."); 
  opt("parametrisation-filename",
	  value<std::string>()->default_value("parameters.root"),
	  "Filename of the parametrisation");
}

ActsExamples::RootNuclearInteractionParametersWriter::Config readMyConfig(
    const boost::program_options::variables_map& variables) {
  
   ActsExamples::RootNuclearInteractionParametersWriter::Config cfg;
   cfg.multiplicityMin = variables["multiplicity-Min"].as<unsigned int>();
   cfg.multiplicityMax = variables["multiplicity-Max"].as<unsigned int>();
   cfg.recordSoft = variables["record-Soft"].as<bool>();
   cfg.writeOptionalHistograms = variables["write-hist"].as<bool>();
   cfg.nSimulatedEvents = variables["num-simulated-events"].as<unsigned int>();
   cfg.outputFilename = variables["parametrisation-filename"].as<std::string>();

	return cfg;
}
}

///
/// Straight forward example of reading a HepMC3 file.
///
int main(int argc, char** argv) {
  // Declare the supported program options.
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addInputOptions(desc);
  ActsExamples::Options::addHepMC3ReaderOptions(desc);
  addMyOptions(desc);

  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Create the reader
  auto hepMC3ReaderConfig = ActsExamples::Options::readHepMC3ReaderOptions(vm);
  hepMC3ReaderConfig.outputEvents = "hepmc-events";

  ActsExamples::HepMC3ProcessExtractor::Config extractionConfig;
  extractionConfig.inputEvents = hepMC3ReaderConfig.outputEvents;
  extractionConfig.extractionProcess = "Inelastic";

  ActsExamples::RootNuclearInteractionParametersWriter::Config writerCfg = readMyConfig(vm);
  writerCfg.inputSimulationProcesses =
      extractionConfig.outputSimulationProcesses;
  hepMC3ReaderConfig.processExtractorCfg = extractionConfig;

  // Add to the sequencer
  sequencer.addReader(std::make_shared<ActsExamples::HepMC3AsciiReader>(
      hepMC3ReaderConfig, logLevel));
  //~ sequencer.addAlgorithm(std::make_shared<ActsExamples::HepMCProcessExtractor>(
      //~ std::move(extractionConfig), logLevel));
  sequencer.addWriter(
      std::make_shared<ActsExamples::RootNuclearInteractionParametersWriter>(
          writerCfg, logLevel));

  // Run
  return sequencer.run();
}
