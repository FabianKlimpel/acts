// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <HepMC3/Units.h>
#include <HepMC3/ReaderRootTree.h>

bool ActsExamples::HepMC3AsciiReader::readEvent(HepMC3::ReaderAscii& reader,
                                                HepMC3::GenEvent& event) {
  // Read event and store it
  return reader.read_event(event);
}

bool ActsExamples::HepMC3AsciiReader::status(HepMC3::ReaderAscii& reader) {
  return !reader.failed();
}

ActsExamples::HepMC3AsciiReader::HepMC3AsciiReader(
    const ActsExamples::HepMC3AsciiReader::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(
          determineEventFilesRange(cfg.inputDir, cfg.inputStem + ".hepmc3")),
      m_logger(Acts::getDefaultLogger("HepMC3AsciiReader", lvl)),
      m_processExtractor(HepMC3ProcessExtractor(cfg.processExtractorCfg, lvl)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputEvents.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
}

std::string ActsExamples::HepMC3AsciiReader::HepMC3AsciiReader::name() const {
  return "HepMC3AsciiReader";
}

std::pair<size_t, size_t> ActsExamples::HepMC3AsciiReader::availableEvents()
    const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::HepMC3AsciiReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  std::vector<HepMC3::GenEvent> events;
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);

  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".hepmc3",
                               ctx.eventNumber);

  ACTS_DEBUG("Attempting to read event from " << path);
  //HepMC3::ReaderAscii reader(path);
  HepMC3::ReaderRootTree reader(path);

  ActsExamples::ExtractedSimulationProcessContainer interactions;

  reader.read_event(event);
  while (!reader.failed()) {
    interactions.push_back(m_processExtractor.execute(ctx, event));
    //events.push_back(std::move(event));
    event.clear();
    reader.read_event(event);
  }

  //if (events.empty())
  //  return ActsExamples::ProcessCode::ABORT;

  ACTS_VERBOSE(interactions.size() << " events read");
  //ctx.eventStore.add(m_cfg.outputEvents, std::move(interactions));
  ctx.eventStore.add("event-fraction", std::move(interactions));  

  reader.close();
  
  //m_processExtractor.execute(ctx, events);  
  return ActsExamples::ProcessCode::SUCCESS;
}
