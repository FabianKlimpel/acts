// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include "ActsExamples/Utilities/Paths.hpp"
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

ActsExamples::HepMC3AsciiWriter::HepMC3AsciiWriter(const Config&& cfg,
                                                   Acts::Logging::Level lvl)
    : WriterT(cfg.inputEvents, "HepMC3EventWriter", lvl), m_cfg(cfg) {
  if (m_cfg.outputStem.empty())
    throw std::invalid_argument("Missing output stem file name");
}

ActsExamples::ProcessCode ActsExamples::HepMC3AsciiWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const std::vector<HepMC3::GenEvent>& events) {
  auto path = perEventFilepath(m_cfg.outputDir, m_cfg.outputStem + ".hepmc3",
                               ctx.eventNumber);

  ACTS_DEBUG("Attempting to write event to " << path);
  HepMC3::WriterAscii writer(path);

std::cout << "Number of events: " << events.size() << std::endl;
for(unsigned int i = 78; i < events.size(); i++)
{
	std::cout << "Summary: " << i << " | " << events[i].particles().size() << " " <<  events[i].vertices().size() << std::endl;
	for(const auto& part : events[i].particles())
	{
		std::cout << "part id: " << part->id() << " | " << part->attribute_names().size() << std::endl;
		for(unsigned int j = 0; j < part->attribute_names().size(); j++)
		{
			std::cout << "Part Attr: " << j << " | " << part->attribute_names()[j] << " " << part->attribute_as_string(part->attribute_names()[j]) << std::endl;
		}
	}
	for(const auto& vtx : events[i].vertices())
	{
		std::cout << "vtx id: " << vtx->id() << " | " << vtx->attribute_names().size() << std::endl;
		for(unsigned int j = 0; j < vtx->attribute_names().size(); j++)
		{
			std::cout << "Vtx Attr: " << j << " | " << vtx->attribute_names()[j] << " " << vtx->attribute_as_string(vtx->attribute_names()[j]) << std::endl;
		}
	}
}
  for (const auto& event : events) {
    writer.write_event(event);
    if (writer.failed())
      return ActsExamples::ProcessCode::ABORT;
  }

  writer.close();
  return ActsExamples::ProcessCode::SUCCESS;
}
