// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include "ActsExamples/Utilities/Paths.hpp"
#include <HepMC3/WriterRootTree.h>
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

  //HepMC3::WriterAscii writer(path);
  HepMC3::WriterRootTree writer(path);

std::cout << "Number of events: " << events.size() << std::endl;
/*
	unsigned int ci = 0, cj = 0;
  //~ for (const auto& event : events) {
  for(unsigned int i = 78; i < events.size(); i++) {
	 
	 const auto attributes = events[i].attributes();
	 for(const auto& attribute : attributes)
	 {
		 std::cout << "First Attribute: " << attribute.first << std::endl;
		 for(const auto& att : attribute.second)
		 {
			 std::cout << "First att: " << att.first << " | " << attributes.size() << " " << attribute.second.size() << " | " << ci << " " << cj << std::endl;
			 std::cout << "shared ptr: " << att.second << std::endl;
			 if(att.second != nullptr)
			 {
				std::string str;
				att.second->to_string(str);
				std::cout << "to String: " << str << std::endl; 
			 }
			 cj++;
		 }
		 ci++;
		 cj = 0;
	 }
	 
	 std::cout << "Loop through map finished" << std::endl;
	  
	std::cout << "Summary: " << i << " | " << events[i].particles().size() << " " <<  events[i].vertices().size() << std::endl;
	for(const auto& part : events[i].particles())
	{
		if(part == nullptr)
		{
			std::cout << "PART IS NULLPTR!!!!!" << std::endl;
			continue;
		}
		//~ std::cout << "part id: " << part->id() << " | " << part->attribute_names().size() << std::endl;
		//~ for(unsigned int j = 0; j < part->attribute_names().size(); j++)
		//~ {
			//~ std::cout << j << " | " << part->attribute_names()[j] << " " << part->attribute_as_string(part->attribute_names()[j]) << " || ";
		//~ }
		//~ std::cout << std::endl; 
	}
	for(const auto& vtx : events[i].vertices())
	{
		if(vtx == nullptr)
		{
			std::cout << "VTX IS NULLPTR!!!!!" << std::endl;
			continue;
		}
		//~ std::cout << "vtx id: " << vtx->id() << " | " << vtx->attribute_names().size() << std::endl;
		//~ for(unsigned int j = 0; j < vtx->attribute_names().size(); j++)
		//~ {
			//~ std::cout << j << " | " << vtx->attribute_names()[j] << " " << vtx->attribute_as_string(vtx->attribute_names()[j]) << " || ";
		//~ }
		//~ std::cout << std::endl;
	}
	*/
	
    writer.write_event(events[i]);
    if (writer.failed())
      return ActsExamples::ProcessCode::ABORT;
      
    std::cout << i << " finished" << std::endl;
  }
std::cout << "Writer ran through" << std::endl;
  writer.close();
std::cout << "Writer closed" << std::endl;
  return ActsExamples::ProcessCode::SUCCESS;
}
