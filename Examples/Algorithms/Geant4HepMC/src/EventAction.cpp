// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "EventAction.hpp"

#include <stdexcept>

#include <G4Event.hh>
#include <G4RunManager.hh>

#include "SteppingAction.hpp"

namespace {
/// @brief This function writes an event to output
///
/// @param [in] evt The event to be written
/// @param [in] outputname The file name (optional)
///
inline void save_event(const HepMC3::GenEvent* evt, const std::string outputname="", const int line=0) {
  std::string filename=outputname;
  if (filename=="") {
    std::string thisfile(__FILE__);
    size_t slash=thisfile.find_last_of("/");
    filename=thisfile.substr(slash)+"_"+std::to_string(line)+"_"+std::to_string(evt->event_number())+".hepmc3";
  }
  HepMC3::WriterAscii writer(filename);
  writer.write_event(*evt);
  writer.close();
}

/// @brief This function tests whether a process is available in the record
///
/// @param [in] vertex The vertex that will be tested
/// @param [in] processFilter List of processes that will be filtered
///
/// @return True if the process was found, false if not
bool findAttribute(HepMC3::ConstGenVertexPtr vertex,
                   const std::vector<std::string>& processFilter) {
  if (!vertex) return false;
  // Consider only 1->1 vertices to keep a correct history
  if ((vertex->particles_in().size() == 1) &&
      (vertex->particles_out().size() == 1)) {
    if (processFilter.empty()) return true;
    // Test for all attributes if one matches the filter pattern
    const std::vector<std::string> vertexAttributes = vertex->attribute_names();
    for (const auto& att : vertexAttributes) {
      const std::string process = vertex->attribute_as_string(att);
      if (std::find(processFilter.begin(), processFilter.end(), process) !=
          processFilter.end()) {
        return true;
      }
    }
  }
  return false;
}

/// @brief This function reduces multiple vertices that should be filtered and
/// combines them in a single dummy vertex
///
/// @param [in, out] event The underlying event
/// @param [in, out] vertex The vertex that will be reduced
/// @param [in] processFilter List of processes that will be filtered
void reduceVertex(HepMC3::GenEvent& event, HepMC3::GenVertexPtr vertex,
                  const std::vector<std::string>& processFilter) {
  if (!vertex) return;
  // Store the particles associated to the vertex
  HepMC3::GenParticlePtr particleIn = vertex->particles_in()[0];
  HepMC3::GenParticlePtr particleOut = vertex->particles_out()[0];

  // Walk backwards to find all vertices in a row that match the pattern
  while (findAttribute(particleIn->production_vertex(), processFilter)) {
    // Take the particles before the vertex and remove particles & vertices in
    // between
    HepMC3::GenVertexPtr currentVertex = particleOut->production_vertex();
    HepMC3::GenParticlePtr nextParticle = currentVertex->particles_in()[0];
    // Cut connections from particle and remove it
    if (particleIn->end_vertex()) {
      particleIn->end_vertex()->remove_particle_in(particleIn);
    }
    currentVertex->remove_particle_out(particleIn);
    event.remove_particle(particleIn);
    // Cut connections from vertext and remove it
    particleIn = nextParticle;
    currentVertex->remove_particle_in(particleIn);
    event.remove_vertex(currentVertex);
  }
  // Walk forwards to find all vertices in a row that match the pattern
  while (findAttribute(particleOut->end_vertex(), processFilter)) {
    // Take the particles after the vertex and remove particles & vertices in
    // between
    HepMC3::GenVertexPtr currentVertex = particleOut->end_vertex();
    HepMC3::GenParticlePtr nextParticle = currentVertex->particles_out()[0];
    // Cut connections from particle and remove it
    if (particleOut->production_vertex()) {
      particleOut->production_vertex()->remove_particle_out(particleOut);
    }
    currentVertex->remove_particle_in(particleOut);
    event.remove_particle(particleOut);
    // Cut connections from vertext and remove it
    particleOut = nextParticle;
    currentVertex->remove_particle_out(particleOut);
    event.remove_vertex(currentVertex);
  }

  // Build the dummy vertex
  auto reducedVertex = std::make_shared<HepMC3::GenVertex>();
  event.add_vertex(reducedVertex);

  reducedVertex->add_particle_in(particleIn);
  reducedVertex->add_particle_out(particleOut);

  // Remove and replace the old vertex
  event.remove_vertex(vertex);
  vertex = reducedVertex;
}

/// @brief This method walks over all vertices and tests whether it should be
/// filtered
///
/// @param [in, out] event The underlying event
/// @param [in, out] The current vertex under investigation
/// @param [in] processFilter List of processes that will be filtered
void followOutgoingParticles(HepMC3::GenEvent& event,
                             HepMC3::GenVertexPtr vertex,
                             const std::vector<std::string>& processFilter) {
  if (!vertex) return;
  // Replace and reduce vertex if it should be filtered
  if (findAttribute(vertex, processFilter)) {
    reduceVertex(event, vertex, processFilter);
  }
  // Move forward to the next vertices
  for (const auto& particle : vertex->particles_out()) {
    followOutgoingParticles(event, particle->end_vertex(), processFilter);
  }
}
  
HepMC3::GenParticlePtr findLastParticleTrivialAncestor( HepMC3::GenParticlePtr particle,const std::vector<std::string>& processFilter) {
  auto endvtx=particle->end_vertex();
  if (!endvtx) return particle;
  if (endvtx->particles_in().size()!=1) return particle;
  if (endvtx->particles_out().size()!=1) return particle;
  if (!findAttribute(endvtx,processFilter)) return particle;
  return  findLastParticleTrivialAncestor(endvtx->particles_out()[0], processFilter);
}

void followOutgoingParticlesTrivialAncestor(HepMC3::GenEvent& event, const std::vector<std::string>& processFilter) {
   bool updated = false;    
   for (;;){
     //To make things simpler, we update not more than one particle per pass.
     //The 'updated' flag signals if any particle was updated.
     updated = false;
     for (auto p: event.particles()) {
       auto plast=findLastParticleTrivialAncestor(p,processFilter);  
       if (p==plast) continue;
       updated = true;
       auto vertex_to_remove=p->end_vertex();
       auto vertex_to_add=plast->end_vertex();   
       if (vertex_to_add) { 
         vertex_to_add->add_particle_in(p);
         vertex_to_add->remove_particle_in(plast);
       }
       if (vertex_to_remove) vertex_to_remove->remove_particle_in(p);
       event.remove_vertex(vertex_to_remove);
       break;
   }
   if (!updated) break;//If there are no particles to update, break
   }
}
  
  
}  // namespace

ActsExamples::EventAction* ActsExamples::EventAction::s_instance = nullptr;

ActsExamples::EventAction* ActsExamples::EventAction::instance() {
  // Static acces function via G4RunManager
  return s_instance;
}

ActsExamples::EventAction::EventAction(std::vector<std::string> processFilter)
    : G4UserEventAction(), m_processFilter(std::move(processFilter)) {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
  m_use_trivial_ancestor=true;//one can add setters/getters
}

ActsExamples::EventAction::~EventAction() {
  s_instance = nullptr;
}

void ActsExamples::EventAction::BeginOfEventAction(const G4Event* e) {
  SteppingAction::instance()->clear();
  m_event = HepMC3::GenEvent(HepMC3::Units::GEV, HepMC3::Units::MM);
  m_event.add_beam_particle(std::make_shared<HepMC3::GenParticle>());
  m_event.set_event_number(e->GetEventID());
}

void ActsExamples::EventAction::EndOfEventAction(const G4Event*) {
  // Fast exit if the event is empty
  if (m_event.vertices().empty()) {
    return;
  }
  std::vector<HepMC3::GenVertexPtr> todelete;
  for (HepMC3::GenVertexPtr v: m_event.vertices()) if (v->particles_out().empty()) todelete.push_back(v);
  for (auto v: todelete) {
    for (auto p: v->particles_in()) p->set_status(1);
    m_event.remove_vertex(v);
  }
  if (m_event.vertices().empty()) {
    return;
  }
  save_event(&m_event,"",__LINE__);
  // Filter irrelevant processes
  if (m_use_trivial_ancestor) {
    followOutgoingParticlesTrivialAncestor(m_event, m_processFilter);
  } else {
    auto currentVertex = m_event.vertices()[0];
    followOutgoingParticles(m_event, currentVertex, m_processFilter);
  }
  save_event(&m_event,"",__LINE__);
}

void ActsExamples::EventAction::clear() {
  SteppingAction::instance()->clear();
}

HepMC3::GenEvent& ActsExamples::EventAction::event() {
  return m_event;
}
