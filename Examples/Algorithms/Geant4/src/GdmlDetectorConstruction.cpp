// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

#include <G4GDMLParser.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
//~ #include <SystemOfUnits.hh>

using namespace ActsExamples;

GdmlDetectorConstruction::GdmlDetectorConstruction(std::string path)
    : G4VUserDetectorConstruction(), m_path(std::move(path)) {}

G4VPhysicalVolume* GdmlDetectorConstruction::Construct() {
  G4GDMLParser parser;
  // TODO how to handle errors
  parser.Read(m_path,false);
  G4VPhysicalVolume* world = parser.GetWorldVolume();
  G4LogicalVolume* logicalWorld = world->GetLogicalVolume();
  
  G4MagneticField* magneticField = new G4UniformMagField(G4ThreeVector(0.,0.,2.*CLHEP::tesla));
  G4FieldManager* fieldManager = new G4FieldManager(magneticField);
  logicalWorld->SetFieldManager(fieldManager, true); 

  return world;
}
