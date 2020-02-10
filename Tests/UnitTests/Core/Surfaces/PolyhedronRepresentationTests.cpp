// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/PolyhedronRepresentation.hpp"

// Cone surface
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"

// Disc Surface
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
// Plane Surface
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "Acts/Utilities/ObjHelper.hpp"
#include "Acts/Utilities/Units.hpp"

#include <fstream>

namespace Acts {

using namespace UnitLiterals;

using IdentifiedPolyderon = std::pair<std::string, PolyhedronRepresentation>;

/// Helper method to be called from sub tests
void writeObj(const std::vector<IdentifiedPolyderon>& iphs) {
  for (const auto& iph : iphs) {
    std::ofstream ostream;
    ostream.open(iph.first + ".obj");
    ObjHelper objH;
    iph.second.draw(objH);
    objH.write(ostream);
    ostream.close();
  }
}

std::vector<std::pair<std::string, unsigned int>> testModes = {
    {"", 72}, {"Extremas", 1}, {"TriangleMesh", 72}};

auto transform = std::make_shared<Transform3D>(Transform3D::Identity());

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit tests for Cone Surfaces
BOOST_AUTO_TEST_CASE(ConeSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  for (const auto& mode : testModes) {
    // For ConeSurface types the polyhedron representation is
    /// always a triangular mesh
    bool triangulate = (mode.first == "TriangleMesh");
    if (triangulate) {
      continue;
    }

    auto cone = std::make_shared<ConeBounds>(0.235, 0_mm, 100_mm);
    auto oneCone = Surface::makeShared<ConeSurface>(transform, cone);
    auto oneConePh = oneCone->polyhedronRepresentation(tgContext, mode.second);
    testTypes.push_back({"ConeOneFull" + mode.first, oneConePh});
  }
  writeObj(testTypes);
}

/// Unit tests for Disc Surfaces
BOOST_AUTO_TEST_CASE(DiscSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  for (const auto& mode : testModes) {
    bool triangulate = (mode.first == "TriangleMesh");

    // Full disc
    auto disc = std::make_shared<RadialBounds>(0_mm, 25_mm);
    auto fullDisc = Surface::makeShared<DiscSurface>(transform, disc);
    auto fullPh =
        fullDisc->polyhedronRepresentation(tgContext, mode.second, triangulate);
    testTypes.push_back({"DiscFull" + mode.first, fullPh});
    /*
        // Ring disc
        auto radial = std::make_shared<RadialBounds>(10_mm, 25_mm);
        auto radialDisc = Surface::makeShared<DiscSurface>(transform, radial);
        auto radialPh =
            radialDisc->polyhedronRepresentation(tgContext, mode.second,
       triangulate); testTypes.push_back({"DiscRing" + mode.first, radialPh});

        // Sectoral disc - around 0.
        auto sector = std::make_shared<RadialBounds>(10_mm, 25_mm, 1.25);
        auto sectorDisc = Surface::makeShared<DiscSurface>(transform, sector);
        auto sectorPh =
            sectorDisc->polyhedronRepresentation(tgContext, mode.second,
       triangulate); testTypes.push_back({"DiscSectorCentered" + mode.first,
       sectorPh});

        // Sectoral disc - shifted
        auto sectorShifted =
            std::make_shared<RadialBounds>(10_mm, 25_mm, 0.25, 0.25);
        auto sectorDiscShifted =
            Surface::makeShared<DiscSurface>(transform, sectorShifted);
        auto sectorPhShifted =
            sectorDiscShifted->polyhedronRepresentation(tgContext, mode.second,
       triangulate); testTypes.push_back({"DiscSectorShifted" + mode.first,
       sectorPhShifted});
    */
  }
  writeObj(testTypes);
}

/// Unit tests for Plane Surfaces
BOOST_AUTO_TEST_CASE(PlaneSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  for (const auto& mode : testModes) {
    bool triangulate = (mode.first == "TriangleMesh");

    /// Rectangular Plane
    auto rectangular = std::make_shared<RectangleBounds>(10_mm, 25_mm);
    auto rectangularPlane =
        Surface::makeShared<PlaneSurface>(transform, rectangular);
    auto rectangularPh = rectangularPlane->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    BOOST_CHECK(rectangularPh.vertices.size() == 4);
    if (not triangulate) {
      BOOST_CHECK(rectangularPh.faces.size() == 1);
      std::vector<size_t> expectedRect = {0, 1, 2, 3};
      BOOST_CHECK(rectangularPh.faces[0] == expectedRect);
    } else {
      BOOST_CHECK(rectangularPh.faces.size() == 2);
    }
    testTypes.push_back({"PlaneRectangle" + mode.first, rectangularPh});

    /// Trapezoidal Plane
    auto trapezoid = std::make_shared<TrapezoidBounds>(10_mm, 25_mm, 35_mm);
    auto trapezoidalPlane =
        Surface::makeShared<PlaneSurface>(transform, trapezoid);
    auto trapezoidalPh = trapezoidalPlane->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    BOOST_CHECK(trapezoidalPh.vertices.size() == 4);
    if (not triangulate) {
      BOOST_CHECK(trapezoidalPh.faces.size() == 1);
      std::vector<size_t> expectedTra = {0, 1, 2, 3};
      BOOST_CHECK(trapezoidalPh.faces[0] == expectedTra);
    } else {
      BOOST_CHECK(trapezoidalPh.faces.size() == 2);
    }
    testTypes.push_back({"PlaneTrapezoid" + mode.first, trapezoidalPh});

    /// Full ellispoidal Plane
    /**
    auto ellipse = std::make_shared<EllipseBounds>(10_mm,20_mm,30_mm,40_mm);
    auto ellipsoidPlane = Surface::makeShared<PlaneSurface>(transform,ellipse);
    auto ellispoidPh = ellipsoidPlane->polyhedronRepresentation(tgContext,72);
    BOOST_CHECK(ellispoidPh.vertices.size()==72);
    testTypes.push_back({"PlaneFullEllipse",ellispoidPh});
    */

    /// Sectoral epplise

    /// Diamond shaped plane
    auto diamond =
        std::make_shared<DiamondBounds>(10_mm, 20_mm, 15_mm, 40_mm, 50_mm);
    auto diamondPlane = Surface::makeShared<PlaneSurface>(transform, diamond);
    auto diamondPh = diamondPlane->polyhedronRepresentation(
        tgContext, mode.second, triangulate);
    BOOST_CHECK(diamondPh.vertices.size() == 6);
    if (not triangulate) {
      BOOST_CHECK(diamondPh.faces.size() == 1);
    } else {
      BOOST_CHECK(diamondPh.faces.size() == 4);
    }
    testTypes.push_back({"PlaneDiamond" + mode.first, diamondPh});
  }
  writeObj(testTypes);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts