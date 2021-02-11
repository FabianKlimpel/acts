// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <stdexcept>
#include <string>

namespace FW {

/// Source link class for simulation in the acts-framework.
///
/// The source link stores the measuremts, surface, and the associated simulated
/// truth hit.
///
/// @todo Allow multiple truth hits e.g. for merged hits.
class EffectiveSourceLink {
 public:
  EffectiveSourceLink(const Acts::GeometryObject& referenceObject, const ActsFatras::Hit& truthHit,
                size_t dim, Acts::FreeVector values, Acts::FreeMatrix cov)
      : m_values(values), m_cov(cov), m_dim(dim),
        m_geometryId(truthHit.geometryId()),
        m_referenceObject(&referenceObject),
        m_truthHit(&truthHit) {m_bound = false;}
  EffectiveSourceLink(const Acts::GeometryObject& referenceObject, const ActsFatras::Hit& truthHit,
                size_t dim, Acts::BoundVector values, Acts::BoundMatrix cov)
      : m_dim(dim),
        m_geometryId(truthHit.geometryId()),
        m_referenceObject(&referenceObject),
        m_truthHit(&truthHit) {m_values.template head<6>() = values;
			m_cov.template topLeftCorner<6, 6>() = cov;
			m_bound = true;}
  
  /// Must be default_constructible to satisfy SourceLinkConcept.
  EffectiveSourceLink() = default;
  EffectiveSourceLink(EffectiveSourceLink&&) = default;
  EffectiveSourceLink(const EffectiveSourceLink&) = default;
  EffectiveSourceLink& operator=(EffectiveSourceLink&&) = default;
  EffectiveSourceLink& operator=(const EffectiveSourceLink&) = default;

  using MeasurementType = std::variant<Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices, Acts::BoundParametersIndices::eLOC_0>,
						  Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,Acts::BoundParametersIndices::eLOC_0, Acts::BoundParametersIndices::eLOC_1>, 
						  Acts::Measurement<EffectiveSourceLink, Acts::FreeParametersIndices, Acts::FreeParametersIndices::eFreePos0, Acts::FreeParametersIndices::eFreePos1, Acts::FreeParametersIndices::eFreePos2>,
						  Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,Acts::BoundParametersIndices::eLOC_0, Acts::BoundParametersIndices::eLOC_1, Acts::BoundParametersIndices::eBoundQOverP>,
						  Acts::Measurement<EffectiveSourceLink, Acts::FreeParametersIndices, Acts::FreeParametersIndices::eFreePos0, Acts::FreeParametersIndices::eFreePos1, 
									Acts::FreeParametersIndices::eFreePos2, Acts::FreeParametersIndices::eFreeQOverP>,
						  Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,Acts::BoundParametersIndices::eLOC_0, Acts::BoundParametersIndices::eLOC_1, Acts::BoundParametersIndices::eBoundPhi, Acts::BoundParametersIndices::eBoundTheta, Acts::BoundParametersIndices::eBoundQOverP>>;
						  
  constexpr Acts::GeometryID geometryId() const { return m_geometryId; }
  constexpr const Acts::GeometryObject& referenceObject() const { return *m_referenceObject; }
  constexpr const ActsFatras::Hit& truthHit() const { return *m_truthHit; }
  MeasurementType operator*() const 
  { 
	  	if (m_dim == 0) {
		  throw std::runtime_error("Cannot create dim 0 measurement");
		} else if (m_dim == 1) {
		  return Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,
								   Acts::ParDef::eLOC_0>(dynamic_cast<const Acts::Surface*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<1, 1>(), m_values[0]);
		} else if (m_dim == 2) {
		  return Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,
								   Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>(
			  dynamic_cast<const Acts::Surface*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<2, 2>(),
			  m_values[0], m_values[1]);
		} else if (m_dim == 5) {
		  return Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,
								   Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1, Acts::ParDef::eBoundPhi, Acts::ParDef::eBoundTheta, Acts::ParDef::eBoundQOverP>(
			  dynamic_cast<const Acts::Surface*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<5, 5>(),
			  m_values[0], m_values[1], m_values[2], m_values[3], m_values[4]);
		} else if (m_dim == 4) {
			Acts::ActsSymMatrixD<4> mat;
			for(unsigned int row = 0; row < 8; row++)
				for(unsigned int col = 0; col < 8; col++)
				{
					if((row >= 3 && row < 7) || (col >= 3 && col < 7))
						continue;
					unsigned int i = (row == 7) ? 4 : row;
					unsigned int j = (col == 7) ? 4 : col;
					mat(i,j) = m_cov(row,col);
				}
		  return Acts::Measurement<EffectiveSourceLink, Acts::FreeParametersIndices,
								    Acts::FreeParametersIndices::eFreePos0, Acts::FreeParametersIndices::eFreePos1, 
									Acts::FreeParametersIndices::eFreePos2, Acts::FreeParametersIndices::eFreeQOverP>(
			  dynamic_cast<const Acts::Volume*>(m_referenceObject)->getSharedPtr(), *this, mat,
			  m_values[0], m_values[1], m_values[2], m_values[7]);
	    } else if (m_dim == 3) {
			/// Looks deprecated
			//~ if(m_bound) {
				//~ Acts::ActsSymMatrixD<3> mat;
				//~ for(unsigned int row = 0; row < 5; row++)
					//~ for(unsigned int col = 0; col < 5; col++)
					//~ {
						//~ if((row == 2 || row == 3) || (col == 2 || col == 3))
							//~ continue;
						//~ unsigned int i = (row == 4) ? 3 : row;
						//~ unsigned int j = (col == 4) ? 3 : col;
						//~ mat(i,j) = m_cov(row,col);
					//~ }
				//~ return Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,Acts::BoundParametersIndices::eLOC_0, 
									//~ Acts::BoundParametersIndices::eLOC_1, Acts::BoundParametersIndices::eBoundQOverP>(
									//~ dynamic_cast<const Acts::Surface*>(m_referenceObject)->getSharedPtr(), *this, mat,
									//~ m_values[0], m_values[1], m_values[4]);
			//~ } else {
			return Acts::Measurement<EffectiveSourceLink, Acts::FreeParametersIndices,
								   Acts::FreeParametersIndices::eFreePos0, Acts::FreeParametersIndices::eFreePos1, Acts::FreeParametersIndices::eFreePos2>(
			  dynamic_cast<const Acts::Volume*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<3, 3>(),
			  m_values[0], m_values[1], m_values[2]);
		  //~ }
		} else {
		  throw std::runtime_error("Dim " + std::to_string(m_dim) +
								   " currently not supported.");
		}
	}
	  
 private:
 bool m_bound;
   Acts::FreeVector m_values;
  Acts::FreeMatrix m_cov;
  size_t m_dim = 0u;
	 Acts::GeometryID m_geometryId;
	 const Acts::GeometryObject* m_referenceObject;
	 const ActsFatras::Hit* m_truthHit;

  friend constexpr bool operator==(const EffectiveSourceLink& lhs,
                                   const EffectiveSourceLink& rhs) {
    return lhs.m_truthHit == rhs.m_truthHit;
  }
};

/// Store source links ordered by geometry identifier.
using EffectiveSourceLinkContainer = GeometryIdMultiset<EffectiveSourceLink>;

}  // end of namespace FW