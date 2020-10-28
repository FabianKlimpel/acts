// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/detail/MomentumParametrisations.hpp"

#include <limits>
#include <TMath.h>
#include <Eigen/Eigenvalues> 

namespace {

/// @brief Evaluate the location in a standard normal distribution for a value from a probability distribution
///
/// @param [in] histo The probability distribution
/// @param [in] mom The abscissa value in @p histo
/// 
/// @return The location in a standard normal distribution
float
gaussianValue(TH1F const* histo, const float mom)
{
	// Get the cumulative probability distribution
	TH1F* normalised = (TH1F*) histo->DrawNormalized();
	TH1F* cumulative = (TH1F*) normalised->GetCumulative();
	// Find the cumulative probability
	const float binContent = cumulative->GetBinContent(cumulative->FindBin(mom));
	// Transform the probability to an entry in a standard normal distribution
	const float value = TMath::ErfInverse(2. * binContent - 1.);
	
	delete(normalised);
	delete(cumulative);
	return value;
}

/// @brief This method transforms a probability distribution into components that represent the cumulative probability distribution
///
/// @note This method is used to avoid Root dependencies afterwards by decomposing the histogram
/// @param [in] hist The probability distribution
///
/// @return Pair containing the bin borders and the bin content
MomentumParametrisations::CumulativeDistribution 
buildMap(TH1F const* hist)
{
	// Retrieve the number of bins & borders
  const int nBins = hist->GetNbinsX();
  std::vector<float> histoBorders(nBins+1);
  std::vector<uint32_t> histoContents(nBins);
  
  // Fill the cumulative histogram
  float integral = 0.;
  std::vector<double> temp_HistoContents(nBins);  
  int iBin;
  for (iBin = 0; iBin < nBins; iBin++){
    float binval = hist->GetBinContent(iBin + 1);
    // Avoid negative bin values
    if(binval < 0) {
      binval = 0.;
    }
    // Store the value
    integral += binval;
    temp_HistoContents[iBin] = integral;
  }
  
  // Ensure that content is available
  if(integral == 0.) {
    histoBorders.clear();
    histoContents.clear();
    return std::make_pair(histoBorders, histoContents);
  }

  // Set the bin borders
  for (iBin = 1; iBin <= nBins; iBin++) 
	histoBorders[iBin - 1]=hist->GetXaxis()->GetBinLowEdge(iBin);
  histoBorders[nBins]=hist->GetXaxis()->GetXmax();
  
  // Set the bin content
  const float invIntegral = 1. / integral;
  for(iBin = 0; iBin < nBins; ++iBin) {
    histoContents[iBin] = MomentumParametrisations::s_MaxValue * (temp_HistoContents[iBin] * invIntegral);
  }
  
  return std::make_pair(histoBorders, histoContents);
}
}

namespace {
float
invariantMass(const ActsExamples::SimParticle::Vector4& fourVector1, const ActsExamples::SimParticle::Vector4& fourVector2)
{
	ActsExamples::SimParticle::Vector4 sum = fourVector1 + fourVector2;
	const ActsExamples::SimParticle::Scalar energy = sum[Acts::eEnergy];
	ActsExamples::SimParticle::Vector3 momentum = sum.template segment<3>(Acts::eMom0);
	return std::sqrt(energy * energy - momentum.norm());
}
}

MomentumParametrisations::EventProperties
MomentumParametrisations::prepateMomenta(const MomentumParametrisations::EventCollection& events, unsigned int multiplicity, bool soft) // TODO: build enum instead of bool
{
	MomentumParametrisations::EventProperties result;
	// Loop over all events
	for(const ActsExamples::EventFraction& event : events)
	{
		// Test the multiplicity and type of the event
		if(event.multiplicity == multiplicity && event.soft == soft)
		{
			const float initialMomentum = event.initialParticle.absMomentum();
			float sum = 0.;
			std::vector<float> momenta;
			momenta.reserve(multiplicity);
			// Fill the vector with the scaled momenta
			for(const ActsExamples::SimParticle& p : event.finalParticles)
			{
				sum += p.absMomentum();
				momenta.push_back(p.absMomentum() / initialMomentum);
			}
			// Add the scaled sum of momenta
			momenta.push_back(sum / initialMomentum);
			result.push_back(std::move(momenta));
		}
	}
	return result;
}

MomentumParametrisations::ProbabilityDistributions
MomentumParametrisations::buildMomPerMult(const MomentumParametrisations::EventProperties& events, unsigned int nBins)
{
	// Fast exit
	if(events.empty())
		return {};
	const unsigned int multMax = events[0].size();
	
	// Find the range of each histogram
	std::vector<float> min(multMax, std::numeric_limits<float>::max());
	std::vector<float> max(multMax, 0);
	for(const std::vector<float>& event : events)
		for(unsigned int i = 0; i < multMax; i++)
		{
			min[i] = std::min(event[i], min[i]);
			max[i] = std::max(event[i], max[i]);
		}
	
	// Evaluate the range of the histograms
	// This is used to avoid entries in over-/underflow bins
	std::vector<float> diff(multMax);
	for(unsigned int i = 0; i < multMax; i++)
		diff[i] = (max[i] - min[i]) * 0.1;

	// Build the histograms
	MomentumParametrisations::ProbabilityDistributions histos(multMax);
	for(unsigned int i = 0; i < multMax; i++)
	{
		histos[i] = new TH1F("", "", nBins, min[i] - diff[i], max[i] + diff[i]);
	}

	// Fill the histograms
	for(const std::vector<float>& event : events)
	{
		for(unsigned int i = 0; i < multMax; i++)
		{
			histos[i]->Fill(event[i]);
		}
	}
	return histos;
}

MomentumParametrisations::EventProperties
MomentumParametrisations::convertEventToGaussian(const MomentumParametrisations::ProbabilityDistributions& histos, const MomentumParametrisations::EventProperties& events)
{
	// Fast exit
	if(events.empty())
		return {};
	const unsigned int multMax = events[0].size();
	
	// Loop over the events
	MomentumParametrisations::EventProperties gaussianEvents;
	for(const std::vector<float>& event : events)
	{
		// Transform the properties in the events
		std::vector<float> gaussianEvent;
		for(unsigned int i = 0; i < multMax; i++)
		{
			gaussianEvent.push_back(gaussianValue(histos[i], event[i]));
		}
		// Store the transformed event
		gaussianEvents.push_back(gaussianEvent);
	}
	return gaussianEvents;
}

std::vector<MomentumParametrisations::CumulativeDistribution>
MomentumParametrisations::buildMaps(const MomentumParametrisations::ProbabilityDistributions& histos)
{
	std::vector<MomentumParametrisations::CumulativeDistribution> maps;
	for(auto& h : histos)
	{
		maps.push_back(buildMap(h));
	}
	return maps;
}

MomentumParametrisations::EventProperties
MomentumParametrisations::prepateInvariantMasses(const MomentumParametrisations::EventCollection& events, unsigned int multiplicity, bool soft)
{	
	MomentumParametrisations::EventProperties result;
	// Loop over all events
	for(const ActsExamples::EventFraction& event : events)
	{
		// Test the multiplicity and type of the event
		if(event.multiplicity == multiplicity && event.soft == soft)
		{
			const auto initialFourVector = event.initialParticle.momentum4();
			std::vector<float> invariantMasses;
			invariantMasses.reserve(multiplicity);
			// Fill the vector with the invariant masses
			for(const ActsExamples::SimParticle& p : event.finalParticles)
			{
				const auto fourVector = p.momentum4();
				invariantMasses.push_back(invariantMass(initialFourVector, fourVector));
			}
			result.push_back(invariantMasses);
		}
	}
	return result;
}