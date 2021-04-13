// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

using namespace Acts::UnitLiterals;

unsigned int itest = 0;

//~ // datasets for Boost::Test data-driven test cases
//~ namespace ds {
//~ // track parameters
//~ auto pT = bdata::xrange(0.5_GeV, 10_GeV, 500_MeV);
//~ auto phi = bdata::xrange(-180_degree, 180_degree, 30_degree);
//~ auto theta = bdata::xrange(15_degree, 90_degree, 15_degree);
//~ auto charge = bdata::make({1_e, -1_e});
//~ // combined track parameters as cartesian product of all possible combinations
//~ auto trackParameters = (pT * phi * theta * charge);
//~ auto propagationFraction = bdata::xrange(0.0, 1.0, 0.25);
//~ auto propagationLimit = bdata::xrange(10_cm, 1_m, 10_cm);
//~ // additional random numbers
//~ auto rand1 = bdata::random(
    //~ (bdata::seed = 1515,
     //~ bdata::distribution = std::uniform_real_distribution<>(-1., 1.)));
//~ auto rand2 = bdata::random(
    //~ (bdata::seed = 1616,
     //~ bdata::distribution = std::uniform_real_distribution<>(-1., 1.)));
//~ auto rand3 = bdata::random(
    //~ (bdata::seed = 1717,
     //~ bdata::distribution = std::uniform_real_distribution<>(-1., 1.)));
//~ auto threeRandom = (rand1 ^ rand2 ^ rand2);
//~ }  // namespace ds

/**
// The constant field test
/// test forward propagation in constant magnetic field
BOOST_AUTO_TEST_CASE(constant_bfieldforward_propagation_) {

  TFile tf("IntegrationtestsMean.root", "RECREATE");
   
  // Random
  std::default_random_engine gen(42);
  std::normal_distribution<double> dx(0., sqrt(50_mm));
  std::normal_distribution<double> dy(0., sqrt(50_mm));
  std::normal_distribution<double> dz(0., sqrt(100_mm));
 
  double pT = 10_GeV;
  double phi = 0;
  double theta = 0.5 * M_PI;
  
	double px = pT * cos(phi);
	double py = pT * sin(phi);
	double pz = pT / tan(theta);
	double q = -1_e;
	Vector3D mom(px, py, pz);
	
  GeometryContext gctx;
  auto surface = Surface::makeShared<const PlaneSurface>(Vector3D::Zero(), mom.normalized());
for(unsigned int dist = 1; dist < 4; dist++)
{  
  TGraph2D* tgCS = new TGraph2D();
  TGraph2D* tgFS = new TGraph2D();
  TGraph2D* tgCCE = new TGraph2D();
  TGraph2D* tgCFE = new TGraph2D();
  TGraph2D* tgFCE = new TGraph2D();
  TGraph2D* tgFFE = new TGraph2D();
  
  for(unsigned int i = 0; i < 10000; i++)
  {
	// Start parameters
	double x = dx(gen);
	double y = dy(gen);
	BoundVector startVector;
	startVector << x, y, phi, theta, q / mom.norm(), 0;
	BoundParameters startC(gctx, std::nullopt, startVector, surface);
  
    tgCS->SetPoint(i, startC.position().x(), startC.position().y(), startC.position().z());
  
  	Vector3D cc = constant_field_propagation<CurvilinearParameters>(
		  epropagator, startC, pT, phi, theta, Bz, 0.3_m * dist);
	Vector3D cf = constant_field_propagation<FreeTrackParameters>(
		  epropagator, startC, pT, phi, theta, Bz, 0.3_m * dist);

	tgCCE->SetPoint(i, cc.x(), cc.y(), cc.z());
	tgCFE->SetPoint(i, cf.x(), cf.y(), cf.z());
  }
   
  for(unsigned int i = 0; i < 10000; i++)
  {
	  // Start parameters
	  double x = dx(gen);
	  double y = dy(gen);
	  double z = dz(gen);

	  Vector3D dir = mom.normalized();
	  FreeVector pars;
	  pars << x, y, z, 0., dir.x(), dir.y(), dir.z(), q / mom.norm();
	  FreeTrackParameters startF(std::nullopt, pars);

	  tgFS->SetPoint(i, x, y, z);
	  // constant field propagation eigen stepper
	  Vector3D fc = constant_field_propagation<CurvilinearParameters>(
		  epropagator, startF, pT, phi, theta, Bz, 0.3_m * dist);
	  Vector3D ff = constant_field_propagation<FreeTrackParameters>(
		  epropagator, startF, pT, phi, theta, Bz, 0.3_m * dist);
	  tgFCE->SetPoint(i, fc.x(), fc.y(), fc.z());
	  tgFFE->SetPoint(i, ff.x(), ff.y(), ff.z());
  }

	gDirectory->cd();
	gDirectory->WriteObject(tgCS, ("BoundStart" + std::to_string(2_m * dist)).c_str());
	gDirectory->WriteObject(tgFS, ("FreeStart" + std::to_string(2_m * dist)).c_str());
	gDirectory->WriteObject(tgCCE, ("CCE" + std::to_string(2_m * dist)).c_str());
	gDirectory->WriteObject(tgCFE, ("CFE" + std::to_string(2_m * dist)).c_str());
	gDirectory->WriteObject(tgFCE, ("FCE" + std::to_string(2_m * dist)).c_str());
	gDirectory->WriteObject(tgFFE, ("FFE" + std::to_string(2_m * dist)).c_str());
	tf.Write();
	
	delete(tgCS);
	delete(tgFS);
	delete(tgCCE);
	delete(tgCFE);
	delete(tgFCE);
	delete(tgFFE);
}

// Test different distances
  TGraph2D* tgCCE = new TGraph2D();
  TGraph2D* tgCFE = new TGraph2D();
  TGraph2D* tgFCE = new TGraph2D();
  TGraph2D* tgFFE = new TGraph2D();
  TGraph* tgPhi = new TGraph();
  TGraph* tgThe = new TGraph();
  TGraph* tgMom = new TGraph();
  TGraph* tgCharge = new TGraph();
  
  std::uniform_real_distribution ds(0., 5_m);
  std::uniform_real_distribution dphi(-M_PI, M_PI);
  std::uniform_real_distribution dthe(10_degree, 170_degree);
  std::uniform_real_distribution dp(50_MeV, 100_GeV);
  std::uniform_int_distribution<> dq(-1, 1);
   
  for(unsigned int i = 0; i < 100000; i++)
  {
	// Start parameters
	double x = dx(gen);
	double y = dy(gen);
	double z = dz(gen);
	double s = ds(gen);
	
	phi = dphi(gen);
	theta = dthe(gen);
	pT = dp(gen);
	q = dq(gen);
	
	px = pT * cos(phi);
	py = pT * sin(phi);
	pz = pT / tan(theta);
	mom << px, py, pz;

	if(q == 0)
	{
		NeutralCurvilinearTrackParameters startC(std::nullopt, Vector3D(x,y,z), mom, 0);
		
		Vector3D cc = constant_field_propagation<CurvilinearParameters>(
			  epropagator, startC, pT, phi, theta, Bz, s);
		Vector3D cf = constant_field_propagation<FreeTrackParameters>(
			  epropagator, startC, pT, phi, theta, Bz, s);

		tgCCE->SetPoint(i, cc.x(), cc.y(), cc.z());
		tgCFE->SetPoint(i, cf.x(), cf.y(), cf.z()); 		
	} else {
		CurvilinearParameters startC(std::nullopt, Vector3D(x,y,z), mom, q, 0);
		
		Vector3D cc = constant_field_propagation<CurvilinearParameters>(
			  epropagator, startC, pT, phi, theta, Bz, s);
		Vector3D cf = constant_field_propagation<FreeTrackParameters>(
			  epropagator, startC, pT, phi, theta, Bz, s);

		tgCCE->SetPoint(i, cc.x(), cc.y(), cc.z());
		tgCFE->SetPoint(i, cf.x(), cf.y(), cf.z()); 
	}
    
	Vector3D dir = mom.normalized();
	FreeVector pars;
	pars << x, y, z, 0., dir.x(), dir.y(), dir.z(), (q == 0 ? 1 : q) / mom.norm();

	if(q == 0)
	{
		NeutralFreeTrackParameters startF(std::nullopt, pars);

		// constant field propagation eigen stepper
		Vector3D fc = constant_field_propagation<CurvilinearParameters>(
		  epropagator, startF, pT, phi, theta, Bz, s);
		Vector3D ff = constant_field_propagation<FreeTrackParameters>(
		  epropagator, startF, pT, phi, theta, Bz, s);
		tgFCE->SetPoint(i, fc.x(), fc.y(), fc.z());
		tgFFE->SetPoint(i, ff.x(), ff.y(), ff.z());
	} else {
		FreeTrackParameters startF(std::nullopt, pars);

		// constant field propagation eigen stepper
		Vector3D fc = constant_field_propagation<CurvilinearParameters>(
		  epropagator, startF, pT, phi, theta, Bz, s);
		Vector3D ff = constant_field_propagation<FreeTrackParameters>(
		  epropagator, startF, pT, phi, theta, Bz, s);
		tgFCE->SetPoint(i, fc.x(), fc.y(), fc.z());
		tgFFE->SetPoint(i, ff.x(), ff.y(), ff.z());
	}
	tgPhi->SetPoint(i, s, phi);
	tgThe->SetPoint(i, s, theta);
	tgCharge->SetPoint(i, s, q);
	tgMom->SetPoint(i, s, mom.norm());
  }

	gDirectory->cd();
	gDirectory->WriteObject(tgCCE, "CCpos");
	gDirectory->WriteObject(tgCFE, "CFpos");
	gDirectory->WriteObject(tgFCE, "FCpos");
	gDirectory->WriteObject(tgFFE, "FFpos");
	gDirectory->WriteObject(tgPhi, "phi");
	gDirectory->WriteObject(tgThe, "theta");
	gDirectory->WriteObject(tgMom, "momentum");
	gDirectory->WriteObject(tgCharge, "charge");
	tf.Write();
	
tf.Close();
}


/// test correct covariance transport for curvilinear parameters
/// this test only works within the
/// s_curvilinearProjTolerance (in: Definitions.hpp)
BOOST_AUTO_TEST_CASE(covariance_transport_to_curvilinear) {
  
	TFile tf("IntegrationtestsCov.root", "RECREATE");
	
	std::default_random_engine gen(42);
	std::normal_distribution<double> dx(0., sqrt(50_mm));
	std::normal_distribution<double> dy(0., sqrt(50_mm));
	std::normal_distribution<double> dz(0., sqrt(100_mm));
	std::uniform_real_distribution ds(0., 5_m);
	std::uniform_real_distribution dphi(-M_PI, M_PI);
	std::uniform_real_distribution dthe(10_degree, 170_degree);
	std::uniform_real_distribution dp(50_MeV, 100_GeV);
	std::uniform_int_distribution<> dq(-1, 1);
  
    TGraph2D* start1 = new TGraph2D();
    TGraph* start2 = new TGraph();
    TGraph2D* end1 = new TGraph2D();
    TGraph2D* end2 = new TGraph2D();
    TGraph2D* rend1 = new TGraph2D();
    TGraph2D* rend2 = new TGraph2D();
    
	// Sample
	for(unsigned int i = 0; i < 100000; i++)
	{		
		/// Random parameters
		double x = dx(gen);
		double y = dy(gen);
		double z = dz(gen);
		double phi = dphi(gen);
		double theta = dthe(gen);
		double pT = dp(gen);
		double s = ds(gen);
		double q = dq(gen);

		start1->SetPoint(i, phi, theta, pT);
		start2->SetPoint(i, s, q);

		// define start parameters
		double px = pT * cos(phi);
		double py = pT * sin(phi);
		double pz = pT / tan(theta);
		double time = 0.;
		Vector3D pos(x, y, z);
		Vector3D mom(px, py, pz);
		Vector3D dir = mom.normalized();
		
		std::optional<Covariance> covOpt = std::nullopt;
		Covariance cov;
		cov <<
		10_mm, 0, 0.123, 0, 0.5, 0,
		0, 10_mm, 0, 0.162, 0, 0,
		0.123, 0, 0.1, 0, 0, 0,
		0, 0.162, 0, 0.1, 0, 0,
		0.5, 0, 0, 0, 1_e / 10_GeV, 0,
		0, 0, 0, 0, 0, 1_ns;
		// clang-format on
		covOpt = cov;
  
		if(q == 0)
		{
			NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, 0);
			auto covCurvSN = covariance_curvilinear<CurvilinearParameters>(spropagator, startN, s);
			auto covCurvRN = covariance_curvilinear<CurvilinearParameters>(rspropagator, startN, s);
			
			end1->SetPoint(i, covCurvSN(0,0), covCurvSN(1,1), covCurvSN(2,2));
			end2->SetPoint(i, covCurvSN(3,3), covCurvSN(4,4), covCurvSN.determinant());
			rend1->SetPoint(i, covCurvRN(0,0), covCurvRN(1,1), covCurvRN(2,2));
			rend2->SetPoint(i, covCurvRN(3,3), covCurvRN(4,4), covCurvRN.determinant());
		} else {
			CurvilinearParameters startC(covOpt, pos, mom, q, 0);
			auto covCurvEC = covariance_curvilinear<CurvilinearParameters>(epropagator, startC, s);
			auto covCurvRC = covariance_curvilinear<CurvilinearParameters>(repropagator, startC, s);

			end1->SetPoint(i, covCurvEC(0,0), covCurvEC(1,1), covCurvEC(2,2));
			end2->SetPoint(i, covCurvEC(3,3), covCurvEC(4,4), covCurvEC.determinant());
			rend1->SetPoint(i, covCurvRC(0,0), covCurvRC(1,1), covCurvRC(2,2));
			rend2->SetPoint(i, covCurvRC(3,3), covCurvRC(4,4), covCurvRC.determinant());
		}
	}
	
	gDirectory->cd();
	gDirectory->WriteObject(start1, "StartGraph1");
	gDirectory->WriteObject(start2, "StartGraph2");
	gDirectory->WriteObject(end1, "EndGraph1");
	gDirectory->WriteObject(end2, "EndGraph2");
	gDirectory->WriteObject(rend1, "REndGraph1");
	gDirectory->WriteObject(rend2, "REndGraph2");
	
	tf.Write();
	tf.Close();
}


BOOST_AUTO_TEST_CASE(covariance_transport_to_curvilinear) {
  
	TFile tf("IntegrationtestsCov.root", "RECREATE");
	
	std::default_random_engine gen(42);
	std::normal_distribution<double> dx(0., sqrt(50_mm));
	std::normal_distribution<double> dy(0., sqrt(50_mm));
	std::normal_distribution<double> dz(0., sqrt(100_mm));
	std::uniform_real_distribution ds(0., 5_m);
	std::uniform_real_distribution dphi(-M_PI, M_PI);
	std::uniform_real_distribution dthe(10_degree, 170_degree);
	std::uniform_real_distribution dp(50_MeV, 100_GeV);
	std::uniform_int_distribution<> dq(-1, 1);
  
    TGraph2D* start1 = new TGraph2D();
    TGraph* start2 = new TGraph();
    TGraph2D* end1 = new TGraph2D();
    TGraph* end2 = new TGraph();
    TGraph* end3 = new TGraph();
    TGraph2D* rend1 = new TGraph2D();
    TGraph* rend2 = new TGraph();
    TGraph* rend3 = new TGraph();
    
	// Sample
	for(unsigned int i = 0; i < 100000; i++)
	{		
		/// Random parameters
		double x = dx(gen);
		double y = dy(gen);
		double z = dz(gen);
		double phi = dphi(gen);
		double theta = dthe(gen);
		double pT = dp(gen);
		double s = ds(gen);
		double q = dq(gen);

		//~ std::vector<float> startParams{x, y, z, phi, theta, pT, s, q};
		start1->SetPoint(i, phi, theta, pT);
		start2->SetPoint(i, s, q);

		// define start parameters
		double px = pT * cos(phi);
		double py = pT * sin(phi);
		double pz = pT / tan(theta);
		double time = 0.;
		Vector3D pos(x, y, z);
		Vector3D mom(px, py, pz);
		Vector3D dir = mom.normalized();
		
		std::optional<Covariance> covOpt = std::nullopt;
		Covariance cov;
		cov <<
		10_mm, 0, 0.123, 0, 0.5, 0,
		0, 10_mm, 0, 0.162, 0, 0,
		0.123, 0, 0.1, 0, 0, 0,
		0, 0.162, 0, 0.1, 0, 0,
		0.5, 0, 0, 0, 1_e / 10_GeV, 0,
		0, 0, 0, 0, 0, 1_ns;
		// clang-format on
		covOpt = cov;
  
		if(q == 0)
		{
			NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, 0);
			auto covCurvSN = covariance_curvilinear<CurvilinearParameters>(spropagator, startN, s);
			auto covCurvRN = covariance_curvilinear<CurvilinearParameters>(rspropagator, startN, s);
			
			end1->SetPoint(i, covCurvSN(0,0), covCurvSN(1,1), covCurvSN(2,2));
			end2->SetPoint(i, covCurvSN(3,3), covCurvSN(4,4));
			end3->SetPoint(i, covCurvSN(5,5), covCurvSN.determinant());
			rend1->SetPoint(i, covCurvRN(0,0), covCurvRN(1,1), covCurvRN(2,2));
			rend2->SetPoint(i, covCurvRN(3,3), covCurvRN(4,4));
			rend3->SetPoint(i, covCurvRN(5,5), covCurvRN.determinant());
		} else {
			CurvilinearParameters startC(covOpt, pos, mom, q, 0);
			auto covCurvEC = covariance_curvilinear<CurvilinearParameters>(epropagator, startC, s);
			auto covCurvRC = covariance_curvilinear<CurvilinearParameters>(repropagator, startC, s);

			end1->SetPoint(i, covCurvEC(0,0), covCurvEC(1,1), covCurvEC(2,2));
			end2->SetPoint(i, covCurvEC(3,3), covCurvEC(4,4));
			end3->SetPoint(i, covCurvEC(5,5), covCurvEC.determinant());
			rend1->SetPoint(i, covCurvRC(0,0), covCurvRC(1,1), covCurvRC(2,2));
			rend2->SetPoint(i, covCurvRC(3,3), covCurvRC(4,4));
			rend3->SetPoint(i, covCurvRC(5,5), covCurvRC.determinant());
		}
	}
	
	gDirectory->cd();
	gDirectory->WriteObject(start1, "StartGraph1");
	gDirectory->WriteObject(start2, "StartGraph2");
	gDirectory->WriteObject(end1, "EndGraph1");
	gDirectory->WriteObject(end2, "EndGraph2");
	gDirectory->WriteObject(end3, "EndGraph3");
	gDirectory->WriteObject(rend1, "REndGraph1");
	gDirectory->WriteObject(rend2, "REndGraph2");
	gDirectory->WriteObject(rend3, "REndGraph3");
	
	tf.Write();
	tf.Close();
}
*/
BOOST_AUTO_TEST_CASE(dense_covariance_transport) {

	TFile tf("IntegrationtestsCovDense.root", "RECREATE");
	
	std::default_random_engine gen(42);
	std::normal_distribution<double> dx(0., sqrt(50_mm));
	std::normal_distribution<double> dy(0., sqrt(50_mm));
	std::normal_distribution<double> dz(0., sqrt(100_mm));
	std::uniform_real_distribution ds(0., 5_m);
	std::uniform_real_distribution dphi(-M_PI, M_PI);
	std::uniform_real_distribution dthe(10_degree, 170_degree);
	std::uniform_real_distribution dp(50_MeV, 100_GeV);
	std::uniform_int_distribution<> dq(-1, 1);
  
    TGraph2D* start1 = new TGraph2D();
    TGraph* start2 = new TGraph();
    TGraph2D* end1 = new TGraph2D();
    TGraph* end2 = new TGraph();
    TGraph* end3 = new TGraph();
    TGraph2D* rend1 = new TGraph2D();
    TGraph* rend2 = new TGraph();
    TGraph* rend3 = new TGraph();
    
    DensePropagatorType dpropagator = setupDensePropagator();
	RiddersPropagator rdpropagator(dpropagator);
  
    // Sample
	for(unsigned int i = 0; i < 100000; i++)
	{		
std::cout << "Running iteration " << i << std::endl;
		/// Random parameters
		double x = dx(gen);
		double y = dy(gen);
		double z = dz(gen);
		double phi = dphi(gen);
		double theta = dthe(gen);
		double pT = dp(gen);
		double s = ds(gen);
		double q = dq(gen);

		start1->SetPoint(i, phi, theta, pT);
		start2->SetPoint(i, s, q);

		// define start parameters
		double px = pT * cos(phi);
		double py = pT * sin(phi);
		double pz = pT / tan(theta);
		double time = 0.;
		Vector3D pos(x, y, z);
		Vector3D mom(px, py, pz);
		Vector3D dir = mom.normalized();
		
		std::optional<Covariance> covOpt = std::nullopt;
		Covariance cov;
		cov <<
		10_mm, 0, 0.123, 0, 0.5, 0,
		0, 10_mm, 0, 0.162, 0, 0,
		0.123, 0, 0.1, 0, 0, 0,
		0, 0.162, 0, 0.1, 0, 0,
		0.5, 0, 0, 0, 1_e / 10_GeV, 0,
		0, 0, 0, 0, 0, 1_ns;
		// clang-format on
		covOpt = cov;
  
		if(q == 0)
		{
			NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, 0);
			auto covCurvSN = covariance_curvilinear<CurvilinearParameters>(dpropagator, startN, s);
			auto covCurvRN = covariance_curvilinear<CurvilinearParameters>(rdpropagator, startN, s);
			
			end1->SetPoint(i, covCurvSN(0,0), covCurvSN(1,1), covCurvSN(2,2));
			end2->SetPoint(i, covCurvSN(3,3), covCurvSN(4,4));
			end3->SetPoint(i, covCurvSN(5,5), covCurvSN.determinant());
			rend1->SetPoint(i, covCurvRN(0,0), covCurvRN(1,1), covCurvRN(2,2));
			rend2->SetPoint(i, covCurvRN(3,3), covCurvRN(4,4));
			rend3->SetPoint(i, covCurvRN(5,5), covCurvRN.determinant());
		} else {
			CurvilinearParameters startC(covOpt, pos, mom, q, 0);
			auto covCurvEC = covariance_curvilinear<CurvilinearParameters>(dpropagator, startC, s);
			auto covCurvRC = covariance_curvilinear<CurvilinearParameters>(rdpropagator, startC, s);

			end1->SetPoint(i, covCurvEC(0,0), covCurvEC(1,1), covCurvEC(2,2));
			end2->SetPoint(i, covCurvEC(3,3), covCurvEC(4,4));
			end3->SetPoint(i, covCurvEC(5,5), covCurvEC.determinant());
			rend1->SetPoint(i, covCurvRC(0,0), covCurvRC(1,1), covCurvRC(2,2));
			rend2->SetPoint(i, covCurvRC(3,3), covCurvRC(4,4));
			rend3->SetPoint(i, covCurvRC(5,5), covCurvRC.determinant());
		}		
	}
	
	gDirectory->cd();
	gDirectory->WriteObject(start1, "StartGraph1");
	gDirectory->WriteObject(start2, "StartGraph2");
	gDirectory->WriteObject(end1, "EndGraph1");
	gDirectory->WriteObject(end2, "EndGraph2");
	gDirectory->WriteObject(end3, "EndGraph3");
	gDirectory->WriteObject(rend1, "REndGraph1");
	gDirectory->WriteObject(rend2, "REndGraph2");
	gDirectory->WriteObject(rend3, "REndGraph3");
	
	tf.Write();
	tf.Close();



}

//~ /// test correct covariance transport for curvilinear parameters
//~ /// this test only works within the
//~ /// s_curvilinearProjTolerance (in: Definitions.hpp)
//~ BOOST_AUTO_TEST_CASE(covariance_transport_to_curvilinear) {

	//~ TFile tf("IntegrationtestsCov.root", "RECREATE");
	
	//~ std::default_random_engine gen(42);
	//~ std::normal_distribution<double> dx(0., sqrt(50_mm));
	//~ std::normal_distribution<double> dy(0., sqrt(50_mm));
	//~ std::normal_distribution<double> dz(0., sqrt(100_mm));
	//~ std::uniform_real_distribution ds(0., 5_m);
	//~ std::uniform_real_distribution dphi(-M_PI, M_PI);
	//~ std::uniform_real_distribution dthe(10_degree, 170_degree);
	//~ std::uniform_real_distribution dp(50_MeV, 100_GeV);
	//~ std::uniform_int_distribution<> dq(-1, 1);
  
    //~ TGraph2D* start1 = new TGraph2D();
    //~ TGraph* start2 = new TGraph();
    //~ TGraph2D* end1 = new TGraph2D();
    //~ TGraph2D* end2 = new TGraph2D();
    //~ TGraph2D* rend1 = new TGraph2D();
    //~ TGraph2D* rend2 = new TGraph2D();
    
	//~ // Sample
	//~ for(unsigned int i = 0; i < 100000; i++)
	//~ {		
		//~ /// Random parameters
		//~ double x = dx(gen);
		//~ double y = dy(gen);
		//~ double z = dz(gen);
		//~ double phi = dphi(gen);
		//~ double theta = dthe(gen);
		//~ double pT = dp(gen);
		//~ double s = ds(gen);
		//~ double q = dq(gen);
		
		//~ start1->SetPoint(i, phi, theta, pT);
		//~ start2->SetPoint(i, s, q);

		//~ // define start parameters
		//~ double px = pT * cos(phi);
		//~ double py = pT * sin(phi);
		//~ double pz = pT / tan(theta);
		//~ double time = 0.;
		//~ Vector3D pos(x, y, z);
		//~ Vector3D mom(px, py, pz);
		//~ Vector3D dir = mom.normalized();
		
		//~ std::optional<FreeCovariance> covOptFree = std::nullopt;
		//~ FreeCovariance covFree;
		//~ // take some major correlations (off-diagonals)
		//~ // clang-format off
		//~ covFree <<
		 //~ 10_mm, 0, 0.123, 0, 0.01, 0.01, 0.01, 0,
		 //~ 0, 10_mm, 0, 0, 0.01, 0.01, 0.01, 0,
		 //~ 0.123, 0, 10_mm, 0, 0.01, 0.01, 0.1, 0,
		 //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
		 //~ 0.01, 0.01, 0.01, 0, 0.0123, 0, 0, 0,
		 //~ 0.01, 0.01, 0.01, 0, 0, 0.0123, 0, 0,
		 //~ 0.01, 0.01, 0.01, 0, 0, 0, 0.0123, 0,
		 //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
		//~ // clang-format on
		//~ covOptFree = covFree;

		//~ if(q == 0)
		//~ {
			//~ FreeVector parsN;
			//~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
			//~ NeutralFreeTrackParameters startNF(covOptFree, parsN);
			
			//~ auto covCurvSN = covariance_curvilinear<CurvilinearParameters>(spropagator, startNF, s);
			//~ auto covCurvRN = covariance_curvilinear<CurvilinearParameters>(rspropagator, startNF, s);
			
			//~ end1->SetPoint(i, covCurvSN(0,0), covCurvSN(1,1), covCurvSN(2,2));
			//~ end2->SetPoint(i, covCurvSN(3,3), covCurvSN(4,4), covCurvSN.determinant());
			//~ rend1->SetPoint(i, covCurvRN(0,0), covCurvRN(1,1), covCurvRN(2,2));
			//~ rend2->SetPoint(i, covCurvRN(3,3), covCurvRN(4,4), covCurvRN.determinant());
		//~ } else {
			//~ FreeVector parsC;
			//~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
			//~ FreeTrackParameters startCF(covOptFree, parsC);
			//~ auto covCurvEC = covariance_curvilinear<CurvilinearParameters>(epropagator, startCF, s);
			//~ auto covCurvRC = covariance_curvilinear<CurvilinearParameters>(repropagator, startCF, s);

			//~ end1->SetPoint(i, covCurvEC(0,0), covCurvEC(1,1), covCurvEC(2,2));
			//~ end2->SetPoint(i, covCurvEC(3,3), covCurvEC(4,4), covCurvEC.determinant());
			//~ rend1->SetPoint(i, covCurvRC(0,0), covCurvRC(1,1), covCurvRC(2,2));
			//~ rend2->SetPoint(i, covCurvRC(3,3), covCurvRC(4,4), covCurvRC.determinant());
		//~ }
	//~ }
	
	//~ gDirectory->cd();
	//~ gDirectory->WriteObject(start1, "StartGraph1");
	//~ gDirectory->WriteObject(start2, "StartGraph2");
	//~ gDirectory->WriteObject(end1, "EndGraph1");
	//~ gDirectory->WriteObject(end2, "EndGraph2");
	//~ gDirectory->WriteObject(rend1, "REndGraph1");
	//~ gDirectory->WriteObject(rend2, "REndGraph2");
	
	//~ tf.Write();
	//~ tf.Close();
//~ }


//~ BOOST_AUTO_TEST_CASE(covariance_transport_to_free) {

	//~ TFile tf("IntegrationtestsCov2.root", "RECREATE");
	
	//~ std::default_random_engine gen(42);
	//~ std::normal_distribution<double> dx(0., sqrt(50_mm));
	//~ std::normal_distribution<double> dy(0., sqrt(50_mm));
	//~ std::normal_distribution<double> dz(0., sqrt(100_mm));
	//~ std::uniform_real_distribution ds(0., 5_m);
	//~ std::uniform_real_distribution dphi(-M_PI, M_PI);
	//~ std::uniform_real_distribution dthe(10_degree, 170_degree);
	//~ std::uniform_real_distribution dp(50_MeV, 100_GeV);
	//~ std::uniform_int_distribution<> dq(-1, 1);
  
	//~ TGraph2D* start1 = new TGraph2D();
    //~ TGraph* start2 = new TGraph();
    //~ TGraph2D* cend1 = new TGraph2D();
    //~ TGraph2D* cend2 = new TGraph2D();
    //~ TGraph2D* cend3 = new TGraph2D();
    //~ TGraph2D* crend1 = new TGraph2D();
    //~ TGraph2D* crend2 = new TGraph2D();
    //~ TGraph2D* crend3 = new TGraph2D();
    //~ TGraph2D* fend1 = new TGraph2D();
    //~ TGraph2D* fend2 = new TGraph2D();
    //~ TGraph2D* fend3 = new TGraph2D();
    //~ TGraph2D* frend1 = new TGraph2D();
    //~ TGraph2D* frend2 = new TGraph2D();
    //~ TGraph2D* frend3 = new TGraph2D();
    
	//~ // Sample
	//~ for(unsigned int i = 0; i < 100000; i++)
	//~ {		
		//~ /// Random parameters
		//~ double x = dx(gen);
		//~ double y = dy(gen);
		//~ double z = dz(gen);
		//~ double phi = dphi(gen);
		//~ double theta = dthe(gen);
		//~ double pT = dp(gen);
		//~ double s = ds(gen);
		//~ double q = dq(gen);

		//~ start1->SetPoint(i, phi, theta, pT);
		//~ start2->SetPoint(i, s, q);

		//~ // define start parameters
		//~ double px = pT * cos(phi);
		//~ double py = pT * sin(phi);
		//~ double pz = pT / tan(theta);
		//~ double time = 0.;
		//~ Vector3D pos(x, y, z);
		//~ Vector3D mom(px, py, pz);
		//~ Vector3D dir = mom.normalized();
		
		//~ std::optional<Covariance> covOpt = std::nullopt;
		//~ std::optional<FreeCovariance> covOptFree = std::nullopt;
		//~ Covariance cov;
		//~ // take some major correlations (off-diagonals)
		//~ // clang-format off
		//~ cov <<
		 //~ 10_um, 0, 0.123, 0, 0.5, 0,
		 //~ 0, 10_um, 0, 0.162, 0, 0,
		 //~ 0.123, 0, 0.1, 0, 0, 0,
		 //~ 0, 0.162, 0, 0.1, 0, 0,
		 //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
		 //~ 0, 0, 0, 0, 0, 1_us;
		//~ // clang-format on
		//~ covOpt = cov;

		//~ FreeCovariance covFree;
		//~ // take some major correlations (off-diagonals)
		//~ // clang-format off
		//~ covFree <<
		 //~ 10_mm, 0, 0.123, 0, 0.1, 0.1, 0.1, 0,
		 //~ 0, 10_mm, 0, 0, 0.1, 0.1, 0.1, 0,
		 //~ 0.123, 0, 10_mm, 0, 0.1, 0.1, 0.1, 0,
		 //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
		 //~ 0.1, 0.1, 0.1, 0, 0.123, 0, 0, 0,
		 //~ 0.1, 0.1, 0.1, 0, 0, 0.123, 0, 0,
		 //~ 0.1, 0.1, 0.1, 0, 0, 0, 0.123, 0,
		 //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
		//~ // clang-format on
		//~ covOptFree = covFree;

		//~ ///
		//~ /// Curvilinear to Free Tests
		//~ ///
		//~ if(q == 0) 
		//~ {
			//~ NeutralCurvilinearTrackParameters startNC(covOpt, pos, mom, 0);
			
			//~ // covariance check for straight line stepper
			//~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
				//~ spropagator, startNC, s);
			//~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
				//~ rspropagator, startNC, s);
			
			//~ cend1->SetPoint(i, covObtained(0,0), covObtained(1,1), covObtained(2,2));
			//~ cend2->SetPoint(i, covObtained(3,3), covObtained(4,4), covObtained(5,5));
			//~ cend3->SetPoint(i, covObtained(6,6), covObtained(7,7), covObtained.determinant());
			//~ crend1->SetPoint(i, covCalculated(0,0), covCalculated(1,1), covCalculated(2,2));
			//~ crend2->SetPoint(i, covCalculated(3,3), covCalculated(4,4), covCalculated(5,5));
			//~ crend3->SetPoint(i, covCalculated(6,6), covCalculated(7,7), covCalculated.determinant());
		//~ } else {
			//~ CurvilinearParameters startCC(covOpt, pos, mom, q, 0);
				
			//~ // covariance check for eigen stepper
			//~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
				//~ repropagator, startCC, s);
				
			//~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
																	  //~ startCC, s);

			//~ cend1->SetPoint(i, covObtained(0,0), covObtained(1,1), covObtained(2,2));
			//~ cend2->SetPoint(i, covObtained(3,3), covObtained(4,4), covObtained(5,5));
			//~ cend3->SetPoint(i, covObtained(6,6), covObtained(7,7), covObtained.determinant());
			//~ crend1->SetPoint(i, covCalculated(0,0), covCalculated(1,1), covCalculated(2,2));
			//~ crend2->SetPoint(i, covCalculated(3,3), covCalculated(4,4), covCalculated(5,5));
			//~ crend3->SetPoint(i, covCalculated(6,6), covCalculated(7,7), covCalculated.determinant());
		//~ }
		
		//~ ///
		//~ /// Free to Free Tests
		//~ ///
		//~ if(q == 0) {
			//~ FreeVector parsN;
			//~ parsN << x, y, z, 0, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
			//~ NeutralFreeTrackParameters startNF(covOptFree, parsN);
			
			//~ // covariance check for straight line stepper
			//~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
				//~ spropagator, startNF, s);
			//~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
				//~ rspropagator, startNF, s);

			//~ fend1->SetPoint(i, covObtained(0,0), covObtained(1,1), covObtained(2,2));
			//~ fend2->SetPoint(i, covObtained(3,3), covObtained(4,4), covObtained(5,5));
			//~ fend3->SetPoint(i, covObtained(6,6), covObtained(7,7), covObtained.determinant());
			//~ frend1->SetPoint(i, covCalculated(0,0), covCalculated(1,1), covCalculated(2,2));
			//~ frend2->SetPoint(i, covCalculated(3,3), covCalculated(4,4), covCalculated(5,5));
			//~ frend3->SetPoint(i, covCalculated(6,6), covCalculated(7,7), covCalculated.determinant());
		//~ } else {
			//~ FreeVector parsC;
			//~ parsC << x, y, z, 0, dir.x(), dir.y(), dir.z(), q / mom.norm();
			//~ FreeTrackParameters startCF(covOptFree, parsC);
			
			//~ // covariance check for eigen stepper
			//~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
																	  //~ startCF, s);
			//~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
				//~ repropagator, startCF, s);

			//~ fend1->SetPoint(i, covObtained(0,0), covObtained(1,1), covObtained(2,2));
			//~ fend2->SetPoint(i, covObtained(3,3), covObtained(4,4), covObtained(5,5));
			//~ fend3->SetPoint(i, covObtained(6,6), covObtained(7,7), covObtained.determinant());
			//~ frend1->SetPoint(i, covCalculated(0,0), covCalculated(1,1), covCalculated(2,2));
			//~ frend2->SetPoint(i, covCalculated(3,3), covCalculated(4,4), covCalculated(5,5));
			//~ frend3->SetPoint(i, covCalculated(6,6), covCalculated(7,7), covCalculated.determinant());
		//~ }
	//~ }
	
	//~ gDirectory->cd();
	//~ gDirectory->WriteObject(start1, "StartGraph1");
	//~ gDirectory->WriteObject(start2, "StartGraph2");
	//~ gDirectory->WriteObject(cend1, "CEndGraph1");
	//~ gDirectory->WriteObject(cend2, "CEndGraph2");
	//~ gDirectory->WriteObject(cend2, "CEndGraph3");
	//~ gDirectory->WriteObject(crend1, "CREndGraph1");
	//~ gDirectory->WriteObject(crend2, "CREndGraph2");
	//~ gDirectory->WriteObject(crend3, "CREndGraph3");
	//~ gDirectory->WriteObject(fend1, "FEndGraph1");
	//~ gDirectory->WriteObject(fend2, "FEndGraph2");
	//~ gDirectory->WriteObject(fend2, "FEndGraph3");
	//~ gDirectory->WriteObject(frend1, "FREndGraph1");
	//~ gDirectory->WriteObject(frend2, "FREndGraph2");
	//~ gDirectory->WriteObject(frend3, "FREndGraph3");
	
	//~ tf.Write();
	//~ tf.Close();
//~ }

/**
BOOST_DATA_TEST_CASE(covariance_transport_to_free,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  (void)rand3;
  // define start parameters
  double x = 1;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  std::optional<FreeCovariance> covOptFree = std::nullopt;
  if (covtpr) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_um, 0, 0.123, 0, 0.5, 0,
     0, 10_um, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;

    FreeCovariance covFree;
    // take some major correlations (off-diagonals)
    // clang-format off
    covFree <<
     10_mm, 0, 0.123, 0, 0.1, 0.1, 0.1, 0,
     0, 10_mm, 0, 0, 0.1, 0.1, 0.1, 0,
     0.123, 0, 10_mm, 0, 0.1, 0.1, 0.1, 0,
     0, 0, 0, 1_ns, 0, 0, 0, 0,
     0.1, 0.1, 0.1, 0, 0.123, 0, 0, 0,
     0.1, 0.1, 0.1, 0, 0, 0.123, 0, 0,
     0.1, 0.1, 0.1, 0, 0, 0, 0.123, 0,
     0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
    // clang-format on
    covOptFree = covFree;
  }

  ///
  /// Curvilinear to Free Tests
  ///
  {
    CurvilinearParameters startCC(covOpt, pos, mom, q, time);
    NeutralCurvilinearTrackParameters startNC(covOpt, pos, mom, time);

    // covariance check for straight line stepper
    auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        spropagator, startNC, plimit);
    auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        rspropagator, startNC, plimit);
    // Numerical fluctuations in the covariances cause errors in relative
    // comparison. This needs to be tested and avoided by setting both entries
    // to 1
    for (unsigned int i = 0; i < covObtained.rows(); i++)
      for (unsigned int j = 0; j < covObtained.cols(); j++) {
        if (std::abs(covObtained(i, j)) <
                std::numeric_limits<double>::epsilon() ||
            std::abs(covCalculated(i, j)) <
                std::numeric_limits<double>::epsilon()) {
          covObtained(i, j) = 1.;
          covCalculated(i, j) = 1.;
        }
      }
    CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
    // covariance check for eigen stepper
    covCalculated = covariance_curvilinear<FreeTrackParameters>(
        repropagator, startCC, plimit);
        
    covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              startCC, plimit);


    // Numerical fluctuations in the covariances cause errors in relative
    // comparison. This needs to be tested and avoided by setting both entries
    // to 1
    for (unsigned int i = 0; i < covObtained.rows(); i++)
      for (unsigned int j = 0; j < covObtained.cols(); j++) {
        if (std::abs(covObtained(i, j)) <
                std::numeric_limits<double>::epsilon() ||
            std::abs(covCalculated(i, j)) <
                std::numeric_limits<double>::epsilon()) {
          covObtained(i, j) = 1.;
          covCalculated(i, j) = 1.;
        }
      }
    CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
  }

  ///
  /// Free to Free Tests
  ///
  {
    Vector3D dir = mom.normalized();
    FreeVector parsC, parsN;
    parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
    parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
    FreeTrackParameters startCF(covOptFree, parsC);
    NeutralFreeTrackParameters startNF(covOptFree, parsN);

    // covariance check for straight line stepper
    auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        spropagator, startNF, plimit);
    auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        rspropagator, startNF, plimit);
    // Numerical fluctuations in the covariances cause errors in relative
    // comparison. This needs to be tested and avoided by setting both entries
    // to 1
    for (unsigned int i = 0; i < covObtained.rows(); i++)
      for (unsigned int j = 0; j < covObtained.cols(); j++) {
        if (std::abs(covObtained(i, j)) <
                std::numeric_limits<double>::epsilon() ||
            std::abs(covCalculated(i, j)) <
                std::numeric_limits<double>::epsilon()) {
          covObtained(i, j) = 1.;
          covCalculated(i, j) = 1.;
        }
      }
    CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);

    // covariance check for eigen stepper
    covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              startCF, plimit);
    covCalculated = covariance_curvilinear<FreeTrackParameters>(
        repropagator, startCF, plimit);
    // Numerical fluctuations in the covariances cause errors in relative
    // comparison. This needs to be tested and avoided by setting both entries
    // to 1
    for (unsigned int i = 0; i < covObtained.rows(); i++)
      for (unsigned int j = 0; j < covObtained.cols(); j++) {
        if (std::abs(covObtained(i, j)) <
                std::numeric_limits<double>::epsilon() ||
            std::abs(covCalculated(i, j)) <
                std::numeric_limits<double>::epsilon()) {
          covObtained(i, j) = 1.;
          covCalculated(i, j) = 1.;
        }
      }
    CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
  }
}
*/

//~ // The constant field test
//~ /// test forward propagation in constant magnetic field
//~ BOOST_DATA_TEST_CASE(
    //~ constant_bfieldforward_propagation_,
    //~ bdata::random((bdata::seed = 0,
                   //~ bdata::distribution =
                       //~ std::uniform_real_distribution<>(0.4_GeV, 10_GeV))) ^
        //~ bdata::random((bdata::seed = 1,
                       //~ bdata::distribution =
                           //~ std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        //~ bdata::random((bdata::seed = 2,
                       //~ bdata::distribution =
                           //~ std::uniform_real_distribution<>(0.1, M_PI - 0.1))) ^
        //~ bdata::random(
            //~ (bdata::seed = 3,
             //~ bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        //~ bdata::random(
            //~ (bdata::seed = 4,
             //~ bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        //~ bdata::xrange(ntests),
    //~ pT, phi, theta, charge, time, index) {
  //~ if (index < skip) {
    //~ return;
  //~ }

  //~ // define start parameters
  //~ double x = 0;
  //~ double y = 0;
  //~ double z = 0;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = -1 + 2 * charge;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);
  //~ CurvilinearParameters startC(std::nullopt, pos, mom, q, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector pars;
  //~ pars << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ FreeTrackParameters startF(std::nullopt, pars);

  //~ // constant field propagation atlas stepper
  //~ auto aposition = constant_field_propagation<CurvilinearParameters>(
      //~ apropagator, startC, pT, phi, theta, Bz);
  //~ // constant field propagation eigen stepper
  //~ auto epositionCC = constant_field_propagation<CurvilinearParameters>(
      //~ epropagator, startC, pT, phi, theta, Bz);
  //~ auto epositionFC = constant_field_propagation<FreeTrackParameters>(
      //~ epropagator, startC, pT, phi, theta, Bz);
  //~ auto epositionCF = constant_field_propagation<CurvilinearParameters>(
      //~ epropagator, startF, pT, phi, theta, Bz);
  //~ auto epositionFF = constant_field_propagation<FreeTrackParameters>(
      //~ epropagator, startF, pT, phi, theta, Bz);
  //~ // check consistency
  //~ CHECK_CLOSE_REL(epositionCC, aposition, 1e-6);
  //~ CHECK_CLOSE_REL(epositionCC, epositionFC, 1e-6);
  //~ CHECK_CLOSE_REL(epositionCC, epositionCF, 1e-6);
  //~ CHECK_CLOSE_REL(epositionCC, epositionFF, 1e-6);
//~ }

//~ /// test consistency of forward-backward propagation
//~ BOOST_DATA_TEST_CASE(forward_backward_propagation_,
                     //~ ds::trackParameters* ds::propagationLimit, pT, phi, theta,
                     //~ charge, plimit) {
  //~ ++itest;

  //~ // define start parameters
  //~ double x = 0;
  //~ double y = 0;
  //~ double z = 0;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);
  //~ CurvilinearParameters startC(std::nullopt, pos, mom, q, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector pars;
  //~ pars << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ FreeTrackParameters startF(std::nullopt, pars);

  //~ // foward backward check atlas stepper
  //~ foward_backward<CurvilinearParameters, CurvilinearParameters>(
      //~ apropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ // foward backward check eigen stepper
  //~ foward_backward<CurvilinearParameters, CurvilinearParameters>(
      //~ epropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, CurvilinearParameters>(
      //~ epropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<CurvilinearParameters, FreeTrackParameters>(
      //~ epropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, FreeTrackParameters>(
      //~ epropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<CurvilinearParameters, CurvilinearParameters>(
      //~ epropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, CurvilinearParameters>(
      //~ epropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ foward_backward<CurvilinearParameters, FreeTrackParameters>(
      //~ epropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, FreeTrackParameters>(
      //~ epropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ // foward backward check straight line stepper
  //~ foward_backward<CurvilinearParameters, CurvilinearParameters>(
      //~ spropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, CurvilinearParameters>(
      //~ spropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<CurvilinearParameters, FreeTrackParameters>(
      //~ spropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, FreeTrackParameters>(
      //~ spropagator, plimit, startC, 1_um, 1_eV, debug);
  //~ foward_backward<CurvilinearParameters, CurvilinearParameters>(
      //~ spropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, CurvilinearParameters>(
      //~ spropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ foward_backward<CurvilinearParameters, FreeTrackParameters>(
      //~ spropagator, plimit, startF, 1_um, 1_eV, debug);
  //~ foward_backward<FreeTrackParameters, FreeTrackParameters>(
      //~ spropagator, plimit, startF, 1_um, 1_eV, debug);
//~ }

//~ /// test consistency of propagators when approaching a cylinder
//~ BOOST_DATA_TEST_CASE(propagation_to_cylinder_,
                     //~ ds::trackParameters* ds::propagationFraction ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, pfrac, rand1, rand2, rand3) {
  //~ // just make sure we can reach it
  //~ double r = (pfrac > 0. ? pfrac : 0.01) * std::abs(pT / Bz);
  //~ r = (r > 2.5_m) ? 2.5_m : r;

  //~ // define start parameters
  //~ double x = 0;
  //~ double y = 0;
  //~ double z = 0;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;
  //~ }
  //~ CurvilinearParameters startC(covOpt, pos, mom, q, time);
  //~ NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(std::nullopt, parsC);
  //~ NeutralFreeTrackParameters startNF(std::nullopt, parsN);

  //~ // check atlas stepper
  //~ auto a_at_cylinder =
      //~ to_cylinder(apropagator, startC, r, rand1, rand2, rand3, debug);
  //~ // check eigen stepper
  //~ auto e_at_cylinder =
      //~ to_cylinder(epropagator, startC, r, rand1, rand2, rand3, debug);
  //~ auto e_free_at_cylinder =
      //~ to_cylinder(epropagator, startCF, r, rand1, rand2, rand3, debug);
  //~ CHECK_CLOSE_ABS(e_at_cylinder.first, a_at_cylinder.first, 10_um);
  //~ CHECK_CLOSE_ABS(e_at_cylinder.first, e_free_at_cylinder.first, 10_um);

  //~ // check without charge
  //~ auto s_at_cylinder =
      //~ to_cylinder(spropagator, startN, r, rand1, rand2, rand3, debug);
  //~ auto s_free_at_cylinder =
      //~ to_cylinder(spropagator, startNF, r, rand1, rand2, rand3, debug);
  //~ e_at_cylinder =
      //~ to_cylinder(epropagator, startN, r, rand1, rand2, rand3, debug);
  //~ e_free_at_cylinder =
      //~ to_cylinder(epropagator, startNF, r, rand1, rand2, rand3, debug);
  //~ CHECK_CLOSE_ABS(s_at_cylinder.first, s_free_at_cylinder.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_cylinder.first, e_at_cylinder.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_cylinder.first, e_free_at_cylinder.first, 1_um);
//~ }

//~ /// test consistency of propagators to a plane
//~ BOOST_DATA_TEST_CASE(propagation_to_plane_,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ // define start parameters
  //~ double x = 0;
  //~ double y = 0;
  //~ double z = 0;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;
  //~ }
  //~ CurvilinearParameters startC(covOpt, pos, mom, q, time);
  //~ NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(std::nullopt, parsC);
  //~ NeutralFreeTrackParameters startNF(std::nullopt, parsN);

  //~ // to a plane with the atlas stepper
  //~ auto a_at_plane = to_surface<AtlasPropagatorType, PlaneSurface>(
      //~ apropagator, startC, plimit, rand1, rand2, rand3, true);
  //~ // to a plane with the eigen stepper
  //~ auto e_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      //~ epropagator, startC, plimit, rand1, rand2, rand3, true);
  //~ auto e_free_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      //~ epropagator, startCF, plimit, rand1, rand2, rand3, true);
  //~ CHECK_CLOSE_ABS(e_at_plane.first, a_at_plane.first, 1_um);
  //~ CHECK_CLOSE_ABS(e_at_plane.first, e_free_at_plane.first, 1_um);

  //~ // to a plane with the straight line stepper
  //~ auto s_at_plane = to_surface<StraightPropagatorType, PlaneSurface>(
      //~ spropagator, startN, plimit, rand1, rand2, rand3, true);
  //~ auto s_free_at_plane = to_surface<StraightPropagatorType, PlaneSurface>(
      //~ spropagator, startNF, plimit, rand1, rand2, rand3, true);
  //~ // to a plane with the eigen stepper without charge
  //~ e_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      //~ epropagator, startN, plimit, rand1, rand2, rand3, true);
  //~ e_free_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      //~ epropagator, startNF, plimit, rand1, rand2, rand3, true);
  //~ CHECK_CLOSE_ABS(s_at_plane.first, s_free_at_plane.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_plane.first, e_at_plane.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_plane.first, e_free_at_plane.first, 1_um);
//~ }

//~ /// test consistency of propagators to a disc
//~ BOOST_DATA_TEST_CASE(propagation_to_disc_,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ // define start parameters
  //~ double x = 0;
  //~ double y = 0;
  //~ double z = 0;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;
  //~ }
  //~ CurvilinearParameters startC(covOpt, pos, mom, q, time);
  //~ NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(std::nullopt, parsC);
  //~ NeutralFreeTrackParameters startNF(std::nullopt, parsN);

  //~ // to a disc with the  atlas stepper
  //~ auto a_at_disc = to_surface<AtlasPropagatorType, DiscSurface>(
      //~ apropagator, startC, plimit, rand1, rand2, rand3, true);
  //~ // to a disc with the eigen stepper
  //~ auto e_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      //~ epropagator, startC, plimit, rand1, rand2, rand3, true);
  //~ auto e_free_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      //~ epropagator, startCF, plimit, rand1, rand2, rand3, true);
  //~ CHECK_CLOSE_ABS(a_at_disc.first, e_at_disc.first, 1_um);
  //~ CHECK_CLOSE_ABS(a_at_disc.first, e_free_at_disc.first, 1_um);

  //~ // to a disc with the straight line stepper
  //~ auto s_at_disc = to_surface<StraightPropagatorType, DiscSurface>(
      //~ spropagator, startN, plimit, rand1, rand2, rand3, true);
  //~ auto s_free_at_disc = to_surface<StraightPropagatorType, DiscSurface>(
      //~ spropagator, startNF, plimit, rand1, rand2, rand3, true);
  //~ // to a disc with the eigen stepper without charge
  //~ e_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      //~ epropagator, startN, plimit, rand1, rand2, rand3, true);
  //~ e_free_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      //~ epropagator, startNF, plimit, rand1, rand2, rand3, true);
  //~ CHECK_CLOSE_ABS(s_at_disc.first, s_free_at_disc.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_disc.first, e_at_disc.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_disc.first, e_free_at_disc.first, 1_um);
//~ }

//~ /// test consistency of propagators to a line
//~ BOOST_DATA_TEST_CASE(propagation_to_line_,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ std::cout << "Running (first patch) tests : " << itest << std::endl;
  //~ ++itest;

  //~ // define start parameters
  //~ double x = 1.;
  //~ double y = 0.;
  //~ double z = 0.;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;
  //~ }
  //~ CurvilinearParameters startC(covOpt, pos, mom, q, time);
  //~ NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(std::nullopt, parsC);
  //~ NeutralFreeTrackParameters startNF(std::nullopt, parsN);

  //~ // to a line with the atlas stepper
  //~ if (debug) {
    //~ std::cout << "[ >>>> Testing Atlas Propagator <<<< ]" << std::endl;
  //~ }
  //~ auto a_at_line = to_surface<AtlasPropagatorType, StrawSurface>(
      //~ apropagator, startC, plimit, rand1, rand2, rand3, false, debug);
  //~ // to a line with the eigen stepper
  //~ if (debug) {
    //~ std::cout << "[ >>>> Testing Eigen Propagator <<<< ]" << std::endl;
  //~ }
  //~ auto e_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      //~ epropagator, startC, plimit, rand1, rand2, rand3, false, debug);
  //~ auto e_free_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      //~ epropagator, startCF, plimit, rand1, rand2, rand3, false, debug);
  //~ CHECK_CLOSE_ABS(a_at_line.first, e_at_line.first, 10_um);
  //~ CHECK_CLOSE_ABS(a_at_line.first, e_free_at_line.first, 10_um);

  //~ if (debug) {
    //~ std::cout << "[ >>>> Testing Neutral Propagators <<<< ]" << std::endl;
  //~ }
  //~ // to a straw with the straight line stepper
  //~ auto s_at_line = to_surface<StraightPropagatorType, StrawSurface>(
      //~ spropagator, startN, plimit, rand1, rand2, rand3, false, debug);
  //~ auto s_free_at_line = to_surface<StraightPropagatorType, StrawSurface>(
      //~ spropagator, startNF, plimit, rand1, rand2, rand3, false, debug);
  //~ // to a straw with the eigen stepper without charge
  //~ e_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      //~ epropagator, startN, plimit, rand1, rand2, rand3, false, debug);
  //~ e_free_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      //~ epropagator, startNF, plimit, rand1, rand2, rand3, false, debug);
  //~ CHECK_CLOSE_ABS(s_at_line.first, s_free_at_line.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_line.first, e_at_line.first, 1_um);
  //~ CHECK_CLOSE_ABS(s_at_line.first, e_free_at_line.first, 1_um);
//~ }

//~ /// test correct covariance transport for curvilinear parameters
//~ /// this test only works within the
//~ /// s_curvilinearProjTolerance (in: Definitions.hpp)
//~ BOOST_DATA_TEST_CASE(covariance_transport_to_curvilinear,
                     //~ ds::trackParameters* ds::propagationLimit, pT, phi, theta,
                     //~ charge, plimit) {
  //~ // define start parameters
  //~ double x = 1.;
  //~ double y = 0.;
  //~ double z = 0.;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ std::optional<FreeCovariance> covOptFree = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;

    //~ FreeCovariance covFree;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ covFree <<
     //~ 10_mm, 0, 0.123, 0, 0.01, 0.01, 0.01, 0,
     //~ 0, 10_mm, 0, 0, 0.01, 0.01, 0.01, 0,
     //~ 0.123, 0, 10_mm, 0, 0.01, 0.01, 0.1, 0,
     //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
     //~ 0.01, 0.01, 0.01, 0, 0.0123, 0, 0, 0,
     //~ 0.01, 0.01, 0.01, 0, 0, 0.0123, 0, 0,
     //~ 0.01, 0.01, 0.01, 0, 0, 0, 0.0123, 0,
     //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
    //~ // clang-format on
    //~ covOptFree = covFree;
  //~ }
  //~ CurvilinearParameters startC(covOpt, pos, mom, q, time);
  //~ NeutralCurvilinearTrackParameters startN(covOpt, pos, mom, time);

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(covOptFree, parsC);
  //~ NeutralFreeTrackParameters startNF(covOptFree, parsN);

  //~ ///
  //~ /// Curvilinear to Curvilinear Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ CHECK_CLOSE_COVARIANCE(covariance_curvilinear<CurvilinearParameters>(
                             //~ rspropagator, startN, plimit),
                         //~ covariance_curvilinear<CurvilinearParameters>(
                             //~ spropagator, startN, plimit),
                         //~ 1e-3);
  //~ // covariance check for eigen stepper
  //~ CHECK_CLOSE_COVARIANCE(covariance_curvilinear<CurvilinearParameters>(
                             //~ repropagator, startC, plimit),
                         //~ covariance_curvilinear<CurvilinearParameters>(
                             //~ epropagator, startC, plimit),
                         //~ 1e-3);
  //~ // covariance check fo atlas stepper
  //~ CHECK_CLOSE_COVARIANCE(covariance_curvilinear<CurvilinearParameters>(
                             //~ rapropagator, startC, plimit),
                         //~ covariance_curvilinear<CurvilinearParameters>(
                             //~ apropagator, startC, plimit),
                         //~ 1e-3);

  //~ ///
  //~ /// Free to Curvilinear Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ CHECK_CLOSE_COVARIANCE(covariance_curvilinear<CurvilinearParameters>(
                             //~ spropagator, startNF, plimit),
                         //~ covariance_curvilinear<CurvilinearParameters>(
                             //~ rspropagator, startNF, plimit),
                         //~ 1e-3);
  //~ // covariance check for eigen stepper
  //~ CHECK_CLOSE_COVARIANCE(covariance_curvilinear<CurvilinearParameters>(
                             //~ epropagator, startCF, plimit),
                         //~ covariance_curvilinear<CurvilinearParameters>(
                             //~ repropagator, startCF, plimit),
                         //~ 2e-3);
//~ }

//~ // test correct covariance transport from disc to disc
//~ BOOST_DATA_TEST_CASE(covariance_transport_to_disc,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ // define start parameters
  //~ double x = 1.;
  //~ double y = 0.;
  //~ double z = 0.;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ std::optional<FreeCovariance> covOptFree = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;

    //~ FreeCovariance covFree;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ covFree <<
     //~ 10_mm, 0, 0.123, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 10_mm, 0, 0, 0.1, 0.1, 0.1, 0,
     //~ 0.123, 0, 10_mm, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0.123, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0.123, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0, 0.123, 0,
     //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
    //~ // clang-format on
    //~ covOptFree = covFree;
  //~ }

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(covOptFree, parsC);
  //~ NeutralFreeTrackParameters startNF(covOptFree, parsN);

  //~ auto ssTransform =
      //~ createPlanarTransform(pos, mom.normalized(), 0.05 * rand1, 0.05 * rand2);
  //~ auto startSurface = Surface::makeShared<DiscSurface>(ssTransform, nullptr);
  //~ BoundParameters startC(tgContext, covOpt, pos, mom, q, time, startSurface);
  //~ NeutralBoundTrackParameters startN(tgContext, covOpt, pos, mom, time,
                                     //~ startSurface);

  //~ ///
  //~ /// Disc to Disc Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ auto covCalculated =
      //~ covariance_bound<RiddersStraightPropagatorType, DiscSurface>(
          //~ rspropagator, plimit, rand1, rand2, rand3, startN);
  //~ auto covObtained = covariance_bound<StraightPropagatorType, DiscSurface>(
      //~ spropagator, plimit, rand1, rand2, rand3, startN);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  //~ }

  //~ // covariance check for eigen stepper
  //~ covCalculated = covariance_bound<RiddersEigenPropagatorType, DiscSurface>(
      //~ repropagator, plimit, rand1, rand2, rand3, startC);
  //~ covObtained = covariance_bound<EigenPropagatorType, DiscSurface>(
      //~ epropagator, plimit, rand1, rand2, rand3, startC);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  //~ }

  //~ // covariance check for atlas stepper
  //~ covCalculated = covariance_bound<RiddersAtlasPropagatorType, DiscSurface>(
      //~ rapropagator, plimit, rand1, rand2, rand3, startC);
  //~ covObtained = covariance_bound<AtlasPropagatorType, DiscSurface>(
      //~ apropagator, plimit, rand1, rand2, rand3, startC);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  //~ }

  //~ ///
  //~ /// Free to Disc Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ covCalculated = covariance_bound<RiddersStraightPropagatorType, DiscSurface>(
      //~ rspropagator, plimit, rand1, rand2, rand3, startNF);
  //~ covObtained = covariance_bound<StraightPropagatorType, DiscSurface>(
      //~ spropagator, plimit, rand1, rand2, rand3, startNF);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  //~ }

  //~ // covariance check for eigen stepper
  //~ covCalculated = covariance_bound<RiddersEigenPropagatorType, DiscSurface>(
      //~ repropagator, plimit, rand1, rand2, rand3, startCF);
  //~ covObtained = covariance_bound<EigenPropagatorType, DiscSurface>(
      //~ epropagator, plimit, rand1, rand2, rand3, startCF);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  //~ }
//~ }

//~ // test correct covariance transport from plane to plane
//~ BOOST_DATA_TEST_CASE(covariance_transport_to_plane,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ // define start parameters
  //~ double x = 1.;
  //~ double y = 0.;
  //~ double z = 0.;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ std::optional<FreeCovariance> covOptFree = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;

    //~ FreeCovariance covFree;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ covFree <<
     //~ 10_mm, 0, 0.123, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 10_mm, 0, 0, 0.1, 0.1, 0.1, 0,
     //~ 0.123, 0, 10_mm, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0.123, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0.123, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0, 0.123, 0,
     //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
    //~ // clang-format on
    //~ covOptFree = covFree;
  //~ }

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(covOptFree, parsC);
  //~ NeutralFreeTrackParameters startNF(covOptFree, parsN);

  //~ auto ssTransform =
      //~ createPlanarTransform(pos, mom.normalized(), 0.05 * rand1, 0.05 * rand2);

  //~ auto startSurface = Surface::makeShared<PlaneSurface>(ssTransform, nullptr);
  //~ BoundParameters startC(tgContext, covOpt, pos, mom, q, time, startSurface);
  //~ NeutralBoundTrackParameters startN(tgContext, covOpt, pos, mom, time,
                                     //~ startSurface);

  //~ ///
  //~ /// Plane to Plane Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ auto covCalculated =
      //~ covariance_bound<RiddersStraightPropagatorType, PlaneSurface>(
          //~ rspropagator, plimit, rand1, rand2, rand3, startN);
  //~ auto covObtained = covariance_bound<StraightPropagatorType, PlaneSurface>(
      //~ spropagator, plimit, rand1, rand2, rand3, startN);
  //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-3);

  //~ // covariance check for eigen stepper
  //~ covCalculated = covariance_bound<RiddersEigenPropagatorType, PlaneSurface>(
      //~ repropagator, plimit, rand1, rand2, rand3, startC);
  //~ covObtained = covariance_bound<EigenPropagatorType, PlaneSurface>(
      //~ epropagator, plimit, rand1, rand2, rand3, startC);
  //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-2);

  //~ // covariance check for atlas stepper
  //~ covCalculated = covariance_bound<RiddersAtlasPropagatorType, PlaneSurface>(
      //~ rapropagator, plimit, rand1, rand2, rand3, startC);
  //~ covObtained = covariance_bound<AtlasPropagatorType, PlaneSurface>(
      //~ apropagator, plimit, rand1, rand2, rand3, startC);
  //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-2);

  //~ ///
  //~ /// Free to Plane Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ covCalculated = covariance_bound<RiddersStraightPropagatorType, PlaneSurface>(
      //~ rspropagator, plimit, rand1, rand2, rand3, startNF);
  //~ covObtained = covariance_bound<StraightPropagatorType, PlaneSurface>(
      //~ spropagator, plimit, rand1, rand2, rand3, startNF);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 2e-2);
  //~ }
  //~ // covariance check for eigen stepper
  //~ covCalculated = covariance_bound<RiddersEigenPropagatorType, PlaneSurface>(
      //~ repropagator, plimit, rand1, rand2, rand3, startCF);
  //~ covObtained = covariance_bound<EigenPropagatorType, PlaneSurface>(
      //~ epropagator, plimit, rand1, rand2, rand3, startCF);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);
  //~ }
//~ }

//~ // test correct covariance transport from straw to straw
//~ // for straw surfaces the numerical fixture is actually more difficult
//~ // to calculate
//~ BOOST_DATA_TEST_CASE(covariance_transport_to_line,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ // define start parameters
  //~ double x = 1.;
  //~ double y = 0.;
  //~ double z = 0.;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ std::optional<FreeCovariance> covOptFree = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;

    //~ FreeCovariance covFree;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ covFree <<
     //~ 10_mm, 0, 0.123, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 10_mm, 0, 0, 0.1, 0.1, 0.1, 0,
     //~ 0.123, 0, 10_mm, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0.123, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0.123, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0, 0.123, 0,
     //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
    //~ // clang-format on
    //~ covOptFree = covFree;
  //~ }

  //~ Vector3D dir = mom.normalized();
  //~ FreeVector parsC, parsN;
  //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  //~ FreeTrackParameters startCF(covOptFree, parsC);
  //~ NeutralFreeTrackParameters startNF(covOptFree, parsN);

  //~ auto ssTransform = createCylindricTransform(pos, 0.01 * rand1, 0.01 * rand2);

  //~ auto startSurface = Surface::makeShared<StrawSurface>(ssTransform, nullptr);
  //~ BoundParameters startC(tgContext, covOpt, pos, mom, q, time, startSurface);
  //~ NeutralBoundTrackParameters startN(tgContext, covOpt, pos, mom, time,
                                     //~ startSurface);

  //~ ///
  //~ /// Line to Line Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ auto covCalculated =
      //~ covariance_bound<RiddersStraightPropagatorType, StrawSurface>(
          //~ rspropagator, plimit, rand1, rand2, rand3, startN, false);
  //~ auto covObtained = covariance_bound<StraightPropagatorType, StrawSurface>(
      //~ spropagator, plimit, rand1, rand2, rand3, startN, false);
  //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);

  //~ // covariance check for eigen stepper
  //~ covCalculated = covariance_bound<RiddersEigenPropagatorType, StrawSurface>(
      //~ repropagator, plimit, rand1, rand2, rand3, startC, false);
  //~ covObtained = covariance_bound<EigenPropagatorType, StrawSurface>(
      //~ epropagator, plimit, rand1, rand2, rand3, startC, false);
  //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);

  //~ // covariance check for atlas stepper
  //~ covCalculated = covariance_bound<RiddersAtlasPropagatorType, StrawSurface>(
      //~ rapropagator, plimit, rand1, rand2, rand3, startC, false);
  //~ covObtained = covariance_bound<AtlasPropagatorType, StrawSurface>(
      //~ apropagator, plimit, rand1, rand2, rand3, startC, false);
  //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);

  //~ ///
  //~ /// Free to Line Tests
  //~ ///
  //~ // covariance check for straight line stepper
  //~ covCalculated = covariance_bound<RiddersStraightPropagatorType, StrawSurface>(
      //~ rspropagator, plimit, rand1, rand2, rand3, startNF, false);
  //~ covObtained = covariance_bound<StraightPropagatorType, StrawSurface>(
      //~ spropagator, plimit, rand1, rand2, rand3, startNF, false);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  //~ }

  //~ // covariance check for eigen stepper
  //~ covCalculated = covariance_bound<RiddersEigenPropagatorType, StrawSurface>(
      //~ repropagator, plimit, rand1, rand2, rand3, startCF, false);
  //~ covObtained = covariance_bound<EigenPropagatorType, StrawSurface>(
      //~ epropagator, plimit, rand1, rand2, rand3, startCF, false);
  //~ if (covCalculated != Covariance::Zero()) {
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);
  //~ }
//~ }

//~ BOOST_DATA_TEST_CASE(covariance_transport_to_free,
                     //~ ds::trackParameters* ds::propagationLimit ^
                         //~ ds::threeRandom,
                     //~ pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  //~ (void)rand3;
  //~ // define start parameters
  //~ double x = 1;
  //~ double y = 0;
  //~ double z = 0;
  //~ double px = pT * cos(phi);
  //~ double py = pT * sin(phi);
  //~ double pz = pT / tan(theta);
  //~ double q = charge;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ std::optional<FreeCovariance> covOptFree = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_um, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_um, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;

    //~ FreeCovariance covFree;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ covFree <<
     //~ 10_mm, 0, 0.123, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 10_mm, 0, 0, 0.1, 0.1, 0.1, 0,
     //~ 0.123, 0, 10_mm, 0, 0.1, 0.1, 0.1, 0,
     //~ 0, 0, 0, 1_ns, 0, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0.123, 0, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0.123, 0, 0,
     //~ 0.1, 0.1, 0.1, 0, 0, 0, 0.123, 0,
     //~ 0, 0, 0, 0, 0, 0, 0, 1_e / 10_GeV;
    //~ // clang-format on
    //~ covOptFree = covFree;
  //~ }

  //~ ///
  //~ /// Curvilinear to Free Tests
  //~ ///
  //~ {
    //~ CurvilinearParameters startCC(covOpt, pos, mom, q, time);
    //~ NeutralCurvilinearTrackParameters startNC(covOpt, pos, mom, time);

    //~ // covariance check for straight line stepper
    //~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        //~ spropagator, startNC, plimit);
    //~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ rspropagator, startNC, plimit);
    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
    //~ // covariance check for eigen stepper
    //~ covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ repropagator, startCC, plimit);
        
    //~ covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              //~ startCC, plimit);


    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
  //~ }

  //~ ///
  //~ /// Disc to Free Tests
  //~ ///
  //~ {
    //~ auto ssTransform = createPlanarTransform(pos, mom.normalized(),
                                             //~ 0.05 * rand1, 0.05 * rand2);
    //~ auto startSurface = Surface::makeShared<DiscSurface>(ssTransform, nullptr);
    //~ BoundParameters startCB(tgContext, covOpt, pos, mom, q, time, startSurface);
    //~ NeutralBoundTrackParameters startNB(tgContext, covOpt, pos, mom, time,
                                        //~ startSurface);

    //~ // covariance check for straight line stepper
    //~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        //~ spropagator, startNB, plimit);
    //~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ rspropagator, startNB, plimit);
    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);

    //~ // covariance check for eigen stepper
    //~ covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              //~ startCB, plimit);
    //~ covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ repropagator, startCB, plimit);

    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
  //~ }

  //~ ///
  //~ /// Plane to Free Tests
  //~ ///
  //~ {
    //~ auto ssTransform = createPlanarTransform(pos, mom.normalized(),
                                             //~ 0.05 * rand1, 0.05 * rand2);
    //~ auto startSurface = Surface::makeShared<PlaneSurface>(ssTransform, nullptr);
    //~ BoundParameters startCB(tgContext, covOpt, pos, mom, q, time, startSurface);
    //~ NeutralBoundTrackParameters startNB(tgContext, covOpt, pos, mom, time,
                                        //~ startSurface);

    //~ // covariance check for straight line stepper
    //~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        //~ spropagator, startNB, plimit);
    //~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ rspropagator, startNB, plimit);
    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);

    //~ // covariance check for eigen stepper
    //~ covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              //~ startCB, plimit);
    //~ covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ repropagator, startCB, plimit);

    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
  //~ }

  //~ ///
  //~ /// Line to Free Tests
  //~ ///
  //~ {
    //~ auto ssTransform =
        //~ createCylindricTransform(pos, 0.005 * rand1, 0.005 * rand2);
    //~ auto startSurface = Surface::makeShared<StrawSurface>(ssTransform, nullptr);
    //~ BoundParameters startCB(tgContext, covOpt, pos, mom, q, time, startSurface);
    //~ NeutralBoundTrackParameters startNB(tgContext, covOpt, pos, mom, time,
                                        //~ startSurface);

    //~ // covariance check for straight line stepper
    //~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        //~ spropagator, startNB, plimit);
    //~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ rspropagator, startNB, plimit);
    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 5e-1);

    //~ // covariance check for eigen stepper
    //~ covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              //~ startCB, plimit);
    //~ covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ repropagator, startCB, plimit);

    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-1);
  //~ }

  //~ ///
  //~ /// Free to Free Tests
  //~ ///
  //~ {
    //~ Vector3D dir = mom.normalized();
    //~ FreeVector parsC, parsN;
    //~ parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
    //~ parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
    //~ FreeTrackParameters startCF(covOptFree, parsC);
    //~ NeutralFreeTrackParameters startNF(covOptFree, parsN);

    //~ // covariance check for straight line stepper
    //~ auto covObtained = covariance_curvilinear<FreeTrackParameters>(
        //~ spropagator, startNF, plimit);
    //~ auto covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ rspropagator, startNF, plimit);
    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);

    //~ // covariance check for eigen stepper
    //~ covObtained = covariance_curvilinear<FreeTrackParameters>(epropagator,
                                                              //~ startCF, plimit);
    //~ covCalculated = covariance_curvilinear<FreeTrackParameters>(
        //~ repropagator, startCF, plimit);
    //~ // Numerical fluctuations in the covariances cause errors in relative
    //~ // comparison. This needs to be tested and avoided by setting both entries
    //~ // to 1
    //~ for (unsigned int i = 0; i < covObtained.rows(); i++)
      //~ for (unsigned int j = 0; j < covObtained.cols(); j++) {
        //~ if (std::abs(covObtained(i, j)) <
                //~ std::numeric_limits<double>::epsilon() ||
            //~ std::abs(covCalculated(i, j)) <
                //~ std::numeric_limits<double>::epsilon()) {
          //~ covObtained(i, j) = 1.;
          //~ covCalculated(i, j) = 1.;
        //~ }
      //~ }
    //~ CHECK_CLOSE_COVARIANCE(covObtained, covCalculated, 1e-3);
  //~ }
//~ }

//~ /// test correct covariance transport for curvilinear parameters in dense
//~ /// environment
//~ /// this test only works within the
//~ /// s_curvilinearProjTolerance (in: Definitions.hpp)
//~ BOOST_DATA_TEST_CASE(dense_covariance_transport,
                     //~ ds::pT* ds::propagationLimit ^ ds::threeRandom, pT, plimit,
                     //~ rand1, rand2, rand3) {
  //~ // define start parameters
  //~ double x = 1.;
  //~ double y = 0.;
  //~ double z = 0.;
  //~ double px = pT * cos(0_degree);
  //~ double py = pT * sin(0_degree);
  //~ double pz = pT / tan(45_degree);
  //~ double q = 1_e;
  //~ double time = 0.;
  //~ Vector3D pos(x, y, z);
  //~ Vector3D mom(px, py, pz);

  //~ std::optional<Covariance> covOpt = std::nullopt;
  //~ if (covtpr) {
    //~ Covariance cov;
    //~ // take some major correlations (off-diagonals)
    //~ // clang-format off
    //~ cov <<
     //~ 10_mm, 0, 0.123, 0, 0.5, 0,
     //~ 0, 10_mm, 0, 0.162, 0, 0,
     //~ 0.123, 0, 0.1, 0, 0, 0,
     //~ 0, 0.162, 0, 0.1, 0, 0,
     //~ 0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     //~ 0, 0, 0, 0, 0, 1_us;
    //~ // clang-format on
    //~ covOpt = cov;
  //~ }

  //~ auto ssTransform =
      //~ createPlanarTransform(pos, mom.normalized(), 0.05 * rand1, 0.05 * rand2);

  //~ // covariance check for eigen stepper in dense environment
  //~ DensePropagatorType dpropagator = setupDensePropagator();
  //~ RiddersPropagator rdpropagator(dpropagator);
  //~ {
    //~ CurvilinearParameters startCC(covOpt, pos, mom, q, time);

    //~ CHECK_CLOSE_COVARIANCE(covariance_curvilinear<CurvilinearParameters>(
                               //~ rdpropagator, startCC, plimit),
                           //~ covariance_curvilinear<CurvilinearParameters>(
                               //~ dpropagator, startCC, plimit),
                           //~ 3e-1);
  //~ }
  //~ {
    //~ auto startSurface = Surface::makeShared<DiscSurface>(ssTransform, nullptr);
    //~ BoundParameters startC(tgContext, covOpt, pos, mom, q, time, startSurface);
    //~ auto covCalculated =
        //~ covariance_bound<RiddersPropagator<DensePropagatorType>, DiscSurface>(
            //~ rdpropagator, plimit, rand1, rand2, rand3, startC);
    //~ auto covObtained = covariance_bound<DensePropagatorType, DiscSurface>(
        //~ dpropagator, plimit, rand1, rand2, rand3, startC);
    //~ if (covCalculated != Covariance::Zero()) {
      //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 8e-1);
    //~ }
  //~ }
  //~ {
    //~ auto startSurface = Surface::makeShared<PlaneSurface>(ssTransform, nullptr);
    //~ BoundParameters startC(tgContext, covOpt, pos, mom, q, time, startSurface);
    //~ auto covCalculated =
        //~ covariance_bound<RiddersPropagator<DensePropagatorType>, PlaneSurface>(
            //~ rdpropagator, plimit, rand1, rand2, rand3, startC);
    //~ auto covObtained = covariance_bound<DensePropagatorType, PlaneSurface>(
        //~ dpropagator, plimit, rand1, rand2, rand3, startC);
    //~ CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 8e-1);
  //~ }
//~ }