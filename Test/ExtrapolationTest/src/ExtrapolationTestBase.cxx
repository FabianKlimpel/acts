///////////////////////////////////////////////////////////////////
// ExtrapolationTestBase.cxx, ATS project
//////////////////////////////////////////////////////////////////

// Test module
#include "ExtrapolationTest/ExtrapolationTestBase.h"

Ats::ExtrapolationTestBase::ExtrapolationTestBase(const std::string& name, ISvcLocator* pSvcLocator):
    AlgorithmBase(name,pSvcLocator),
    m_gaussDist(0),
    m_flatDist(0),
    m_numTests(100),
    m_scanMode(false)
{
  declareProperty("NumberOfTestsPerEvent",   m_numTests);
  declareProperty("EnableScanMode",          m_scanMode);
}

Ats::ExtrapolationTestBase::~ExtrapolationTestBase()
{
    delete m_gaussDist;
    delete m_flatDist;
    delete m_landauDist;
}

StatusCode Ats::ExtrapolationTestBase::initialize()
{
    MSG_INFO( "Creating random number services, call bookTree() and initializeTest()" );

    // intialize the random number generators
    m_gaussDist  = new Rndm::Numbers(randSvc(), Rndm::Gauss(0.,1.));
    m_flatDist   = new Rndm::Numbers(randSvc(), Rndm::Flat(0.,1.));
    m_landauDist = new Rndm::Numbers(randSvc(), Rndm::Landau(0.,1.)); 

    if  (bookTree().isFailure()){
        MSG_FATAL( "Could not book the TTree object" );
        return StatusCode::FAILURE;
    }
    
    if (initializeTest().isFailure()){
        MSG_FATAL( "Could not initialize the test" );
        return StatusCode::FAILURE;
    }
    
    return StatusCode::SUCCESS;
}
 
StatusCode Ats::ExtrapolationTestBase::execute()
{
  return runTest();
}

StatusCode Ats::ExtrapolationTestBase::finalize()
{

   MSG_INFO( "finalize()" );
   return StatusCode::SUCCESS;
}

StatusCode Ats::ExtrapolationTestBase::bookTree() { return StatusCode::SUCCESS; }

StatusCode Ats::ExtrapolationTestBase::initializeTest() { return StatusCode::SUCCESS; }

