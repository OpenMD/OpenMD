#include "SquareMatrixTestCase.hpp"
#include "io/basic_ifstrstream.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( IfstrstreamTestCase );

#ifndef IS_MPI
void IfstrstreamTestCase::setUp() {

}

void IfstrstreamTestCase::tearDown() {

}


void IfstrstreamTestCase::testConstrutor() {

}

void IfstrstreamTestCase::testOpen() {

}

void IfstrstreamTestCase::testIs_open() {

}

void IfstrstreamTestCase::testClose() {

}

#else
void IfstrstreamTestCase::testMasterConstrutor() {

}    

void IfstrstreamTestCase::testMasterOpen() {

}

void IfstrstreamTestCase::testMasterIs_open() {

}

void IfstrstreamTestCase::testMasterClose() {


}

void IfstrstreamTestCase::testSlaveConstrutor() {

}    

void IfstrstreamTestCase::testSlaveOpen() {

}

void IfstrstreamTestCase::testSlaveIs_open() {

}

void IfstrstreamTestCase::testSlaveClose() {


}
#endif