#include "io/basic_ifstrstream.hpp"
#include "IfstrstreamTestCase.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( IfstrstreamTestCase );

#ifndef IS_MPI
void IfstrstreamTestCase::setUp() {

}

void IfstrstreamTestCase::tearDown() {

}


void IfstrstreamTestCase::testConstructor() {

}

void IfstrstreamTestCase::testOpen() {

}

void IfstrstreamTestCase::testIs_open() {

}

void IfstrstreamTestCase::testClose() {

}

#else
void IfstrstreamTestCase::setUp() {

}

void IfstrstreamTestCase::tearDown() {

}

void IfstrstreamTestCase::testMasterConstructor() {
        CPPUNIT_ASSERT(1 == 0);
}    

void IfstrstreamTestCase::testMasterOpen() {

}

void IfstrstreamTestCase::testMasterIs_open() {

}

void IfstrstreamTestCase::testMasterClose() {


}

void IfstrstreamTestCase::testSlaveConstructor() {

}    

void IfstrstreamTestCase::testSlaveOpen() {

}

void IfstrstreamTestCase::testSlaveIs_open() {

}

void IfstrstreamTestCase::testSlaveClose() {


}
#endif
