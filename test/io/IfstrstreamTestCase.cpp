#include <iostream>
#include <string>

#include "io/basic_ifstrstream.hpp"
#include "IfstrstreamTestCase.hpp"

// Registers the fixture into the 'registry'

using namespace std;
using namespace oopse;
CPPUNIT_TEST_SUITE_REGISTRATION( IfstrstreamTestCase );

#ifndef IS_MPI
void IfstrstreamTestCase::setUp() {

    
}

void IfstrstreamTestCase::tearDown() {

}


void IfstrstreamTestCase::testConstructor() {
    ifstrstream fin1;
    fin1.open("DUFF.frc");
    CPPUNIT_ASSERT(fin1.is_open());

    ifstrstream fin2;
    fin2.open("NonExistFile");
    CPPUNIT_ASSERT(!fin2.is_open());

    ifstrstream fin3("DUFF.frc");
    CPPUNIT_ASSERT(fin3.is_open());

    ifstrstream fin4("NonExistFile");
    CPPUNIT_ASSERT(!fin4.is_open());
}

void IfstrstreamTestCase::testOpen() {
    const int MAXLEN = 1024;
    char buffer[MAXLEN];

    string firstLine = "! This is the forcefield file for the Dipolar Unified-atom Force Field (DUFF).";

    ifstrstream fin1;
    fin1.open("DUFF.frc");

    fin1.getline(buffer, MAXLEN);
    CPPUNIT_ASSERT(buffer == firstLine);
    cout << buffer;
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
    const int MAXLEN = 1024;
    char buffer[MAXLEN];

    string firstLine = "! This is the forcefield file for the Dipolar Unified-atom Force Field (DUFF).";

    ifstrstream fin1;
    fin1.open("DUFF.frc");

    fin1.getline(buffer, MAXLEN);
    CPPUNIT_ASSERT(buffer == firstLine);
    cout << buffer;
    
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
