#ifndef TEST_OOPSERANDNUMGENTESTCASE_HPP
#define TEST_OOPSERANDNUMGENTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/OOPSERandNumGen.hpp"
 


class OOPSERandNumGenTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( OOPSERandNumGenTestCase );
    CPPUNIT_TEST(testUniform);
    CPPUNIT_TEST(testGaussian);
    CPPUNIT_TEST(testMPIRNG);    

    CPPUNIT_TEST_SUITE_END();

    public:

        void testUniform();

        void testGaussian();

        void testMPIRNG();

        
};


#endif //TEST_POLYNOMIALTESTCASE_HPP



