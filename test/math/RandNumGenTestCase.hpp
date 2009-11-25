#ifndef TEST_RANDNUMGENTESTCASE_HPP
#define TEST_RANDNUMGENTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/ParallelRandNumGen.hpp"
 


class RandNumGenTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( RandNumGenTestCase );
    CPPUNIT_TEST(testUniform);
    CPPUNIT_TEST(testGaussian);
    CPPUNIT_TEST(testMPIRNG);    

    CPPUNIT_TEST_SUITE_END();

    public:

        void testUniform();

        void testGaussian();

        void testMPIRNG();

        
};


#endif
