#ifndef TEST_OOPSERANDNUMGENTESTCASE_HPP
#define TEST_OOPSERANDNUMGENTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/ParallelandNumGenTestCase.hpp"
 


class ParallelandNumGenTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( ParallelandNumGenTestCase );
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



