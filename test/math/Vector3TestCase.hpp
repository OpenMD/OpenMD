#ifndef TEST_VECTOR3TESTCASE_HPP
#define TEST_VECTOR3TESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/Vector3.hpp"
 
using namespace oopse;

class Vector3TestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( Vector3TestCase );
    CPPUNIT_TEST(testConstructors);
    CPPUNIT_TEST(testArithmetic);
    CPPUNIT_TEST(testAccessEntries);
    CPPUNIT_TEST(testOtherTemplateFunctions);
    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();
        virtual void tearDown();
        void testConstructors();
        void testArithmetic();
        void testOperators();
        void testAccessEntries();
        void testOtherTemplateFunctions();        
};

#endif //TEST_VECTO3RTESTCASE_HPP
