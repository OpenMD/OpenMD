#ifndef TEST_QUATERNIONTESTCASE_HPP
#define TEST_QUATERNIONTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/Quaternion.hpp"
 
using namespace oopse;

class QuaternionTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( QuaternionTestCase );
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

#endif //TEST_QUATERNIONTESTCASE_HPP

