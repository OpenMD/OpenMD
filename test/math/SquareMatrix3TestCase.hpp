#ifndef TEST_SQUAREMATRIX3TESTCASE_HPP
#define TEST_SQUAREMATRIX3TESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/SquareMatrix3.hpp"
 
using namespace oopse;

class SquareMatrix3TestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( SquareMatrix3TestCase );
    CPPUNIT_TEST(testSetupRotationMatrix);
    CPPUNIT_TEST(testOtherMemberFunctions);
    CPPUNIT_TEST(testTransformation);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testSetupRotationMatrix();
        void testOtherMemberFunctions();
        void testTransformation();      
};

#endif //TEST_SQUAREMATRIX3TESTCASE_HPP
