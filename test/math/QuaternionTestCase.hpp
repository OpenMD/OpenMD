#ifndef TEST_QUATERNIONTESTCASE_HPP
#define TEST_QUATERNIONTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/Quaternion.hpp"
 
using namespace oopse;

typedef Vector<double, 4> Vec4;
class QuaternionTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( QuaternionTestCase );
    CPPUNIT_TEST(testConstructors);
    CPPUNIT_TEST(testArithmetic);
    CPPUNIT_TEST(testAccessEntries);
    CPPUNIT_TEST(testOtherMemberFunctions);
    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();

        void testConstructors();
        void testArithmetic();
        void testOperators();
        void testAccessEntries();
        void testOtherMemberFunctions();
        
    private:
        Quat4d q1;
        Quat4d q2;
        Quat4d q3;
        Quat4d q4;
        
};


#endif //TEST_QUATERNIONTESTCASE_HPP

