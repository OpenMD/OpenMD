#ifndef TEST_VECTORTESTCASE_HPP
#define TEST_VECTORTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/Vector.hpp"
 
using namespace oopse;


typedef Vector<double, 3> Vec3; 
typedef Vector<double, 4> Vec4;

class RectMatrixTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( RectMatrixTestCase );
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testNegate);
    CPPUNIT_TEST(testAdd);
    CPPUNIT_TEST(testSub);
    CPPUNIT_TEST(testMul);
    CPPUNIT_TEST(testDiv);
    CPPUNIT_TEST(testAccessEntries);
    CPPUNIT_TEST(testTranspose);
    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();
        virtual void tearDown();

        
};

#endif //TEST_VECTORTESTCASE_HPP
