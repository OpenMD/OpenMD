#ifndef TEST_RECTMATRIXTEST_HPP
#define TEST_RECTMATRIXTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/RectMatrix.hpp"

/**
 * @namespace oopse
 */
using namespace oopse;

typedef RectMatrix<double, 2, 2> RMat2x2;
typedef RectMatrix<double, 2, 3> RMat2x3;
typedef RectMatrix<double, 3, 2> RMat3x2;
typedef RectMatrix<double, 3, 3> RMat3x3;
typedef RectMatrix<double, 3, 4> RMat3x4;

typedef Vector<double, 3> Vec3;

class RectMatrixTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( RectMatrixTestCase );
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testNegate);
    CPPUNIT_TEST(testAdd);
    CPPUNIT_TEST(testSub);
    CPPUNIT_TEST(testMul);
    CPPUNIT_TEST(testDiv);
    CPPUNIT_TEST(testAccessEntries);
    CPPUNIT_TEST(testRowColOperations);
    CPPUNIT_TEST(testOtherMemberFunctions);
    CPPUNIT_TEST_SUITE_END();

    public:
        
        virtual void setUp();
        
        void testConstructor();
        void testEqual();
        void testNegate();
        void testAdd();
        void testSub();
        void testMul();
        void testDiv();
        void testAccessEntries();
        void testRowColOperations();
        void testOtherMemberFunctions();

    private:
        RMat2x2 m1;
        RMat2x2 m2;
        RMat2x2 m3;
        RMat2x2 m4;
        RMat2x2 zero;
        RMat2x2 one;
        RMat2x2 two;
        
        RMat2x3 a;
        RMat3x2 b;
        RMat2x2 c;
        
        RMat3x3 d;
        RMat3x3 e;
        RMat3x3 f;

        RMat3x3 g;
        RMat3x4 h;
        RMat3x4 i;

        Vec3 v1;
        Vec3 v2;

        double s1;

        double s2;
                
};
#endif //TEST_RECTMATRIXTEST_HPP
