#ifndef TEST_SQUAREMATRIXTESTCASE_HPP
#define TEST_SQUAREMATRIXTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/SquareMatrix.hpp"

using namespace oopse;

typedef SquareMatrix<double, 3> SMat3x3;

class SquareMatrixTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( SquareMatrixTestCase );
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testIdentity);
    CPPUNIT_TEST(testInverse);
    CPPUNIT_TEST(testDeterminant);
    CPPUNIT_TEST(testTrace);
    CPPUNIT_TEST(testIsSymmertric);
    CPPUNIT_TEST(testIsOrthogonal);
    CPPUNIT_TEST(testIsDiagonal);
    CPPUNIT_TEST(testIsUnitMatrix);
    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();

        void testConstructor();
        void testIdentity();
        void testInverse();
        void testDeterminant();
        void testTrace();
        void testIsSymmertric();
        void testIsOrthogonal();
        void testIsDiagonal();
        void testIsUnitMatrix();

    private:

        SMat3x3 identMat;
        SMat3x3 invMat;
        SMat3x3 symMat;
        SMat3x3 ortMat;
        SMat3x3 diagMat;        
        SMat3x3 unitMat;        
        
};

#endif // TEST_SQUAREMATRIXTESTCASE_HPP
