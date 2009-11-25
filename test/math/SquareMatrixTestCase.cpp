#include "SquareMatrixTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SquareMatrixTestCase );

void SquareMatrixTestCase::setUp() {

    identMat(0, 0) = 1.0;
    identMat(0, 1) = 0.0;
    identMat(0, 2) = 0.0;
    identMat(1, 0) = 0.0;
    identMat(1, 1) = 1.0;
    identMat(1, 2) = 0.0;    
    identMat(2, 0) = 0.0;
    identMat(2, 1) = 0.0;
    identMat(2, 2) = 1.0;

    symMat(0, 0) = 2.0;
    symMat(0, 1) = 4.0;
    symMat(0, 2) = 0.5;
    symMat(1, 0) = 4.0;
    symMat(1, 1) = 1.0;
    symMat(1, 2) = 3.0;    
    symMat(2, 0) = 0.5;
    symMat(2, 1) = 3.0;
    symMat(2, 2) = 1.0;    

    ortMat(0, 0) = 1.0;
    ortMat(0, 1) = 0.0;
    ortMat(0, 2) = 0.0;
    ortMat(1, 0) = 0.0;
    ortMat(1, 1) = cos(0.5);
    ortMat(1, 2) = -sin(0.5);    
    ortMat(2, 0) = 0.0;
    ortMat(2, 1) = sin(0.5);
    ortMat(2, 2) = cos(0.5); 
    
    diagMat(0, 0) = 8.0;
    diagMat(0, 1) = 0.0;
    diagMat(0, 2) = 0.0;
    diagMat(1, 0) = 0.0;
    diagMat(1, 1) = 1.0;
    diagMat(1, 2) = 0.0;    
    diagMat(2, 0) = 0.0;
    diagMat(2, 1) = 0.0;
    diagMat(2, 2) = 3.0;    

}

void SquareMatrixTestCase::testConstructor() {

}

void SquareMatrixTestCase::testIdentity() {
    CPPUNIT_ASSERT(SMat3x3::identity() == identMat);
}


void SquareMatrixTestCase::testJacobi() {
    SMat3x3 a;
    Vector<double, 3> w1L;
    Vector<double, 3> w1R;    
    SMat3x3 v;
    a(0, 0) = 3.0;
    a(0, 1) = 4.0;
    a(0, 2) = 5.0;
    a(1, 0) = 4.0;
    a(1, 1) = 5.0;
    a(1, 2) = 6.0;    
    a(2, 0) = 5.0;
    a(2, 1) = 6.0;
    a(2, 2) = 7.0;   

    w1R[0] = 15.3899;
    w1R[1] = 0.0;
    w1R[2] = -0.389867;
    
    SMat3x3::jacobi(a, w1L, v);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(w1L[0], w1R[0], 0.0001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(w1L[1], w1R[1], OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(w1L[2], w1R[2], OpenMD::NumericConstant::epsilon);

}

void SquareMatrixTestCase::testTrace() {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(identMat.trace(), 3.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(symMat.trace(), 4.0, OpenMD::NumericConstant::epsilon);
}
void SquareMatrixTestCase::testIsSymmertric(){
    CPPUNIT_ASSERT(identMat.isSymmetric());
    CPPUNIT_ASSERT(symMat.isSymmetric());
}

void SquareMatrixTestCase::testIsOrthogonal(){
    CPPUNIT_ASSERT(ortMat.isOrthogonal());
}

void SquareMatrixTestCase::testIsDiagonal() {
    CPPUNIT_ASSERT(identMat.isDiagonal());
    CPPUNIT_ASSERT(diagMat.isDiagonal());
    CPPUNIT_ASSERT(!symMat.isDiagonal());
}
void SquareMatrixTestCase::testIsUnitMatrix() {
    CPPUNIT_ASSERT(identMat.isUnitMatrix());
    CPPUNIT_ASSERT(!symMat.isUnitMatrix());
}
