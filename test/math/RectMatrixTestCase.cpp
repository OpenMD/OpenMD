#include "RectMatrixTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RectMatrixTestCase );

void RectMatrixTestCase::setUp(){

    m1(0, 0) = 1.0;
    m1(0, 1) = 1.0;
    m1(1, 0) = 3.0;
    m1(1, 1) = 1.0;
    
    m2(0, 0) = 1.0;
    m2(0, 1) = 0.0;
    m2(1, 0) = 0.0;
    m2(1, 1) = 1.0;
 
    
    m3(0, 0) = 1.0;
    m3(0, 1) = 0.0;
    m3(1, 0) = 0.0;
    m3(1, 1) = 1.0;

    m4(0, 0) = -1.0;
    m4(0, 1) = -1.0;
    m4(1, 0) = -3.0;
    m4(1, 1) = -1.0;

    zero(0, 0) = 0.0;
    zero(0, 1) = 0.0;
    zero(1, 0) = 0.0;
    zero(1, 1) = 0.0;    

    one(0, 0) = 1.0;
    one(0, 1) = 1.0;
    one(1, 0) = 1.0;
    one(1, 1) = 1.0;    

    two(0, 0) = 2.0;
    two(0, 1) = 2.0;
    two(1, 0) = 2.0;
    two(1, 1) = 2.0;  
    
    a(0, 0) = 1.0;
    a(0, 1) = 0.0;
    a(0, 2) = 0.0;
    a(1, 0) = 0.0;
    a(1, 1) = 1.0;
    a(1, 2) = 0.0;

    b(0, 0) = 1.0;
    b(0, 1) = 0.0;
    b(1, 0) = 0.0;
    b(1, 1) = 1.0;
    b(2, 0) = 0.0;
    b(2, 1) = 0.0;

    
}

void RectMatrixTestCase::testConstructor(){

    //test default constructor
    RMat2x2 tmp1;
    CPPUNIT_ASSERT(tmp1 == zero);

    //test RectMatrix(Real s)
    RMat2x2 tmp2;
    CPPUNIT_ASSERT(tmp2 == zero);

    //test copy constructor
    RMat2x2 tmp3(m1);
    CPPUNIT_ASSERT(tmp3 == m1);

    //test copy assignment
    RMat2x2 tmp4 = m2;
    tmp4 = tmp4;
    CPPUNIT_ASSERT(tmp4 == m2);
    
    
}

void RectMatrixTestCase::testEqual() {
    CPPUNIT_ASSERT(m2 == m3);
}

void RectMatrixTestCase::testNegate() {

    CPPUNIT_ASSERT(m1 == -m4);
}

void RectMatrixTestCase::testAdd() {
    RMat2x2 tmp;
    
    CPPUNIT_ASSERT(m1 + m4 == zero);

    tmp.add(m1, m1);

    CPPUNIT_ASSERT(m1 * 2.0 == tmp);

    tmp = one;
    tmp *= 2.0;

    CPPUNIT_ASSERT(tmp == two);
    
}

void RectMatrixTestCase::testSub() {
    
    RMat2x2 tmp(m2);
    tmp.sub(m2);
    CPPUNIT_ASSERT(tmp == zero);

    tmp = m2;
    tmp -= m2;
    CPPUNIT_ASSERT(tmp == zero);    

    tmp = m1;
    tmp.sub(m1 , m4);
    CPPUNIT_ASSERT(tmp == m1 * 2.0);
    
    CPPUNIT_ASSERT(m1 -m4 == m1 * 2.0);

    
}

void RectMatrixTestCase::testMul() {

    CPPUNIT_ASSERT(m1 * 1.0 == m1);
    CPPUNIT_ASSERT(m1 * 0.0 == zero);
    CPPUNIT_ASSERT(2.0 *m1 == m1 + m1);
    
}

void RectMatrixTestCase::testDiv() {

    CPPUNIT_ASSERT(m1 / 1.0 == m1);

    CPPUNIT_ASSERT(m1 / 2.0 * 2.0 == m1);

}

void RectMatrixTestCase::testAccessEntries(){
    CPPUNIT_ASSERT(m1(1, 0) == 3.0);
}

void RectMatrixTestCase::testTranspose(){

    CPPUNIT_ASSERT((a.transpose()).transpose() == a);
    
    CPPUNIT_ASSERT(a.transpose() == b);
}
