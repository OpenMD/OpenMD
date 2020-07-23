/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

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

    c(0, 0) = 1.0;
    c(0, 1) = 0.0;
    c(1, 0) = 0.0;
    c(1, 1) = 1.0;       

    d(0, 0) = 1.0;
    d(0, 1) = 0.0;
    d(0, 2) = 0.0;
    d(1, 0) = 0.0;
    d(1, 1) = 0.0;
    d(1, 2) = 1.0;    
    d(2, 0) = 0.0;
    d(2, 1) = 1.0;
    d(2, 2) = 0.0;  

    e(0, 0) = 2.0;
    e(0, 1) = 4.0;
    e(0, 2) = 1.0;
    e(1, 0) = 0.0;
    e(1, 1) = 0.0;
    e(1, 2) = 3.0;    
    e(2, 0) = 0.0;
    e(2, 1) = 6.0;
    e(2, 2) = 5.0;  

    f(0, 0) = 2.0;
    f(0, 1) = 4.0;
    f(0, 2) = 1.0;
    f(1, 0) = 0.0;
    f(1, 1) = 6.0;
    f(1, 2) = 5.0;    
    f(2, 0) = 0.0;
    f(2, 1) = 0.0;
    f(2, 2) = 3.0;  

    f(0, 0) = 2.0;
    f(0, 1) = 4.0;
    f(0, 2) = 1.0;
    f(1, 0) = 0.0;
    f(1, 1) = 6.0;
    f(1, 2) = 5.0;    
    f(2, 0) = 0.0;
    f(2, 1) = 0.0;
    f(2, 2) = 3.0;  

    g(0, 0) = 1.0;
    g(0, 1) = 0.0;
    g(0, 2) = 0.0;
    g(1, 0) = -2.0;
    g(1, 1) = 1.0;
    g(1, 2) = 0.0;    
    g(2, 0) = 0.0;
    g(2, 1) = 0.0;
    g(2, 2) = 1.0;  

    h(0, 0) = 2.0;
    h(0, 1) = 4.0;
    h(0, 2) = -2.0;
    h(0, 3) = 2.0;
    h(1, 0) = 4.0;
    h(1, 1) = 9.0;
    h(1, 2) = -3.0;
    h(1, 3) = 8.0;
    h(2, 0) = -2.0;
    h(2, 1) = -3.0;
    h(2, 2) = 7.0;
    h(2, 3) = 10.0;

    i(0, 0) = 2.0;
    i(0, 1) = 4.0;
    i(0, 2) = -2.0;
    i(0, 3) = 2.0;
    i(1, 0) = 0.0;
    i(1, 1) = 1.0;
    i(1, 2) = 1.0;
    i(1, 3) = 4.0;
    i(2, 0) = -2.0;
    i(2, 1) = -3.0;
    i(2, 2) = 7.0;
    i(2, 3) = 10.0;    


    v1(0) = 2.0;
    v1(1) = 4.0;
    v1(2) = -2.0;

    v2(0) = 2.0;
    v2(1) = 0.0;
    v2(2) = -2.0;

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
    
    double tmp5[4];
    tmp5[0] = 1.0;
    tmp5[1] = 1.0;
    tmp5[2] = 3.0;
    tmp5[3] = 1.0; 

    RMat2x2 tmp6(tmp5);
    CPPUNIT_ASSERT(tmp6 == m1);
   
    
}

void RectMatrixTestCase::testEqual() {
    CPPUNIT_ASSERT(m2 == m3);
    CPPUNIT_ASSERT(m2 != m3);
    
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

    //test matrix multiplication
    CPPUNIT_ASSERT(a * b == c);
    CPPUNIT_ASSERT(d  * e == f);
    CPPUNIT_ASSERT(g  * h == i);

    //test matrix vector multiplication
    CPPUNIT_ASSERT(g  * v1 == v2);

}

void RectMatrixTestCase::testDiv() {

    CPPUNIT_ASSERT(m1 / 1.0 == m1);

    CPPUNIT_ASSERT(m1 / 2.0 * 2.0 == m1);


}

void RectMatrixTestCase::testAccessEntries(){
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m1(1, 0), 3.0, OpenMD::NumericConstant::epsilon);
}

void RectMatrixTestCase::testRowColOperations() {
    Vec3 row;
    Vec3 col;
    RMat3x3 m;
    
    //test getRow
    row = e.getRow(0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(row[0], 2.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(row[1], 4.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(row[2], 1.0, OpenMD::NumericConstant::epsilon);
    //test setRow
    row[0] = 2.0;    
    row[1] = 4.0;    
    row[2] = 1.0;    
    m.setRow(0, row);
    row[0] = 0.0;    
    row[1] = 0.0;    
    row[2] = 3.0;    
    m.setRow(1, row);
    row[0] = 0.0;    
    row[1] = 6.0;    
    row[2] = 5.0;    
    m.setRow(2, row);
    CPPUNIT_ASSERT(m == e);
    
    //test getCol
    col = e.getColumn(1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(col[0], 4.0, OpenMD::NumericConstant::epsilon);    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(col[1], 0.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(col[2], 6.0, OpenMD::NumericConstant::epsilon);
    //test setCol
    col[0] = 2.0;    
    col[1] = 0.0;    
    col[2] = 0.0;    
    m.setColumn(0, col);
    col[0] = 4.0;    
    col[1] = 0.0;    
    col[2] = 6.0;    
    m.setColumn(1, col);
    col[0] = 1.0;    
    col[1] = 3.0;    
    col[2] = 5.0;    
    m.setColumn(2, col);
    CPPUNIT_ASSERT(m == e);

    //test swapRow
    RMat2x3 r;
    r(0, 0) = 0.0;
    r(0, 1) = 1.0;
    r(0, 2) = 0.0;
    r(1, 0) = 1.0;
    r(1, 1) = 0.0;
    r(1, 2) = 0.0;
    r.swapRow(0, 1);
    CPPUNIT_ASSERT(r == a);

    //test swapCol
    RMat3x3 s;
    s(0, 0) = 4.0;
    s(0, 1) = 2.0;
    s(0, 2) = 1.0;
    s(1, 0) = 0.0;    
    s(1, 1) = 0.0;
    s(1, 2) = 3.0;    
    s(2, 0) = 6.0;
    s(2, 1) = 0.0;
    s(2, 2) = 5.0;

    s.swapColumn(0, 1);
    CPPUNIT_ASSERT(s == e);

    double* p = s.getArrayPointer();

    p[0] = 2.0;
    p[1] = 4.0;
    p[2] = 1.0;
    p[3] = 0.0;
    p[4] = 6.0;
    p[5] = 5.0;    
    p[6] = 0.0;
    p[7] = 0.0;
    p[8] = 3.0;  

    CPPUNIT_ASSERT(s == f);    
}    

void RectMatrixTestCase::testOtherMemberFunctions(){
    //test transpose
    CPPUNIT_ASSERT((a.transpose()).transpose() == a);
    
    CPPUNIT_ASSERT(a.transpose() == b);

    //test getArray

    double tmp[4];
    m4.getArray(tmp);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp[0], -1.0, OpenMD::NumericConstant::epsilon);    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp[1], -1.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp[2], -3.0, OpenMD::NumericConstant::epsilon);    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp[3], -1.0, OpenMD::NumericConstant::epsilon);    
    
}
