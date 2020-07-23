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
