/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "SquareMatrixTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(SquareMatrixTestCase);

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

void SquareMatrixTestCase::testConstructor() {}

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
  CPPUNIT_ASSERT_DOUBLES_EQUAL(w1L[1], w1R[1],
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(w1L[2], w1R[2],
                               OpenMD::NumericConstant::epsilon);
}

void SquareMatrixTestCase::testTrace() {
  CPPUNIT_ASSERT_DOUBLES_EQUAL(identMat.trace(), 3.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(symMat.trace(), 4.0,
                               OpenMD::NumericConstant::epsilon);
}
void SquareMatrixTestCase::testIsSymmertric() {
  CPPUNIT_ASSERT(identMat.isSymmetric());
  CPPUNIT_ASSERT(symMat.isSymmetric());
}

void SquareMatrixTestCase::testIsOrthogonal() {
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
