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

#include "math/VectorTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(VectorTestCase);

void VectorTestCase::setUp() {
  zero[0] = 0.0;
  zero[1] = 0.0;
  zero[2] = 0.0;
  zero[3] = 0.0;

  one[0] = 1.0;
  one[1] = 1.0;
  one[2] = 1.0;
  one[3] = 1.0;

  two[0] = 2.0;
  two[1] = 2.0;
  two[2] = 2.0;
  two[3] = 2.0;

  v1[0] = 3.0;
  v1[1] = 0.0;
  v1[2] = 2.0;
  v1[3] = 1.0;

  v2[0] = -3.0;
  v2[1] = 0.0;
  v2[2] = -2.0;
  v2[3] = -1.0;

  v3[0] = 4.0;
  v3[1] = 1.0;
  v3[2] = 3.0;
  v3[3] = 2.0;

  s1 = 1.0;
  s2 = 2.0;
}

void VectorTestCase::testConstructors() {
  Vec4 a0;
  Vec4 a1(1);
  Vec4 a2(2);

  CPPUNIT_ASSERT(a0 == zero);
  CPPUNIT_ASSERT(a1 == one);
  CPPUNIT_ASSERT(a2 == two);

  CPPUNIT_ASSERT(a1 != two);

  // test copy constructor
  Vec4 b1(v1);
  CPPUNIT_ASSERT(b1 == v1);

  // test operator =
  b1 = v2;
  CPPUNIT_ASSERT(b1 == v2);

  // test constructor from an array
  double tempArray[] = {1.0, 2.0, 5.0, 8.0};
  Vec4 tempV;
  tempV[0] = 1.0;
  tempV[1] = 2.0;
  tempV[2] = 5.0;
  tempV[3] = 8.0;

  Vec4 b2(tempV);
  CPPUNIT_ASSERT(b2 == tempArray);
}

void VectorTestCase::testArithmetic() {
  // test negate
  Vec4 a0 = v2;
  a0.negate();
  CPPUNIT_ASSERT(a0 == v1);

  Vec4 a1;
  a1.negate(v2);
  CPPUNIT_ASSERT(a1 == v1);

  // test add
  Vec4 a2;
  a2 = v1;
  a2.add(v2);
  CPPUNIT_ASSERT(a2 == zero);

  Vec4 a3;
  a3.add(v2, v3);
  CPPUNIT_ASSERT(a3 == one);

  // test sub
  Vec4 a4;
  a4 = two;
  a4.sub(one);
  CPPUNIT_ASSERT(a4 == one);

  Vec4 a5;
  a5.sub(two, one);
  CPPUNIT_ASSERT(a5 == one);

  // test mul
  Vec4 a6;
  a6 = one;
  a6.mul(2.0);
  CPPUNIT_ASSERT(a6 == two);

  Vec4 a7;
  a7.mul(one, 2.0);
  CPPUNIT_ASSERT(a7 == two);
  Vec4 a8;
  a7.mul(zero, 2.0);
  CPPUNIT_ASSERT(a7 == zero);

  // test div
  Vec4 a9;
  a9 = two;
  a9.div(2.0);
  CPPUNIT_ASSERT(a9 == one);

  Vec4 a10;
  a10.div(two, 2.0);
  CPPUNIT_ASSERT(a10 == one);
  Vec4 a11;
  a11.mul(zero, 2.0);
  CPPUNIT_ASSERT(a11 == zero);
}

void VectorTestCase::testOperators() {
  // test unary minus
  Vec4 a0 = v2;
  a0      = -a0;
  CPPUNIT_ASSERT(a0 == v1);

  // test add
  Vec4 a1;
  a1 = v1;
  a1 += v2;
  CPPUNIT_ASSERT(a1 == zero);

  Vec4 a2;
  a2 = v2 + v3;
  CPPUNIT_ASSERT(a2 == one);

  // test sub
  Vec4 a3;
  a3 = two;
  a3 -= one;
  CPPUNIT_ASSERT(a3 == one);

  Vec4 a4;
  a4 = two - one;
  CPPUNIT_ASSERT(a4 == one);

  Vec4 a5;
  a5.sub(two, one);
  CPPUNIT_ASSERT(a5 == one);

  // test mul
  Vec4 a6;
  a6 = one;
  a6 *= 2.0;
  CPPUNIT_ASSERT(a6 == two);

  Vec4 a7;
  a7 = one * 2.0;
  CPPUNIT_ASSERT(a7 == two);
  a7 = 2.0 * one;
  CPPUNIT_ASSERT(a7 == two);

  // test div
  Vec4 a8;
  a8 = two;
  a8 /= 2.0;
  CPPUNIT_ASSERT(a8 == one);
  a8 = two / 2.0;
  CPPUNIT_ASSERT(a8 == one);
}

void VectorTestCase::testAccessEntries() {
  // test [] operator

  CPPUNIT_ASSERT_DOUBLES_EQUAL(zero[0], 0.0, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(one[0], 1.0, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v3[0], 4.0, OpenMD::NumericConstant::epsilon);
  // test () operator
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v3(0), 4.0, OpenMD::NumericConstant::epsilon);

  Vec4 a1;
  double* pa1 = a1.getArrayPointer();

  pa1[0] = 4.0;
  pa1[1] = 1.0;
  pa1[2] = 3.0;
  pa1[3] = 2.0;

  CPPUNIT_ASSERT(a1 == v3);
}

void VectorTestCase::testOtherMemberFunctions() {
  // test length()
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zero.length(), 0.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(one.length(), 2.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v2.length(), sqrt(14.0),
                               OpenMD::NumericConstant::epsilon);

  // test lengthSquare()
  CPPUNIT_ASSERT_DOUBLES_EQUAL(zero.lengthSquare(), 0.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(one.lengthSquare(), 4.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v2.lengthSquare(), 14.0,
                               OpenMD::NumericConstant::epsilon);

  // test normalize()
  Vec4 a1 = one;
  Vec4 a2 = two;

  a1.normalize();
  a2.normalize();
  CPPUNIT_ASSERT(a1 == a2);

  // test isNormalized();
  CPPUNIT_ASSERT(a1.isNormalized());
  CPPUNIT_ASSERT(!one.isNormalized());

  // test getArray
  double tempV[4];
  v3.getArray(tempV);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tempV[0], v3[0],
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tempV[1], v3[1],
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tempV[2], v3[2],
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tempV[3], v3[3],
                               OpenMD::NumericConstant::epsilon);
}
void VectorTestCase::testOtherTemplateFunctions() {
  // test dot
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dot(one, two), 8.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dot(v1, v3), 20.0,
                               OpenMD::NumericConstant::epsilon);

  // test distance
  CPPUNIT_ASSERT_DOUBLES_EQUAL(distance(one, two), 2.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(distance(v1, v2), sqrt(56.0),
                               OpenMD::NumericConstant::epsilon);

  // test distanceSquare
  CPPUNIT_ASSERT_DOUBLES_EQUAL(distanceSquare(one, two), 4.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(distanceSquare(v1, v2), 56,
                               OpenMD::NumericConstant::epsilon);
}
