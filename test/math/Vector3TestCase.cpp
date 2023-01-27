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

#include "math/Vector3TestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(Vector3TestCase);

void Vector3TestCase::setUp() {
  zero[0] = 0.0;
  zero[1] = 0.0;
  zero[2] = 0.0;

  one[0] = 1.0;
  one[1] = 1.0;
  one[2] = 1.0;

  two[0] = 2.0;
  two[1] = 2.0;
  two[2] = 2.0;

  v1[0] = 1.0;
  v1[1] = 2.0;
  v1[2] = 3.0;

  v2[0] = 4.0;
  v2[1] = 1.0;
  v2[2] = 2.0;

  v3[0] = 1.0;
  v3[1] = 10.0;
  v3[2] = -7.0;
}

void Vector3TestCase::tearDown() {}

void Vector3TestCase::testConstructors() {
  double b[] = {2.9, 3.2, 1.2};
  Vector3d v(b);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v.x(), 2.9, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v.y(), 3.2, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v.z(), 1.2, OpenMD::NumericConstant::epsilon);
}

void Vector3TestCase::testArithmetic() {}

void Vector3TestCase::testOperators() {
  Vector3d tmp;
  // test +=

  // test +
  CPPUNIT_ASSERT(one + one == two);

  // test -=
  tmp = two;
  tmp -= one;
  CPPUNIT_ASSERT(tmp == one);

  // test -
  CPPUNIT_ASSERT(two - one == one);

  // test *=
  tmp = two;
  tmp *= 0.5;
  CPPUNIT_ASSERT(tmp == one);

  // test *
  CPPUNIT_ASSERT(two * 0.5 == one);
  CPPUNIT_ASSERT(0.5 * two == one);

  // test /=
  tmp = two;
  tmp *= 2.0;
  CPPUNIT_ASSERT(tmp == one * 4.0);

  // test /
  CPPUNIT_ASSERT(two / 2.0 == one);
  CPPUNIT_ASSERT(two / 4.0 == one * 0.5);
}

void Vector3TestCase::testAccessEntries() {
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v1.z(), 3.0, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v2.x(), 4.0, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v3.y(), 10.0, OpenMD::NumericConstant::epsilon);

  Vector3d tmp;
  tmp.x() = 78.01;
  tmp.y() = 21.0;
  tmp.z() = 133.12;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp.x(), 78.01,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp.y(), 21.0, OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp.z(), 133.12,
                               OpenMD::NumericConstant::epsilon);
}

void Vector3TestCase::testOtherTemplateFunctions() {
  // test cross
  CPPUNIT_ASSERT(cross(v1, v2) == v3);
  CPPUNIT_ASSERT(cross(one, two) == zero);
}
