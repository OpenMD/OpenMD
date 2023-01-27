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

#include "utils/GenericDataTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(GenericDataTestCase);

void GenericDataTestCase::testGenericData() {
  // test constructor
  GenericData gd1;

  GenericData gd2("OpenMD_Generic");

  CPPUNIT_ASSERT(gd1.getID() == "UndefinedGenericData");
  CPPUNIT_ASSERT(gd2.getID() == "OpenMD_Generic");

  gd1.setID("Dummy");
  CPPUNIT_ASSERT(gd1.getID() == "Dummy");
}

void GenericDataTestCase::testSimpleTypeData() {
  // test IntGenericData
  BoolGenericData b1("Dummy_Bool");
  b1.setData(true);
  CPPUNIT_ASSERT(b1.getData());
  CPPUNIT_ASSERT(b1.getID() == "Dummy_Bool");

  b1.setData(10.0);
  CPPUNIT_ASSERT(b1.getData());

  b1.setData(0.0);
  CPPUNIT_ASSERT(!b1.getData());

  // test IntGenericData
  IntGenericData i1("Dummy_Int");
  i1.setData(10);
  CPPUNIT_ASSERT_EQUAL(i1.getData(), 10);
  CPPUNIT_ASSERT(i1.getID() == "Dummy_Int");

  IntGenericData i2("Dummy_Int");
  i2.setData(10.000);
  CPPUNIT_ASSERT_EQUAL(i1.getData(), 10);

  // test FloatGenericData
  FloatGenericData f1("Dummy_Float");
  f1.setData(10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(f1.getData(), 10, 0.000001);
  CPPUNIT_ASSERT(f1.getID() == "Dummy_Float");

  // test DoubleGenericData
  DoubleGenericData d1("Dummy_Double");
  d1.setData(232);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(d1.getData(), 232, 0.000001);
  CPPUNIT_ASSERT(d1.getID() == "Dummy_Double");

  // test StringGenericData
  StringGenericData s1("Dummy_String");
  s1.setData("Hello World");
  CPPUNIT_ASSERT(s1.getData() == "Hello World");
  CPPUNIT_ASSERT(s1.getID() == "Dummy_String");

  DoubleGenericData d2("Dummy_Double");
  d2 = f1;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(d2.getData(), 10, 0.000001);

  FloatGenericData f2("Dummy_Float");
  f2 = d1;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(d1.getData(), 232, 0.000001);

  // test getData (return reference)
  f2.getData() = 0.004;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(f2.getData(), 0.004, 0.000001);
}

void GenericDataTestCase::testSTLContainerTypeData() {
  // test IntVectorGenericData
  IntVectorGenericData iv1("IntVector1");
  CPPUNIT_ASSERT(iv1.getID() == "IntVector1");
  iv1.setID("Dummy_Test");
  CPPUNIT_ASSERT(iv1.getID() == "Dummy_Test");

  iv1.push_back(2.30);
  iv1.push_back(1);
  iv1.push_back(324);

  CPPUNIT_ASSERT_EQUAL(iv1[0], 2);
  CPPUNIT_ASSERT_EQUAL(iv1[1], 1);
  CPPUNIT_ASSERT_EQUAL(iv1[2], 324);

  IntVectorGenericData iv2("IntVector2");

  iv2.push_back(1);
  iv2.push_back(3);
  iv2.push_back(5);
  iv2.push_back(7);

  iv1 = iv2;
  CPPUNIT_ASSERT(iv1.getID() == "Dummy_Test");
  CPPUNIT_ASSERT(iv1.size() == 4);
  CPPUNIT_ASSERT_EQUAL(iv1[0], 1);
  CPPUNIT_ASSERT_EQUAL(iv1[1], 3);
  CPPUNIT_ASSERT_EQUAL(iv1[2], 5);
  CPPUNIT_ASSERT_EQUAL(iv1[3], 7);

  // test FloatVectorGenericData
  FloatVectorGenericData fv2("FloatVector2");

  fv2.push_back(251.21);
  fv2.push_back(42.90);
}
