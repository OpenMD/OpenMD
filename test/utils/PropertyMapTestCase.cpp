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

#include "utils/PropertyMapTestCase.hpp"

#include <algorithm>
#include <iostream>
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(PropertyMapTestCase);

void PropertyMapTestCase::testPropertyMap() {
  PropertyMap props;

  // test addProperty
  BoolGenericData* b0 = new BoolGenericData("BoolData");
  b0->setData(false);
  props.addProperty(b0);
  CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b0);

  BoolGenericData* b1 = new BoolGenericData("BoolData");
  b1->setData(true);
  props.addProperty(b1);
  CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b1);

  IntGenericData* i1 = new IntGenericData("IntData");
  i1->setData(89);
  props.addProperty(i1);

  FloatGenericData* f1 = new FloatGenericData("FloatData");
  f1->setData(49.328);
  props.addProperty(f1);

  DoubleGenericData* d1 = new DoubleGenericData("DoubleData");
  d1->setData(95.1933432);
  props.addProperty(d1);

  StringGenericData* s1 = new StringGenericData("StringData");
  s1->setData("Hello");
  props.addProperty(s1);

  IntVectorGenericData* iv1 = new IntVectorGenericData("IntVector");
  iv1->push_back(2);
  iv1->push_back(1);
  iv1->push_back(324);
  props.addProperty(iv1);

  // test getPropertyNames
  std::vector<std::string> propNames = props.getPropertyNames();

  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") !=
                 propNames.end());

  // test getProperties
  std::vector<GenericData*> propPointers = props.getProperties();
  CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), b1) !=
                 propPointers.end());
  CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), i1) !=
                 propPointers.end());
  CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), f1) !=
                 propPointers.end());
  CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), d1) !=
                 propPointers.end());
  CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), s1) !=
                 propPointers.end());
  CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), iv1) !=
                 propPointers.end());

  // test getPropertyByName
  CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b1);
  CPPUNIT_ASSERT(props.getPropertyByName("IntData") == i1);
  CPPUNIT_ASSERT(props.getPropertyByName("FloatData") == f1);
  CPPUNIT_ASSERT(props.getPropertyByName("DoubleData") == d1);
  CPPUNIT_ASSERT(props.getPropertyByName("StringData") == s1);
  CPPUNIT_ASSERT(props.getPropertyByName("IntVector") == iv1);

  CPPUNIT_ASSERT(b1->getData() == true);
  CPPUNIT_ASSERT(i1->getData() == 89);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getData(), 49.328, 0.000001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(d1->getData(), 95.1933432, 0.000001);
  CPPUNIT_ASSERT(s1->getData() == "Hello");
  CPPUNIT_ASSERT_EQUAL((*iv1)[0], 2);
  CPPUNIT_ASSERT_EQUAL((*iv1)[1], 1);
  CPPUNIT_ASSERT_EQUAL((*iv1)[2], 324);

  // test removeProperty
  props.removeProperty("DoubleData");
  props.removeProperty("FloatData");
  props.removeProperty("IntVector");
  propNames = props.getPropertyNames();
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") !=
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") ==
                 propNames.end());

  // test clearProperties
  props.clearProperties();
  propNames = props.getPropertyNames();
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") ==
                 propNames.end());
  CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") ==
                 propNames.end());
}
