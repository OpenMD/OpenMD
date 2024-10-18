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

#include "utils/StringTokenizerTestCase.hpp"

#include <algorithm>
#include <iostream>

#include "utils/StringTokenizer.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(StringTokenizerTestCase);

void StringTokenizerTestCase::testStringTokenizer() {
  std::string str1 = "    \t  Hello \r World \n";
  StringTokenizer tokenizer1(str1);

  CPPUNIT_ASSERT(tokenizer1.getOriginal() == str1);

  CPPUNIT_ASSERT(tokenizer1.countTokens() == 2);

  std::string token1 = tokenizer1.nextToken();
  std::string token2 = tokenizer1.nextToken();

  CPPUNIT_ASSERT(token1 == "Hello");
  CPPUNIT_ASSERT(token2 == "World");

  // test reading and converting tokens to other data type
  std::string str2 = "1991.2\t129\t1e2 1 OpenMD\n";
  StringTokenizer tokenizer2(str2);

  CPPUNIT_ASSERT(tokenizer2.countTokens() == 5);

  float floatVal        = tokenizer2.nextTokenAsFloat();
  int intVal            = tokenizer2.nextTokenAsInt();
  double doubleVal      = tokenizer2.nextTokenAsDouble();
  bool boolVal          = tokenizer2.nextTokenAsBool();
  std::string stringVal = tokenizer2.nextToken();

  CPPUNIT_ASSERT_DOUBLES_EQUAL(floatVal, 1991.2, 0.0001);
  CPPUNIT_ASSERT(intVal == 129);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(doubleVal, 100, 0.0001);
  CPPUNIT_ASSERT(boolVal);
  CPPUNIT_ASSERT(stringVal == "OpenMD");

  CPPUNIT_ASSERT(!tokenizer2.hasMoreTokens());

  // test peekNextToken and using different delimeters
  StringTokenizer tokenizer3(str2, " \n");

  CPPUNIT_ASSERT(tokenizer3.getDelimiters() == " \n");

  CPPUNIT_ASSERT(tokenizer3.countTokens() == 3);

  CPPUNIT_ASSERT(tokenizer3.peekNextToken() == "1991.2\t129\t1e2");

  CPPUNIT_ASSERT(tokenizer3.countTokens() == 3);

  CPPUNIT_ASSERT(tokenizer3.nextToken() == "1991.2\t129\t1e2");

  CPPUNIT_ASSERT(tokenizer3.countTokens() == 2);

  // test return tokens
  StringTokenizer tokenizer4(str2, " \n", true);
  CPPUNIT_ASSERT(tokenizer4.countTokens() == 6);
  CPPUNIT_ASSERT(tokenizer4.nextToken() == "1991.2\t129\t1e2");
  CPPUNIT_ASSERT(tokenizer4.nextToken() == " ");
  CPPUNIT_ASSERT(tokenizer4.nextToken() == "1");
  CPPUNIT_ASSERT(tokenizer4.nextToken() == " ");
  CPPUNIT_ASSERT(tokenizer4.nextToken() == "OpenMD");
  CPPUNIT_ASSERT(tokenizer4.nextToken() == "\n");
}
