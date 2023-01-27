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

#include "utils/NextCombinationTestCase.hpp"

#include <algorithm>
#include <iostream>

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(NextCombinationTestCase);

void NextCombinationTestCase::testNextCombination() {
  std::vector<int> iv;
  std::vector<std::vector<int>::iterator> ic;
  std::vector<std::vector<int>::iterator>::iterator i;
  iv.push_back(0);
  iv.push_back(1);
  iv.push_back(4);
  std::cout << std::endl;

  std::vector<std::vector<int>> results;
  while (next_combination(ic, iv.begin(), iv.end())) {
    std::vector<int> v;
    for (i = ic.begin(); i != ic.end(); ++i) {
      v.push_back(**i);
    }
    results.push_back(v);
  }

  CPPUNIT_ASSERT(results.size() == 7);
  CPPUNIT_ASSERT(results[0][0] == 0 && results[0].size() == 1);
  CPPUNIT_ASSERT(results[1][0] == 1 && results[1].size() == 1);
  CPPUNIT_ASSERT(results[2][0] == 4 && results[2].size() == 1);
  CPPUNIT_ASSERT(results[3][0] == 0 && results[3][1] == 1 &&
                 results[3].size() == 2);
  CPPUNIT_ASSERT(results[4][0] == 0 && results[4][1] == 4 &&
                 results[4].size() == 2);
  CPPUNIT_ASSERT(results[5][0] == 1 && results[5][1] == 4 &&
                 results[5].size() == 2);
  CPPUNIT_ASSERT(results[6][0] == 0 && results[6][1] == 1 &&
                 results[6][2] == 4 && results[6].size() == 3);

  std::vector<std::string> sv;
  std::vector<std::vector<std::string>::iterator> sic;
  std::vector<std::vector<std::string>::iterator>::iterator j;
  std::vector<std::vector<std::string>> resultStrings;
  std::vector<std::string> resultString;
  sv.push_back("H");
  sv.push_back("C");
  sv.push_back("N");

  while (replaceWithWildCard(sic, sv, resultString)) {
    resultStrings.push_back(resultString);
  }

  CPPUNIT_ASSERT(resultStrings.size() == 7);
  CPPUNIT_ASSERT(resultStrings[0][0] == "H" && resultStrings[0][1] == "X" &&
                 resultStrings[0][2] == "X" && resultStrings[0].size() == 3);
  CPPUNIT_ASSERT(resultStrings[1][0] == "X" && resultStrings[1][1] == "C" &&
                 resultStrings[1][2] == "X" && resultStrings[1].size() == 3);
  CPPUNIT_ASSERT(resultStrings[2][0] == "X" && resultStrings[2][1] == "X" &&
                 resultStrings[2][2] == "N" && resultStrings[2].size() == 3);
  CPPUNIT_ASSERT(resultStrings[3][0] == "H" && resultStrings[3][1] == "C" &&
                 resultStrings[3][2] == "X" && resultStrings[3].size() == 3);
  CPPUNIT_ASSERT(resultStrings[4][0] == "H" && resultStrings[4][1] == "X" &&
                 resultStrings[4][2] == "N" && resultStrings[4].size() == 3);
  CPPUNIT_ASSERT(resultStrings[5][0] == "X" && resultStrings[5][1] == "C" &&
                 resultStrings[5][2] == "N" && resultStrings[5].size() == 3);
  CPPUNIT_ASSERT(resultStrings[6][0] == "H" && resultStrings[6][1] == "C" &&
                 resultStrings[6][2] == "N" && resultStrings[6].size() == 3);
}
