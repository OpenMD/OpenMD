#include "utils/NextCombinationTestCase.hpp"
#include <iostream>
#include <algorithm>

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( NextCombinationTestCase);


void NextCombinationTestCase::testNextCombination() {
    std::vector<int> iv;
    std::vector<std::vector<int>::iterator> ic;
    std::vector<std::vector<int>::iterator>::iterator i; 
    iv.push_back(0);
    iv.push_back(1);
    iv.push_back(4);
    std::cout << std::endl;

    std::vector<std::vector<int> > results;
    while (next_combination(ic, iv.begin(), iv.end())) {
      std::vector<int> v;  
      for(i = ic.begin();i != ic.end(); ++i) {
         v.push_back(**i);  
      }       
      results.push_back(v);  
    }

    CPPUNIT_ASSERT(results.size() == 7);
    CPPUNIT_ASSERT(results[0][0] == 0  && results[0].size() == 1);
    CPPUNIT_ASSERT(results[1][0] == 1 && results[1].size() == 1);
    CPPUNIT_ASSERT(results[2][0] == 4  && results[2].size() == 1);
    CPPUNIT_ASSERT(results[3][0] == 0  && results[3][1] == 1  && results[3].size() == 2);
    CPPUNIT_ASSERT(results[4][0] == 0 && results[4][1] == 4 && results[4].size() == 2);
    CPPUNIT_ASSERT(results[5][0] == 1 && results[5][1] == 4  && results[5].size() == 2);
    CPPUNIT_ASSERT(results[6][0] == 0  && results[6][1] == 1 && results[6][2] == 4 && results[6].size() == 3);

    std::vector<std::string> sv;
    std::vector<std::vector<std::string>::iterator> sic;
    std::vector<std::vector<std::string>::iterator>::iterator j; 
    std::vector<std::vector<std::string> > resultStrings;
    std::vector<std::string> resultString;
    sv.push_back("H");
    sv.push_back("C");
    sv.push_back("N");

    while (replaceWildCard(sic, sv, resultString)) {   
      resultStrings.push_back(resultString);  
    }

    CPPUNIT_ASSERT(resultStrings.size() == 7);
    CPPUNIT_ASSERT(resultStrings[0][0] == "H"  && resultStrings[0][1] == "X" && resultStrings[0][2] == "X" && resultStrings[0].size() == 3);
    CPPUNIT_ASSERT(resultStrings[1][0] == "X"  && resultStrings[1][1] == "C" && resultStrings[1][2] == "X" && resultStrings[1].size() == 3);
    CPPUNIT_ASSERT(resultStrings[2][0] == "X"  && resultStrings[2][1] == "X" && resultStrings[2][2] == "N" && resultStrings[2].size() == 3);
    CPPUNIT_ASSERT(resultStrings[3][0] == "H"  && resultStrings[3][1] == "C" && resultStrings[3][2] == "X" && resultStrings[3].size() == 3);
    CPPUNIT_ASSERT(resultStrings[4][0] == "H"  && resultStrings[4][1] == "X" && resultStrings[4][2] == "N" && resultStrings[4].size() == 3);
    CPPUNIT_ASSERT(resultStrings[5][0] == "X"  && resultStrings[5][1] == "C" && resultStrings[5][2] == "N" && resultStrings[5].size() == 3);
    CPPUNIT_ASSERT(resultStrings[6][0] == "H"  && resultStrings[6][1] == "C" && resultStrings[6][2] == "N" && resultStrings[6].size() == 3);
    
}

