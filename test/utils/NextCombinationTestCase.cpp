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
    
}
