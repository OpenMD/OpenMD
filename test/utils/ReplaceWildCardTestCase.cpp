#include "utils/ReplaceWildCardTestCase.hpp"
#include <iostream>
#include <algorithm>
#include "utils/ReplaceWildCard.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ReplaceWildCardTestCase );

void ReplaceWildCardTestCase::testReplaceWildCard(){
    std::vector<std::vector<int> > sequences;

    sequences = ReplaceAll(3);   
    CPPUNIT_ASSERT(sequences.size() == 8);
    
    CPPUNIT_ASSERT(sequences[0][0] == 0  && sequences[0][1] == 1  && sequences[0][2] == 2);
    CPPUNIT_ASSERT(sequences[1][0] == -1 && sequences[1][1] == 1  && sequences[1][2] == 2);
    CPPUNIT_ASSERT(sequences[2][0] == 0  && sequences[2][1] == -1 && sequences[2][2] == 2);
    CPPUNIT_ASSERT(sequences[3][0] == 0  && sequences[3][1] == 1  && sequences[3][2] == -1);
    CPPUNIT_ASSERT(sequences[4][0] == -1 && sequences[4][1] == -1 && sequences[4][2] == 2);
    CPPUNIT_ASSERT(sequences[5][0] == -1 && sequences[5][1] == 1  && sequences[5][2] == -1);
    CPPUNIT_ASSERT(sequences[6][0] == 0  && sequences[6][1] == -1 && sequences[6][2] == -1);
    CPPUNIT_ASSERT(sequences[7][0] == -1 && sequences[7][1] == -1 && sequences[7][2] == -1);

}    
