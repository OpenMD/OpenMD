#include "utils/ReplaceWildCardTestCase.hpp"
#include <iostream>
#include <algorithm>
#include "utils/ReplaceWildCard.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ReplaceWildCardTestCase );

void ReplaceWildCardTestCase::testReplaceWildCard(){
    std::vector<std::vector<int> > sequences;

    sequences = ReplaceWildCardTestCase(3);
    CPPUNIT_ASSERT(sequences.size() == 7);
}    
