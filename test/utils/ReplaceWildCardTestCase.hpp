#ifndef TEST_REPLACEWILDCARDTESTCASE_HPP
#define TEST_REPLACEWILDCARDTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "utils/ReplaceWildCard.hpp"

using namespace oopse;

class ReplaceWildCardTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( ReplaceWildCardTestCase );
    CPPUNIT_TEST(testReplaceWildCard);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testReplaceWildCard(); 
};


#endif //TEST_REPLACEWILDCARDTESTCASE_HPP

