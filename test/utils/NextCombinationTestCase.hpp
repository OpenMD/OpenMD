#ifndef TEST_NEXTCOMBINATIONTESTCASE_HPP
#define TEST_NEXTCOMBINATIONTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>

using namespace oopse;

class NextCombinationTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( NextCombinationTestCase );
    CPPUNIT_TEST(testNextCombination);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testNextCombination(); 
};


#endif //TEST_NEXTCOMBINATIONTESTCASE_HPP