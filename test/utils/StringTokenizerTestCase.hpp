#ifndef TEST_STRINGTOKENIZERTESTCASE_HPP
#define TEST_STRINGTOKENIZERTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "utils/StringTokenizer.hpp"

using namespace OpenMD;

class StringTokenizerTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( StringTokenizerTestCase );
    CPPUNIT_TEST(testStringTokenizer);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testStringTokenizer(); 
};


#endif //TEST_PROPERTYMAPTESTCASE_HPP

