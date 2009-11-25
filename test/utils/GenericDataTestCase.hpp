#ifndef TEST_GENERICDATATESTCASE_HPP
#define TEST_GENERICDATATESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "utils/GenericData.hpp"

using namespace OpenMD;

class GenericDataTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( GenericDataTestCase );
    CPPUNIT_TEST(testGenericData );
    CPPUNIT_TEST(testSimpleTypeData);
    CPPUNIT_TEST(testSTLContainerTypeData);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testGenericData(); 
        void testSimpleTypeData();
        void testSTLContainerTypeData();
};


#endif //TEST_GENERICDATATESTCASE_HPP

