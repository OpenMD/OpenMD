#ifndef TEST_PROPERTYMAPTESTCASE_HPP
#define TEST_PROPERTYMAPTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "utils/GenericFactory.hpp"

using namespace oopse;

class PropertyMapTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( GenericFactoryTestCase );
    CPPUNIT_TEST(testPropertyMap);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testPropertyMap(); 
};


#endif //TEST_PROPERTYMAPTESTCASE_HPP
