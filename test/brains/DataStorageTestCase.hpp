#ifndef BRAINS_DATASTORAGETESTCASE_HPP
#define BRAINS_DATASTORAGETESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "brains/DataStorage.hpp"
 #include <vector>
using namespace OpenMD;

class DataStorageTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( DataStorageTestCase );
    CPPUNIT_TEST(testDataStorage);
    CPPUNIT_TEST_SUITE_END();

    public:
        void testDataStorage();
    private:

};


#endif //BRAINS_DATASTORAGETESTCASE_HPP

