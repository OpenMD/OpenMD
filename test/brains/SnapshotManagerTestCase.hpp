#ifndef BRAINS_SNAPSHOTMANAGERTESTCASE_HPP
#define BRAINS_SNAPSHOTMANAGERTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "brains/SnapshotManager.hpp"
 
using namespace oopse;

class SnapshotManagerTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( SnapshotManagerTestCase );
    CPPUNIT_TEST(testConstructors);

    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();
        virtual void tearDown();
        void testConstructors();

    private:

};


#endif //BRAINS_SNAPSHOTMANAGERTESTCASE_HPP