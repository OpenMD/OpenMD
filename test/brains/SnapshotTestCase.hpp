#ifndef BRAINS_SNAPSHOTTESTCASE_HPP
#define BRAINS_SNAPSHOTTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "brains/Snapshot.hpp"
 
using namespace oopse;

class SnapshotTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( SnapshotTestCase );
    CPPUNIT_TEST(testConstructors);

    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();
        virtual void tearDown();
        void testConstructors();

    private:

};


#endif
