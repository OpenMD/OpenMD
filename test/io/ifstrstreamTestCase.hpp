#ifndef TEST_IFSTRSTREAMTESTCASE_HPP
#define TEST_IFSTRSTREAMTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/SquareMatrix.hpp"

using namespace oopse;

typedef SquareMatrix<double, 3> SMat3x3;

class IfstrstreamTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( IfstrstreamTestCase );
#ifndef IS_MPI    
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testOpen);
    CPPUNIT_TEST(testIs_open);
    CPPUNIT_TEST(testClose);
#else
    CPPUNIT_TEST(testMasterConstructor);
    CPPUNIT_TEST(testMasterOpen);
    CPPUNIT_TEST(testMasterIs_open);
    CPPUNIT_TEST(testMasterClose);
    CPPUNIT_TEST(testSlaveConstructor);
    CPPUNIT_TEST(testSlaveOpen);
    CPPUNIT_TEST(testSlaveIs_open);
    CPPUNIT_TEST(testSlaveClose);
    CPPUNIT_TEST_SUITE_END();
#endif

    public:
        virtual void setUp();
        virtual void tearDown();

#ifndef IS_MPI
        void testConstrutor();
        void testOpen();
        void testIs_open();
        void testClose();
        
#else
        void testMasterConstrutor();
        void testMasterOpen();
        void testMasterIs_open();
        void testMasterClose();
        void testSlaveConstrutor();
        void testSlaveOpen();
        void testSlaveIs_open();
        void testSlaveClose();                
#endif
};

#endif //TEST_IFSTRSTREAMTESTCASE_HPP
