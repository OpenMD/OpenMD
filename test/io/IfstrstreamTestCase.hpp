#ifndef TEST_IFSTRSTREAMTESTCASE_HPP
#define TEST_IFSTRSTREAMTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>

//using namespace oopse;

class IfstrstreamTestCase : public CPPUNIT_NS::TestFixture {

    CPPUNIT_TEST_SUITE( IfstrstreamTestCase);    

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
#endif

    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();
        virtual void tearDown();

#ifndef IS_MPI
        void testConstructor();
        void testOpen();
        void testIs_open();
        void testClose();
        
#else
        void testMasterConstructor();
        void testMasterOpen();
        void testMasterIs_open();
        void testMasterClose();
        void testSlaveConstructor();
        void testSlaveOpen();
        void testSlaveIs_open();
        void testSlaveClose();                
#endif
};

#endif //TEST_IFSTRSTREAMTESTCASE_HPP
