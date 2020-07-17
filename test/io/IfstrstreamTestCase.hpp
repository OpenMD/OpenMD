#ifndef TEST_IFSTRSTREAMTESTCASE_HPP
#define TEST_IFSTRSTREAMTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>



class IfstrstreamTestCase : public CPPUNIT_NS::TestFixture {

    CPPUNIT_TEST_SUITE( IfstrstreamTestCase);    

#ifndef IS_MPI    
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testOpen);
    CPPUNIT_TEST(testIs_open);
    CPPUNIT_TEST(testClose);
#else
    CPPUNIT_TEST(testPrimaryConstructor);
    CPPUNIT_TEST(testPrimaryOpen);
    CPPUNIT_TEST(testPrimaryIs_open);
    CPPUNIT_TEST(testPrimaryClose);
    CPPUNIT_TEST(testSecondaryConstructor);
    CPPUNIT_TEST(testSecondaryOpen);
    CPPUNIT_TEST(testSecondaryIs_open);
    CPPUNIT_TEST(testSecondaryClose);
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
        void testPrimaryConstructor();
        void testPrimaryOpen();
        void testPrimaryIs_open();
        void testPrimaryClose();
        void testSecondaryConstructor();
        void testSecondaryOpen();
        void testSecondaryIs_open();
        void testSecondaryClose();                
#endif
};

#endif //TEST_IFSTRSTREAMTESTCASE_HPP
