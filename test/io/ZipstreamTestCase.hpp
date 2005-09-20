#ifndef TEST_IFSTRSTREAMTESTCASE_HPP
#define TEST_IFSTRSTREAMTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>



class ZipstreamTeseCase : public CPPUNIT_NS::TestFixture {

    CPPUNIT_TEST_SUITE( ZipstreamTeseCase);    


    CPPUNIT_TEST(test_buffer_to_buffer);
    CPPUNIT_TEST(test_wbuffer_to_wbuffer);
    CPPUNIT_TEST(test_string_string);
    CPPUNIT_TEST(test_wstring_wstring);
    CPPUNIT_TEST(test_file_file);
    CPPUNIT_TEST(test_crc);
    CPPUNIT_TEST_SUITE_END();

    public:


        void test_buffer_to_buffer();
        void test_wbuffer_to_wbuffer();
        void test_string_string();
        void test_wstring_wstring();
        void test_file_file();
        void test_crc();

};

#endif //TEST_IFSTRSTREAMTESTCASE_HPP
