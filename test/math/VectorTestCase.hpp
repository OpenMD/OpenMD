#ifndef TEST_VECTORTESTCASE_HPP
#define TEST_VECTORTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/Vector.hpp"
 
using namespace oopse;

typedef Vector<double, 4> Vec4;

class VectorTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( VectorTestCase );
    CPPUNIT_TEST(testConstructors);
    CPPUNIT_TEST(testArithmetic);
    CPPUNIT_TEST(testAccessEntries);
    CPPUNIT_TEST(testOtherTemplateFunctions);
    CPPUNIT_TEST_SUITE_END();

    public:
        virtual void setUp();

        void testConstructors();
        void testArithmetic();
        void testOperators();
        void testAccessEntries();
        void testOtherTemplateFunctions();

    private:
        Vec4 zero;
        Vec4 one;
        Vec4 two;
        Vec4 v1;
        Vec4 v2;
        Vec4 v3;

        double s1;
        double s2;
        
   
};

#endif //TEST_VECTORTESTCASE_HPP
