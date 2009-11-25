#ifndef TEST_POLYNOMIALTESTCASE_HPP
#define TEST_POLYNOMIALTESTCASE_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "math/Polynomial.hpp"
 
using namespace OpenMD;

class PolynomialTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( PolynomialTestCase );
    CPPUNIT_TEST(testPolynomial);

    CPPUNIT_TEST_SUITE_END();

    public:

        void testPolynomial();
        
    private:

        
};


#endif //TEST_POLYNOMIALTESTCASE_HPP


