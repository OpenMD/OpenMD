#include "math/PolynomialTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PolynomialTestCase );


void PolynomialTestCase::testPolynomial(){
    DoublePolynomial dp1;

    dp1.setCoefficient(1, 3.0);
    dp1.setCoefficient(2, 1.0);
    dp1.setCoefficient(3, 2.0);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp1.getCoefficient(1), 3.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp1.getCoefficient(2), 1.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp1.getCoefficient(3), 2.0, 0.000001);

    CPPUNIT_ASSERT(dp1.size() == 3);
    CPPUNIT_ASSERT(dp1.begin()->first == 1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp1.begin()->second ,3.0, 0.000001);
    

    dp1.addCoefficient(2, 1.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp1.getCoefficient(2), 2.5, 0.000001);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dp1.evaluate(0.0), 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.5, dp1.evaluate(1.0), 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(32.0, dp1.evaluate(2.0), 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.5, dp1.evaluate(-1.0), 0.000001);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, dp1.evaluateDerivative(1.0), 0.000001);


    //dp2 = x^2 + 3x +2
    DoublePolynomial dp2;
    dp2.setCoefficient(0, 2.0);
    dp2.setCoefficient(1, 3.0);
    dp2.setCoefficient(2, 1.0);

    //dp3 = x^3 + 2x
    DoublePolynomial dp3;
    dp3.setCoefficient(1, 2.0);
    dp3.setCoefficient(3, 1.0);

    //dp4 = x^3 + x^2 +5x + 2
    DoublePolynomial dp4 = dp2 + dp3;
    CPPUNIT_ASSERT(dp4.size() == 4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp4.getCoefficient(0), 2.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp4.getCoefficient(1), 5.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp4.getCoefficient(2), 1.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp4.getCoefficient(3), 1.0, 0.000001);
    

    
    //dp5 = -x^3 + x^2 +x + 2
    DoublePolynomial dp5 = dp2 - dp3;
    CPPUNIT_ASSERT(dp5.size() == 4);    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp5.getCoefficient(0), 2.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp5.getCoefficient(1), 1.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp5.getCoefficient(2), 1.0, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp5.getCoefficient(3), -1.0, 0.000001);


    //dp6 = x^5 + 3x^4 + 4x^3 + 6x^2 + 4x
    DoublePolynomial dp6 = dp2 * dp3;
    CPPUNIT_ASSERT(dp6.size() == 5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp6.getCoefficient(1), 4, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp6.getCoefficient(2), 6, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp6.getCoefficient(3), 4, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp6.getCoefficient(4), 3, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dp6.getCoefficient(5), 1, 0.000001);
    
}

