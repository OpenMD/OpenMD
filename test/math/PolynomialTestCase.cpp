/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

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

