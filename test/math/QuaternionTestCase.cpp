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

#include "math/QuaternionTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( QuaternionTestCase );


void QuaternionTestCase::setUp(){
    q1[0] = 1.0;
    q1[1] = 0.0;
    q1[2] = 0.0;
    q1[3] = 0.0;

    q2[0] = 0.0;
    q2[1] = 0.6;
    q2[2] = 0.0;
    q2[3] = 0.8;

    q3[0] = 0.5;
    q3[1] = 0.5;
    q3[2] = 0.5;
    q3[3] = 0.5;

    q4[0] = 0.0;
    q4[1] = -0.6;
    q4[2] = 0.0;
    q4[3] = -0.8;
        
}

void QuaternionTestCase::testConstructors(){
    Vec4 v;
    v[0] = 0.0;
    v[1] = 0.6;
    v[2] = 0.0;
    v[3] = 0.8;

    //construct from a vector
    Quat4d tmp1(v);

    CPPUNIT_ASSERT(tmp1 == q2);

    //copy constructor
    Quat4d tmp2(q1);
    CPPUNIT_ASSERT(tmp2 == q1);

    //copy assignment
    Quat4d tmp3;
    tmp3 = q4;
    CPPUNIT_ASSERT(tmp3 == q4);

    //construct from w, x, y, z
    Quat4d tmp4(0.0, -0.6, 0.0, -0.8);

    CPPUNIT_ASSERT(tmp4 == q4);
    
}

void QuaternionTestCase::testArithmetic(){
    Quat4d qd3(1, 1, 2, 1);
    Quat4d qd4(1, 1, 1, 2);
    
    //test mul
    Quat4d qda3L = qd3;
    qda3L.mul(2.0);
    CPPUNIT_ASSERT_EQUAL( Quat4d(2,2,4,2) , qda3L );
    
    Quat4d qda5R = 2.0  * qd4;
    qda5R.div(qd4);
    CPPUNIT_ASSERT_EQUAL( Quat4d(2.0, 0.0, 0.0, 0.0) , qda5R );    

}

void QuaternionTestCase::testOperators(){

    Quat4d qd0(1, 1, 1, 1);
    Quat4d qd1(2, 1, 1, 1);
    Quat4d qd2(1, 2, 1, 1);
    Quat4d qd3(1, 1, 2, 1);
    Quat4d qd4(1, 1, 1, 2);


    Quat4d qda2Q = qd2 - qd0;
    Quat4d qda3L = qd3 * 2.0;
    Quat4d qda3R = 2.0 * qd3;
    Quat4d qda3Q = qd3 * qd0;
    Quat4d qda4L = qd4 / 2.0;
    Quat4d qda4R = 2.0 / qd4;
    Quat4d qda4Q = qd4 / qd0;
    Quat4d qda5L = qda4L * 2.0;
    Quat4d qda5R = qda4R * qd4;
    Quat4d qda5Q = qda4Q * qd0;    
    
    CPPUNIT_ASSERT_EQUAL( Quat4d(0,1,0,0) , qda2Q );
    CPPUNIT_ASSERT_EQUAL( Quat4d(2,2,4,2) , qda3L );
    CPPUNIT_ASSERT_EQUAL( Quat4d(2,2,4,2) , qda3R );
    CPPUNIT_ASSERT_EQUAL( Quat4d(-3,3,3,1) , qda3Q );
    CPPUNIT_ASSERT_EQUAL( Quat4d(0.5,0.5,0.5,1) , qda4L );
    //CPPUNIT_ASSERT_EQUAL( 2.0 * qd4.inverse() , qda4R );
    CPPUNIT_ASSERT_EQUAL( Quat4d(1.25,0.25,-0.25,0.25) , qda4Q );
    CPPUNIT_ASSERT_EQUAL( qd4 , qda5L );
    CPPUNIT_ASSERT_EQUAL( Quat4d(2.0, 0.0, 0.0, 0.0) , qda5R );
    CPPUNIT_ASSERT_EQUAL( qd4 , qda5Q );

    Quat4d quat0(1, 1, 1, 1);

    Quat4d quat1a = quat0;
    Quat4d quat1s = quat0;
    Quat4d quat1p = quat0;
    Quat4d quat1d = quat0;
    
    quat1a += 2.0;
    quat1s -= 2.0;
    quat1p *= 2.0;
    quat1d /= 2.0;
    
    CPPUNIT_ASSERT_EQUAL( Quat4d(3,1,1,1) , quat1a );
    CPPUNIT_ASSERT_EQUAL( Quat4d(-1,1,1,1) , quat1s );
    CPPUNIT_ASSERT_EQUAL( Quat4d(2,2,2,2) , quat1p );
    CPPUNIT_ASSERT_EQUAL( Quat4d(0.5,0.5,0.5,0.5) , quat1d );

}

void QuaternionTestCase::testAccessEntries(){
    CPPUNIT_ASSERT_DOUBLES_EQUAL(q1.w(), 1.0, OpenMD::NumericConstants::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(q2.x(), 0.6, OpenMD::NumericConstants::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(q3.y(), 0.5, OpenMD::NumericConstants::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(q4.z(), -0.8, OpenMD::NumericConstants::epsilon);

    Quat4d tmp;

    tmp.w() = 0.0;
    tmp.x() = -0.6;
    tmp.y() = 0.0;
    tmp.z() = -0.8;
    
    CPPUNIT_ASSERT(tmp == q4);


}

void QuaternionTestCase::testOtherMemberFunctions(){        
    //test inverse
    CPPUNIT_ASSERT(q2.inverse() == q4);

    //test conjugate
    CPPUNIT_ASSERT(q2 == q4.conjugate());

    //test toRotationMatrix3
    SquareMatrix<double, 3> rotMat;

    rotMat(0,0) = -0.28;
    rotMat(0,1) =  0;
    rotMat(0,2) =  0.96;

    rotMat(1,0) =  0.0;
    rotMat(1,1) = -1.0;
    rotMat(1,2) =  0.0;

    rotMat(2,0) = 0.96;
    rotMat(2,1) = 0;
    rotMat(2,2) = 0.28;

    CPPUNIT_ASSERT(rotMat == q4.toRotationMatrix3());

}
