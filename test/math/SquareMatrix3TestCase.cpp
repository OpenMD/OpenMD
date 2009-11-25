#include "math/SquareMatrix3TestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SquareMatrix3TestCase );


void SquareMatrix3TestCase::testSetupRotationMatrix(){
    //test setupRotationMatrix by quaternion

    RotMat3x3d m1L(0.0, -0.6, 0.0, -0.8);
    RotMat3x3d m1R;
    m1R(0,0) = -0.28;
    m1R(0,1) =  0;
    m1R(0,2) =  0.96;
    m1R(1,0) =  0.0;
    m1R(1,1) = -1.0;
    m1R(1,2) =  0.0;
    m1R(2,0) = 0.96;
    m1R(2,1) = 0;
    m1R(2,2) = 0.28;

    CPPUNIT_ASSERT(m1L == m1R);

    Quat4d q1(0.0, -0.6, 0.0, -0.8);
    RotMat3x3d m2L(q1);
    RotMat3x3d m2R(m1R);
    CPPUNIT_ASSERT(m2L == m2R);
    
    //test setupRotationMatrix by euler angles
    Vector3d v1(0.0, M_PI/2.0, 0.0);
    RotMat3x3d m3L(v1);
    RotMat3x3d m3R;
    m3R(0,0) = 1.0;
    m3R(0,1) = 0;
    m3R(0,2) = 0.0;
    m3R(1,0) = 0.0;
    m3R(1,1) = 0.0;
    m3R(1,2) = 1.0;
    m3R(2,0) = 0.0;
    m3R(2,1) = -1.0;
    m3R(2,2) = 0.0;    
    CPPUNIT_ASSERT( m3L == m3R);

    RotMat3x3d m4L(0.0, M_PI/2.0, 0.0);
    RotMat3x3d m4R = m3R;    
    CPPUNIT_ASSERT( m4L == m4R);    

    Vector3d v2(M_PI/4.0, M_PI/4.0, M_PI/4.0);
    RotMat3x3d m5L(v2);
    RotMat3x3d m5R;
    double root2Over4 = sqrt(2)/4.0;
    m5R(0,0) = 0.5 - root2Over4;
    m5R(0,1) = 0.5 + root2Over4;
    m5R(0,2) = 0.5;
    m5R(1,0) = -0.5 -root2Over4;
    m5R(1,1) = -0.5 + root2Over4;
    m5R(1,2) = 0.5;
    m5R(2,0) = 0.5;
    m5R(2,1) = -0.5;
    m5R(2,2) = sqrt(2)/2.0;    
    CPPUNIT_ASSERT( m5L == m5R);

    
}

void SquareMatrix3TestCase::testOtherMemberFunctions() {
    //test inverse
    RotMat3x3d ident = RotMat3x3d::identity();
    CPPUNIT_ASSERT(ident == ident.inverse());

    RotMat3x3d m1;
    m1(0,0) = 1.0;
    m1(0,1) = 5.0;
    m1(0,2) = 3.0;
    m1(1,0) = 3.0;
    m1(1,1) = 1.0;
    m1(1,2) = 2.0;
    m1(2,0) = 0.0;
    m1(2,1) = -21.0;
    m1(2,2) = -81.0; 
    
    CPPUNIT_ASSERT(m1 == (m1.inverse()).inverse());
    
    //test determinant
    RotMat3x3d m2;
    m2(0,0) = 1.0;
    m2(0,1) = 5.0;
    m2(0,2) = 3.0;
    m2(1,0) = 6.0;
    m2(1,1) = 0.0;
    m2(1,2) = 2.0;
    m2(2,0) = 0.0;
    m2(2,1) = -1.0;
    m2(2,2) = 1.0; 
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m2.determinant(), -46.0, OpenMD::NumericConstant::epsilon);
}
void SquareMatrix3TestCase::testTransformation(){

    //test toQuaternion
    RotMat3x3d m1;
    Quat4d q1L;
    Quat4d q1R(0.0, -0.6, 0.0, -0.8);
    m1(0,0) = -0.28;
    m1(0,1) =  0;
    m1(0,2) =  0.96;
    m1(1,0) =  0.0;
    m1(1,1) = -1.0;
    m1(1,2) =  0.0;
    m1(2,0) = 0.96;
    m1(2,1) = 0;
    m1(2,2) = 0.28;
    q1L = m1.toQuaternion();
    //CPPUNIT_ASSERT( q1L == q1R);    

    RotMat3x3d m2;
    Quat4d q2L(0.4, -0.6, 0.3, -0.8);
    Quat4d q2R;
    q2L.normalize();

    m2 = q2L.toRotationMatrix3();    
    q2R = m2.toQuaternion();
    CPPUNIT_ASSERT( q2L == q2R);    
    
    //test toEuler
    Vector3d v1L;
    Vector3d v1R(M_PI/4.0, M_PI/4.0, M_PI/4.0);
    RotMat3x3d m3;
    double root2Over4 = sqrt(2)/4.0;
    m3(0,0) = 0.5 - root2Over4;
    m3(0,1) = 0.5 + root2Over4;
    m3(0,2) = 0.5;
    m3(1,0) = -0.5 -root2Over4;
    m3(1,1) = -0.5 + root2Over4;
    m3(1,2) = 0.5;
    m3(2,0) = 0.5;
    m3(2,1) = -0.5;
    m3(2,2) = sqrt(2)/2.0;    
    v1L = m3.toEulerAngles();
    CPPUNIT_ASSERT( v1L == v1R);

    //test diagonalize

    RotMat3x3d m4;    
    RotMat3x3d a;
    Vector3d w;
    RotMat3x3d m5L;
    RotMat3x3d m5R;
    m4(0, 0) = 3.0;
    m4(0, 1) = 4.0;
    m4(0, 2) = 5.0;
    m4(1, 0) = 4.0;
    m4(1, 1) = 5.0;
    m4(1, 2) = 6.0;    
    m4(2, 0) = 5.0;
    m4(2, 1) = 6.0;
    m4(2, 2) = 7.0; 
    a = m4;
    
    RotMat3x3d::diagonalize(a, w, m5L);

    m5R(0, 0) = 0.789067 ;
    m5R(0, 1) = -0.408248;
    m5R(0, 2) = 0.459028;
    m5R(1, 0) = 0.090750;
    m5R(1, 1) = 0.816497;
    m5R(1, 2) = 0.570173;    
    m5R(2, 0) = -0.607567;
    m5R(2, 1) = -0.408248 ;
    m5R(2, 2) = 0.681319; 

    CPPUNIT_ASSERT(m5L == m5R);
}
