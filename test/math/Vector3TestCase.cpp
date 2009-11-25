#include "math/Vector3TestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( Vector3TestCase );


void Vector3TestCase::setUp(){
    zero[0] = 0.0;
    zero[1] = 0.0;
    zero[2] = 0.0;

    one[0] = 1.0;
    one[1] = 1.0;
    one[2] = 1.0;

    two[0] = 2.0;
    two[1] = 2.0;
    two[2] = 2.0;

    v1[0] = 1.0;
    v1[1] = 2.0;
    v1[2] = 3.0;

    v2[0] = 4.0;
    v2[1] = 1.0;
    v2[2] = 2.0;

    v3[0] = 1.0;
    v3[1] = 10.0;
    v3[2] = -7.0;

    
}

void Vector3TestCase::tearDown(){
}

void Vector3TestCase::testConstructors(){
    double b[] = {2.9, 3.2, 1.2};
    Vector3d v(b);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.x(), 2.9, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.y(), 3.2, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.z(), 1.2, OpenMD::NumericConstant::epsilon);    
}

void Vector3TestCase::testArithmetic(){

    
}

void Vector3TestCase::testOperators(){
    Vector3d tmp;
    //test +=

    //test +
    CPPUNIT_ASSERT(one + one == two);
    
    //test -=
    tmp = two;
    tmp -= one;
    CPPUNIT_ASSERT(tmp == one);
    
    //test -
    CPPUNIT_ASSERT(two -one == one);
    
    //test *=
    tmp = two;
    tmp *= 0.5;
    CPPUNIT_ASSERT(tmp == one);
    
    //test *
    CPPUNIT_ASSERT(two * 0.5 == one);
    CPPUNIT_ASSERT(0.5 * two == one);

    //test /=
    tmp = two;
    tmp *= 2.0;
    CPPUNIT_ASSERT(tmp == one * 4.0);    
    
    //test /
    CPPUNIT_ASSERT( two /2.0 == one);
    CPPUNIT_ASSERT( two /4.0 == one * 0.5);

}

void Vector3TestCase::testAccessEntries(){

    CPPUNIT_ASSERT_DOUBLES_EQUAL(v1.z(), 3.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v2.x(), 4.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3.y(), 10.0, OpenMD::NumericConstant::epsilon);

    Vector3d tmp;
    tmp.x() = 78.01;
    tmp.y() = 21.0;
    tmp.z() =133.12;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp.x(), 78.01, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp.y(), 21.0, OpenMD::NumericConstant::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp.z(), 133.12, OpenMD::NumericConstant::epsilon);

}

void Vector3TestCase::testOtherTemplateFunctions(){      
    //test cross
    CPPUNIT_ASSERT(cross(v1, v2) == v3);
    CPPUNIT_ASSERT(cross(one, two) == zero);
    
}
