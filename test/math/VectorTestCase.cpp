#include "math/VectorTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( VectorTestCase );


void VectorTestCase::setUp(){
    zero[0] = 0.0;
    zero[1] = 0.0;
    zero[2] = 0.0;
    zero[3] = 0.0;

    one[0] = 1.0;
    one[1] = 1.0;
    one[2] = 1.0;
    one[3] = 1.0;

    two[0] = 2.0;
    two[1] = 2.0;
    two[2] = 2.0;
    two[3] = 2.0;


    v1[0] = 3.0;
    v1[1] = 0.0;
    v1[2] = 2.0;
    v1[3] = 1.0;        

    v2[0] = -3.0;
    v2[1] = 0.0;
    v2[2] = -2.0;
    v2[3] = -1.0;        

    v3[0] = 4.0;
    v3[1] = 1.0;
    v3[2] = 3.0;
    v3[3] = 2.0;        

    s1 = 1.0;
    s2 = 2.0;
    
}

void VectorTestCase::testConstructors(){
    Vec4 a0;
    Vec4 a1(1);
    Vec4 a2(2);
    
    CPPUNIT_ASSERT( a0 == zero);
    CPPUNIT_ASSERT( a1 == one);
    CPPUNIT_ASSERT( a2 == two);

    CPPUNIT_ASSERT( a1 != two);

    //test copy constructor
    Vec4 b1(v1);
    CPPUNIT_ASSERT( b1 == v1);

    //test operator =
    b1 = v2;
    CPPUNIT_ASSERT( b1 == v2);

    //test constructor from an array    
    double tempArray[] = {1.0, 2.0, 5.0, 8.0};
    Vec4 tempV;
    tempV[0] = 1.0;
    tempV[1] = 2.0;
    tempV[2] = 5.0;
    tempV[3] = 8.0;
    
    Vec4 b2(tempV);
    CPPUNIT_ASSERT( b2 == tempArray);
    
}

void VectorTestCase::testArithmetic(){
    //test negate
    Vec4 a0 = v2;
    a0.negate();
    CPPUNIT_ASSERT (a0 == v1);
    
    Vec4 a1;
    a1.negate(v2);
    CPPUNIT_ASSERT(a1 == v1);

    //test add
    Vec4 a2;
    a2 = v1;
    a2.add(v2);    
    CPPUNIT_ASSERT( a2 == zero);

    Vec4 a3;
    a3.add(v2, v3);       
    CPPUNIT_ASSERT( a3 == one);
    

    //test sub
    Vec4 a4;
    a4 = two;
    a4.sub(one);
    CPPUNIT_ASSERT( a4 == one);    

    Vec4 a5;
    a5.sub(two, one);
    CPPUNIT_ASSERT( a5 == one);      
    
    //test mul
    Vec4 a6;
    a6 = one;
    a6.mul(2.0);
    CPPUNIT_ASSERT( a6 == two);      

    Vec4 a7;
    a7.mul(one, 2.0);
    CPPUNIT_ASSERT( a7 == two);      
    Vec4 a8;
    a7.mul(zero, 2.0);
    CPPUNIT_ASSERT( a7 == zero);      

    //test div
    Vec4 a9;
    a9 = two;
    a9.div(2.0);
    CPPUNIT_ASSERT( a9 == one);      

    Vec4 a10;
    a10.div(two, 2.0);
    CPPUNIT_ASSERT( a10 == one);      
    Vec4 a11;
    a11.mul(zero, 2.0);
    CPPUNIT_ASSERT( a11 == zero);      
    
     
}

void VectorTestCase::testOperators(){
    //test unary minus
    Vec4 a0 = v2;
    a0 = - a0;
    CPPUNIT_ASSERT (a0 == v1);
   

    //test add
    Vec4 a1;
    a1 = v1;
    a1 += v2;
    CPPUNIT_ASSERT( a1 == zero);

    Vec4 a2;
    a2 = v2 + v3;
    CPPUNIT_ASSERT( a2 == one);

    //test sub
    Vec4 a3;
    a3 = two;
    a3 -= one;
    CPPUNIT_ASSERT( a3 == one);    

    Vec4 a4;
    a4 = two - one;
    CPPUNIT_ASSERT( a4 == one);

    Vec4 a5;
    a5.sub(two, one);
    CPPUNIT_ASSERT( a5 == one);      
    
    //test mul
    Vec4 a6;
    a6 = one;
    a6 *= 2.0;
    CPPUNIT_ASSERT( a6 == two);      

    Vec4 a7;
    a7 =  one * 2.0;
    CPPUNIT_ASSERT( a7 == two);      
    a7 =  2.0 * one;
    CPPUNIT_ASSERT( a7 == two);      

    //test div
    Vec4 a8;
    a8 = two;
    a8 /= 2.0;
    CPPUNIT_ASSERT( a8 == one);      
    a8 = two /2.0;
    CPPUNIT_ASSERT( a8 == one);      
}

void VectorTestCase::testAccessEntries(){
    //test [] operator

    CPPUNIT_ASSERT_DOUBLES_EQUAL(zero[0], 0.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(one[0] , 1.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3[0] , 4.0, oopse::epsilon);
    //test () operator
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3(0) , 4.0, oopse::epsilon);

    Vec4 a1;
    double *pa1 = a1.getArrayPointer();
    
    pa1[0] = 4.0;
    pa1[1] = 1.0;
    pa1[2] = 3.0;
    pa1[3] = 2.0;        

    CPPUNIT_ASSERT(a1 == v3);    
}

void VectorTestCase::testOtherMemberFunctions(){
    //test length()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(zero.length(), 0.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(one.length(), 2.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v2.length(), sqrt(14.0), oopse::epsilon);
    
    //test lengthSquare()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(zero.lengthSquare(), 0.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(one.lengthSquare(), 4.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v2.lengthSquare(), 14.0, oopse::epsilon);

    //test normalize()
    Vec4 a1 = one;
    Vec4 a2 = two;

    a1.normalize();
    a2.normalize();
    CPPUNIT_ASSERT(a1 == a2);

    //test isNormalized();
    CPPUNIT_ASSERT(a1.isNormalized());
    CPPUNIT_ASSERT(!one.isNormalized());
    

}
void VectorTestCase::testOtherTemplateFunctions(){        
    //test dot
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dot(one, two), 8.0, oopse::epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dot(v1, v3), 20.0, oopse::epsilon);

    //test distance
    CPPUNIT_ASSERT_DOUBLES_EQUAL(distance(one, two), 2.0, oopse::epsilon);    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(distance(v1, v2), sqrt(56.0), oopse::epsilon);
    
    //test distanceSquare
    CPPUNIT_ASSERT_DOUBLES_EQUAL(distanceSquare(one, two), 4.0, oopse::epsilon);    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(distanceSquare(v1, v2), 56, oopse::epsilon);

}
