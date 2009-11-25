#include "utils/GenericDataTestCase.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataTestCase );

void GenericDataTestCase::testGenericData(){

    //test constructor
    GenericData gd1;
    
    GenericData gd2("OpenMD_Generic");

    CPPUNIT_ASSERT( gd1.getID() == "UndefinedGenericData");
    CPPUNIT_ASSERT( gd2.getID() == "OpenMD_Generic");

    gd1.setID("Dummy");
    CPPUNIT_ASSERT( gd1.getID() == "Dummy");
}

void GenericDataTestCase::testSimpleTypeData(){

    //test IntGenericData
    BoolGenericData b1("Dummy_Bool");
    b1.setData(true);
    CPPUNIT_ASSERT(b1.getData());
    CPPUNIT_ASSERT(b1.getID() == "Dummy_Bool");

    b1.setData(10.0);
    CPPUNIT_ASSERT(b1.getData());

    b1.setData(0.0);
    CPPUNIT_ASSERT(!b1.getData());
    
    //test IntGenericData
    IntGenericData i1("Dummy_Int");
    i1.setData(10);
    CPPUNIT_ASSERT_EQUAL(i1.getData(), 10);
    CPPUNIT_ASSERT(i1.getID() == "Dummy_Int");

    IntGenericData i2("Dummy_Int");
    i2.setData(10.000);
    CPPUNIT_ASSERT_EQUAL(i1.getData(), 10);
    
    //test FloatGenericData
    FloatGenericData f1("Dummy_Float");
    f1.setData(10);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1.getData(), 10, 0.000001);
    CPPUNIT_ASSERT(f1.getID() == "Dummy_Float");

    
    //test DoubleGenericData
    DoubleGenericData d1("Dummy_Double");
    d1.setData(232);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(d1.getData(), 232, 0.000001);
    CPPUNIT_ASSERT(d1.getID() == "Dummy_Double");

    //test StringGenericData
    StringGenericData s1("Dummy_String");
    s1.setData("Hello World");
    CPPUNIT_ASSERT(s1.getData() == "Hello World");
    CPPUNIT_ASSERT(s1.getID() == "Dummy_String");

    DoubleGenericData d2("Dummy_Double");
    d2 = f1;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(d2.getData(), 10, 0.000001);
    
     FloatGenericData f2("Dummy_Float");
     f2 = d1;
     CPPUNIT_ASSERT_DOUBLES_EQUAL(d1.getData(), 232, 0.000001);

     //test getData (return reference)
     f2.getData() = 0.004;
     CPPUNIT_ASSERT_DOUBLES_EQUAL(f2.getData(), 0.004, 0.000001);
     
    
}

void GenericDataTestCase::testSTLContainerTypeData(){

    //test IntVectorGenericData
    IntVectorGenericData iv1("IntVector1");
    CPPUNIT_ASSERT(iv1.getID() == "IntVector1");
    iv1.setID("Dummy_Test");
    CPPUNIT_ASSERT(iv1.getID() == "Dummy_Test");

    iv1.push_back(2.30);
    iv1.push_back(1);
    iv1.push_back(324);

    CPPUNIT_ASSERT_EQUAL(iv1[0], 2);
    CPPUNIT_ASSERT_EQUAL(iv1[1], 1);
    CPPUNIT_ASSERT_EQUAL(iv1[2], 324);

    IntVectorGenericData iv2("IntVector2");

    iv2.push_back(1);
    iv2.push_back(3);
    iv2.push_back(5);
    iv2.push_back(7);

    iv1 = iv2;
    CPPUNIT_ASSERT(iv1.getID() == "Dummy_Test");
    CPPUNIT_ASSERT(iv1.size() == 4);
    CPPUNIT_ASSERT_EQUAL(iv1[0], 1);
    CPPUNIT_ASSERT_EQUAL(iv1[1], 3);
    CPPUNIT_ASSERT_EQUAL(iv1[2], 5);
    CPPUNIT_ASSERT_EQUAL(iv1[3], 7);

    //test FloatVectorGenericData
    FloatVectorGenericData fv2("FloatVector2");

    fv2.push_back(251.21);
    fv2.push_back(42.90);

}
