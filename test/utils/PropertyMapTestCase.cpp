#include "utils/PropertyMapTestCase.hpp"
#include <iostream>
#include <algorithm>
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PropertyMapTestCase );


void PropertyMapTestCase::testPropertyMap(){
    PropertyMap props;

    //test addProperty
    BoolGenericData* b0 = new BoolGenericData("BoolData");
    b0->setData(false);
    props.addProperty(b0);
    CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b0);
    
    BoolGenericData* b1 = new BoolGenericData("BoolData");
    b1->setData(true);
    props.addProperty(b1);
    CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b1);


    IntGenericData* i1 = new IntGenericData("IntData");
    i1->setData(89);
    props.addProperty(i1);
    
    FloatGenericData* f1 = new FloatGenericData("FloatData");
    f1->setData(49.328);
    props.addProperty(f1);

    DoubleGenericData* d1 = new DoubleGenericData("DoubleData");
    d1->setData(95.1933432);
    props.addProperty(d1);

    StringGenericData* s1 = new StringGenericData("StringData");
    s1->setData("Hello");
    props.addProperty(s1);


    IntVectorGenericData* iv1 = new IntVectorGenericData("IntVector");
    iv1->push_back(2);
    iv1->push_back(1);
    iv1->push_back(324);
    props.addProperty(iv1);

    //test getPropertyNames
    std::vector<std::string> propNames = props.getPropertyNames();


    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") != propNames.end());

    //test getProperties    
    std::vector<GenericData*> propPointers = props.getProperties();    
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), b1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), i1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), f1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), d1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), s1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), iv1) != propPointers.end());

    //test getPropertyByName
    CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b1);
    CPPUNIT_ASSERT(props.getPropertyByName("IntData") == i1);
    CPPUNIT_ASSERT(props.getPropertyByName("FloatData") == f1);
    CPPUNIT_ASSERT(props.getPropertyByName("DoubleData") == d1);
    CPPUNIT_ASSERT(props.getPropertyByName("StringData") == s1);
    CPPUNIT_ASSERT(props.getPropertyByName("IntVector") == iv1);

    CPPUNIT_ASSERT(b1->getData() == true);
    CPPUNIT_ASSERT(i1->getData() == 89);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getData(), 49.328, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(d1->getData(), 95.1933432, 0.000001);   
    CPPUNIT_ASSERT(s1->getData() == "Hello");
    CPPUNIT_ASSERT_EQUAL((*iv1)[0], 2);
    CPPUNIT_ASSERT_EQUAL((*iv1)[1], 1);
    CPPUNIT_ASSERT_EQUAL((*iv1)[2], 324);        

    //test removeProperty
    props.removeProperty("DoubleData");
    props.removeProperty("FloatData");
    props.removeProperty("IntVector");
    propNames = props.getPropertyNames();
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") == propNames.end());

    //test clearProperties
    props.clearProperties();
    propNames = props.getPropertyNames();
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") == propNames.end());    

}
