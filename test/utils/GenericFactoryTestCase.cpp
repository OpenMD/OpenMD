#include "utils/GenericFactoryTestCase.hpp"
#include <iostream>
#include <algorithm>
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( GenericFactoryTestCase );
void GenericFactoryTestCase::testGenericFactory() {
    //test getIdents
    std::vector<std::string> idents;
    idents = ShapeFactory::getInstance()->getIdents();
    CPPUNIT_ASSERT(std::find(idents.begin(), idents.end(), "MyLine") != idents.end());
    CPPUNIT_ASSERT(std::find(idents.begin(), idents.end(), "MyCircle") != idents.end());

    //test createObject
    Shape* line = ShapeFactory::getInstance()->createObject("MyLine");
    CPPUNIT_ASSERT(line != NULL && line->getType() == "Line");
    delete line;

    Shape* circle = ShapeFactory::getInstance()->createObject("MyCircle");
    CPPUNIT_ASSERT(circle != NULL && circle->getType() == "Circle");
    delete circle;

    //test registerCreator
    bool registeredCreateAnotherCircle = 
      ShapeFactory::getInstance()->registerCreator("MyCircle", createAnotherCircle);
    CPPUNIT_ASSERT(!registeredCreateAnotherCircle);

    //test unregisterCreator
    ShapeFactory::getInstance()->unregisterCreator("MyCircle");
    idents = ShapeFactory::getInstance()->getIdents();
    CPPUNIT_ASSERT(idents.size() == 1 && std::find(idents.begin(), idents.end(), "MyLine") != idents.end());


    //test registerCreator
    registeredCreateAnotherCircle = 
      ShapeFactory::getInstance()->registerCreator("MyCircle", createAnotherCircle);
    CPPUNIT_ASSERT(registeredCreateAnotherCircle);
    idents = ShapeFactory::getInstance()->getIdents();
    CPPUNIT_ASSERT(std::find(idents.begin(), idents.end(), "MyLine") != idents.end());
    CPPUNIT_ASSERT(std::find(idents.begin(), idents.end(), "MyCircle") != idents.end());
    
    //expect createAnotherCircle will replace createCircle
    Shape* anotherCircle = ShapeFactory::getInstance()->createObject("MyCircle");
    CPPUNIT_ASSERT(anotherCircle != NULL && anotherCircle->getType() == "AnotherCircle");
    delete anotherCircle;

    //test macro REGISTER_CREATOR
    REGISTER_CREATOR(ShapeFactory, "MyCubic", Cubic);
    Shape* cubic = ShapeFactory::getInstance()->createObject("MyCubic");
    CPPUNIT_ASSERT(circle != NULL && circle->getType() == "Cubic");
}
