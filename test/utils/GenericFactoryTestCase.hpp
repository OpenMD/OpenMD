#ifndef TEST_PROPERTYMAPTESTCASE_HPP
#define TEST_PROPERTYMAPTESTCASE_HPP

#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include "utils/GenericFactory.hpp"

using namespace oopse;


class Shape {
    public:
        std::string getType() {return type_;}
        void setType(const std::string& type) { type_ = type;}
    protected:
        Shape(){};
    private:
        std::string type_;
    
};

typedef GenericFactory<Shape> ShapeFactory;

class Line : public Shape {
    public:
        Line() {
            setType("Line");
        }
};

Shape* createLine(){
    return new Line();
}

//register createLine
const bool registeredCreateLine = ShapeFactory::getInstance()->registerCreator("MyLine", createLine);


class Circle : public Shape {
    public:
        Circle() {
            setType("Circle");
        }

};

Shape* createCircle() {
    return new Circle();
}

//register createLine
const bool registeredCreateCircle = ShapeFactory::getInstance()->registerCreator("MyCircle", createCircle);


class AnotherCircle : public Shape {
    public:
        AnotherCircle() {
            setType("AnotherCircle");
        }

};

Shape* createAnotherCircle() {
    return new AnotherCircle();
}    

class Cubic : public Shape {
    public:
        Cubic() {
            setType("Cubic");
        }
};

DECLARE_CREATOR(Shape, Cubic)
    
class GenericFactoryTestCase : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE( GenericFactoryTestCase );
    CPPUNIT_TEST(testGenericFactory);
    CPPUNIT_TEST_SUITE_END();

    public:

        void testGenericFactory(); 
};


#endif //TEST_PROPERTYMAPTESTCASE_HPP

