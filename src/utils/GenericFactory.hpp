/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/**
 * @file GenericFactory.hpp
 * @author Teng Lin
 * @date 10/24/2004
 * @version 1.0
 */
#ifndef UTIL_GENERICFACTORY_HPP
#define UTIL_GENERICFACTORY_HPP
#include <cassert>
#include <map>
#include <string>
#include <vector>

namespace oopse {

/**
 * @class GenericFactory GenericFactory.hpp "utils/GenericFactory.hpp"
 * @brief GenericFactory is a template based Object Factory 
 * Factory pattern is used to define an interface for creating an object.
 * 
 * @param Object the base class of the hierarchy for which you provide the object factory.
 * @param IdentType the object that identifies the type of the concrete object. Default type is string
 * @param Creator  the callable entity that creates objects. This type must support operator(),
 * taking no parameters and returning a pointer to Object. Default type is function pointer.
 *
 * Usage:
 * @code
 * //Shape class
 * class Shape {
 * ...
 * };
 * 
 * //instantiating a new object factory
 * typedef GenericFactory<Shape> ShapeFactory;
 *
 * //Line class
 * class Line : public Shape{
 * ...
 * };
 *
 * //declare function to create Line
 * Shape* createLine() {
 *   return new Line;
 * }
 *
 * //register createLine
 * //note: must put ShapeFactory::getInstance()->registerCreator("Line", createLine) on the right
 * //hand side, otherwise the compiler will consider it as a function declaration
 * const bool registeredLine = ShapeFactory::getInstance()->registerCreator("Line", createLine);
 *
 * //Circle class
 * class Circle : public Shape{
 * ...
 * };
 *
 * //declare function to create Circle
 * Shape* createCircle() { 
 *   return new Circle;
 * }
 *
 * //register createCircle
 * const bool registeredCircle = ShapeFactory::getInstance()->registerCreator("Circle", createCircle); 
 *
 * //create object by ident
 * Line* line = ShapeFactory::getInstance()->createObject("Line");
 * Circle* circle = ShapeFactory::getInstance()->createObject("Circle"); 
 * @endcode
 *
 * Or the user can use predefined macro DECLARE_CREATOR and REGISTER_CREATOR
 * @code
 * //Shape class
 * class Shape {
 * ...
 * };
 * 
 * //instantiating a new object factory
 * typedef GenericFactory<Shape> ShapeFactory;
 *
 * //Line class
 * class Line : public Shape{
 * ...
 * };
 *
 * //declare function using macro 
 * DECLARE_CREATOR(Shape, Line)
 * 
 * //register using macro
 * REGISTER_CREATOR(ShapeFactory, "Line", Line);
 
 * //Circle class
 * class Circle : public Shape{
 * ...
 * };
 *
 * //declare function using macro 
 * DECLARE_CREATOR(Shape, Circle)
 * 
 * //register using macro
 * REGISTER_CREATOR(ShapeFactory, "Circle", Circle);
 * @endcode
 */
template<class Object, typename IdentType = std::string, typename Creator = Object* (*)()>
class GenericFactory {
    public:
        typedef GenericFactory<Object, IdentType, Creator> FactoryType;
        typedef std::map<IdentType, Creator> CreatorMapType;
        
        /**
         * Returns an instance of object factory
         * @return an instance of object factory
         */        
        static FactoryType* getInstance(){
            if (instance_ == NULL)
                instance_ = new FactoryType;
            return instance_;
        }

        /**
         * Registers a creator with a type identifier
         * @return true if registration is succeed, otherwise return false
         * @id the identification of the concrete object
         * @creator the object responsible to create the concrete object 
         */
        bool registerCreator(const IdentType& id, Creator creator) {
            return creatorMap_.insert(
                CreatorMapType::value_type(id, creator)).second;
        }

        /**
         * Unregisters the creator for the given type identifier. If the type identifier 
         * was previously registered, the function returns true.
         * @return truethe type identifier was previously registered and the creator is removed,
         * otherwise return false
         * @id the identification of the concrete object
         */
        bool unregisterCreator(const IdentType& id) {
            return creatorMap_.erase(id) == 1;
        }

        /**
         * Looks up the type identifier in the internal map. If it is found, it invokes the
         * corresponding creator for the type identifier and returns its result. 
         * @return a pointer of the concrete object, return NULL if no creator is registed for 
         * creating this concrete object
         * @param id the identification of the concrete object
         */
        Object* createObject(const IdentType& id) {
            typename CreatorMapType::iterator i = creatorMap_.find(id);
            if (i != creatorMap_.end()) {
                //invoke functor to create object
                return (i->second)();
            } else {
                return NULL;
            }
        }

        /** 
         *  Returns all of the registed  type identifiers
         * @return all of the registed  type identifiers
         */
        std::vector<IdentType> getIdents() {
            std::vector<IdentType> idents;
            typename CreatorMapType::iterator i;

            for (i = creatorMap_.begin(); i != creatorMap_.end(); ++i) {
                idents.push_back(i->first);
            }
            
            return idents;
        }

    public:
        static FactoryType* instance_;
        CreatorMapType creatorMap_;
};

/** write out all of the type identifier to a output stream */
template<typename O, typename I, typename C>
std::ostream& operator <<(std::ostream& o, GenericFactory<O, I, C>& factory) {
    std::vector<I> idents;
    std::vector<I>::iterator i;

    idents = factory.getIdents();

    o << "Avaliable type identifiers in this factory: " << std::endl;
    for (i = idents.begin(); i != idents.end(); ++i) {
        o << *i << std::endl;
    }

    return o;
}

//static template class member
template<class Object, typename IdentType,typename Creator>
GenericFactory<Object,IdentType,Creator>* GenericFactory<Object,IdentType,Creator>::instance_ ; 


#define DECLARE_CREATOR(abstractObject, concreteObject) \
    inline abstractObject* create##concreteObject(){\
        return new concreteObject;\
    }

#define REGISTER_CREATOR(factory, ident, concreteObject) \
        const bool registered##concreteObject = factory::getInstance()->registerCreator(ident, create##concreteObject); 


}//namespace oopse
#endif //UTIL_GENERICFACTORY_HPP

