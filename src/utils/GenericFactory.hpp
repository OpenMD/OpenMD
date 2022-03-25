/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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

namespace OpenMD {

  /**
   * @class GenericFactory GenericFactory.hpp "utils/GenericFactory.hpp"
   * @brief GenericFactory is a template based Object Factory
   * Factory pattern is used to define an interface for creating an object.
   *
   * @param Object the base class of the hierarchy for which you provide the
   object factory.
   * @param IdentType the object that identifies the type of the concrete
   object. Default type is  std::string * @param Creator  the callable entity
   that creates objects. This type must support operator(),
   * taking no parameters and returning a pointer to Object. Default type is
   function pointer.
   *
   * Usage:
   * @code
   * //Shape class
   * class Shape {
   * ...
   * };
   *
   * //instantiating a new object factory
   * using ShapeFactory = GenericFactory<Shape>;
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
   * //note: must put ShapeFactory::getInstance()->registerCreator("Line",
   createLine) on the right
   * //hand side, otherwise the compiler will consider it as a function
   declaration
   * const bool registeredLine =
   ShapeFactory::getInstance()->registerCreator("Line", createLine);
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
   * const bool registeredCircle =
   ShapeFactory::getInstance()->registerCreator("Circle", createCircle);
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
   * using ShapeFactory = GenericFactory<Shape>;
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
  template<class Object, typename IdentType = std::string,
           typename Creator = Object* (*)()>
  class GenericFactory {
  public:
    using FactoryType    = GenericFactory<Object, IdentType, Creator>;
    using CreatorMapType = std::map<IdentType, Creator>;

    /**
     * Returns an instance of object factory
     * @return an instance of object factory
     */
    static FactoryType* getInstance() {
      if (instance_ == NULL) instance_ = new FactoryType;
      return instance_;
    }

    /**
     * Registers a creator with a type identifier
     * @return true if registration is succeed, otherwise return false
     * @param id the identification of the concrete object
     * @param creator the object responsible to create the concrete object
     */
    bool registerCreator(const IdentType& id, Creator creator) {
      return creatorMap_.insert(CreatorMapType::value_type(id, creator)).second;
    }

    /**
     * Unregisters the creator for the given type identifier. If the type
     * identifier was previously registered, the function returns true.
     * @return truethe type identifier was previously registered and the creator
     * is removed, otherwise return false
     * @param id the identification of the concrete object
     */
    bool unregisterCreator(const IdentType& id) {
      return creatorMap_.erase(id) == 1;
    }

    /**
     * Looks up the type identifier in the internal map. If it is found, it
     * invokes the corresponding creator for the type identifier and returns its
     * result.
     * @return a pointer of the concrete object, return NULL if no creator is
     * registed for creating this concrete object
     * @param id the identification of the concrete object
     */
    Object* createObject(const IdentType& id) {
      typename CreatorMapType::iterator i = creatorMap_.find(id);
      if (i != creatorMap_.end()) {
        // invoke functor to create object
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

  /** write out all of the type identifiers to an output stream */
  template<typename O, typename I, typename C>
  std::ostream& operator<<(std::ostream& o, GenericFactory<O, I, C>& factory) {
    std::vector<I> idents;
    std::vector<I>::iterator i;

    idents = factory.getIdents();

    o << "Avaliable type identifiers in this factory: " << std::endl;
    for (i = idents.begin(); i != idents.end(); ++i) {
      o << *i << std::endl;
    }

    return o;
  }

  // static template class member
  template<class Object, typename IdentType, typename Creator>
  GenericFactory<Object, IdentType, Creator>*
      GenericFactory<Object, IdentType, Creator>::instance_;

#define DECLARE_CREATOR(abstractObject, concreteObject) \
  inline abstractObject* create##concreteObject() { return new concreteObject; }

#define REGISTER_CREATOR(factory, ident, concreteObject) \
  const bool registered##concreteObject =                \
      factory::getInstance()->registerCreator(ident, create##concreteObject);
}  // namespace OpenMD

#endif  // UTIL_GENERICFACTORY_HPP
