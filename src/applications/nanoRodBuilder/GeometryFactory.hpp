/*
 *  GeometryFactory.hpp
 *  OOPSE-2.0
 *
 *  Created by Charles F. Vardeman II on 5/3/05.
 *  Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
#ifndef NANORODBUILDER_GEOMETRYFACTORY_HPP
#define NANORODBUILDER_GEOMETRYFACTORY_HPP
#include <cassert>
#include <map>
#include <string>
#include <vector>
#include <iostream>
namespace oopse {
   
   //forward declaration
   class Geometry;
   class GeometryCreator;
   /**
      * @class GeometryFactory GeometryFactory.hpp 
    * Factory pattern and Singleton Pattern are used to define an interface for creating an Lattice.
    */
   class GeometryFactory {
public:
      
      typedef std::map<std::string, LatticeCreator*> CreatorMapType;
      typedef std::vector<std::string> IdentVectorType;
      typedef std::vector<std::string>::iterator IdentVectorIterator;
      
      ~GeometryFactory();
      
      /**
         * Returns an instance of Geometry factory
       * @return an instance of Geometry factory
       */        
      static GeometryFactory* getInstance() {
         
         if (instance_ == NULL) {
            instance_ = new GeometryFactory();
         }
         return instance_;
         
      }
      
      /**
         * Registers a creator with a type identifier
       * @return true if registration is succeed, otherwise return false
       * @id the identification of the concrete object
       * @creator the object responsible to create the concrete object 
       */
      bool registerGeometry(GeometryCreator* creator);
      
      /**
         * Unregisters the creator for the given type identifier. If the type identifier 
       * was previously registered, the function returns true.
       * @return truethe type identifier was previously registered and the creator is removed,
       * otherwise return false
       * @id the identification of the concrete object
       */
      bool unregisterGeometry(const std::string& id);
      /**
         * Looks up the type identifier in the internal map. If it is found, it invokes the
       * corresponding creator for the type identifier and returns its result. 
       * @return a pointer of the concrete object, return NULL if no creator is registed for 
       * creating this concrete object
       * @param id the identification of the concrete object
       */
      Lattice* createGeometry(const std::string& id);
      
      /** 
         *  Returns all of the registed  type identifiers
         * @return all of the registed  type identifiers
         */
      IdentVectorType getIdents();
      
private:
         GeometryFactory() {}
      
      static GeometryFactory* instance_;
      CreatorMapType creatorMap_;
   };
   
   /** write out all of the type identifiers to an output stream */
   std::ostream& operator <<(std::ostream& o, GeometryFactory& factory);
   
}//namespace oopse
#endif //NANORODBUILDER_GEOMETRYFACTORY_HPP
