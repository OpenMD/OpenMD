/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
#include "applications/hydrodynamics/ShapeBuilder.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"
namespace oopse {
  
  Shape* ShapeBuilder::createShape(StuntDouble* sd) {
    Shape* currShape = NULL;
    if (sd->isDirectionalAtom()) {
      currShape = internalCreateShape(static_cast<DirectionalAtom*>(sd));
    } else if (sd->isAtom()) {
      currShape = internalCreateShape(static_cast<Atom*>(sd));
    } else if (sd->isRigidBody()) {
      currShape = internalCreateShape(static_cast<RigidBody*>(sd));
    }       
    return currShape;   
  }
  
  Shape* ShapeBuilder::internalCreateShape(Atom* atom) {
    AtomType* atomType = atom->getAtomType();
    Shape* currShape = NULL;
    if (atomType->isLennardJones()){
      GenericData* data = atomType->getPropertyByName("LennardJones");
      if (data != NULL) {
        LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
        
        if (ljData != NULL) {
          LJParam ljParam = ljData->getData();
          currShape = new Sphere(atom->getPos(), ljParam.sigma/2.0);
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to LJParam\n");
          painCave.severity = OOPSE_ERROR;
          painCave.isFatal = 1;
          simError();          
        }       
      }
    } else {
      int obanum = etab.GetAtomicNum((atom->getType()).c_str());
      if (obanum != 0) {
        currShape = new Sphere(atom->getPos(), etab.GetVdwRad(obanum));
      } else {
        sprintf( painCave.errMsg,
                 "Could not find atom type in default element.txt\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
    return currShape;
  }
  
  Shape* ShapeBuilder::internalCreateShape(DirectionalAtom* datom) {
    AtomType* atomType = datom->getAtomType();
    Shape* currShape = NULL;
    if (atomType->isGayBerne()) {
      DirectionalAtomType* dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);
      
      GenericData* data = dAtomType->getPropertyByName("GayBerne");
      if (data != NULL) {
        GayBerneParamGenericData* gayBerneData = dynamic_cast<GayBerneParamGenericData*>(data);
        
        if (gayBerneData != NULL) {  
          GayBerneParam gayBerneParam = gayBerneData->getData();
          currShape = new Ellipsoid(datom->getPos(), gayBerneParam.GB_d/2.0, gayBerneParam.GB_l/2.0, datom->getA());
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to GayBerneParam\n");
          painCave.severity = OOPSE_ERROR;
          painCave.isFatal = 1;
          simError();   
        }
      } else {
        sprintf( painCave.errMsg, "Can not find Parameters for GayBerne\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();    
      }            
    } else if (atomType->isLennardJones()) {
      GenericData* data = atomType->getPropertyByName("LennardJones");
      if (data != NULL) {
        LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
        
        if (ljData != NULL) {
          LJParam ljParam = ljData->getData();
          currShape = new Sphere(datom->getPos(), ljParam.sigma/2.0);
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to LJParam\n");
          painCave.severity = OOPSE_ERROR;
          painCave.isFatal = 1;
          simError();          
        }       
      } else {
        int obanum = etab.GetAtomicNum((datom->getType()).c_str());
        if (obanum != 0) {
          currShape = new Sphere(datom->getPos(), etab.GetVdwRad(obanum));
        } else {
          sprintf( painCave.errMsg,
                   "Could not find atom type in default element.txt\n");
          painCave.severity = OOPSE_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
      }      
    }
    return currShape;
  }

  Shape* ShapeBuilder::internalCreateShape(RigidBody* rb) {
    
    std::vector<Atom*>::iterator ai; 
    CompositeShape* compositeShape = new CompositeShape;
    Atom* atom;
    for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
      Shape* currShape = NULL;
      if (atom->isDirectionalAtom()){
        currShape = internalCreateShape(static_cast<DirectionalAtom*>(atom));
      }else if (atom->isAtom()){
        currShape =  internalCreateShape(static_cast<Atom*>(atom));
      }
      if (currShape != NULL)
        compositeShape->addShape(currShape);
    }
    
    return compositeShape;
  }  
}
