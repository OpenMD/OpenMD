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

#include "applications/hydrodynamics/CompositeShape.hpp" 
#include "utils/MemoryUtils.hpp"
#include "applications/hydrodynamics/HydrodynamicsModel.hpp"
namespace oopse {
  
  CompositeShape::~CompositeShape() {
    MemoryUtils::deletePointers(shapes_);
  }
  bool CompositeShape::isInterior(Vector3d pos) {
    bool result = false;
    std::vector<Shape*>::iterator iter;
    for (iter = shapes_.begin(); iter != shapes_.end(); ++ iter) {
      if ((*iter)->isInterior(pos)) {
        result = true;
        break;
      }
    }
    
    return result;
  }
  
  template<class Cont, class Predict>
  void swap_if(Cont& b1, Cont& b2, Predict predict) {
    unsigned int size = b1.size();
    assert(size == b2.size());
    for (unsigned int i = 0 ; i < size; ++i) {
      if (predict(b1[i], b2[i]))
        std::swap(b1[i], b2[i]);
    }
    
  }
  
  std::pair<Vector3d, Vector3d> CompositeShape::getBoundingBox() {
    std::vector<Shape*>::iterator iter = shapes_.begin();
    std::pair<Vector3d, Vector3d>  boundary = (*iter)->getBoundingBox();
    for (++iter; iter != shapes_.end(); ++iter) {
      std::pair<Vector3d, Vector3d> currBoundary = (*iter)->getBoundingBox();
      swap_if(boundary.first, currBoundary.first, std::greater<RealType>());
      swap_if(boundary.second, currBoundary.second, std::less<RealType>());        
    }
    
    return boundary;
  }
  
  HydroProp* CompositeShape::getHydroProp(RealType viscosity, RealType temperature) {
    HydroProp* props;
    props->setCOR(V3Zero);
    sprintf( painCave.errMsg,
             "CompositeShape was asked to return an analytic HydroProps.\n");
    painCave.severity = OOPSE_ERROR;
    painCave.isFatal = 1;
    simError();          
    return props;
  }   
}
