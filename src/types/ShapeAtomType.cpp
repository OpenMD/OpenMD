/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#include "types/ShapeAtomType.hpp"
#include "UseTheForce/DarkSide/shapes_interface.h"
#include <cstdio>

namespace OpenMD {
  
  ShapeAtomType::~ShapeAtomType() {
    std::vector<RealSphericalHarmonic*>::iterator iter;
    for (iter = contactFuncs.begin(); iter != contactFuncs.end(); ++iter) 
      delete (*iter);
    for (iter = rangeFuncs.begin(); iter != rangeFuncs.end(); ++iter) 
      delete (*iter);
    for (iter = strengthFuncs.begin(); iter != strengthFuncs.end(); ++iter) 
      delete (*iter);
    contactFuncs.clear();
    rangeFuncs.clear();
    strengthFuncs.clear();
  }
  
  RealType ShapeAtomType::getContactValueAt(RealType costheta, RealType phi) {
    
    std::vector<RealSphericalHarmonic*>::iterator contactIter;
    RealType contactVal;
    
    contactVal = 0.0;
    
    for(contactIter = contactFuncs.begin();  contactIter != contactFuncs.end();
        ++contactIter) 
      contactVal += (*contactIter)->getValueAt(costheta, phi);
    
    return contactVal;
  }
  
  RealType ShapeAtomType::getRangeValueAt(RealType costheta, RealType phi) {
    
    std::vector<RealSphericalHarmonic*>::iterator rangeIter;
    RealType rangeVal;
    
    rangeVal = 0.0;
    
    for(rangeIter = rangeFuncs.begin();  rangeIter != rangeFuncs.end(); 
        ++rangeIter)     
      rangeVal += (*rangeIter)->getValueAt(costheta, phi);
    
    return rangeVal;
  }
  
  RealType ShapeAtomType::getStrengthValueAt(RealType costheta, RealType phi) {
    
    std::vector<RealSphericalHarmonic*>::iterator strengthIter;
    RealType strengthVal;
    
    strengthVal = 0.0;
    
    for(strengthIter = strengthFuncs.begin();  
        strengthIter != strengthFuncs.end(); 
        ++strengthIter)     
      strengthVal += (*strengthIter)->getValueAt(costheta, phi);
    
    return strengthVal;
  }
  
  void ShapeAtomType::complete() {
    
    // first complete all the non-shape atomTypes
    DirectionalAtomType::complete();
    
    int isError = 0;
    
    //setup dipole atom  type in fortran side
    if (isShape()) {
      // vectors for shape transfer to fortran
      std::vector<RealSphericalHarmonic*> tempSHVector;
      std::vector<int> contactL;
      std::vector<int> contactM;
      std::vector<int> contactFunc;
      std::vector<RealType> contactCoeff;
      std::vector<int> rangeL;
      std::vector<int> rangeM;
      std::vector<int> rangeFunc;
      std::vector<RealType> rangeCoeff;
      std::vector<int> strengthL;
      std::vector<int> strengthM;
      std::vector<int> strengthFunc;
      std::vector<RealType> strengthCoeff;
      
      tempSHVector.clear();
      contactL.clear();
      contactM.clear();
      contactFunc.clear();
      contactCoeff.clear();
      
      tempSHVector = getContactFuncs();
      
      int nContact = tempSHVector.size();
      for (int i=0; i<nContact; i++){
        contactL.push_back(tempSHVector[i]->getL());
        contactM.push_back(tempSHVector[i]->getM());
        contactFunc.push_back(tempSHVector[i]->getFunctionType());
        contactCoeff.push_back(tempSHVector[i]->getCoefficient());
      }
      
      tempSHVector.clear();
      rangeL.clear();
      rangeM.clear();
      rangeFunc.clear();
      rangeCoeff.clear();
      
      tempSHVector = getRangeFuncs();
      
      int nRange = tempSHVector.size();
      for (int i=0; i<nRange; i++){
        rangeL.push_back(tempSHVector[i]->getL());
        rangeM.push_back(tempSHVector[i]->getM());
        rangeFunc.push_back(tempSHVector[i]->getFunctionType());
        rangeCoeff.push_back(tempSHVector[i]->getCoefficient());
      }
      
      tempSHVector.clear();
      strengthL.clear();
      strengthM.clear();
      strengthFunc.clear();
      strengthCoeff.clear();
      
      tempSHVector = getStrengthFuncs();
      
      int nStrength = tempSHVector.size();
      for (int i=0; i<nStrength; i++){
        strengthL.push_back(tempSHVector[i]->getL());
        strengthM.push_back(tempSHVector[i]->getM());
        strengthFunc.push_back(tempSHVector[i]->getFunctionType());
        strengthCoeff.push_back(tempSHVector[i]->getCoefficient());
      }
      
      int myATID = getIdent();
      
      makeShape( &nContact, &contactL[0], &contactM[0], &contactFunc[0], 
                 &contactCoeff[0],
                 &nRange, &rangeL[0], &rangeM[0], &rangeFunc[0], 
		 &rangeCoeff[0],
                 &nStrength, &strengthL[0], &strengthM[0], &strengthFunc[0], 
                 &strengthCoeff[0],
                 &myATID, 
                 &isError);
      
      if( isError ){
        sprintf( painCave.errMsg,
                 "Error initializing the \"%s\" shape in fortran\n",
                 (getName()).c_str() );
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }
    }
  }  
}
