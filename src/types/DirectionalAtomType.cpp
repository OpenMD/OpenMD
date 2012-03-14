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

#include <cstdio> 
#include "types/DirectionalAtomType.hpp"
#include "UseTheForce/DarkSide/electrostatic_interface.h"
#include "UseTheForce/DarkSide/sticky_interface.h"
#include "UseTheForce/DarkSide/gb_interface.h"
#include "utils/simError.h"
namespace OpenMD {
  
  DirectionalAtomType::DirectionalAtomType() : AtomType() { 
    myResponsibilities_["is_Directional"] = true; 
    atp.is_Directional = 1; 

    myResponsibilities_["Dipole"] = false;
    myResponsibilities_["SplitDipole"] = false;
    myResponsibilities_["Quadrupole"] = false;
    myResponsibilities_["GayBerne"] = false;
    myResponsibilities_["Sticky"] = false;
    myResponsibilities_["StickyPower"] = false;
    myResponsibilities_["Shape"] = false;

  }   
  
  void DirectionalAtomType::copyAllData(AtomType* at) {
    AtomType::copyAllData(at);
    myResponsibilities_["is_Directional"] = true; 
    atp.is_Directional = 1; 
  }

  void DirectionalAtomType::setDipole() { 
    myResponsibilities_["is_Dipole"] = true;
    atp.is_Dipole = 1; 
  }
  void DirectionalAtomType::setSplitDipole() { 
    myResponsibilities_["is_SplitDipole"] = true;
    atp.is_SplitDipole = 1; 
  }
  void DirectionalAtomType::setQuadrupole() { 
    myResponsibilities_["is_Quadrupole"] = true;
    atp.is_Quadrupole = 1; 
  }
  void DirectionalAtomType::setGayBerne() { 
    myResponsibilities_["is_GayBerne"] = true;
    atp.is_GayBerne = 1; 
  }
  void DirectionalAtomType::setSticky() { 
    myResponsibilities_["is_Sticky"] = true;
    atp.is_Sticky = 1; 
  }
  void DirectionalAtomType::setStickyPower() { 
    myResponsibilities_["is_StickyPower"] = true;
    atp.is_StickyPower = 1; 
  }
  void DirectionalAtomType::setShape() { 
    myResponsibilities_["is_Shape"] = true;
    atp.is_Shape = 1; 
  }

  void DirectionalAtomType::complete() {
    
    AtomType::complete();
    
    int isError  = 0;
    GenericData* data;
    
    //setup dipole atom type in fortran side
    if (isDipole()) {
      data = getPropertyByName("Dipole");
      if (data != NULL) {
        DoubleGenericData* doubleData= dynamic_cast<DoubleGenericData*>(data);
        
        if (doubleData != NULL) {
          RealType dipole = doubleData->getData();
          
          setDipoleMoment(&atp.ident, &dipole, &isError);
          if (isError != 0) {
            sprintf( painCave.errMsg,
                     "Fortran rejected setDipoleMoment\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();          
          }
          
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to DoubleGenericData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
      } else {
        sprintf( painCave.errMsg, "Can not find Dipole Parameters\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
    
    if (isSplitDipole()) {
      data = getPropertyByName("SplitDipoleDistance");
      if (data != NULL) {
        DoubleGenericData* doubleData= dynamic_cast<DoubleGenericData*>(data);
        
        if (doubleData != NULL) {
          RealType splitDipoleDistance = doubleData->getData();
          
          setSplitDipoleDistance(&atp.ident, &splitDipoleDistance, &isError);
          if (isError != 0) {
            sprintf( painCave.errMsg,
                     "Fortran rejected setSplitDipoleDistance\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();          
          }
          
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to DoubleGenericData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
      } else {
        sprintf( painCave.errMsg, "Can not find SplitDipole distance parameter\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
    
    //setup quadrupole atom type in fortran side
    if (isQuadrupole()) {
      data = getPropertyByName("QuadrupoleMoments");
      if (data != NULL) {
        Vector3dGenericData* vector3dData= dynamic_cast<Vector3dGenericData*>(data);
        
        // Quadrupoles in OpenMD are set as the diagonal elements
        // of the diagonalized traceless quadrupole moment tensor.
        // The column vectors of the unitary matrix that diagonalizes 
        // the quadrupole moment tensor become the eFrame (or the
        // electrostatic version of the body-fixed frame.
        
        if (vector3dData != NULL) {
          Vector3d diagElem= vector3dData->getData();
          
          setQuadrupoleMoments(&atp.ident, diagElem.getArrayPointer(), 
                               &isError);
          if (isError != 0) {
            sprintf( painCave.errMsg,
                     "Fortran rejected setQuadrupoleMoments\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();          
          }
          
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to Vector3dGenericData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
      } else {
        sprintf( painCave.errMsg, "Can not find QuadrupoleMoments\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
      
    }
    
    //setup sticky atom type in fortran side
    if (isSticky() || isStickyPower()) {
      data = getPropertyByName("Sticky");
      if (data != NULL) {
        StickyParamGenericData* stickyData = dynamic_cast<StickyParamGenericData*>(data);
        
        if (stickyData != NULL) {
          StickyParam stickyParam = stickyData->getData();
          
          newStickyType(&atp.ident,&stickyParam.w0, &stickyParam.v0, 
                        &stickyParam.v0p, &stickyParam.rl, &stickyParam.ru, 
                        &stickyParam.rlp, &stickyParam.rup, &isError);
          if (isError != 0) {
            sprintf( painCave.errMsg,
                     "Fortran rejected newLJtype\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();          
          }
          
        } else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to StickyParam\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }            
      } else {
        sprintf( painCave.errMsg, "Can not find Parameters for Sticky\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
    
    //setup GayBerne type in fortran side
    if (isGayBerne()) {
      data = getPropertyByName("GayBerne");
      if (data != NULL) {
        GayBerneParamGenericData* gayBerneData = dynamic_cast<GayBerneParamGenericData*>(data);
        
        if (gayBerneData != NULL) {
          GayBerneParam gayBerneParam = gayBerneData->getData();
          
          newGayBerneType(&atp.ident, 
                          &gayBerneParam.GB_d, 
                          &gayBerneParam.GB_l, 
                          &gayBerneParam.GB_eps_X,
                          &gayBerneParam.GB_eps_S, 
                          &gayBerneParam.GB_eps_E, 
                          &gayBerneParam.GB_dw, 
                          &isError);
          
          if (isError != 0) {
            sprintf( painCave.errMsg,
                     "Fortran rejected newGayBerneType\n");
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();          
          }
          
        } 
	
        else {
          sprintf( painCave.errMsg,
                   "Can not cast GenericData to GayBerneParam\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }            
      } 
      else {
        sprintf( painCave.errMsg, "Can not find Parameters for GayBerne\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
  }
} //end namespace OpenMD
