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
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/AtomType.hpp"
#include "utils/simError.h"
#define __C
#include "UseTheForce/DarkSide/atype_interface.h"
#include "UseTheForce/DarkSide/lj_interface.h"
#include "UseTheForce/DarkSide/eam_interface.h"
#include "UseTheForce/DarkSide/electrostatic_interface.h"
namespace oopse {
  AtomType::AtomType(){
    
    // initialize to an error:
    atp.ident = -1;

    // and mass_less:
    mass_ = 0.0;
    
    // atom type is a Tabula Rasa:
    atp.is_Directional = 0;
    atp.is_LennardJones = 0;
    atp.is_Charge = 0;
    atp.is_Dipole = 0;
    atp.is_Quadrupole = 0;
    atp.is_Sticky = 0;
    atp.is_GayBerne = 0;
    atp.is_EAM = 0;
    atp.is_Shape = 0;
    atp.is_FLARB = 0;  
  }
    
  void AtomType::makeFortranAtomType() {
    
    int status;

    if (name_.empty()) {
      sprintf( painCave.errMsg,
               "Attempting to complete an AtomType without giving "
               "it a name_!\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    if (atp.ident == -1) {
      sprintf( painCave.errMsg,
               "Attempting to complete AtomType %s without setting the"
               " ident!/n", name_.c_str());
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
 
    status = 0;

    makeAtype(&atp, &status);   
    
    if (status != 0) {
      sprintf( painCave.errMsg,
               "Fortran rejected AtomType %s!\n", name_.c_str());
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
  }


void AtomType::complete() {
    int isError;
    GenericData* data;
    
    //notify a new LJtype atom type is created
    if (isLennardJones()) {
        data = getPropertyByName("LennardJones");
        if (data != NULL) {
            LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);

            if (ljData != NULL) {
                LJParam ljParam = ljData->getData();
                
                newLJtype(&atp.ident, &ljParam.sigma, &ljParam.epsilon, &ljParam.soft_pot, &isError);

                if (isError != 0) {
                    sprintf( painCave.errMsg,
                           "Fortran rejected newLJtype\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
                }
                
            } else {
                    sprintf( painCave.errMsg,
                           "Can not cast GenericData to LJParam\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
            }            
        } else {
            sprintf( painCave.errMsg, "Can not find Parameters for LennardJones\n");
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();          
        }
    }

    if (isElectrostatic()) {
      newElectrostaticType(&atp, &isError);
      if (isError != 0) {
        sprintf( painCave.errMsg,
                 "Fortran rejected newElectrostaticType\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
      
    if (isCharge()) {
        data = getPropertyByName("Charge");
        if (data != NULL) {
            DoubleGenericData* doubleData= dynamic_cast<DoubleGenericData*>(data);

            if (doubleData != NULL) {
                double charge = doubleData->getData();
                setCharge(&atp.ident, &charge, &isError);
                
                if (isError != 0) {
                    sprintf( painCave.errMsg,
                           "Fortran rejected setCharge\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
                }
            } else {
                    sprintf( painCave.errMsg,
                           "Can not cast GenericData to DoubleGenericData\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
            }
        } else {
            sprintf( painCave.errMsg, "Can not find Charge Parameters\n");
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();          
        }
    }

    if (isEAM()) {
        data = getPropertyByName("EAM");
        if (data != NULL) {
            EAMParamGenericData* eamData = dynamic_cast<EAMParamGenericData*>(data);

            if (eamData != NULL) {

                EAMParam eamParam = eamData->getData();
                

                newEAMtype(&eamParam.latticeConstant, &eamParam.nrho, &eamParam.drho,  &eamParam.nr, &eamParam.dr, &eamParam.rcut,
                                &eamParam.rvals[0], &eamParam.rhovals[0], &eamParam.Frhovals[0], &atp.ident, &isError );

                if (isError != 0) {
                    sprintf( painCave.errMsg,
                           "Fortran rejected newEAMtype\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
                }
            } else {
                    sprintf( painCave.errMsg,
                           "Can not cast GenericData to EAMParam\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
            }
        } else {
            sprintf( painCave.errMsg, "Can not find EAM Parameters\n");
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();          
        }
    }

    
}

void AtomType::addProperty(GenericData* genData) {
    properties_.addProperty(genData);  
}

void AtomType::removeProperty(const std::string& propName) {
    properties_.removeProperty(propName);  
}

void AtomType::clearProperties() {
    properties_.clearProperties(); 
}

std::vector<std::string> AtomType::getPropertyNames() {
    return properties_.getPropertyNames();  
}
      
std::vector<GenericData*> AtomType::getProperties() { 
    return properties_.getProperties(); 
}

GenericData* AtomType::getPropertyByName(const std::string& propName) {
    return properties_.getPropertyByName(propName); 
}  
}
