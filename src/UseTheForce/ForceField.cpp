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
 
/**
 * @file ForceField.cpp
 * @author tlin
 * @date 11/04/2004
 * @time 22:51am
 * @version 1.0
 */
  
#include "UseTheForce/ForceField.hpp"
#include "utils/simError.h"
#include "UseTheForce/DarkSide/atype_interface.h"
#include "UseTheForce/DarkSide/fForceOptions_interface.h"
#include "UseTheForce/DarkSide/switcheroo_interface.h"
namespace oopse {

  ForceField::ForceField() { 
    char* tempPath; 
    tempPath = getenv("FORCE_PARAM_PATH");

    if (tempPath == NULL) {
      //convert a macro from compiler to a string in c++
      STR_DEFINE(ffPath_, FRC_PATH );
    } else {
      ffPath_ = tempPath;
    }
  }


  ForceField::~ForceField() {
    deleteAtypes();
    deleteSwitch();
  }

  AtomType* ForceField::getAtomType(const std::string &at) {
    std::vector<std::string> keys;
    keys.push_back(at);
    return atomTypeCont_.find(keys);
  }

  BondType* ForceField::getBondType(const std::string &at1, const std::string &at2) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    

    //try exact match first
    BondType* bondType = bondTypeCont_.find(keys);
    if (bondType) {
      return bondType;
    } else {
      //if no exact match found, try wild card match
      return bondTypeCont_.find(keys, wildCardAtomTypeName_);
    }

  }

  BendType* ForceField::getBendType(const std::string &at1, const std::string &at2,
				    const std::string &at3) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    

    //try exact match first
    BendType* bendType = bendTypeCont_.find(keys);
    if (bendType) {
      return bendType;
    } else {
      //if no exact match found, try wild card match
      return bendTypeCont_.find(keys, wildCardAtomTypeName_);
    }
  }

  TorsionType* ForceField::getTorsionType(const std::string &at1, const std::string &at2,
					  const std::string &at3, const std::string &at4) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);    

    TorsionType* torsionType = torsionTypeCont_.find(keys);
    if (torsionType) {
      return torsionType;
    } else {
      //if no exact match found, try wild card match
      return torsionTypeCont_.find(keys, wildCardAtomTypeName_);
    }
    
    return torsionTypeCont_.find(keys, wildCardAtomTypeName_);

  }

  BondType* ForceField::getExactBondType(const std::string &at1, const std::string &at2){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    return bondTypeCont_.find(keys);
  }

  BendType* ForceField::getExactBendType(const std::string &at1, const std::string &at2,
					 const std::string &at3){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    return bendTypeCont_.find(keys);
  }

  TorsionType* ForceField::getExactTorsionType(const std::string &at1, const std::string &at2,
					       const std::string &at3, const std::string &at4){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);   
    return torsionTypeCont_.find(keys);
  }
  bool ForceField::addAtomType(const std::string &at, AtomType* atomType) {
    std::vector<std::string> keys;
    keys.push_back(at);
    return atomTypeCont_.add(keys, atomType);
  }

  bool ForceField::addBondType(const std::string &at1, const std::string &at2, BondType* bondType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    return bondTypeCont_.add(keys, bondType);

  }

  bool ForceField::addBendType(const std::string &at1, const std::string &at2,
			       const std::string &at3, BendType* bendType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    return bendTypeCont_.add(keys, bendType);
  }

  bool ForceField::addTorsionType(const std::string &at1, const std::string &at2,
				  const std::string &at3, const std::string &at4, TorsionType* torsionType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);    
    return torsionTypeCont_.add(keys, torsionType);
  }

  double ForceField::getRcutFromAtomType(AtomType* at) {
    /**@todo */
    GenericData* data;
    double rcut = 0.0;

    if (at->isLennardJones()) {
      data = at->getPropertyByName("LennardJones");
      if (data != NULL) {
	LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);

	if (ljData != NULL) {
	  LJParam ljParam = ljData->getData();

	  //by default use 2.5*sigma as cutoff radius
	  rcut = 2.5 * ljParam.sigma;
                
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

    return rcut;    
  }


  ifstrstream* ForceField::openForceFieldFile(const std::string& filename) {
    std::string forceFieldFilename(filename);
    ifstrstream* ffStream = new ifstrstream();
    
    //try to open the force filed file in current directory first    
    ffStream->open(forceFieldFilename.c_str());
    if(!ffStream->is_open()){

      forceFieldFilename = ffPath_ + "/" + forceFieldFilename;
      ffStream->open( forceFieldFilename.c_str() );

      //if current directory does not contain the force field file,
      //try to open it in the path        
      if(!ffStream->is_open()){

	sprintf( painCave.errMsg,
		 "Error opening the force field parameter file:\n"
		 "\t%s\n"
		 "\tHave you tried setting the FORCE_PARAM_PATH environment "
		 "variable?\n",
		 forceFieldFilename.c_str() );
	painCave.severity = OOPSE_ERROR;
	painCave.isFatal = 1;
	simError();
      }
    }  

    return ffStream;

  }

  void ForceField::setFortranForceOptions(){
    ForceOptions theseFortranOptions;
    forceFieldOptions_.makeFortranOptions(theseFortranOptions);
    setfForceOptions(&theseFortranOptions);
  }
} //end namespace oopse
