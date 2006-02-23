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
 
#include <iostream>
#include <fstream>
#include <string>

#include "applications/hydrodynamics/HydroCmd.h"
#include "applications/hydrodynamics/HydrodynamicsModel.hpp"
#include "applications/hydrodynamics/HydrodynamicsModelCreator.hpp"
#include "applications/hydrodynamics/HydrodynamicsModelFactory.hpp"
#include "applications/hydrodynamics/BeadModel.hpp"
#include "applications/hydrodynamics/RoughShell.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"

using namespace oopse;

/** Register different hydrodynamics models */
void registerHydrodynamicsModels();

bool calcHydrodynamicsProp(const std::string& modelType, StuntDouble* sd, const DynamicProperty& param, std::ostream& os,  const std::string& prefix);

int main(int argc, char* argv[]){
  //register force fields
  registerForceFields();    
  registerHydrodynamicsModels();
  
  gengetopt_args_info args_info;
  std::string dumpFileName;
  std::string mdFileName;
  std::string prefix;
  
  //parse the command line option
  if (cmdline_parser (argc, argv, &args_info) != 0) {
    exit(1) ;
  }
  
  //get the dumpfile name and meta-data file name
  if (args_info.input_given){
    dumpFileName = args_info.input_arg;
  } else {
    std::cerr << "Does not have input file name" << std::endl;
    exit(1);
  }
  
  mdFileName = dumpFileName;
  mdFileName = mdFileName.substr(0, mdFileName.rfind(".")) + ".md";

  if (args_info.output_given){
    prefix = args_info.output_arg;
  } else {
    prefix = "hydro";    
  }
  std::string outputFilename = prefix + ".diff";

  DynamicProperty param;
  param.insert(DynamicProperty::value_type("Viscosity", args_info.viscosity_arg));
  param.insert(DynamicProperty::value_type("Temperature", args_info.temperature_arg));
  
  if (args_info.sigma_given) {
    param.insert(DynamicProperty::value_type("Sigma", args_info.sigma_arg));
  }
  
  
  //parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(mdFileName, true);
    
  SimInfo::MoleculeIterator mi;
  Molecule* mol;
  Molecule::IntegrableObjectIterator  ii;
  StuntDouble* integrableObject;
  Mat3x3d identMat;
  identMat(0,0) = 1.0;
  identMat(1,1) = 1.0;
  identMat(2,2) = 1.0;
  
 
  std::map<std::string, StuntDouble*> uniqueStuntDoubles;
  
  for (mol = info->beginMolecule(mi); mol != NULL; mol = info->nextMolecule(mi)) {
      for (integrableObject = mol->beginIntegrableObject(ii); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(ii)) {
          if (uniqueStuntDoubles.find(integrableObject->getType()) ==  uniqueStuntDoubles.end()) {
            uniqueStuntDoubles.insert(std::map<std::string, StuntDouble*>::value_type(integrableObject->getType(), integrableObject));
            integrableObject->setPos(V3Zero);
            integrableObject->setA(identMat);
            if (integrableObject->isRigidBody()) {
                RigidBody* rb = static_cast<RigidBody*>(integrableObject);
                rb->updateAtoms();
            }
          }
	}
  }

  std::map<std::string, StuntDouble*>::iterator iter;
  std::ofstream outputDiff(outputFilename.c_str());
  for (iter = uniqueStuntDoubles.begin(); iter != uniqueStuntDoubles.end(); ++iter) {
    calcHydrodynamicsProp(args_info.model_arg, iter->second, param, outputDiff, prefix);
  }
  
  delete info;
   
}

void registerHydrodynamicsModels() {
    HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(new HydrodynamicsModelBuilder<RoughShell>("RoughShell"));
    HydrodynamicsModelFactory::getInstance()->registerHydrodynamicsModel(new HydrodynamicsModelBuilder<BeadModel>("BeadModel"));

}

bool calcHydrodynamicsProp(const std::string& modelType, StuntDouble* sd, const DynamicProperty& param, std::ostream& os, const std::string& prefix) {
    HydrodynamicsModel* hydroModel = HydrodynamicsModelFactory::getInstance()->createHydrodynamicsModel(modelType, sd, param);
    bool ret = false;
    if (hydroModel == NULL) {
        std::cout << "Integrator Factory can not create " << modelType <<std::endl;
    }    

    if (hydroModel->calcHydrodyanmicsProps()) {
        ret = true;        
        hydroModel->writeDiffCenterAndDiffTensor(os);

        std::ofstream ofs;
        std::stringstream outputBeads;
        outputBeads << prefix << "_" << sd->getType() << ".xyz";
        ofs.open(outputBeads.str().c_str());
        hydroModel->writeBeads(ofs);
        ofs.close();
    }

    delete hydroModel;

    return ret;    
}
