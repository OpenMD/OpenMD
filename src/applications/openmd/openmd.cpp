/*
 * Copyright (c) 2010 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2010).                        
 */
 
#ifdef IS_MPI
#include <mpi.h>
#endif

#include <fstream>
#include "utils/simError.h"
#include "utils/CaseConversion.hpp"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "constraints/ZconstraintForceManager.hpp"
#include "restraints/RestraintForceManager.hpp"
#include "integrators/IntegratorFactory.hpp"
#include "integrators/Integrator.hpp"
#include "minimizers/MinimizerFactory.hpp"
#include "minimizers/Minimizer.hpp"
#include "restraints/ThermoIntegrationForceManager.hpp"
#include "utils/Report.hpp"
#include "utils/Exception.hpp"

using namespace OpenMD;

int main(int argc,char* argv[]){
  
  // first things first, all of the initializations

#ifdef IS_MPI
  MPI_Init( &argc, &argv ); // the MPI communicators
#endif
   
  initSimError();           // the error handler
  srand48( 1337 );          // the random number generator.

  std::string svnrev;
  //convert a macro from compiler to a string in c++
  STR_DEFINE(svnrev, SVN_REV );

  std::string revision;

  if (!svnrev.empty()) {
     revision.assign("  Revision: " + svnrev);
  }

  revision.resize(19,' ');

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif
    std::cerr << 
      "  +-----------------------------------------------------------------------+\n"<<
      "  |    ____                    __  ___ ____                               |\n"<<
      "  |   / __ \\____  ___  ____   /  |/  // __ \\  The Open Molecular Dynamics |\n"<<
      "  |  / / / / __ \\/ _ \\/ __ \\ / /|_/ // / / /  Engine (formerly OOPSE).    |\n"<<
      "  | / /_/ / /_/ /  __/ / / // /  / // /_/ /                               |\n"<<
      "  | \\____/ .___/\\___/_/ /_//_/  /_//_____/    Copyright 2004-2010 by the  |\n"<<
      "  |     /_/                                   University of Notre Dame.   |\n"<<
      "  |                                                                       |\n"<<
      "  |        version " << 
      OPENMD_VERSION_MAJOR << "." << OPENMD_VERSION_MINOR << revision << 
      "     http://www.openmd.net       |\n"<<
      "  |                                                                       |\n"<<
      "  | OpenMD is an OpenScience project.  All source code is available for   |\n"<<
      "  | any use whatsoever under a BSD-style license.                         |\n"<<
      "  |                                                                       |\n"<<
      "  | Support OpenScience!  If you use OpenMD or its source code in your    |\n"<<
      "  | research, please cite the appropriate papers when you publish your    |\n"<<
      "  | work.  Good starting points are:                                      |\n"<<
      "  |                                                                       |\n"<<
      "  | [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).              |\n"<<
      "  | [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).           |\n"<<
      "  | [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          |\n"<<
      "  | [4]  Vardeman & Gezelter, in progress (2010).                         |\n"<<
      "  +-----------------------------------------------------------------------+\n"<<
      "\n";
    
    if( argc < 2 ){
      strcpy( painCave.errMsg, "Error, a meta-data file is needed to run.\n" );
      painCave.isFatal = 1;
      simError();
    }
#ifdef IS_MPI
  }
#endif
  
  strcpy( checkPointMsg, "Successful number of arguments" );
  errorCheckPoint();

  //register forcefields, integrators and minimizers
  registerAll();

  //create simulation model
  SimCreator creator;
  SimInfo* info = creator.createSim(argv[1]);

  Globals* simParams = info->getSimParams();

  if (simParams->haveMinimizer() && simParams->haveEnsemble()) {
    sprintf(painCave.errMsg, "Minimizer keyword and Ensemble keyword can not exist together\n");
    painCave.isFatal = 1;
    simError();        
  }
  

  if (simParams->haveMinimizer()) {
    //create minimizer
    Minimizer* myMinimizer = MinimizerFactory::getInstance()->createMinimizer(toUpperCopy(simParams->getMinimizer()), info);

    if (myMinimizer == NULL) {
      sprintf(painCave.errMsg, "Minimizer Factory can not create %s Minimizer\n",
	      simParams->getMinimizer().c_str());
      painCave.isFatal = 1;
      simError();
    }

    myMinimizer->minimize();
    delete myMinimizer;
  } else if (simParams->haveEnsemble()) {
    //create Integrator

    Integrator* myIntegrator = IntegratorFactory::getInstance()->createIntegrator(toUpperCopy(simParams->getEnsemble()), info);
 
    if (myIntegrator == NULL) {
      sprintf(painCave.errMsg, "Integrator Factory can not create %s Integrator\n",
	      simParams->getEnsemble().c_str());
      painCave.isFatal = 1;
      simError();
    }
                
    //Thermodynamic Integration Method
    //set the force manager for thermodynamic integration if specified
    if (simParams->getUseThermodynamicIntegration()){
      ForceManager* fman = new ThermoIntegrationForceManager(info);
      myIntegrator->setForceManager(fman);
    }

    // Restraints
    if (simParams->getUseRestraints() && !simParams->getUseThermodynamicIntegration()) {
      ForceManager* fman = new RestraintForceManager(info);
      myIntegrator->setForceManager(fman);
    }

    //Zconstraint-Method
    if (simParams->getNZconsStamps() > 0) {
      info->setNZconstraint(simParams->getNZconsStamps());
      ForceManager* fman = new ZconstraintForceManager(info);
      myIntegrator->setForceManager(fman);
    }
        
    myIntegrator->integrate();
    delete myIntegrator;
  }else {
    sprintf(painCave.errMsg, "Integrator Factory can not create %s Integrator\n",
            simParams->getEnsemble().c_str());
    painCave.isFatal = 1;
    simError();
  }
    
  delete info;


  strcpy( checkPointMsg, "Great googly moogly!  It worked!" );
  errorCheckPoint();

#ifdef IS_MPI  
  MPI_Finalize();
#endif

  return 0 ;
}
