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

#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "utils/simError.h"

#include "applications/dynamicProps/DynamicPropsCmd.h"
#include "applications/dynamicProps/DipoleCorrFunc.hpp"
#include "applications/dynamicProps/RCorrFunc.hpp"
#include "applications/dynamicProps/VCorrFunc.hpp"
#include "applications/dynamicProps/LegendreCorrFunc.hpp"
#include "applications/dynamicProps/RadialRCorrFunc.hpp"
#include "applications/dynamicProps/ThetaCorrFunc.hpp"
#include "applications/dynamicProps/DirectionalRCorrFunc.hpp"
#include "applications/dynamicProps/EnergyCorrFunc.hpp"


using namespace oopse;

int main(int argc, char* argv[]){
  
  //register force fields
  registerForceFields();

  gengetopt_args_info args_info;

  //parse the command line option
  if (cmdline_parser (argc, argv, &args_info) != 0) {
    exit(1) ;
  }


  //get the dumpfile name and meta-data file name
  std::string dumpFileName = args_info.input_arg;
    
  std::string sele1;
  std::string sele2;

  if (args_info.sele1_given) {
    sele1 = args_info.sele1_arg;
  }else {
    char*  sele1Env= getenv("OOPSE_SELE1");
    if (sele1Env) {
      sele1 = sele1Env;
    }else {
      sprintf( painCave.errMsg,
               "neither --sele1 option nor $OOPSE_SELE1 is set");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }
    
  if (args_info.sele2_given) {
    sele2 = args_info.sele2_arg;
  }else {
    char* sele2Env = getenv("OOPSE_SELE2");
    if (sele2Env) {
      sele2 = sele2Env;            
    } else {
      sele2 = sele1;
    }
  }

  //parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName, false);


  TimeCorrFunc* corrFunc;
  if (args_info.dcorr_given){
    corrFunc = new DipoleCorrFunc(info, dumpFileName, sele1, sele2);
  } else if (args_info.rcorr_given) {
    corrFunc = new RCorrFunc(info, dumpFileName, sele1, sele2);
  } else if (args_info.r_rcorr_given) {
    corrFunc = new RadialRCorrFunc(info, dumpFileName, sele1, sele2);
  } else if (args_info.thetacorr_given) {
    corrFunc = new ThetaCorrFunc(info, dumpFileName, sele1, sele2);
  } else if (args_info.drcorr_given) {
    corrFunc = new DirectionalRCorrFunc(info, dumpFileName, sele1, sele2);
  } else if (args_info.vcorr_given) {
    corrFunc = new VCorrFunc(info, dumpFileName, sele1, sele2); 
  } else if (args_info.helfandEcorr_given){
    corrFunc = new EnergyCorrFunc(info, dumpFileName, sele1, sele2);  
  } else if (args_info.lcorr_given) {
    int order;
    if (args_info.order_given)
        order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --lcoor is set\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
        
    corrFunc = new LegendreCorrFunc(info, dumpFileName, sele1, sele2, order); 
  }

  if (args_info.output_given) {
    corrFunc->setOutputName(args_info.output_arg);
  }


  corrFunc->doCorrelate();

  delete corrFunc;    
  delete info;

  return 0;   
}


