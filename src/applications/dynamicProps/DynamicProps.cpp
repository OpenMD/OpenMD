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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <iostream>
#include <fstream>
#include <string>

#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

#include "applications/dynamicProps/DynamicPropsCmd.h"
#include "applications/dynamicProps/SelectionCorrFunc.hpp"
#include "applications/dynamicProps/DipoleCorrFunc.hpp"
#include "applications/dynamicProps/RCorrFunc.hpp"
#include "applications/dynamicProps/VCorrFunc.hpp"
#include "applications/dynamicProps/LegendreCorrFunc.hpp"
#include "applications/dynamicProps/LegendreCorrFuncZ.hpp"
#include "applications/dynamicProps/RadialRCorrFunc.hpp"
#include "applications/dynamicProps/ThetaCorrFunc.hpp"
#include "applications/dynamicProps/DirectionalRCorrFunc.hpp"
#include "applications/dynamicProps/EnergyCorrFunc.hpp"
#include "applications/dynamicProps/StressCorrFunc.hpp"
#include "applications/dynamicProps/SystemDipoleCorrFunc.hpp"
#include "applications/dynamicProps/MomentumCorrFunc.hpp"
#include "applications/dynamicProps/cOHz.hpp"
#include "applications/dynamicProps/BondCorrFunc.hpp"
#include "applications/dynamicProps/FreqFlucCorrFunc.hpp"

using namespace OpenMD;

int main(int argc, char* argv[]){
  
  gengetopt_args_info args_info;

  //parse the command line option
  if (cmdline_parser (argc, argv, &args_info) != 0) {
    exit(1) ;
  }

  //get the dumpfile name and meta-data file name
  std::string dumpFileName = args_info.input_arg;
    
  std::string sele1;
  std::string sele2;
  
  // check the first selection argument, or set it to the environment
  // variable, or failing that, set it to "select all"
  
  if (args_info.sele1_given) {
    sele1 = args_info.sele1_arg;
  } else {
    char*  sele1Env= getenv("SELECTION1");
    if (sele1Env) {
      sele1 = sele1Env;
    } else {
      sele1 = "select all";
    }
  }
  
  // check the second selection argument, or set it to the environment
  // variable, or failing that, set it to the first selection
  
  if (args_info.sele2_given) {
    sele2 = args_info.sele2_arg;
  } else {
    char* sele2Env = getenv("SELECTION2");
    if (sele2Env) {
      sele2 = sele2Env;            
    } else { 
      //If sele2 is not specified, then the default behavior
      //should be what is already intended for sele1
      sele2 = sele1;
    }
  }

  // use the memory string to figure out how much memory we can use:
  char *end;
  long long int memSize = memparse(args_info.memory_arg, &end);
  sprintf( painCave.errMsg,
           "Amount of memory being used: %llu bytes\n", memSize);
  painCave.severity = OPENMD_INFO;
  painCave.isFatal = 0;
  simError();     

  //parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName, false);

  TimeCorrFunc* corrFunc;
  if(args_info.sdcorr_given){
    corrFunc = new SystemDipoleCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.selecorr_given){
    corrFunc = new SelectionCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.dcorr_given){
    corrFunc = new DipoleCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.rcorr_given) {
    corrFunc = new RCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.r_rcorr_given) {
    corrFunc = new RadialRCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.thetacorr_given) {
    corrFunc = new ThetaCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.drcorr_given) {
    corrFunc = new DirectionalRCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.vcorr_given) {
    corrFunc = new VCorrFunc(info, dumpFileName, sele1, sele2, memSize); 
  } else if (args_info.vcorrZ_given) {
    corrFunc = new VCorrFuncZ(info, dumpFileName, sele1, sele2, memSize); 
  } else if (args_info.vcorrR_given) {
    corrFunc = new VCorrFuncR(info, dumpFileName, sele1, sele2, memSize); 
  } else if (args_info.bondcorr_given) {
    corrFunc = new BondCorrFunc(info, dumpFileName, sele1, memSize); 
  } else if (args_info.helfandEcorr_given){
    corrFunc = new EnergyCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.stresscorr_given){
    corrFunc = new StressCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.momentum_given){
    corrFunc = new MomentumCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.freqfluccorr_given){
    corrFunc = new FreqFlucCorrFunc(info, dumpFileName, sele1, sele2, memSize);
  } else if (args_info.lcorr_given) {
    int order;
    if (args_info.order_given)
        order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --lcorr is set\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
        
    corrFunc = new LegendreCorrFunc(info, dumpFileName, sele1, sele2, order, memSize); 
  } else if (args_info.lcorrZ_given) {
    int order;
    if (args_info.order_given)
        order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --lcorrZ is set\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
        
    corrFunc = new LegendreCorrFuncZ(info, dumpFileName, sele1, sele2, order, args_info.nzbins_arg, memSize); 

  } else if (args_info.cohZ_given) {
    int order;
    if (args_info.order_given)
        order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --cohZ is set\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
        
    corrFunc = new COHZ(info, dumpFileName, sele1, sele2, order, args_info.nzbins_arg, memSize); 

  }


  if (args_info.output_given) {
    corrFunc->setOutputName(args_info.output_arg);
  }


  corrFunc->doCorrelate();

  delete corrFunc;    
  delete info;

  return 0;   
}


