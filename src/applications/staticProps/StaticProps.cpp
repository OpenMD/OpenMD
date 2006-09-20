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
#include "io/DumpReader.hpp"
#include "utils/simError.h"

#include "applications/staticProps/StaticPropsCmd.h"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "applications/staticProps/GofR.hpp"
#include "applications/staticProps/GofRAngle.hpp"
#include "applications/staticProps/GofAngle2.hpp"
#include "applications/staticProps/GofXyz.hpp"
#include "applications/staticProps/P2OrderParameter.hpp"
#include "applications/staticProps/BondOrderParameter.hpp"
#include "applications/staticProps/RippleOP.hpp"
#include "applications/staticProps/SCDOrderParameter.hpp"
#include "applications/staticProps/DensityPlot.hpp"
#include "applications/staticProps/RhoZ.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "applications/staticProps/Hxy.hpp"
#endif

using namespace oopse;

int main(int argc, char* argv[]){
  
  //register force fields
  registerForceFields();

  gengetopt_args_info args_info;

  //parse the command line option
  if (cmdline_parser (argc, argv, &args_info) != 0) {
    exit(1) ;
  }

  //get the dumpfile name
  std::string dumpFileName = args_info.input_arg;
  std::string sele1;
  std::string sele2;
  bool userSpecifiedSelect1;
  bool userSpecifiedSelect2;

  // check the first selection argument, or set it to the environment
  // variable, or failing that, set it to "select all"

  if (args_info.sele1_given) {
    sele1 = args_info.sele1_arg;
  } else {
    char*  sele1Env= getenv("OOPSE_SELE1");
    if (sele1Env) {
      sele1 = sele1Env;
    } else {
      sele1 = "select all";
    }
  }

  // check the second selection argument, or set it to the environment
  // variable, or failing that, set it to "select all"
  
  if (args_info.sele2_given) {
    sele2 = args_info.sele2_arg;
  } else {
    char* sele2Env = getenv("OOPSE_SELE2");
    if (sele2Env) {
      sele2 = sele2Env;            
    } else { 
      sele2 = "select all";
    }
  }


  // Problems if sele1 wasn't specified, but 
// if (!args_info.scd_given) {
//       sprintf( painCave.errMsg,
//                "neither --sele1 option nor $OOPSE_SELE1 is set");
//       painCave.severity = OOPSE_ERROR;
//       painCave.isFatal = 1;
//       simError();
//     }
//   }

  // Problems if sele1 wasn't specified

//     if(!args_info.scd_given && !args_info.density_given && !args_info.slab_density_given)  {
//       sprintf( painCave.errMsg,
//                "neither --sele2 option nor $OOPSE_SELE2 is set");
//       painCave.severity = OOPSE_ERROR;
//       painCave.isFatal = 1;
//       simError();        
//     }
//   }

  bool batchMode;
  if (args_info.scd_given){
    if (args_info.sele1_given && args_info.sele2_given && args_info.sele3_given) {
      batchMode = false;
    } else if (args_info.molname_given && args_info.begin_given && args_info.end_given) {
      if (args_info.begin_arg < 0 || args_info.end_arg < 0 || args_info.begin_arg > args_info.end_arg-2) {
        sprintf( painCave.errMsg,
                 "below conditions are not satisfied:\n"
                 "0 <= begin && 0<= end && begin <= end-2\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();                    
      }
      batchMode = true;        
    } else{
      sprintf( painCave.errMsg,
               "either --sele1, --sele2, --sele3 are specified,"
               " or --molname, --begin, --end are specified\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();        
    
    }
  }

  //parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName);

  RealType maxLen;
  if (args_info.length_given) {
    maxLen = args_info.length_arg;
  } else {
    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    maxLen = std::min(std::min(hmat(0, 0), hmat(1, 1)), hmat(2, 2)) /2.0;        
  }    

  StaticAnalyser* analyser;
  if (args_info.gofr_given){
    analyser= new GofR(info, dumpFileName, sele1, sele2, maxLen, 
                       args_info.nbins_arg);        
  } else if (args_info.r_theta_given) {
    analyser  = new GofRTheta(info, dumpFileName, sele1, sele2, maxLen, 
                              args_info.nbins_arg, args_info.nanglebins_arg);
  } else if (args_info.r_omega_given) {
    analyser  = new GofROmega(info, dumpFileName, sele1, sele2, maxLen, 
                              args_info.nbins_arg, args_info.nanglebins_arg);
  } else if (args_info.theta_omega_given) {
    analyser  = new GofAngle2(info, dumpFileName, sele1, sele2, 
                              args_info.nanglebins_arg);
  } else if (args_info.gxyz_given) {
    if (args_info.refsele_given) {
      analyser= new GofXyz(info, dumpFileName, sele1, sele2,args_info.refsele_arg, 
                           maxLen, args_info.nbins_arg);        
    } else {
      sprintf( painCave.errMsg,
               "--refsele must set when --gxyz is used");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();  
    }
  } else if (args_info.p2_given) {
    analyser  = new P2OrderParameter(info, dumpFileName, sele1, sele2);
  } else if (args_info.rp2_given){
    analyser = new RippleOP(info, dumpFileName, sele1, sele2);
  } else if (args_info.bo_given){
    if (args_info.rcut_given && args_info.LegendreL_given) {
      analyser = new BondOrderParameter(info, dumpFileName, sele1, 
                                        args_info.rcut_arg, 
                                        args_info.LegendreL_arg,
                                        args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
               "Both the cutoff radius (rcut) and LegendreL must be specified when calculating Bond Order Parameters");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.scd_given) {
    if (batchMode) {
      analyser  = new SCDOrderParameter(info, dumpFileName, args_info.molname_arg, 
                                        args_info.begin_arg, args_info.end_arg);
    } else{
      std::string sele3 = args_info.sele3_arg;
      analyser  = new SCDOrderParameter(info, dumpFileName, sele1, sele2, sele3);
    }
  }else if (args_info.density_given) {
    analyser= new DensityPlot(info, dumpFileName, sele1, sele2, maxLen, 
                              args_info.nbins_arg);  
  } else if (args_info.slab_density_given) {
    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    analyser = new RhoZ(info, dumpFileName, sele1, hmat(2, 2), args_info.nbins_arg);
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
  }else if (args_info.hxy_given) {
    analyser = new Hxy(info, dumpFileName, sele1, args_info.nbins_x_arg, 
                       args_info.nbins_y_arg, args_info.nbins_arg);
#endif
  }
  
  if (args_info.output_given) {
    analyser->setOutputName(args_info.output_arg);
  }
  if (args_info.step_given) {
    analyser->setStep(args_info.step_arg);
  }

  analyser->process();

  delete analyser;    
  delete info;

  return 0;   
}

