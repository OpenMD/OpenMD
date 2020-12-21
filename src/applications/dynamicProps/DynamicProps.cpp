/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"

#include "DynamicPropsCmd.hpp"
#include "applications/dynamicProps/SelectionCorrFunc.hpp"
#include "applications/dynamicProps/DipoleCorrFunc.hpp"
#include "applications/dynamicProps/RCorrFunc.hpp"
#include "applications/dynamicProps/VCorrFunc.hpp"
#include "applications/dynamicProps/LegendreCorrFunc.hpp"
#include "applications/dynamicProps/LegendreCorrFuncZ.hpp"
#include "applications/dynamicProps/ThetaCorrFunc.hpp"
#include "applications/dynamicProps/DirectionalRCorrFunc.hpp"
#include "applications/dynamicProps/StressCorrFunc.hpp"
#include "applications/dynamicProps/SystemDipoleCorrFunc.hpp"
#include "applications/dynamicProps/MomAngMomCorrFunc.hpp"
#include "applications/dynamicProps/ForTorCorrFunc.hpp"
#include "applications/dynamicProps/ForceAutoCorrFunc.hpp"
#include "applications/dynamicProps/TorForCorrFunc.hpp"
#include "applications/dynamicProps/TorqueAutoCorrFunc.hpp"
#include "applications/dynamicProps/cOHz.hpp"
#include "applications/dynamicProps/BondCorrFunc.hpp"
#include "applications/dynamicProps/FreqFlucCorrFunc.hpp"
#include "applications/dynamicProps/HBondJump.hpp"
#include "applications/dynamicProps/HBondPersistence.hpp"
#include "applications/dynamicProps/Displacement.hpp"
#include "applications/dynamicProps/WCorrFunc.hpp"
#include "applications/dynamicProps/ChargeKineticCorrFunc.hpp"
#include "applications/dynamicProps/ChargeOrientationCorrFunc.hpp"
#include "applications/dynamicProps/CurrentDensityAutoCorrFunc.hpp"
#include "applications/dynamicProps/CollectiveDipoleDisplacement.hpp"
#include "applications/dynamicProps/VelocityAutoOutProductCorrFunc.hpp"
#include "applications/dynamicProps/AngularVelocityAutoOutProductCorrFunc.hpp"
#include "applications/dynamicProps/VelAngularVelOutProdCorrFunc.hpp"
#include "applications/dynamicProps/AngularVelVelOutProdCorrFunc.hpp"


using namespace OpenMD;

int main(int argc, char* argv[]){

  struct gengetopt_args_info args_info;

  //parse the command line option

  if (cmdline_parser (argc, argv, &args_info) != 0) {
    cmdline_parser_print_help();
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

  // convert privilegedAxis to corresponding integer
  // x axis -> 0
  // y axis -> 1
  // z axis -> 2 (default)

  int privilegedAxis;
  switch (args_info.privilegedAxis_arg) {
  case privilegedAxis_arg_x:
    privilegedAxis = 0;
    break;
  case privilegedAxis_arg_y:
    privilegedAxis = 1;
    break;
  case privilegedAxis_arg_z:
  default:
    privilegedAxis = 2;
    break;
  }

  // use the memory string to figure out how much memory we can use:
  //char *end;
  //long long int memSize = memparse(args_info.memory_arg, &end);

  // We don't really need to print this out anymore:
  // sprintf( painCave.errMsg,
  //          "Amount of memory being used: %llu bytes\n", memSize);
  // painCave.severity = OPENMD_INFO;
  // painCave.isFatal = 0;
  // simError();

  //parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName, false);

//   DynamicProperty* corrFunc = NULL;
  std::unique_ptr<DynamicProperty> corrFunc {nullptr};

  if(args_info.sdcorr_given){
    corrFunc = MemoryUtils::make_unique<SystemDipoleCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.selecorr_given){
    corrFunc = MemoryUtils::make_unique<SelectionCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.dcorr_given){
    corrFunc = MemoryUtils::make_unique<DipoleCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.rcorr_given) {
    corrFunc = MemoryUtils::make_unique<RCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.r_rcorr_given) {
    corrFunc = MemoryUtils::make_unique<RCorrFuncR>(info, dumpFileName, sele1, sele2);
  } else if (args_info.thetacorr_given) {
    corrFunc = MemoryUtils::make_unique<ThetaCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.drcorr_given) {
    corrFunc = MemoryUtils::make_unique<DirectionalRCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.rcorrZ_given) {
    corrFunc = MemoryUtils::make_unique<RCorrFuncZ>(info, dumpFileName, sele1, sele2,
                              args_info.nzbins_arg, privilegedAxis);
  } else if (args_info.vcorr_given) {
    corrFunc = MemoryUtils::make_unique<VCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.vcorrZ_given) {
    corrFunc = MemoryUtils::make_unique<VCorrFuncZ>(info, dumpFileName, sele1, sele2);
  } else if (args_info.vcorrR_given) {
    corrFunc = MemoryUtils::make_unique<VCorrFuncR>(info, dumpFileName, sele1, sele2);
  } else if (args_info.wcorr_given){
    corrFunc = MemoryUtils::make_unique<WCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.pjcorr_given){
    corrFunc = MemoryUtils::make_unique<MomAngMomCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.ftcorr_given){
    corrFunc = MemoryUtils::make_unique<ForTorCorrFunc>(info, dumpFileName, sele1, sele2);
  }else if (args_info.ckcorr_given){
    corrFunc = MemoryUtils::make_unique<ChargeKineticCorrFunc>(info, dumpFileName, sele1, sele2,
                                         args_info.rcut_arg);
  }else if (args_info.cscorr_given){
    if(args_info.dipoleX_given && args_info.dipoleY_given && args_info.dipoleZ_given){
      corrFunc = MemoryUtils::make_unique<ChargeOrientationCorrFunc>(info, dumpFileName, sele1,
                                               sele2, args_info.dipoleX_arg,
                                               args_info.dipoleY_arg,
                                               args_info.dipoleZ_arg,
                                               args_info.rcut_arg);
    }
  } else if (args_info.facorr_given){
    corrFunc = MemoryUtils::make_unique<ForceAutoCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.tfcorr_given){
    corrFunc = MemoryUtils::make_unique<TorForCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.tacorr_given){
    corrFunc = MemoryUtils::make_unique<TorqueAutoCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.bondcorr_given) {
    corrFunc = MemoryUtils::make_unique<BondCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.stresscorr_given){
    corrFunc = MemoryUtils::make_unique<StressCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.freqfluccorr_given){
    corrFunc = MemoryUtils::make_unique<FreqFlucCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.lcorr_given) {
    int order(0);
    if (args_info.order_given)
      order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --lcorr is set\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    corrFunc = MemoryUtils::make_unique<LegendreCorrFunc>(info, dumpFileName, sele1, sele2, order);
  } else if (args_info.lcorrZ_given) {
    int order(0);
    if (args_info.order_given)
      order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --lcorrZ is set\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    corrFunc = MemoryUtils::make_unique<LegendreCorrFuncZ>(info, dumpFileName, sele1, sele2, order,
				     args_info.nzbins_arg, privilegedAxis);

  } else if (args_info.cohZ_given) {
    int order(0);
    if (args_info.order_given)
      order = args_info.order_arg;
    else {
      sprintf( painCave.errMsg,
               "--order must be set if --cohZ is set\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    corrFunc = MemoryUtils::make_unique<COHZ>(info, dumpFileName, sele1, sele2, order,
			args_info.nzbins_arg, privilegedAxis);

  } else if (args_info.jumptime_given) {
    corrFunc = MemoryUtils::make_unique<HBondJump>(info, dumpFileName, sele1, sele2,
                             args_info.OOcut_arg,
                             args_info.thetacut_arg,
                             args_info.OHcut_arg);
  } else if (args_info.jumptimeZ_given) {
    corrFunc = MemoryUtils::make_unique<HBondJumpZ>(info, dumpFileName, sele1, sele2,
                              args_info.OOcut_arg,
                              args_info.thetacut_arg,
                              args_info.OHcut_arg, args_info.nzbins_arg, privilegedAxis);
  } else if (args_info.persistence_given) {
    corrFunc = MemoryUtils::make_unique<HBondPersistence>(info, dumpFileName, sele1, sele2,
                                    args_info.OOcut_arg,
                                    args_info.thetacut_arg,
                                    args_info.OHcut_arg);
  } else if (args_info.disp_given) {
    corrFunc = MemoryUtils::make_unique<Displacement>(info, dumpFileName, sele1, sele2);
  } else if (args_info.dispZ_given) {
    corrFunc = MemoryUtils::make_unique<DisplacementZ>(info, dumpFileName, sele1, sele2,
                                 args_info.nzbins_arg, privilegedAxis);
  } else if (args_info.current_given) {
    corrFunc = MemoryUtils::make_unique<CurrentDensityAutoCorrFunc>(info, dumpFileName, sele1, sele2);
  } else if (args_info.ddisp_given) {
    corrFunc = MemoryUtils::make_unique<CollectiveDipoleDisplacement>(info, dumpFileName, sele1, sele2);
  }
    else if (args_info.vaOutProdcorr_given){
    corrFunc = MemoryUtils::make_unique<VelocityAutoOutProductCorrFunc>(info, dumpFileName, sele1, sele2);
  }
    else if (args_info.waOutProdcorr_given){
    corrFunc = MemoryUtils::make_unique<AngularVelocityAutoOutProductCorrFunc>(info, dumpFileName, sele1, sele2);
  }
    else if (args_info.vwOutProdcorr_given){
    corrFunc = MemoryUtils::make_unique<VelAngularVelOutProdCorrFunc>(info, dumpFileName, sele1, sele2);
  }
  else if (args_info.wvOutProdcorr_given){
  corrFunc = MemoryUtils::make_unique<AngularVelVelOutProdCorrFunc>(info, dumpFileName, sele1, sele2);
  }

  if (args_info.output_given) {
    corrFunc->setOutputName(args_info.output_arg);
  }

  corrFunc->doCorrelate();

  delete info;
  return 0;
}
