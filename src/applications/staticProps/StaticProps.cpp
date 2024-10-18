/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "StaticPropsCmd.hpp"
#include "applications/staticProps/BOPofR.hpp"
#include "applications/staticProps/BondAngleDistribution.hpp"
#include "applications/staticProps/BondOrderParameter.hpp"
#include "applications/staticProps/DensityPlot.hpp"
#include "applications/staticProps/GofAngle2.hpp"
#include "applications/staticProps/GofR.hpp"
#include "applications/staticProps/GofRAngle.hpp"
#include "applications/staticProps/GofRAngle2.hpp"
#include "applications/staticProps/GofRZ.hpp"
#include "applications/staticProps/GofXyz.hpp"
#include "applications/staticProps/GofZ.hpp"
#include "applications/staticProps/KirkwoodBuff.hpp"
#include "applications/staticProps/NanoLength.hpp"
#include "applications/staticProps/NanoVolume.hpp"
#include "applications/staticProps/ObjectCount.hpp"
#include "applications/staticProps/P2OrderParameter.hpp"
#include "applications/staticProps/P2R.hpp"
#include "applications/staticProps/PipeDensity.hpp"
#include "applications/staticProps/RhoZ.hpp"
#include "applications/staticProps/RippleOP.hpp"
#include "applications/staticProps/SCDOrderParameter.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "applications/staticProps/TwoDGofR.hpp"
#include "applications/staticProps/pAngle.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "utils/simError.h"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "applications/staticProps/Hxy.hpp"
#endif
#include "applications/staticProps/AngleR.hpp"
#include "applications/staticProps/ChargeDensityZ.hpp"
#include "applications/staticProps/ChargeHistogram.hpp"
#include "applications/staticProps/ChargeR.hpp"
#include "applications/staticProps/ChargeZ.hpp"
#include "applications/staticProps/CoordinationNumber.hpp"
#include "applications/staticProps/CurrentDensity.hpp"
#include "applications/staticProps/DensityHistogram.hpp"
#include "applications/staticProps/DipoleOrientation.hpp"
#include "applications/staticProps/Field.hpp"
#include "applications/staticProps/HBondGeometric.hpp"
#include "applications/staticProps/HBondR.hpp"
#include "applications/staticProps/HBondRvol.hpp"
#include "applications/staticProps/HBondZ.hpp"
#include "applications/staticProps/HBondZvol.hpp"
#include "applications/staticProps/Kirkwood.hpp"
#include "applications/staticProps/MassDensityR.hpp"
#include "applications/staticProps/MassDensityZ.hpp"
#include "applications/staticProps/MomentumHistogram.hpp"
#include "applications/staticProps/MultipoleSum.hpp"
#include "applications/staticProps/NitrileFrequencyMap.hpp"
#include "applications/staticProps/NumberR.hpp"
#include "applications/staticProps/NumberZ.hpp"
#include "applications/staticProps/OrderParameterProbZ.hpp"
#include "applications/staticProps/PositionZ.hpp"
#include "applications/staticProps/PotDiff.hpp"
#include "applications/staticProps/RNEMDStats.hpp"
#include "applications/staticProps/RhoR.hpp"
#include "applications/staticProps/SurfaceDiffusion.hpp"
#include "applications/staticProps/TetrahedralityHBMatrix.hpp"
#include "applications/staticProps/TetrahedralityParam.hpp"
#include "applications/staticProps/TetrahedralityParamDens.hpp"
#include "applications/staticProps/TetrahedralityParamR.hpp"
#include "applications/staticProps/TetrahedralityParamXYZ.hpp"
#include "applications/staticProps/TetrahedralityParamZ.hpp"
#include "applications/staticProps/TranslationalOrderParamZ.hpp"
#include "applications/staticProps/VelocityZ.hpp"

using namespace OpenMD;

int main(int argc, char* argv[]) {
  gengetopt_args_info args_info;

  // parse the command line option
  if (cmdline_parser(argc, argv, &args_info) != 0) { exit(1); }

  // get the dumpfile name
  std::string dumpFileName = args_info.input_arg;
  std::string sele1;
  std::string sele2;
  std::string sele3;
  std::string comsele;

  // check the first selection argument, or set it to the environment
  // variable, or failing that, set it to "select all"

  if (args_info.sele1_given) {
    sele1 = args_info.sele1_arg;
  } else {
    char* sele1Env = getenv("SELECTION1");
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
      // If sele2 is not specified, then the default behavior
      // should be what is already intended for sele1
      sele2 = sele1;
    }
  }

  // check the third selection argument, which is only set if
  // requested by the user

  if (args_info.sele3_given) sele3 = args_info.sele3_arg;

  // check the comsele selection argument, which is only set if
  // requested by the user

  if (args_info.comsele_given) comsele = args_info.comsele_arg;

  bool batchMode(false);
  if (args_info.scd_given) {
    if (args_info.sele1_given && args_info.sele2_given &&
        args_info.sele3_given) {
      batchMode = false;
    } else if (args_info.molname_given && args_info.begin_given &&
               args_info.end_given) {
      if (args_info.begin_arg < 0 || args_info.end_arg < 0 ||
          args_info.begin_arg > args_info.end_arg - 2) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "below conditions are not satisfied:\n"
                 "0 <= begin && 0<= end && begin <= end-2\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
      batchMode = true;
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "either --sele1, --sele2, --sele3 are specified,"
               " or --molname, --begin, --end are specified\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  }

  // parse md file and set up the system
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName);

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

  int privilegedAxis2;
  switch (args_info.privilegedAxis2_arg) {
  case privilegedAxis2_arg_x:
    privilegedAxis2 = 0;
    break;
  case privilegedAxis2_arg_y:
    privilegedAxis2 = 1;
    break;
  case privilegedAxis2_arg_z:
  default:
    privilegedAxis2 = 2;
    break;
  }

  RealType maxLen;
  RealType zmaxLen(0.0);
  if (args_info.length_given) {
    maxLen = args_info.length_arg;
    if (args_info.zlength_given) { zmaxLen = args_info.zlength_arg; }
  } else {
    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    // The maximum length for radial distribution functions is actually half
    // the smallest box length
    maxLen = std::min(std::min(hmat(0, 0), hmat(1, 1)), hmat(2, 2)) / 2.0;
    // This should be the extent of the privileged axis:
    zmaxLen = hmat(privilegedAxis, privilegedAxis);
  }

  int nanglebins, nrbins;
  // in case we override nbins with nrbins:
  if (args_info.nrbins_given) {
    nrbins = args_info.nrbins_arg;
  } else {
    nrbins = args_info.nbins_arg;
  }
  // in case we override nbins with nanglebins:
  if (args_info.nanglebins_given) {
    nanglebins = args_info.nanglebins_arg;
  } else {
    nanglebins = args_info.nbins_arg;
  }

  RealType binWidth = args_info.binWidth_arg;

  // override default vander waals radius for fictious atoms in a model
  RealType vRadius;
  if (args_info.v_radius_given) {
    vRadius = args_info.v_radius_arg;
  } else {
    vRadius = 1.52;
  }

  int momentum_type;
  switch (args_info.momentum_arg) {
  case momentum_arg_P:
    momentum_type = 0;
    break;
  case momentum_arg_J:
  default:
    momentum_type = 1;
    break;
  }

  int momentum_comp;
  switch (args_info.component_arg) {
  case component_arg_x:
    momentum_comp = 0;
    break;
  case component_arg_y:
    momentum_comp = 1;
    break;
  case component_arg_z:
  default:
    momentum_comp = 2;
    break;
  }

  std::unique_ptr<StaticAnalyser> analyser {nullptr};

  if (args_info.gofr_given) {
    analyser = std::make_unique<GofR>(info, dumpFileName, sele1, sele2, maxLen,
                                      nrbins);
  } else if (args_info.gofz_given) {
    analyser =
        std::make_unique<GofZ>(info, dumpFileName, sele1, sele2, maxLen,
                               zmaxLen, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.r_z_given) {
    analyser = std::make_unique<GofRZ>(info, dumpFileName, sele1, sele2, maxLen,
                                       zmaxLen, nrbins, args_info.nbins_z_arg,
                                       privilegedAxis);
  } else if (args_info.r_theta_given) {
    if (args_info.sele3_given)
      analyser = std::make_unique<GofRTheta>(info, dumpFileName, sele1, sele2,
                                             sele3, maxLen, nrbins, nanglebins);
    else
      analyser = std::make_unique<GofRTheta>(info, dumpFileName, sele1, sele2,
                                             maxLen, nrbins, nanglebins);
  } else if (args_info.r_omega_given) {
    if (args_info.sele3_given)
      analyser = std::make_unique<GofROmega>(info, dumpFileName, sele1, sele2,
                                             sele3, maxLen, nrbins, nanglebins);
    else
      analyser = std::make_unique<GofROmega>(info, dumpFileName, sele1, sele2,
                                             maxLen, nrbins, nanglebins);
  } else if (args_info.theta_omega_given) {
    if (args_info.sele3_given)
      analyser = std::make_unique<GofAngle2>(info, dumpFileName, sele1, sele2,
                                             sele3, nanglebins);
    else
      analyser = std::make_unique<GofAngle2>(info, dumpFileName, sele1, sele2,
                                             nanglebins);
  } else if (args_info.r_theta_omega_given) {
    if (args_info.sele3_given)
      analyser = std::make_unique<GofRAngle2>(
          info, dumpFileName, sele1, sele2, sele3, maxLen, nrbins, nanglebins);
    else
      analyser = std::make_unique<GofRAngle2>(info, dumpFileName, sele1, sele2,
                                              maxLen, nrbins, nanglebins);
  } else if (args_info.gxyz_given) {
    if (args_info.refsele_given) {
      analyser = std::make_unique<GofXyz>(info, dumpFileName, sele1, sele2,
                                          args_info.refsele_arg, maxLen,
                                          args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--refsele must set when --gxyz is used");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.twodgofr_given) {
    if (args_info.dz_given) {
      analyser = std::make_unique<TwoDGofR>(info, dumpFileName, sele1, sele2,
                                            maxLen, args_info.dz_arg, nrbins);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A slab width (dz) must be specified when calculating TwoDGofR");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.kirkwood_buff_given) {
    if (args_info.sele1_given && args_info.sele2_given) {
      analyser = std::make_unique<KirkwoodBuff>(info, dumpFileName, sele1,
                                                sele2, maxLen, nrbins);
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "Two selection scripts (--sele1 and --sele2) must be specified when "
          "calculating Kirkwood Buff integrals");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.p2_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given) {
        analyser = std::make_unique<P2OrderParameter>(info, dumpFileName, sele1,
                                                      sele2);
      } else if (args_info.seleoffset_given) {
        analyser = std::make_unique<P2OrderParameter>(info, dumpFileName, sele1,
                                                      args_info.seleoffset_arg);
      } else {
        analyser =
            std::make_unique<P2OrderParameter>(info, dumpFileName, sele1);
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "At least one selection script (--sele1) must be specified when "
               "calculating P2 order parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.rp2_given) {
    analyser = std::make_unique<RippleOP>(info, dumpFileName, sele1, sele2);
  } else if (args_info.bo_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<BondOrderParameter>(
          info, dumpFileName, sele1, args_info.rcut_arg, args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating Bond "
               "Order "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.multipole_given) {
    analyser = std::make_unique<MultipoleSum>(info, dumpFileName, sele1, maxLen,
                                              args_info.nbins_arg);

  } else if (args_info.tet_param_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<TetrahedralityParam>(
          info, dumpFileName, sele1, args_info.rcut_arg, args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating "
               "Tetrahedrality "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.tet_param_z_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<TetrahedralityParamZ>(
          info, dumpFileName, sele1, sele2, args_info.rcut_arg,
          args_info.nbins_arg, privilegedAxis);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating "
               "Tetrahedrality "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.tet_param_r_given) {
    if (args_info.rcut_given) {
      if (args_info.sele3_given) {
        analyser = std::make_unique<TetrahedralityParamR>(
            info, dumpFileName, sele1, sele2, sele3, args_info.rcut_arg, maxLen,
            nrbins);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Selection3 (--sele3) must be given when calculating "
                 "Tetrahedrality Parameter Qk(r)");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating "
               "Tetrahedrality "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.tet_param_dens_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<TetrahedralityParamDens>(
          info, dumpFileName, sele1, sele2, args_info.rcut_arg,
          args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating "
               "Tetrahedrality "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.tet_hb_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<TetrahedralityHBMatrix>(
          info, dumpFileName, sele1, args_info.rcut_arg, args_info.OOcut_arg,
          args_info.thetacut_arg, args_info.OHcut_arg, args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating "
               " Tetrahedrality Hydrogen Bonding Matrix");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.tet_param_xyz_given) {
    if (!args_info.rcut_given) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating"
               " Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    if (!args_info.voxelSize_given) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A voxel size must be specified when calculating"
               " volume-resolved Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    if (!args_info.gaussWidth_given) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A gaussian width must be specified when calculating"
               " volume-resolved Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
    analyser = std::make_unique<TetrahedralityParamXYZ>(
        info, dumpFileName, sele1, sele2, args_info.rcut_arg,
        args_info.voxelSize_arg, args_info.gaussWidth_arg);
  } else if (args_info.ior_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<IcosahedralOfR>(
          info, dumpFileName, sele1, args_info.rcut_arg, nrbins, maxLen);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating Bond "
               "Order "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.for_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<FCCOfR>(info, dumpFileName, sele1,
                                          args_info.rcut_arg, nrbins, maxLen);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating Bond "
               "Order "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.bad_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<BondAngleDistribution>(
          info, dumpFileName, sele1, args_info.rcut_arg, args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating Bond "
               "Angle "
               "Distributions");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.scd_given) {
    if (batchMode) {
      analyser = std::make_unique<SCDOrderParameter>(
          info, dumpFileName, args_info.molname_arg, args_info.begin_arg,
          args_info.end_arg);
    } else {
      analyser = std::make_unique<SCDOrderParameter>(info, dumpFileName, sele1,
                                                     sele2, sele3);
    }
  } else if (args_info.density_given) {
    analyser = std::make_unique<DensityPlot>(info, dumpFileName, sele1, sele2,
                                             maxLen, args_info.nbins_arg);
  } else if (args_info.count_given) {
    analyser = std::make_unique<ObjectCount>(info, dumpFileName, sele1);
  } else if (args_info.mcount_given) {
    analyser = std::make_unique<MoleculeCount>(info, dumpFileName, sele1);
  } else if (args_info.slab_density_given) {
    analyser = std::make_unique<RhoZ>(info, dumpFileName, sele1,
                                      args_info.nbins_arg, privilegedAxis);
  } else if (args_info.eam_density_given) {
    analyser = std::make_unique<DensityHistogram>(info, dumpFileName, sele1,
                                                  args_info.nbins_arg);
  } else if (args_info.momentum_distribution_given) {
    analyser = std::make_unique<MomentumHistogram>(
        info, dumpFileName, sele1, args_info.nbins_arg, momentum_type,
        momentum_comp);
  } else if (args_info.net_charge_given) {
    analyser = std::make_unique<ChargeHistogram>(info, dumpFileName, sele1,
                                                 args_info.nbins_arg);
  } else if (args_info.current_density_given) {
    analyser = std::make_unique<CurrentDensity>(
        info, dumpFileName, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.chargez_given) {
    analyser = std::make_unique<ChargeZ>(info, dumpFileName, sele1,
                                         args_info.nbins_arg, privilegedAxis);
  } else if (args_info.charger_given) {
    analyser = std::make_unique<ChargeR>(info, dumpFileName, sele1, maxLen,
                                         args_info.nbins_arg);
  } else if (args_info.numberz_given) {
    analyser = std::make_unique<NumberZ>(info, dumpFileName, sele1,
                                         args_info.nbins_arg, privilegedAxis);
  } else if (args_info.numberr_given) {
    analyser = std::make_unique<NumberR>(info, dumpFileName, sele1, maxLen,
                                         args_info.nbins_arg);
  } else if (args_info.massdensityz_given) {
    analyser = std::make_unique<MassDensityZ>(
        info, dumpFileName, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.massdensityr_given) {
    analyser = std::make_unique<MassDensityR>(info, dumpFileName, sele1, maxLen,
                                              args_info.nbins_arg);
  } else if (args_info.charge_density_z_given) {
    analyser = std::make_unique<ChargeDensityZ>(
        info, dumpFileName, sele1, args_info.nbins_arg, vRadius,
        args_info.atom_name_arg, args_info.gen_xyz_flag, privilegedAxis);
  } else if (args_info.countz_given) {
    analyser = std::make_unique<PositionZ>(info, dumpFileName, sele1,
                                           args_info.nbins_arg, privilegedAxis);
  } else if (args_info.pipe_density_given) {
    switch (privilegedAxis) {
    case 0:
      analyser = std::make_unique<PipeDensity>(
          info, dumpFileName, sele1, args_info.nbins_y_arg,
          args_info.nbins_z_arg, privilegedAxis);
      break;
    case 1:
      analyser = std::make_unique<PipeDensity>(
          info, dumpFileName, sele1, args_info.nbins_z_arg,
          args_info.nbins_x_arg, privilegedAxis);
      break;
    case 2:
    default:
      analyser = std::make_unique<PipeDensity>(
          info, dumpFileName, sele1, args_info.nbins_x_arg,
          args_info.nbins_y_arg, privilegedAxis);
      break;
    }
  } else if (args_info.rnemdz_given) {
    analyser = std::make_unique<RNEMDZ>(info, dumpFileName, sele1,
                                        args_info.nbins_arg, privilegedAxis);
  } else if (args_info.rnemdr_given) {
    analyser = std::make_unique<RNEMDR>(info, dumpFileName, sele1, comsele,
                                        nrbins, binWidth);
  } else if (args_info.rnemdrt_given) {
    analyser = std::make_unique<RNEMDRTheta>(info, dumpFileName, sele1, comsele,
                                             nrbins, binWidth, nanglebins);
  } else if (args_info.nitrile_given) {
    analyser = std::make_unique<NitrileFrequencyMap>(info, dumpFileName, sele1,
                                                     args_info.nbins_arg);
  } else if (args_info.p_angle_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given)
        analyser = std::make_unique<pAngle>(info, dumpFileName, sele1, sele2,
                                            args_info.nbins_arg);
      else if (args_info.seleoffset_given) {
        if (args_info.seleoffset2_given) {
          analyser = std::make_unique<pAngle>(
              info, dumpFileName, sele1, args_info.seleoffset_arg,
              args_info.seleoffset2_arg, args_info.nbins_arg);
        } else {
          analyser = std::make_unique<pAngle>(info, dumpFileName, sele1,
                                              args_info.seleoffset_arg,
                                              args_info.nbins_arg);
        }
      } else
        analyser = std::make_unique<pAngle>(info, dumpFileName, sele1,
                                            args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "At least one selection script (--sele1) must be specified when "
               "calculating P(angle) distributions");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
  } else if (args_info.hxy_given) {
    analyser = std::make_unique<Hxy>(
        info, dumpFileName, sele1, args_info.nbins_x_arg, args_info.nbins_y_arg,
        args_info.nbins_z_arg, args_info.nbins_arg);
#endif
  } else if (args_info.cn_given || args_info.scn_given || args_info.gcn_given) {
    if (args_info.rcut_given) {
      if (args_info.cn_given) {
        analyser = std::make_unique<CoordinationNumber>(
            info, dumpFileName, sele1, sele2, args_info.rcut_arg,
            args_info.nbins_arg);
      } else if (args_info.scn_given) {
        analyser =
            std::make_unique<SCN>(info, dumpFileName, sele1, sele2,
                                  args_info.rcut_arg, args_info.nbins_arg);
      } else if (args_info.gcn_given) {
        analyser =
            std::make_unique<GCN>(info, dumpFileName, sele1, sele2,
                                  args_info.rcut_arg, args_info.nbins_arg);
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating\n"
               "\t Coordination Numbers");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.surfDiffusion_given) {
    analyser =
        std::make_unique<SurfaceDiffusion>(info, dumpFileName, sele1, maxLen);
  } else if (args_info.rho_r_given) {
    if (args_info.radius_given) {
      analyser = std::make_unique<RhoR>(info, dumpFileName, sele1, maxLen,
                                        nrbins, args_info.radius_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A particle radius (radius) must be specified when calculating "
               "Rho(r)");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.hullvol_given) {
    analyser = std::make_unique<NanoVolume>(info, dumpFileName, sele1);
  } else if (args_info.rodlength_given) {
    analyser = std::make_unique<NanoLength>(info, dumpFileName, sele1);
  } else if (args_info.angle_r_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given)
        analyser = std::make_unique<AngleR>(info, dumpFileName, sele1, sele2,
                                            maxLen, nrbins);
      else if (args_info.seleoffset_given) {
        if (args_info.seleoffset2_given) {
          analyser = std::make_unique<AngleR>(
              info, dumpFileName, sele1, args_info.seleoffset_arg,
              args_info.seleoffset2_arg, maxLen, nrbins);
        } else {
          analyser = std::make_unique<AngleR>(info, dumpFileName, sele1,
                                              args_info.seleoffset_arg, maxLen,
                                              nrbins);
        }
      } else
        analyser =
            std::make_unique<AngleR>(info, dumpFileName, sele1, maxLen, nrbins);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "At least one selection script (--sele1) must be specified when "
               "calculating Angle(r) values");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.p2r_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given)
        analyser = std::make_unique<P2R>(info, dumpFileName, sele1, sele2,
                                         args_info.nbins_arg);
      else if (args_info.seleoffset_given) {
        if (args_info.seleoffset2_given) {
          analyser = std::make_unique<P2R>(
              info, dumpFileName, sele1, args_info.seleoffset_arg,
              args_info.seleoffset2_arg, args_info.nbins_arg);
        } else {
          analyser = std::make_unique<P2R>(info, dumpFileName, sele1,
                                           args_info.seleoffset_arg,
                                           args_info.nbins_arg);
        }
      } else
        analyser = std::make_unique<P2R>(info, dumpFileName, sele1,
                                         args_info.nbins_arg);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "At least one selection script (--sele1) must be specified when "
               "calculating P2R values");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.p2z_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given)
        analyser = std::make_unique<P2Z>(info, dumpFileName, sele1, sele2,
                                         args_info.nbins_arg, privilegedAxis);
      else if (args_info.seleoffset_given) {
        if (args_info.seleoffset2_given) {
          analyser = std::make_unique<P2Z>(
              info, dumpFileName, sele1, args_info.seleoffset_arg,
              args_info.seleoffset2_arg, args_info.nbins_arg, privilegedAxis);
        } else {
          analyser = std::make_unique<P2Z>(info, dumpFileName, sele1,
                                           args_info.seleoffset_arg,
                                           args_info.nbins_arg, privilegedAxis);
        }
      } else
        analyser = std::make_unique<P2Z>(info, dumpFileName, sele1,
                                         args_info.nbins_arg, privilegedAxis);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "At least one selection script (--sele1) must be specified when "
               "calculating P2Z values");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.hbond_given) {
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {
        analyser = std::make_unique<HBondGeometric>(
            info, dumpFileName, sele1, sele2, args_info.rcut_arg,
            args_info.thetacut_arg, args_info.nbins_arg);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "A cutoff angle (thetacut) must be specified when calculating "
                 "Hydrogen "
                 "Bonding Statistics");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "A cutoff radius (rcut) must be specified when calculating Hydrogen "
          "Bonding Statistics");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

  } else if (args_info.hbondz_given) {
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {
        analyser = std::make_unique<HBondZ>(
            info, dumpFileName, sele1, sele2, args_info.rcut_arg,
            args_info.thetacut_arg, args_info.nbins_arg);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "A cutoff angle (thetacut) must be specified when calculating "
                 "Hydrogen "
                 "Bonding Statistics");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "A cutoff radius (rcut) must be specified when calculating Hydrogen "
          "Bonding Statistics");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.hbondzvol_given) {
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {
        analyser = std::make_unique<HBondZvol>(
            info, dumpFileName, sele1, sele2, args_info.rcut_arg,
            args_info.thetacut_arg, args_info.nbins_arg);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "A cutoff angle (thetacut) must be specified when calculating "
                 "Hydrogen "
                 "Bonding Statistics");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "A cutoff radius (rcut) must be specified when calculating Hydrogen "
          "Bonding Statistics");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.hbondr_given) {
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {
        analyser = std::make_unique<HBondR>(
            info, dumpFileName, sele1, sele2, sele3, args_info.rcut_arg, maxLen,
            args_info.thetacut_arg, args_info.nrbins_arg);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "A cutoff angle (thetacut) must be specified when calculating "
                 "Hydrogen "
                 "Bonding Statistics");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "A cutoff radius (rcut) must be specified when calculating Hydrogen "
          "Bonding Statistics");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.hbondrvol_given) {
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {
        analyser = std::make_unique<HBondRvol>(
            info, dumpFileName, sele1, sele2, sele3, args_info.rcut_arg, maxLen,
            args_info.thetacut_arg, args_info.nrbins_arg);
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "A cutoff angle (thetacut) must be specified when calculating "
                 "Hydrogen "
                 "Bonding Statistics");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "A cutoff radius (rcut) must be specified when calculating Hydrogen "
          "Bonding Statistics");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.potDiff_given) {
    analyser = std::make_unique<PotDiff>(info, dumpFileName, sele1);
  } else if (args_info.kirkwood_given) {
    analyser = std::make_unique<Kirkwood>(info, dumpFileName, sele1, sele2,
                                          maxLen, nrbins);
  } else if (args_info.kirkwoodQ_given) {
    analyser = std::make_unique<KirkwoodQuadrupoles>(info, dumpFileName, sele1,
                                                     sele2, maxLen, nrbins);
  } else if (args_info.densityfield_given) {
    analyser = std::make_unique<DensityField>(info, dumpFileName, sele1,
                                              args_info.voxelSize_arg);
  } else if (args_info.velocityfield_given) {
    analyser = std::make_unique<VelocityField>(info, dumpFileName, sele1,
                                               args_info.voxelSize_arg);
  } else if (args_info.velocityZ_given) {
    switch (privilegedAxis) {
    case 0:
      if (privilegedAxis2 == 1) {
        analyser = std::make_unique<VelocityZ>(
            info, dumpFileName, sele1, args_info.nbins_x_arg,
            args_info.nbins_y_arg, privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 2) {
        analyser = std::make_unique<VelocityZ>(
            info, dumpFileName, sele1, args_info.nbins_x_arg,
            args_info.nbins_z_arg, privilegedAxis, privilegedAxis2);
      }
      break;
    case 1:
      if (privilegedAxis2 == 0) {
        analyser = std::make_unique<VelocityZ>(
            info, dumpFileName, sele1, args_info.nbins_y_arg,
            args_info.nbins_x_arg, privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 2) {
        analyser = std::make_unique<VelocityZ>(
            info, dumpFileName, sele1, args_info.nbins_y_arg,
            args_info.nbins_z_arg, privilegedAxis, privilegedAxis2);
      }
      break;
    case 2:
    default:
      if (privilegedAxis2 == 0) {
        analyser = std::make_unique<VelocityZ>(
            info, dumpFileName, sele1, args_info.nbins_z_arg,
            args_info.nbins_x_arg, privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 1) {
        analyser = std::make_unique<VelocityZ>(
            info, dumpFileName, sele1, args_info.nbins_z_arg,
            args_info.nbins_y_arg, privilegedAxis, privilegedAxis2);
      }
      break;
    }
  } else if (args_info.dipole_orientation_given) {
    if (args_info.dipoleX_given && args_info.dipoleY_given &&
        args_info.dipoleZ_given)
      analyser = std::make_unique<DipoleOrientation>(
          info, dumpFileName, sele1, args_info.dipoleX_arg,
          args_info.dipoleY_arg, args_info.dipoleZ_arg, args_info.nbins_arg,
          privilegedAxis);
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Dipole components must be provided.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

  } else if (args_info.order_prob_given) {
    if (args_info.dipoleX_given && args_info.dipoleY_given &&
        args_info.dipoleZ_given)
      analyser = std::make_unique<OrderParameterProbZ>(
          info, dumpFileName, sele1, args_info.dipoleX_arg,
          args_info.dipoleY_arg, args_info.dipoleZ_arg, args_info.nbins_arg,
          privilegedAxis);
    else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Dipole components must be provided.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  } else if (args_info.trans_param_z_given) {
    if (args_info.rcut_given) {
      analyser = std::make_unique<TranslationalOrderParamZ>(
          info, dumpFileName, sele1, sele2, args_info.rcut_arg,
          args_info.nbins_arg, args_info.nbins_z_arg, maxLen, zmaxLen,
          privilegedAxis);
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "A cutoff radius (rcut) must be specified when calculating "
               "Translational Order "
               "Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  }

  if (analyser != NULL) {
    if (args_info.output_given) {
      analyser->setOutputName(args_info.output_arg);
    }
    if (args_info.step_given) { analyser->setStep(args_info.step_arg); }

    analyser->process();
  } else {
    snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
             "StaticProps: No Analyser was created, nothing to do!");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
  }

  delete info;

  return 0;
}
