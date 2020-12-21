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
#include "io/DumpReader.hpp"
#include "utils/simError.h"
#include "utils/MemoryUtils.hpp"

#include "StaticPropsCmd.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "applications/staticProps/GofR.hpp"
#include "applications/staticProps/GofZ.hpp"
#include "applications/staticProps/GofRZ.hpp"
#include "applications/staticProps/GofRAngle.hpp"
#include "applications/staticProps/GofAngle2.hpp"
#include "applications/staticProps/GofRAngle2.hpp"
#include "applications/staticProps/GofXyz.hpp"
#include "applications/staticProps/TwoDGofR.hpp"
#include "applications/staticProps/P2OrderParameter.hpp"
#include "applications/staticProps/BondOrderParameter.hpp"
#include "applications/staticProps/BOPofR.hpp"
#include "applications/staticProps/RippleOP.hpp"
#include "applications/staticProps/SCDOrderParameter.hpp"
#include "applications/staticProps/DensityPlot.hpp"
#include "applications/staticProps/ObjectCount.hpp"
#include "applications/staticProps/RhoZ.hpp"
#include "applications/staticProps/PipeDensity.hpp"
#include "applications/staticProps/pAngle.hpp"
#include "applications/staticProps/BondAngleDistribution.hpp"
#include "applications/staticProps/NanoVolume.hpp"
#include "applications/staticProps/NanoLength.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "applications/staticProps/Hxy.hpp"
#endif
#include "applications/staticProps/RhoR.hpp"
#include "applications/staticProps/AngleR.hpp"
#include "applications/staticProps/TetrahedralityParam.hpp"
#include "applications/staticProps/TetrahedralityParamZ.hpp"
#include "applications/staticProps/TetrahedralityParamXYZ.hpp"
#include "applications/staticProps/TetrahedralityParamDens.hpp"
#include "applications/staticProps/RNEMDStats.hpp"
#include "applications/staticProps/NitrileFrequencyMap.hpp"
#include "applications/staticProps/MultipoleSum.hpp"
#include "applications/staticProps/SurfaceDiffusion.hpp"
#include "applications/staticProps/CoordinationNumber.hpp"
#include "applications/staticProps/HBondGeometric.hpp"
#include "applications/staticProps/PotDiff.hpp"
#include "applications/staticProps/TetrahedralityHBMatrix.hpp"
#include "applications/staticProps/Kirkwood.hpp"
#include "applications/staticProps/Field.hpp"
#include "applications/staticProps/VelocityZ.hpp"
#include "applications/staticProps/DensityHistogram.hpp"
#include "applications/staticProps/MomentumHistogram.hpp"
#include "applications/staticProps/ChargeHistogram.hpp"
#include "applications/staticProps/CurrentDensity.hpp"
#include "applications/staticProps/ChargeZ.hpp"
#include "applications/staticProps/PositionZ.hpp"
#include "applications/staticProps/DipoleOrientation.hpp"
#include "applications/staticProps/ChargeDensityZ.hpp"
#include "applications/staticProps/OrderParameterProbZ.hpp"

using namespace OpenMD;

int main(int argc, char* argv[]){


  gengetopt_args_info args_info;

  //parse the command line option
  if (cmdline_parser (argc, argv, &args_info) != 0) {
    exit(1) ;
  }

  //get the dumpfile name
  std::string dumpFileName = args_info.input_arg;
  std::string sele1;
  std::string sele2;
  std::string sele3;

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

  // check the third selection argument, which is only set if
  // requested by the user

  if (args_info.sele3_given) sele3 = args_info.sele3_arg;

  bool batchMode(false);
  if (args_info.scd_given){
    if (args_info.sele1_given &&
        args_info.sele2_given && args_info.sele3_given) {
      batchMode = false;
    } else if (args_info.molname_given &&
               args_info.begin_given && args_info.end_given) {
      if (args_info.begin_arg < 0 ||
          args_info.end_arg < 0 || args_info.begin_arg > args_info.end_arg-2) {
        sprintf( painCave.errMsg,
                 "below conditions are not satisfied:\n"
                 "0 <= begin && 0<= end && begin <= end-2\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }
      batchMode = true;
    } else{
      sprintf( painCave.errMsg,
               "either --sele1, --sele2, --sele3 are specified,"
               " or --molname, --begin, --end are specified\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }

  //parse md file and set up the system
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
    if (args_info.zlength_given){
      zmaxLen = args_info.zlength_arg;
    }
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

//override default vander waals radius for fictious atoms in a model
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

  if (args_info.gofr_given){
    analyser= MemoryUtils::make_unique<GofR>(info, dumpFileName, sele1, sele2, maxLen,
		       nrbins);
  } else if (args_info.gofz_given) {
    analyser= MemoryUtils::make_unique<GofZ>(info, dumpFileName, sele1, sele2, maxLen, zmaxLen,
		       args_info.nbins_arg, privilegedAxis);
  } else if (args_info.r_z_given) {
    analyser  = MemoryUtils::make_unique<GofRZ>(info, dumpFileName, sele1, sele2, maxLen, zmaxLen,
			  nrbins, args_info.nbins_z_arg, privilegedAxis);
  } else if (args_info.r_theta_given) {
    if (args_info.sele3_given)
      analyser  = MemoryUtils::make_unique<GofRTheta>(info, dumpFileName, sele1, sele2, sele3, maxLen,
                                nrbins, nanglebins);
    else
      analyser  = MemoryUtils::make_unique<GofRTheta>(info, dumpFileName, sele1, sele2, maxLen,
                                nrbins, nanglebins);
  } else if (args_info.r_omega_given) {
    if (args_info.sele3_given)
      analyser  = MemoryUtils::make_unique<GofROmega>(info, dumpFileName, sele1, sele2, sele3, maxLen,
                                nrbins, nanglebins);
    else
      analyser  = MemoryUtils::make_unique<GofROmega>(info, dumpFileName, sele1, sele2, maxLen,
                                nrbins, nanglebins);

  } else if (args_info.theta_omega_given) {
    if (args_info.sele3_given)
      analyser  = MemoryUtils::make_unique<GofAngle2>(info, dumpFileName, sele1, sele2, sele3,
                                nanglebins);
    else
      analyser  = MemoryUtils::make_unique<GofAngle2>(info, dumpFileName, sele1, sele2,
                                nanglebins);
  } else if (args_info.r_theta_omega_given) {
    if (args_info.sele3_given)
      analyser  = MemoryUtils::make_unique<GofRAngle2>(info, dumpFileName, sele1, sele2, sele3,
                                 maxLen, nrbins, nanglebins);
    else
      analyser  = MemoryUtils::make_unique<GofRAngle2>(info, dumpFileName, sele1, sele2,
                                 maxLen, nrbins, nanglebins);
  } else if (args_info.gxyz_given) {
    if (args_info.refsele_given) {
      analyser= MemoryUtils::make_unique<GofXyz>(info, dumpFileName, sele1, sele2,
                           args_info.refsele_arg, maxLen, args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
	       "--refsele must set when --gxyz is used");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.twodgofr_given){
    if (args_info.dz_given) {
      analyser= MemoryUtils::make_unique<TwoDGofR>(info, dumpFileName, sele1, sele2, maxLen,
			     args_info.dz_arg, nrbins);
    } else {
      sprintf( painCave.errMsg,
	       "A slab width (dz) must be specified when calculating TwoDGofR");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.p2_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given)
        analyser  = MemoryUtils::make_unique<P2OrderParameter>(info, dumpFileName, sele1, sele2);
      else
        if (args_info.seleoffset_given)
          analyser  = MemoryUtils::make_unique<P2OrderParameter>(info, dumpFileName, sele1,
                                           args_info.seleoffset_arg);
        else
          analyser  = MemoryUtils::make_unique<P2OrderParameter>(info, dumpFileName, sele1);
    } else {
      sprintf( painCave.errMsg,
	       "At least one selection script (--sele1) must be specified when calculating P2 order parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.rp2_given){
    analyser = MemoryUtils::make_unique<RippleOP>(info, dumpFileName, sele1, sele2);
  } else if (args_info.bo_given){
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<BondOrderParameter>(info, dumpFileName, sele1,
					args_info.rcut_arg,
					args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating Bond Order Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.multipole_given){
    analyser = MemoryUtils::make_unique<MultipoleSum>(info, dumpFileName, sele1,
                                maxLen, args_info.nbins_arg);

  } else if (args_info.tet_param_given) {
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<TetrahedralityParam>(info, dumpFileName, sele1,
					 args_info.rcut_arg,
					 args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

  } else if (args_info.tet_param_z_given) {
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<TetrahedralityParamZ>(info, dumpFileName, sele1, sele2,
                                          args_info.rcut_arg,
                                          args_info.nbins_arg,
					  privilegedAxis);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

  } else if (args_info.tet_param_dens_given) {
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<TetrahedralityParamDens>(info, dumpFileName, sele1, sele2,
					     args_info.rcut_arg,
					     args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
               "A cutoff radius (rcut) must be specified when calculating Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.tet_hb_given) {
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<TetrahedralityHBMatrix>(info, dumpFileName, sele1,
                                            args_info.rcut_arg,
                                            args_info.OOcut_arg,
                                            args_info.thetacut_arg,
                                            args_info.OHcut_arg,
                                            args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating "
               " Tetrahedrality Hydrogen Bonding Matrix");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.tet_param_xyz_given) {
    if (!args_info.rcut_given) {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating"
               " Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    if (!args_info.voxelSize_given) {
      sprintf( painCave.errMsg,
	       "A voxel size must be specified when calculating"
               " volume-resolved Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    if (!args_info.gaussWidth_given) {
      sprintf( painCave.errMsg,
	       "A gaussian width must be specified when calculating"
               " volume-resolved Tetrahedrality Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    analyser = MemoryUtils::make_unique<TetrahedralityParamXYZ>(info, dumpFileName, sele1, sele2,
                                          args_info.rcut_arg,
                                          args_info.voxelSize_arg,
                                          args_info.gaussWidth_arg);
  } else if (args_info.ior_given){
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<IcosahedralOfR>(info, dumpFileName, sele1,
                                    args_info.rcut_arg,
                                    nrbins, maxLen);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating Bond Order Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.for_given){
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<FCCOfR>(info, dumpFileName, sele1, args_info.rcut_arg,
			    nrbins, maxLen);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating Bond Order Parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.bad_given){
    if (args_info.rcut_given) {
      analyser = MemoryUtils::make_unique<BondAngleDistribution>(info, dumpFileName, sele1,
                                           args_info.rcut_arg,
					   args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
	       "A cutoff radius (rcut) must be specified when calculating Bond Angle Distributions");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.scd_given) {
    if (batchMode) {
      analyser  = MemoryUtils::make_unique<SCDOrderParameter>(info, dumpFileName,
                                        args_info.molname_arg,
					args_info.begin_arg, args_info.end_arg);
    } else{
      analyser  = MemoryUtils::make_unique<SCDOrderParameter>(info, dumpFileName,
                                        sele1, sele2, sele3);
    }
  } else if (args_info.density_given) {
    analyser= MemoryUtils::make_unique<DensityPlot>(info, dumpFileName, sele1, sele2, maxLen,
			      args_info.nbins_arg);
  } else if (args_info.count_given) {
    analyser = MemoryUtils::make_unique<ObjectCount>(info, dumpFileName, sele1 );
  } else if (args_info.slab_density_given) {
    analyser = MemoryUtils::make_unique<RhoZ>(info, dumpFileName, sele1, args_info.nbins_arg,
                        privilegedAxis);
  } else if (args_info.eam_density_given) {
    analyser = MemoryUtils::make_unique<DensityHistogram>(info, dumpFileName, sele1,
                                    args_info.nbins_arg);
  } else if (args_info.momentum_distribution_given) {

    analyser = MemoryUtils::make_unique<MomentumHistogram>(info, dumpFileName, sele1, args_info.nbins_arg, momentum_type, momentum_comp);
  } else if (args_info.net_charge_given) {
      analyser = MemoryUtils::make_unique<ChargeHistogram>(info, dumpFileName, sele1, args_info.nbins_arg);
  } else if (args_info.current_density_given) {
    analyser = MemoryUtils::make_unique<CurrentDensity>(info, dumpFileName, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.chargez_given) {
    analyser = MemoryUtils::make_unique<ChargeZ>(info, dumpFileName, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.charge_density_z_given) {
    analyser = MemoryUtils::make_unique<ChargeDensityZ>(info, dumpFileName, sele1, args_info.nbins_arg, vRadius,args_info.atom_name_arg, args_info.gen_xyz_flag, privilegedAxis);
  } else if (args_info.countz_given) {
    analyser = MemoryUtils::make_unique<PositionZ>(info, dumpFileName, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.pipe_density_given) {

    switch (privilegedAxis) {
    case 0:
      analyser = MemoryUtils::make_unique<PipeDensity>(info, dumpFileName, sele1,
                                 args_info.nbins_y_arg, args_info.nbins_z_arg,
                                 privilegedAxis);
      break;
    case 1:
      analyser = MemoryUtils::make_unique<PipeDensity>(info, dumpFileName, sele1,
                                 args_info.nbins_z_arg, args_info.nbins_x_arg,
                                 privilegedAxis);
      break;
    case 2:
    default:
      analyser = MemoryUtils::make_unique<PipeDensity>(info, dumpFileName, sele1,
                                 args_info.nbins_x_arg, args_info.nbins_y_arg,
                                 privilegedAxis);
      break;
    }
  } else if (args_info.rnemdz_given) {
    analyser = MemoryUtils::make_unique<RNEMDZ>(info, dumpFileName, sele1, args_info.nbins_arg,
                          privilegedAxis);
  } else if (args_info.rnemdr_given) {
    analyser = MemoryUtils::make_unique<RNEMDR>(info, dumpFileName, sele1, nrbins);
  } else if (args_info.rnemdrt_given) {
    analyser = MemoryUtils::make_unique<RNEMDRTheta>(info, dumpFileName, sele1, nrbins, nanglebins);
  } else if (args_info.nitrile_given) {
    analyser = MemoryUtils::make_unique<NitrileFrequencyMap>(info, dumpFileName, sele1,
                                       args_info.nbins_arg);
  } else if (args_info.p_angle_given) {
    if (args_info.sele1_given) {
      if (args_info.sele2_given)
        analyser  = MemoryUtils::make_unique<pAngle>(info, dumpFileName, sele1, sele2,
                               args_info.nbins_arg);
      else
        if (args_info.seleoffset_given) {
          if (args_info.seleoffset2_given) {
            analyser  = MemoryUtils::make_unique<pAngle>(info, dumpFileName, sele1,
                                   args_info.seleoffset_arg,
                                   args_info.seleoffset2_arg,
                                   args_info.nbins_arg);
          } else {
            analyser  = MemoryUtils::make_unique<pAngle>(info, dumpFileName, sele1,
                                   args_info.seleoffset_arg,
                                   args_info.nbins_arg);
          }
        } else
          analyser  = MemoryUtils::make_unique<pAngle>(info, dumpFileName, sele1,
                                 args_info.nbins_arg);
    } else {
      sprintf( painCave.errMsg,
	       "At least one selection script (--sele1) must be specified when "
               "calculating P(angle) distributions");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
  }else if (args_info.hxy_given) {
    analyser = MemoryUtils::make_unique<Hxy>(info, dumpFileName, sele1, args_info.nbins_x_arg,
		       args_info.nbins_y_arg, args_info.nbins_z_arg,
                       args_info.nbins_arg);
#endif
  }else if(args_info.cn_given || args_info.scn_given || args_info.gcn_given){
    if (args_info.rcut_given) {
      if (args_info.cn_given) {
        analyser = MemoryUtils::make_unique<CoordinationNumber>(info, dumpFileName, sele1, sele2,
                                          args_info.rcut_arg,
                                          args_info.nbins_arg);
      } else if (args_info.scn_given) {
        analyser = MemoryUtils::make_unique<SCN>(info, dumpFileName, sele1, sele2,
                           args_info.rcut_arg, args_info.nbins_arg);
      } else if (args_info.gcn_given) {
        analyser = MemoryUtils::make_unique<GCN>(info, dumpFileName, sele1, sele2,
                           args_info.rcut_arg, args_info.nbins_arg);
      }
    } else {
      sprintf( painCave.errMsg,
               "A cutoff radius (rcut) must be specified when calculating\n"
               "\t Coordination Numbers");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  }
  else if (args_info.surfDiffusion_given){
    analyser = MemoryUtils::make_unique<SurfaceDiffusion>(info, dumpFileName, sele1, maxLen);
  }else if (args_info.rho_r_given) {
    if (args_info.radius_given){
      analyser = MemoryUtils::make_unique<RhoR>(info, dumpFileName, sele1, maxLen, nrbins,
                          args_info.radius_arg);
    }else{
      sprintf( painCave.errMsg,
	       "A particle radius (radius) must be specified when calculating Rho(r)");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.hullvol_given) {
    analyser = MemoryUtils::make_unique<NanoVolume>(info, dumpFileName, sele1);
  } else if (args_info.rodlength_given) {
    analyser = MemoryUtils::make_unique<NanoLength>(info, dumpFileName, sele1);
  } else if (args_info.angle_r_given) {
    analyser = MemoryUtils::make_unique<AngleR>(info, dumpFileName, sele1, maxLen, nrbins);
  } else if (args_info.hbond_given){
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {

        analyser = MemoryUtils::make_unique<HBondGeometric>(info, dumpFileName, sele1, sele2,
                                      args_info.rcut_arg,
                                      args_info.thetacut_arg,
                                      args_info.nbins_arg);
      } else {
        sprintf( painCave.errMsg,
                 "A cutoff angle (thetacut) must be specified when calculating Hydrogen Bonding Statistics");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      }
    } else {
      sprintf( painCave.errMsg,
               "A cutoff radius (rcut) must be specified when calculating Hydrogen Bonding Statistics");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.potDiff_given) {
    analyser = MemoryUtils::make_unique<PotDiff>(info, dumpFileName, sele1);
  } else if (args_info.kirkwood_given) {
    analyser= MemoryUtils::make_unique<Kirkwood>(info, dumpFileName, sele1, sele2, maxLen,
                           nrbins);
  } else if (args_info.kirkwoodQ_given) {
    analyser= MemoryUtils::make_unique<KirkwoodQuadrupoles>(info, dumpFileName, sele1, sele2, maxLen,
                                      nrbins);
  } else if (args_info.densityfield_given) {
    analyser = MemoryUtils::make_unique<DensityField>(info, dumpFileName, sele1,
                                args_info.voxelSize_arg);
  } else if (args_info.velocityfield_given) {
    analyser = MemoryUtils::make_unique<VelocityField>(info, dumpFileName, sele1,
                                 args_info.voxelSize_arg);
  } else if (args_info.velocityZ_given) {

    switch (privilegedAxis) {
    case 0:
      if (privilegedAxis2 == 1) {
	analyser = MemoryUtils::make_unique<VelocityZ>(info, dumpFileName, sele1,
                                 args_info.nbins_x_arg, args_info.nbins_y_arg,
                                 privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 2) {
	analyser = MemoryUtils::make_unique<VelocityZ>(info, dumpFileName, sele1,
                                 args_info.nbins_x_arg, args_info.nbins_z_arg,
                                 privilegedAxis, privilegedAxis2);
      }
      break;
    case 1:
      if (privilegedAxis2 == 0) {
	analyser = MemoryUtils::make_unique<VelocityZ>(info, dumpFileName, sele1,
                                 args_info.nbins_y_arg, args_info.nbins_x_arg,
                                 privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 2) {
	analyser = MemoryUtils::make_unique<VelocityZ>(info, dumpFileName, sele1,
                                 args_info.nbins_y_arg, args_info.nbins_z_arg,
                                 privilegedAxis, privilegedAxis2);
      }
      break;
    case 2:
    default:
      if (privilegedAxis2 == 0) {
	analyser = MemoryUtils::make_unique<VelocityZ>(info, dumpFileName, sele1,
                                 args_info.nbins_z_arg, args_info.nbins_x_arg,
                                 privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 1) {
	analyser = MemoryUtils::make_unique<VelocityZ>(info, dumpFileName, sele1,
                                 args_info.nbins_z_arg, args_info.nbins_y_arg,
                                 privilegedAxis, privilegedAxis2);
      }
      break;
    }
  } else if (args_info.dipole_orientation_given){
    if(args_info.dipoleX_given && args_info.dipoleY_given
       && args_info.dipoleZ_given)
      analyser = MemoryUtils::make_unique<DipoleOrientation>(info, dumpFileName, sele1,
                                       args_info.dipoleX_arg,
                                       args_info.dipoleY_arg,
                                       args_info.dipoleZ_arg,
                                       args_info.nbins_arg, privilegedAxis);
    else{
      sprintf( painCave.errMsg,
      "Dipole components must be provided.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

  } else if (args_info.order_prob_given){
    if(args_info.dipoleX_given && args_info.dipoleY_given
       && args_info.dipoleZ_given)
      analyser = MemoryUtils::make_unique<OrderParameterProbZ>(info, dumpFileName, sele1,
                                       args_info.dipoleX_arg,
                                       args_info.dipoleY_arg,
                                       args_info.dipoleZ_arg,
                                       args_info.nbins_arg, privilegedAxis);
    else{
      sprintf( painCave.errMsg,
      "Dipole components must be provided.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
   }

  }


  if (args_info.output_given) {
    analyser->setOutputName(args_info.output_arg);
  }
  if (args_info.step_given) {
    analyser->setStep(args_info.step_arg);
  }

  analyser->process();

  delete info;

  return 0;
}
