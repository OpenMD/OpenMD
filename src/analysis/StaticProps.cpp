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
#include "io/DumpReader.hpp"
#include "utils/simError.h"

#include "StaticPropsCmd.hpp"
#include "analysis/StaticAnalyser.hpp"
#include "analysis/GofR.hpp"
#include "analysis/GofZ.hpp"
#include "analysis/GofRZ.hpp"
#include "analysis/GofRAngle.hpp"
#include "analysis/GofAngle2.hpp"
#include "analysis/GofRAngle2.hpp"
#include "analysis/GofXyz.hpp"
#include "analysis/TwoDGofR.hpp"
#include "analysis/P2OrderParameter.hpp"
#include "analysis/BondOrderParameter.hpp"
#include "analysis/BOPofR.hpp"
#include "analysis/RippleOP.hpp"
#include "analysis/SCDOrderParameter.hpp"
#include "analysis/DensityPlot.hpp"
#include "analysis/ObjectCount.hpp"
#include "analysis/RhoZ.hpp"
#include "analysis/PipeDensity.hpp"
#include "analysis/pAngle.hpp"
#include "analysis/BondAngleDistribution.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "analysis/Hxy.hpp"
#endif
#include "analysis/RhoR.hpp"
#include "analysis/AngleR.hpp"
#include "analysis/TetrahedralityParam.hpp"
#include "analysis/TetrahedralityParamZ.hpp"
#include "analysis/TetrahedralityParamXYZ.hpp"
#include "analysis/TetrahedralityParamDens.hpp"
#include "analysis/RNEMDStats.hpp"
#include "analysis/NitrileFrequencyMap.hpp"
#include "analysis/MultipoleSum.hpp"
#include "analysis/CoordinationNumber.hpp"
#include "analysis/HBondGeometric.hpp"
#include "analysis/PotDiff.hpp"
#include "analysis/TetrahedralityHBMatrix.hpp"
#include "analysis/Kirkwood.hpp"
#include "analysis/Field.hpp"
#include "analysis/VelocityZ.hpp"

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
  
  /*
    parse md file and set up the system. When we call creator.createSim(),
    the prefix (prefix.omd) and the dump (prefix.dump) file names are 
    established in SimInfo. Therefore, we don't have to pass filenames here.
  */
  SimCreator creator;
  SimInfo* info = creator.createSim(dumpFileName);
  

  RealType maxLen;
  RealType zmaxLen(0.0);
  if (args_info.length_given) {
    maxLen = args_info.length_arg;
    if (args_info.zlength_given){
      zmaxLen = args_info.zlength_arg;
    }
  } else {
    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    maxLen = std::min(std::min(hmat(0, 0), hmat(1, 1)), hmat(2, 2)) /2.0;
    zmaxLen = hmat(2,2);    
  }    

  int nanglebins, nrbins;
  // in case we override nbins with nrbins:
  if (args_info.nrbins_given) {
    nrbins = args_info.nrbins_arg;
  } else {
    nrbins = args_info.nbins_arg;
  }
  // in case we override nbins with nrbins:
  if (args_info.nanglebins_given) {
    nanglebins = args_info.nanglebins_arg;
  } else {
    nanglebins = args_info.nbins_arg;
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
      
  StaticAnalyser* analyser;
  
                                       
  if (args_info.gofr_given){
    analyser= new GofR(info, sele1, sele2, maxLen, 
		       nrbins);        
  } else if (args_info.gofz_given) {
    analyser= new GofZ(info, sele1, sele2, maxLen,
		       args_info.nbins_arg, privilegedAxis);
  } else if (args_info.r_z_given) {
    analyser  = new GofRZ(info, sele1, sele2, maxLen, zmaxLen, 
			  nrbins, args_info.nbins_z_arg, privilegedAxis);
  } else if (args_info.r_theta_given) {
    if (args_info.sele3_given) 
      analyser  = new GofRTheta(info, sele1, sele2, sele3, maxLen,
                                nrbins, nanglebins);
    else 
      analyser  = new GofRTheta(info, sele1, sele2, maxLen, 
                                nrbins, nanglebins);
  } else if (args_info.r_omega_given) {
    if (args_info.sele3_given) 
      analyser  = new GofROmega(info, sele1, sele2, sele3, maxLen,
                                nrbins, nanglebins);
    else 
      analyser  = new GofROmega(info, sele1, sele2, maxLen,
                                nrbins, nanglebins);

  } else if (args_info.theta_omega_given) {
    if (args_info.sele3_given) 
      analyser  = new GofAngle2(info, sele1, sele2, sele3,
                                nanglebins);
    else
      analyser  = new GofAngle2(info, sele1, sele2, 
                                nanglebins);
  } else if (args_info.r_theta_omega_given) {
    if (args_info.sele3_given) 
      analyser  = new GofRAngle2(info, sele1, sele2, sele3,
                                 maxLen, nrbins, nanglebins);
    else
      analyser  = new GofRAngle2(info, sele1, sele2, 
                                 maxLen, nrbins, nanglebins);
  } else if (args_info.gxyz_given) {
    if (args_info.refsele_given) {
      analyser= new GofXyz(info, sele1, sele2,
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
      analyser= new TwoDGofR(info, sele1, sele2, maxLen, 
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
        analyser  = new P2OrderParameter(info, sele1, sele2);
      else 
        if (args_info.seleoffset_given) 
          analyser  = new P2OrderParameter(info, sele1, 
                                           args_info.seleoffset_arg);
        else 
          analyser  = new P2OrderParameter(info, sele1);
    } else {
      sprintf( painCave.errMsg,
	       "At least one selection script (--sele1) must be specified when calculating P2 order parameters");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.rp2_given){
    analyser = new RippleOP(info, sele1, sele2);
  } else if (args_info.bo_given){
    if (args_info.rcut_given) {
      analyser = new BondOrderParameter(info, sele1, 
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
    analyser = new MultipoleSum(info, sele1, 
                                maxLen, args_info.nbins_arg);
    
  } else if (args_info.tet_param_given) {
    if (args_info.rcut_given) {	  
      analyser = new TetrahedralityParam(info, sele1, 
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
      analyser = new TetrahedralityParamZ(info, sele1, sele2,
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
      analyser = new TetrahedralityParamDens(info, sele1, sele2,
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
      analyser = new TetrahedralityHBMatrix(info, sele1,
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
    analyser = new TetrahedralityParamXYZ(info, sele1, sele2,
                                          args_info.rcut_arg, 
                                          args_info.voxelSize_arg,
                                          args_info.gaussWidth_arg);
  } else if (args_info.ior_given){
    if (args_info.rcut_given) {
      analyser = new IcosahedralOfR(info, sele1, 
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
      analyser = new FCCOfR(info, sele1, args_info.rcut_arg,
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
      analyser = new BondAngleDistribution(info, sele1, 
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
      analyser  = new SCDOrderParameter(info, 
                                        args_info.molname_arg, 
					args_info.begin_arg, args_info.end_arg);
    } else{
      analyser  = new SCDOrderParameter(info, 
                                        sele1, sele2, sele3);
    }
  }else if (args_info.density_given) {
    analyser= new DensityPlot(info, sele1, sele2, maxLen,
			      args_info.nbins_arg);  
  } else if (args_info.count_given) {
    analyser = new ObjectCount(info, sele1 );
  } else if (args_info.slab_density_given) {
    analyser = new RhoZ(info, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.pipe_density_given) {

    switch (privilegedAxis) {
    case 0:      
      analyser = new PipeDensity(info, sele1,
                                 args_info.nbins_y_arg, args_info.nbins_z_arg,
                                 privilegedAxis);
      break;
    case 1:
      analyser = new PipeDensity(info, sele1,
                                 args_info.nbins_z_arg, args_info.nbins_x_arg,
                                 privilegedAxis);      
      break;
    case 2:
    default:
      analyser = new PipeDensity(info, sele1,
                                 args_info.nbins_x_arg, args_info.nbins_y_arg,
                                 privilegedAxis);            
      break;
    }
  } else if (args_info.rnemdz_given) {
    analyser = new RNEMDZ(info, sele1, args_info.nbins_arg, privilegedAxis);
  } else if (args_info.rnemdr_given) {
    analyser = new RNEMDR(info, sele1, nrbins);
  } else if (args_info.rnemdrt_given) {
    analyser = new RNEMDRTheta(info, sele1, nrbins, nanglebins);
  } else if (args_info.nitrile_given) {
    analyser = new NitrileFrequencyMap(info, sele1,
                                       args_info.nbins_arg);
  } else if (args_info.p_angle_given) {
    if (args_info.sele1_given) {     
      if (args_info.sele2_given) 
        analyser  = new pAngle(info, sele1, sele2,
                               args_info.nbins_arg);
      else 
        if (args_info.seleoffset_given) {
          if (args_info.seleoffset2_given) {
            analyser  = new pAngle(info, sele1, 
                                   args_info.seleoffset_arg, 
                                   args_info.seleoffset2_arg, 
                                   args_info.nbins_arg);
          } else {
            analyser  = new pAngle(info, sele1, 
                                   args_info.seleoffset_arg, 
                                   args_info.nbins_arg);
          }
        } else 
          analyser  = new pAngle(info, sele1, 
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
    analyser = new Hxy(info, sele1, args_info.nbins_x_arg, 
		       args_info.nbins_y_arg, args_info.nbins_z_arg,
                       args_info.nbins_arg);
#endif
  }else if(args_info.cn_given || args_info.scn_given || args_info.gcn_given){
    if (args_info.rcut_given) {
      if (args_info.cn_given) {
        analyser = new CoordinationNumber(info, sele1, sele2,
                                          args_info.rcut_arg,
                                          args_info.nbins_arg);
      } else if (args_info.scn_given) {
        analyser = new SCN(info, sele1, sele2,
                           args_info.rcut_arg, args_info.nbins_arg);
      } else if (args_info.gcn_given) {
        analyser = new GCN(info, sele1, sele2,
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
  else if (args_info.rho_r_given) {
    if (args_info.radius_given){
      analyser = new RhoR(info, sele1, maxLen, nrbins,
                          args_info.radius_arg);
    }else{
      sprintf( painCave.errMsg,
	       "A particle radius (radius) must be specified when calculating Rho(r)");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
  } else if (args_info.angle_r_given) {
    analyser = new AngleR(info, sele1, maxLen, nrbins);
  } else if (args_info.hbond_given){
    if (args_info.rcut_given) {
      if (args_info.thetacut_given) {
        
        analyser = new HBondGeometric(info, sele1, sele2,
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
    analyser = new PotDiff(info, sele1);
  } else if (args_info.kirkwood_given) {
    analyser= new Kirkwood(info, sele1, sele2, maxLen, 
                           nrbins);
  } else if (args_info.kirkwoodQ_given) {
    analyser= new KirkwoodQuadrupoles(info, sele1, sele2, maxLen, 
                                      nrbins);
  } else if (args_info.densityfield_given) {
    analyser = new DensityField(info, sele1, args_info.voxelSize_arg);
  } else if (args_info.velocityfield_given) {
    analyser = new VelocityField(info, sele1, args_info.voxelSize_arg);
  } else if (args_info.velocityZ_given) {

    switch (privilegedAxis) {
    case 0:
      if (privilegedAxis2 == 1) {
	analyser = new VelocityZ(info, sele1,
                                 args_info.nbins_x_arg, args_info.nbins_y_arg,
                                 privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 2) {
	analyser = new VelocityZ(info, sele1,
                                 args_info.nbins_x_arg, args_info.nbins_z_arg,
                                 privilegedAxis, privilegedAxis2);
      }
      break;
    case 1:
      if (privilegedAxis2 == 0) {
	analyser = new VelocityZ(info, sele1,
				   args_info.nbins_y_arg, args_info.nbins_x_arg,
				   privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 2) {
	analyser = new VelocityZ(info, sele1,
				   args_info.nbins_y_arg, args_info.nbins_z_arg,
				   privilegedAxis, privilegedAxis2);
      }
      break;
    case 2:
    default:
      if (privilegedAxis2 == 0) {
	analyser = new VelocityZ(info, sele1,
				   args_info.nbins_z_arg, args_info.nbins_x_arg,
				   privilegedAxis, privilegedAxis2);
      } else if (privilegedAxis2 == 1) {
	analyser = new VelocityZ(info, sele1,
				   args_info.nbins_z_arg, args_info.nbins_y_arg,
				   privilegedAxis, privilegedAxis2);
      }
      break;
    }
  } 

  if (args_info.output_given) {
    analyser->setOutputName(args_info.output_arg);
  }
  if (args_info.step_given) {
    analyser->setStep(args_info.step_arg);
  }
  
  analyser->processDump();
  
  delete analyser;    
  delete info;
  
  return 0;   
}
