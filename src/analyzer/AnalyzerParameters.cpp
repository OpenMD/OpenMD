/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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
#include <cstdlib>
#include <cstring>

#include "analyzer/AnalyzerParameters.hpp"

namespace OpenMD {
  AnalyzerParameters::AnalyzerParameters() { 
    DefineOptionalParameterWithDefaultValue(UseAnalyzer, "useAnalyzer", false);
    DefineOptionalParameterWithDefaultValue(QueryTime, "queryTime", 100.0);
    DefineOptionalParameter(Method, "method");

    DefineOptionalParameter(Step, "step");
    DefineOptionalParameter(NBins, "nBins");
    DefineOptionalParameter(NBinsX, "nBins_x");
    DefineOptionalParameter(NBinsY, "nBins_y");
    DefineOptionalParameter(NBinsZ, "nBins_z");
    DefineOptionalParameter(NRBins, "nRBins");
    DefineOptionalParameter(NAngleBins, "nAngleBins");
    DefineOptionalParameter(RCut, "rCut");
    DefineOptionalParameter(OOCut, "OOCut");
    DefineOptionalParameter(ThetaCut, "thetaCut");
    DefineOptionalParameter(OHCut, "OHCut");
    DefineOptionalParameter(Dz, "dz");
    DefineOptionalParameter(Length, "length");
    DefineOptionalParameter(ZLength, "zLength");
    DefineOptionalParameter(ZOffSet, "zOffSet");
    DefineOptionalParameter(Sele1, "sele1");
    DefineOptionalParameter(Sele2, "sele2");
    DefineOptionalParameter(Sele3, "sele3");
    DefineOptionalParameter(ComSele, "comSele");
    DefineOptionalParameter(SeleOffSet, "seleOffSet");
    DefineOptionalParameter(SeleOffSet2, "seleOffSet2");
    DefineOptionalParameter(MolName, "molName");
    DefineOptionalParameter(Begin, "begin");
    DefineOptionalParameter(End, "end");
    DefineOptionalParameter(Radius, "radius");
    DefineOptionalParameter(VoxelSize, "voxelSize");
    DefineOptionalParameter(GaussWidth, "gaussWidth");
    DefineOptionalParameter(PrivilegedAxis, "privilegedAxis");
    DefineOptionalParameter(PrivilegedAxis2, "privilegedAxis2");


  }
  
  AnalyzerParameters::~AnalyzerParameters() {    
  }
  
  void AnalyzerParameters::validate() {
    CheckParameter(QueryTime, isPositive());
    CheckParameter(Method,
		   isEqualIgnoreCase("bo") ||
		   isEqualIgnoreCase("ior") ||
		   isEqualIgnoreCase("for") ||
		   isEqualIgnoreCase("bad") ||
		   isEqualIgnoreCase("count") ||
		   isEqualIgnoreCase("gofr") ||
		   isEqualIgnoreCase("gofz") ||
		   isEqualIgnoreCase("r_theta") ||
		   isEqualIgnoreCase("r_omega") ||
		   isEqualIgnoreCase("r_z") ||
		   isEqualIgnoreCase("theta_omega") ||
		   isEqualIgnoreCase("r_theta_omega") ||
		   isEqualIgnoreCase("gxyz") ||
		   isEqualIgnoreCase("twodgofr") ||
		   isEqualIgnoreCase("p2") ||
		   isEqualIgnoreCase("scd") ||
		   isEqualIgnoreCase("density") ||
		   isEqualIgnoreCase("slab_density") ||
		   isEqualIgnoreCase("pipe_density") ||
		   isEqualIgnoreCase("p_angle") ||
		   isEqualIgnoreCase("hxy") ||
		   isEqualIgnoreCase("rho_r") ||
		   isEqualIgnoreCase("angle_r") ||
		   isEqualIgnoreCase("hullvol") ||
		   isEqualIgnoreCase("rodlength") ||
		   isEqualIgnoreCase("tet_param") ||
		   isEqualIgnoreCase("tet_param_z") ||
		   isEqualIgnoreCase("tet_param_dens") ||
		   isEqualIgnoreCase("tet_param_xyz") ||
		   isEqualIgnoreCase("rnemdz") ||
		   isEqualIgnoreCase("rnemdr") ||
		   isEqualIgnoreCase("rnemdrt") ||
		   isEqualIgnoreCase("nitrile") ||
		   isEqualIgnoreCase("multipole") ||
		   isEqualIgnoreCase("surfDiffusion") ||
		   isEqualIgnoreCase("cn") ||
		   isEqualIgnoreCase("scn") ||
		   isEqualIgnoreCase("gcn") ||
		   isEqualIgnoreCase("hbond") ||
		   isEqualIgnoreCase("potDiff") ||
		   isEqualIgnoreCase("tet_hb") ||
		   isEqualIgnoreCase("kirkwood") ||
		   isEqualIgnoreCase("kirkwoodQ") ||
		   isEqualIgnoreCase("densityfield") ||
		   isEqualIgnoreCase("velocityfield") ||
		   isEqualIgnoreCase("velocityZ")
		   );
    CheckParameter(Step, isPositive());
    CheckParameter(NBins, isPositive());
    CheckParameter(NBinsX, isPositive());
    CheckParameter(NBinsY, isPositive());
    CheckParameter(NBinsZ, isPositive());
    CheckParameter(NRBins, isPositive());
    CheckParameter(NAngleBins, isPositive());
    CheckParameter(RCut, isPositive());
    CheckParameter(OOCut, isPositive());
    CheckParameter(ThetaCut, isPositive());
    CheckParameter(OHCut, isPositive());
    CheckParameter(Dz, isPositive());
    CheckParameter(Length, isPositive());
    CheckParameter(ZLength, isPositive());
    CheckParameter(ZOffSet, isPositive());
    CheckParameter(SeleOffSet, isPositive());
    CheckParameter(SeleOffSet2, isPositive());
    CheckParameter(Begin, isPositive());
    CheckParameter(End, isPositive());
    CheckParameter(Radius, isPositive());
    CheckParameter(VoxelSize, isPositive());
    CheckParameter(GaussWidth, isPositive());
    CheckParameter(PrivilegedAxis, 
                   isEqualIgnoreCase("x") ||  
                   isEqualIgnoreCase("y")  ||
                   isEqualIgnoreCase("z"));
    CheckParameter(PrivilegedAxis2, 
                   isEqualIgnoreCase("x") ||  
                   isEqualIgnoreCase("y")  ||
                   isEqualIgnoreCase("z"));
  }
  
}
