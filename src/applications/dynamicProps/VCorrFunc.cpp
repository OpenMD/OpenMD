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

#include "applications/dynamicProps/VCorrFunc.hpp"

namespace OpenMD {
  VCorrFunc::VCorrFunc(SimInfo* info, const std::string& filename, 
                       const std::string& sele1, const std::string& sele2, 
                       long long int memSize)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2, 
                           DataStorage::dslVelocity, memSize){
    
    setCorrFuncType("Velocity Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".vcorr");
    
  }

  VCorrFuncZ::VCorrFuncZ(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2,
                         long long int memSize)
    : VCorrFunc(info, filename, sele1, sele2, memSize){

    setCorrFuncType("Velocity Correlation Function projected along z axis");
    setOutputName(getPrefix(dumpFilename_) + ".vcorrz");

  }
  VCorrFuncR::VCorrFuncR(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2,
                         long long int memSize)
    : VCorrFunc(info, filename, sele1, sele2, memSize){

    // Turn on COM calculation in block snapshot
    bool ncp = true;
    bsMan_->needCOMprops(ncp);
    
    setCorrFuncType("Velocity Correlation Function (radial projection)");
    setOutputName(getPrefix(dumpFilename_) + ".vcorrr");
      
  }

  RealType VCorrFunc::calcCorrVal(int frame1, int frame2, StuntDouble* sd1,  
                                  StuntDouble* sd2) {
    Vector3d v1 = sd1->getVel(frame1);
    Vector3d v2 = sd2->getVel(frame2);
    
    return dot(v1, v2);
  }


  RealType VCorrFuncZ::calcCorrVal(int frame1, int frame2, StuntDouble* sd1,  
                                  StuntDouble* sd2) {
    Vector3d v1 = sd1->getVel(frame1);
    Vector3d v2 = sd2->getVel(frame2);
    
    return v1.z() * v2.z();
  }
  
  RealType VCorrFuncR::calcCorrVal(int frame1, int frame2, StuntDouble* sd1,  
                                  StuntDouble* sd2) {

    Vector3d coord_t0;
    Vector3d coord_t;

    Vector3d r1 = sd1->getPos(frame1);
    Vector3d r2 = sd2->getPos(frame2);

    Vector3d com1 = sd1->getCOM(frame1);
    Vector3d com2 = sd2->getCOM(frame2);

    Vector3d v1 = sd1->getVel(frame1);
    Vector3d v2 = sd2->getVel(frame2);
    
    coord_t0 = r1 - com1;
    coord_t  = r2 - com2;

    coord_t0.normalize();
    coord_t.normalize();

    // project velocity vectors onto the radial vectors:

    RealType v1r = dot(v1, coord_t0);
    RealType v2r = dot(v2, coord_t);

    return v1r*v2r;

  }
}

