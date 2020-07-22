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

#include "applications/dynamicProps/VCorrFunc.hpp"

namespace OpenMD {
  VCorrFunc::VCorrFunc(SimInfo* info, const std::string& filename, 
                       const std::string& sele1, const std::string& sele2)
    : ObjectACF<RealType>(info, filename, sele1, sele2, 
                          DataStorage::dslVelocity | DataStorage::dslAmat | 
                          DataStorage::dslAngularMomentum){
    
    setCorrFuncType("Velocity Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".vcorr");
    setLabelString( "<v(0).v(t)>" );
    velocities_.resize(nFrames_);
  }

  VCorrFuncZ::VCorrFuncZ(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2)
    : ObjectACF<RealType>(info, filename, sele1, sele2, 
                          DataStorage::dslPosition |
                          DataStorage::dslVelocity |
                          DataStorage::dslAmat |
                          DataStorage::dslAngularMomentum){
    
    setCorrFuncType("Velocity Correlation Function projected along z axis");
    setOutputName(getPrefix(dumpFilename_) + ".vcorrz");
    setLabelString( "<vz(0).vz(t)>" );
    velocities_.resize(nFrames_);
  }
  VCorrFuncR::VCorrFuncR(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2)
    : ObjectACF<RealType>(info, filename, sele1, sele2, 
                          DataStorage::dslPosition |
                          DataStorage::dslVelocity |
                          DataStorage::dslAmat |
                          DataStorage::dslAngularMomentum){
    
    // Turn on COM calculation in reader:
    bool ncp = true;
    reader_->setNeedCOMprops(ncp);    
    setCorrFuncType("Velocity Correlation Function (radial projection)");
    setOutputName(getPrefix(dumpFilename_) + ".vcorrr");
    setLabelString( "<vr(0).vr(t)>" );
    velocities_.resize(nFrames_);
  }

  int VCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    velocities_[frame].push_back( sd->getVel() );
    return velocities_[frame].size() - 1;
  }
  RealType VCorrFunc::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType v2 = dot( velocities_[frame1][id1] , velocities_[frame2][id2]);
    return v2;
  }

  int VCorrFuncZ::computeProperty1(int frame, StuntDouble* sd) {
    velocities_[frame].push_back( sd->getVel().z() );
    return velocities_[frame].size() - 1;
  }
  RealType VCorrFuncZ::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType v2 = velocities_[frame1][id1] * velocities_[frame2][id2];
    return v2;
  }

  int VCorrFuncR::computeProperty1(int frame, StuntDouble* sd) {
    // get the radial vector from the frame's center of mass:
    Vector3d coord_t = sd->getPos() - sd->getCOM();
    coord_t.normalize();
      
    // project velocity vectors onto the radial vectors:
    RealType vel = dot(sd->getVel(), coord_t);
    velocities_[frame].push_back( vel );
    return velocities_[frame].size() - 1;
  }
  RealType VCorrFuncR::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType v2;
    v2  = velocities_[frame1][id1] * velocities_[frame2][id2];
    return v2;
  }
}

