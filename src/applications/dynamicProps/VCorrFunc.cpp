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

#include "applications/dynamicProps/VCorrFunc.hpp"

namespace OpenMD {
  VCorrFunc::VCorrFunc(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2) :
      ObjectACF<RealType>(info, filename, sele1, sele2) {
    setCorrFuncType("Velocity Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".vcorr");
    setLabelString("<v(0).v(t)>");
    velocities_.resize(nFrames_);
  }

  VCorrFuncZ::VCorrFuncZ(SimInfo* info, const std::string& filename,
                         const std::string& sele1, const std::string& sele2) :
      ObjectACF<RealType>(info, filename, sele1, sele2) {
    setCorrFuncType("Velocity Correlation Function projected along z axis");
    setOutputName(getPrefix(dumpFilename_) + ".vcorrz");
    setLabelString("<vz(0).vz(t)>");
    velocities_.resize(nFrames_);
  }

  VCorrFuncR::VCorrFuncR(SimInfo* info, const std::string& filename,
                         const std::string& sele1, const std::string& sele2) :
      ObjectACF<RealType>(info, filename, sele1, sele2) {
    // Turn on COM calculation in reader:
    bool ncp = true;
    reader_->setNeedCOMprops(ncp);
    setCorrFuncType("Velocity Correlation Function (radial projection)");
    setOutputName(getPrefix(dumpFilename_) + ".vcorrr");
    setLabelString("<vr(0).vr(t)>");
    velocities_.resize(nFrames_);
  }

  int VCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    velocities_[frame].push_back(sd->getVel());
    return velocities_[frame].size() - 1;
  }

  RealType VCorrFunc::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType v2 = dot(velocities_[frame1][id1], velocities_[frame2][id2]);
    return v2;
  }

  int VCorrFuncZ::computeProperty1(int frame, StuntDouble* sd) {
    velocities_[frame].push_back(sd->getVel().z());
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
    velocities_[frame].push_back(vel);
    return velocities_[frame].size() - 1;
  }

  RealType VCorrFuncR::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType v2;
    v2 = velocities_[frame1][id1] * velocities_[frame2][id2];
    return v2;
  }
}  // namespace OpenMD
