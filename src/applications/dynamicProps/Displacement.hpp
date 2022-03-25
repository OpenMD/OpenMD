/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef APPLICATIONS_DYNAMICPROPS_DISPLACEMENT_HPP
#define APPLICATIONS_DYNAMICPROPS_DISPLACEMENT_HPP

#include "applications/dynamicProps/TimeCorrFunc.hpp"

namespace OpenMD {

  class Displacement : public ObjectACF<Vector3d> {
  public:
    Displacement(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2);

  private:
    virtual int computeProperty1(int frame, StuntDouble* sd);
    virtual Vector3d calcCorrVal(int frame1, int frame2, int id1, int id2);
    std::vector<std::vector<Vector3d>> positions_;
  };

  class DisplacementZ : public ObjectACF<Vector3d> {
  public:
    DisplacementZ(SimInfo* info, const std::string& filename,
                  const std::string& sele1, const std::string& sele2,
                  int nZbins, int axis = 2);

  private:
    virtual void computeFrame(int frame);
    virtual int computeProperty1(int frame, StuntDouble* sd);
    virtual void correlateFrames(int frame1, int frame2, int timeBin);
    Vector3d calcCorrValImpl(int frame1, int frame2, int id1, int id2,
                             int timeBin);
    virtual Vector3d calcCorrVal(int, int, int, int) { return Vector3d(0.0); }

    virtual void postCorrelate();
    virtual void writeCorrelate();

    std::vector<std::vector<Vector3d>> positions_;
    std::vector<std::vector<int>> zBins_;
    std::vector<std::vector<Vector3d>> histograms_;
    std::vector<std::vector<int>> counts_;
    Mat3x3d hmat_;
    RealType halfBoxZ_;
    unsigned int nZBins_;
    int axis_;
    std::string axisLabel_;
  };
}  // namespace OpenMD

#endif
