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

#ifndef APPLICATIONS_STATICPROPS_SPATIALSTATISTICS_HPP
#define APPLICATIONS_STATICPROPS_SPATIALSTATISTICS_HPP

#include <string>
#include <vector>

#include "applications/staticProps/StaticAnalyser.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Snapshot.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  class SpatialStatistics : public StaticAnalyser {
  public:
    SpatialStatistics(SimInfo* info, const std::string& filename,
                      const std::string& sele, int nbins);

    virtual void process();
    virtual void processFrame(int frame);
    virtual int getBin(Vector3d pos)                          = 0;
    virtual void processStuntDouble(StuntDouble* sd, int bin) = 0;

  protected:
    Snapshot* currentSnapshot_;
    int nProcessed_;
    std::string selectionScript_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;
  };

  class SlabStatistics : public SpatialStatistics {
  public:
    SlabStatistics(SimInfo* info, const std::string& filename,
                   const std::string& sele, int nbins, int axis);

    virtual int getBin(Vector3d pos);
    virtual void processFrame(int frame);

  protected:
    Mat3x3d hmat_;
    RealType volume_;
    int axis_;
    std::string axisLabel_;
  };

  class ShellStatistics : public SpatialStatistics {
  public:
    ShellStatistics(SimInfo* info, const std::string& filename,
                    const std::string& sele, const std::string& comsele,
                    int nbins, RealType binWidth);

    virtual void processFrame(int frame);
    virtual int getBin(Vector3d pos);
    void setCoordinateOrigin(Vector3d co) { coordinateOrigin_ = co; }

  protected:
    Vector3d coordinateOrigin_;
    RealType binWidth_;
    std::string comSele_;
    SelectionManager comSeleMan_;
    SelectionEvaluator comEvaluator_;
  };
}  // namespace OpenMD

#endif
