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

#include "applications/staticProps/SpatialStatistics.hpp"

#include <string>
#include <vector>

#include "applications/staticProps/StaticAnalyser.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "io/Globals.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "utils/Accumulator.hpp"
#include "utils/AccumulatorView.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  SpatialStatistics::SpatialStatistics(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele, int nbins) :
      StaticAnalyser(info, filename, nbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    setOutputName(getPrefix(filename) + ".spst");
  }

  void SpatialStatistics::process() {
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      processFrame(istep);
    }

    writeOutput();
  }

  void SpatialStatistics::processFrame(int) {
    StuntDouble* sd;
    int i;

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // loop over the selected atoms:
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {
      // figure out where that object is:
      Vector3d pos = sd->getPos();

      int bin = getBin(pos);

      // forward the work of statistics on to the subclass:
      processStuntDouble(sd, bin);
    }
  }

  SlabStatistics::SlabStatistics(SimInfo* info, const std::string& filename,
                                 const std::string& sele, int nbins, int axis) :
      SpatialStatistics(info, filename, sele, nbins),
      axis_(axis) {
    // Set the axis label for the privileged axis
    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      break;
    case 1:
      axisLabel_ = "y";
      break;
    case 2:
    default:
      axisLabel_ = "z";
      break;
    }
  }

  void SlabStatistics::processFrame(int istep) {
    RealType z;
    hmat_ = currentSnapshot_->getHmat();

    volume_ = currentSnapshot_->getVolume();

    SpatialStatistics::processFrame(istep);
  }

  int SlabStatistics::getBin(Vector3d pos) {
    currentSnapshot_->wrapVector(pos);
    // which bin is this stuntdouble in?
    // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]
    // Shift molecules by half a box to have bins start at 0
    // The modulo operator is used to wrap the case when we are
    // beyond the end of the bins back to the beginning.
    return int(nBins_ * (pos[axis_] / hmat_(axis_, axis_) + 0.5)) % nBins_;
  }

  ShellStatistics::ShellStatistics(SimInfo* info, const std::string& filename,
                                   const std::string& sele,
                                   const std::string& comsele, int nbins,
                                   RealType binWidth) :
      SpatialStatistics(info, filename, sele, nbins),
      coordinateOrigin_(V3Zero), comSele_(comsele), comSeleMan_(info),
      comEvaluator_(info), binWidth_(binWidth) {
    Globals* simParams                  = info->getSimParams();
    RNEMD::RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();
    bool hasCoordinateOrigin            = rnemdParams->haveCoordinateOrigin();

    if (hasCoordinateOrigin) {
      std::vector<RealType> co = rnemdParams->getCoordinateOrigin();
      if (co.size() != 3) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "ShellStatistics: Incorrect number of parameters specified for "
            "coordinateOrigin.\n"
            "\tthere should be 3 parameters, but %zu were specified.\n",
            co.size());
        painCave.isFatal = 1;
        simError();
      }
      coordinateOrigin_.x() = co[0];
      coordinateOrigin_.y() = co[1];
      coordinateOrigin_.z() = co[2];
    } else {
      if (!comSele_.empty()) {
        comEvaluator_.loadScriptString(comSele_);
        if (!comEvaluator_.isDynamic()) {
          comSeleMan_.setSelectionSet(comEvaluator_.evaluate());
          if (comSeleMan_.getSelectionCount() != 1) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "ShellStatistics: More than one selected object in "
                     "comsele.\n"
                     "\tThere are %d selected objects.\n",
                     comSeleMan_.getSelectionCount());
            painCave.isFatal = 1;
            simError();
          }
          int isd;
          StuntDouble* sd   = comSeleMan_.beginSelected(isd);
          coordinateOrigin_ = sd->getPos();
        }
      } else {
        coordinateOrigin_ = V3Zero;
      }
    }
  }

  void ShellStatistics::processFrame(int istep) {
    if (comEvaluator_.isDynamic()) {
      comSeleMan_.setSelectionSet(comEvaluator_.evaluate());
      if (comSeleMan_.getSelectionCount() != 1) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "ShellStatistics: More than one selected object in "
                 "comsele.\n"
                 "\tThere are %d selected objects.\n",
                 comSeleMan_.getSelectionCount());
        painCave.isFatal = 1;
        simError();
      }
      int isd;
      StuntDouble* sd   = comSeleMan_.beginSelected(isd);
      coordinateOrigin_ = sd->getPos();

      SpatialStatistics::processFrame(istep);
    }
  }

  int ShellStatistics::getBin(Vector3d pos) {
    Vector3d rPos = pos - coordinateOrigin_;
    return int(rPos.length() / binWidth_);
  }
}  // namespace OpenMD
