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

#include <string>
#include <vector>

#include "applications/staticProps/SpatialStatistics.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "io/Globals.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "utils/Accumulator.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  SpatialStatistics::SpatialStatistics(SimInfo* info, const std::string& filename,
                                       const std::string& sele, int nbins)
    : StaticAnalyser(info, filename, nbins), selectionScript_(sele),
      evaluator_(info), seleMan_(info) {

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    setOutputName(getPrefix(filename) + ".spst");
  }

  SpatialStatistics::~SpatialStatistics() {

    std::vector<OutputData*>::iterator i;
    OutputData* outputData;

    for(outputData = beginOutputData(i); outputData;
        outputData = nextOutputData(i)) {
      delete outputData;
    }
    data_.clear();
  }

  void SpatialStatistics::process() {

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      processFrame(istep);
    }

    writeOutput();
  }

  void SpatialStatistics::processFrame(int istep) {

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
      processStuntDouble( sd, bin );
    }
  }


  SlabStatistics::SlabStatistics(SimInfo* info, const std::string& filename,
                                 const std::string& sele, int nbins, int axis)
    : SpatialStatistics(info, filename, sele, nbins), axis_(axis) {

    // Set the axis label for the privileged axis
    switch(axis_) {
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

    z_ = new OutputData;
    z_->units =  "Angstroms";
    z_->title =  axisLabel_;
    z_->dataType = odtReal;
    z_->dataHandling = odhAverage;
    z_->accumulator.reserve(nbins);
    for (int i = 0; i < nbins; i++)
      z_->accumulator.push_back( new Accumulator() );
    data_.push_back(z_);
  }

  SlabStatistics::~SlabStatistics() {
  }

  void SlabStatistics::processFrame(int istep) {

    RealType z;
    hmat_ = currentSnapshot_->getHmat();

    for (unsigned int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_,axis_);
      dynamic_cast<Accumulator*>(z_->accumulator[i])->add(z);
    }

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
    return int(nBins_ * (pos[axis_] / hmat_(axis_,axis_) + 0.5)) % nBins_;
  }


  ShellStatistics::ShellStatistics(SimInfo* info, const std::string& filename,
                                   const std::string& sele, int nbins)
    : SpatialStatistics(info, filename, sele, nbins), coordinateOrigin_(V3Zero) {

    binWidth_ = 1.0;

    Globals* simParams = info->getSimParams();
    RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();
    bool hasCoordinateOrigin = rnemdParams->haveCoordinateOrigin();

    if (hasCoordinateOrigin) {
      std::vector<RealType> co = rnemdParams->getCoordinateOrigin();
      if (co.size() != 3) {
        sprintf(painCave.errMsg,
                "RNEMD: Incorrect number of parameters specified for coordinateOrigin.\n"
                "\tthere should be 3 parameters, but %lu were specified.\n",
                co.size());
        painCave.isFatal = 1;
        simError();
      }
      coordinateOrigin_.x() = co[0];
      coordinateOrigin_.y() = co[1];
      coordinateOrigin_.z() = co[2];
    } else {
      coordinateOrigin_ = V3Zero;
    }

    r_ = new OutputData;
    r_->units =  "Angstroms";
    r_->title =  "R";
    r_->dataType = odtReal;
    r_->dataHandling = odhAverage;
    r_->accumulator.reserve(nbins);
    for (int i = 0; i < nbins; i++)
      r_->accumulator.push_back( new Accumulator() );
    data_.push_back(r_);

    for (int i = 0; i < nbins; i++) {
      RealType r = (((RealType)i + 0.5) * binWidth_);
      dynamic_cast<Accumulator*>(r_->accumulator[i])->add(r);
    }
  }

  ShellStatistics::~ShellStatistics() {
  }

  int ShellStatistics::getBin(Vector3d pos) {
    Vector3d rPos = pos - coordinateOrigin_;
    return int(rPos.length() / binWidth_);
  }
}
