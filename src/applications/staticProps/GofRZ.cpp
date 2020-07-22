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

#include <algorithm>
#include <fstream>
#include "applications/staticProps/GofRZ.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofRZ::GofRZ(SimInfo* info, const std::string& filename, const std::string& sele1, 
               const std::string& sele2, RealType len, RealType zlen, int nrbins, int nZBins, int axis)
    : RadialDistrFunc(info, filename, sele1, sele2, nrbins), len_(len), zLen_(zlen), nZBins_(nZBins), axis_(axis){

    setOutputName(getPrefix(filename) + ".gofrz");

    deltaR_ = len_ / (double) nBins_;
    deltaZ_ = zLen_ / (double)nZBins_; 

    histogram_.resize(nBins_);
    avgGofr_.resize(nBins_);
    for (unsigned int i = 0 ; i < nBins_; ++i) {
      histogram_[i].resize(nZBins_);
      avgGofr_[i].resize(nZBins_);
    }

    // Compute complementary axes to the privileged axis
    xaxis_ = (axis_ + 1) % 3;
    yaxis_ = (axis_ + 2) % 3;

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
  }

  void GofRZ::preProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      std::fill(avgGofr_[i].begin(), avgGofr_[i].end(), 0);
    }
  }

  void GofRZ::initializeHistogram() {
    npairs_ = 0;
    for (unsigned int i = 0; i < histogram_.size(); ++i){
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0);
    }
  }
  
  void GofRZ::processHistogram() {
    int nPairs = getNPairs();
    RealType volume = info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = nPairs / volume * 2.0;

    for(unsigned int i = 0 ; i < histogram_.size(); ++i){

      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice = Constants::PI * deltaZ_ * (( rUpper * rUpper ) - ( rLower * rLower ));
      RealType nIdeal = volSlice * pairDensity;

      for (unsigned int j = 0; j < histogram_[i].size(); ++j){
        avgGofr_[i][j] += histogram_[i][j] / nIdeal;    
      }
    }

  }

  void GofRZ::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    if (sd1 == sd2) {
      return;
    }
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = sqrt(pow(r12[xaxis_], 2) + pow(r12[yaxis_], 2));

    int whichRBin = int(distance / deltaR_);

    if (distance <= len_) {
     
      RealType Z = fabs(r12[axis_]);

      if (Z <= zLen_) {
        int whichZBin = int(Z / deltaZ_);
              
        ++histogram_[whichRBin][whichZBin];        
        ++npairs_;
      }
    }
  }

  void GofRZ::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#radial distribution function\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      rdfStream << "#nBins = " << nBins_ << "\t maxLen = " << len_ << "deltaR = " << deltaR_ <<"\n";
      rdfStream << "#n" << axisLabel_ << "Bins =" << nZBins_ << "\t delta" << axisLabel_ << " = " << deltaZ_ << "\n";
      for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
        // RealType r = deltaR_ * (i + 0.5);

        for(unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
          // RealType z = deltaZ_ * (j + 0.5);
          rdfStream << avgGofr_[i][j]/nProcessed_ << "\t";
        }

        rdfStream << "\n";
      }
        
    } else {
      sprintf(painCave.errMsg, "GofRZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();
  }

}


