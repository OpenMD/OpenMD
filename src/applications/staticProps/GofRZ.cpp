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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/GofRZ.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofRZ::GofRZ(SimInfo* info, const std::string& filename, const std::string& sele1, 
               const std::string& sele2, RealType len, RealType zlen, int nrbins, int nZBins)
    : RadialDistrFunc(info, filename, sele1, sele2), len_(len), zLen_(zlen), nRBins_(nrbins), nZBins_(nZBins){

    setOutputName(getPrefix(filename) + ".gofrz");

    deltaR_ = len_ / (double) nRBins_;
    deltaZ_ = zLen_ / (double)nZBins_;    // for solvated_NVT.md4

    histogram_.resize(nRBins_);
    avgGofr_.resize(nRBins_);
    for (int i = 0 ; i < nRBins_; ++i) {
      histogram_[i].resize(nZBins_);
      avgGofr_[i].resize(nZBins_);
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
      RealType volSlice = NumericConstant::PI * deltaZ_ * (( rUpper * rUpper ) - ( rLower * rLower ));
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
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = sqrt(pow(r12.x(), 2) + pow(r12.y(), 2));

    int whichRBin = distance / deltaR_;

    if (distance <= len_) {
     
      RealType Z = fabs(r12.z());

      if (Z <= zLen_) {
        int whichZBin = Z / deltaZ_;
              
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
      rdfStream << "#nRBins = " << nRBins_ << "\t maxLen = " << len_ << "deltaR = " << deltaR_ <<"\n";
      rdfStream << "#nZBins =" << nZBins_ << "\t deltaZ = " << deltaZ_ << "\n";
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


