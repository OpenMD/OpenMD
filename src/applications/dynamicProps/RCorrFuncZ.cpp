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

#include "applications/dynamicProps/RCorrFuncZ.hpp"
#include "utils/simError.h"

namespace OpenMD {
  RCorrFuncZ::RCorrFuncZ(SimInfo* info, 
                         const std::string& filename, 
                         const std::string& sele1, 
                         const std::string& sele2, 
                         int nZbins,
                         long long int memSize)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2, 
                           DataStorage::dslPosition, memSize), nZBins_(nZbins) {
    
      setCorrFuncType("RMSD binned by Z");
      setOutputName(getPrefix(dumpFilename_) + ".rcorrZ");
      histogram_.resize(nTimeBins_);
      counts_.resize(nTimeBins_);
      for (int i = 0; i < nTimeBins_; i++) {
        histogram_[i].resize(nZBins_);
        counts_[i].resize(nZBins_);
      }
    }

  void RCorrFuncZ::correlateFrames(int frame1, int frame2) {
    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    Mat3x3d hmat = snapshot1->getHmat();
    RealType halfBoxZ_ = hmat(2,2) / 2.0;      

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();

    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);

    int i;
    int j;
    int id1, id2;
    StuntDouble* sd1;
    StuntDouble* sd2;

    for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan2_.beginSelected(j); 
         sd1 != NULL && sd2 != NULL;
         sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j) ) {
      
      id1 = sd1->getGlobalIndex();
      id2 = sd2->getGlobalIndex();

      // If the selections are dynamic, they might not have the same
      // objects in both frames, so we need to roll either of the
      // selections until we have the same particle to correlate.

      while (id1 < id2 && sd1 != NULL) {
        sd1 = seleMan1_.nextSelected(i);
        if (sd1 != NULL) id1 = sd1->getGlobalIndex();
      }
      while (id2 < id1 && sd2 != NULL) {
        sd2 = seleMan2_.nextSelected(j);
        if (sd2 != NULL) id2 = sd2->getGlobalIndex();
      }
     
      if (sd1 == NULL || sd2 == NULL) break;

      Vector3d pos1 = sd1->getPos(frame1);
      Vector3d pos2 = sd2->getPos(frame2);
      
      if (info_->getSimParams()->getUsePeriodicBoundaryConditions()) {
        snapshot1->wrapVector(pos1);
        snapshot2->wrapVector(pos2);
      }
      
      int zBin1 = int(nZBins_ * (halfBoxZ_ + pos1.z()) / hmat(2,2));
      int zBin2 = int(nZBins_ * (halfBoxZ_ + pos1.z()) / hmat(2,2));

      if (zBin1 == zBin2) {
        // we need the raw (not wrapped) positions for RMSD:
        pos1 = sd1->getPos(frame1);
        pos2 = sd2->getPos(frame2);
        RealType corrVal = (pos2-pos1).lengthSquare();

        histogram_[timeBin][zBin1] += corrVal; 
        counts_[timeBin][zBin1]++;
      }
    }    
  }

  void RCorrFuncZ::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      for (int j = 0; j < nZBins_; ++j) {
        if (counts_[i][j] > 0) {
          histogram_[i][j] /= counts_[i][j];
        }
      }
    }
  }

  void RCorrFuncZ::preCorrelate() {
    for (int i = 0; i < nTimeBins_; i++) {
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0.0);
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
  }


  void RCorrFuncZ::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#time\t<r^2>(z)\n";

      for (int i = 0; i < nTimeBins_; ++i) {

        ofs << time_[i];

        for (int j = 0; j < nZBins_; ++j) {          
          ofs << "\t" << histogram_[i][j];
        }
        ofs << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "RCorrFuncZ::writeCorrelate Error: fail to open %s\n",
              getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs.close();    
  }
}
