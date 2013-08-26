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

#include "applications/dynamicProps/cOHz.hpp"
#include "math/LegendrePolynomial.hpp"
#include "utils/simError.h"

namespace OpenMD {
  COHZ::COHZ(SimInfo* info, 
             const std::string& filename, 
             const std::string& sele1, 
             const std::string& sele2, 
             int order, int nZbins,
             long long int memSize)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2, 
                           DataStorage::dslAmat, memSize), nZBins_(nZbins) {

      setCorrFuncType("Legendre Correlation Function for OH bond vector of Z");
      setOutputName1(getPrefix(dumpFilename_) + ".cohZ");
      setOutputName2(getPrefix(dumpFilename_) + ".lcorrZ");
      histogram_.resize(nTimeBins_);
      counts_.resize(nTimeBins_);
      for (int i = 0; i < nTimeBins_; i++) {
        histogram_[i].resize(nZBins_);
        counts_[i].resize(nZBins_);
      }
      LegendrePolynomial polynomial(order);
      legendre_ = polynomial.getLegendrePolynomial(order);
    }

  void COHZ::correlateFrames(int frame1, int frame2) {
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
    StuntDouble* sd1;
    StuntDouble* sd2;

    for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan2_.beginSelected(j);
         sd1 != NULL && sd2 != NULL;
         sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j)) {

      Vector3d pos = sd1->getPos();
      if (info_->getSimParams()->getUsePeriodicBoundaryConditions())
        snapshot1->wrapVector(pos);
      int zBin = int(nZBins_ * (halfBoxZ_ + pos.z()) / hmat(2,2));

      Vector3d corrVals = calcCorrVals(frame1, frame2, sd1, sd2);
      histogram_[timeBin][zBin] += corrVals; 
      counts_[timeBin][zBin]++;
    }
    
  }

  void COHZ::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      for (int j = 0; j < nZBins_; ++j) {
        if (counts_[i][j] > 0) {
          histogram_[i][j] /= counts_[i][j];
        }
      }
    }
  }

  void COHZ::preCorrelate() {
    for (int i = 0; i < nTimeBins_; i++) {
      std::fill(histogram_[i].begin(), histogram_[i].end(), Vector3d(0.0));
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
  }

  Vector3d COHZ::calcCorrVals(int frame1, int frame2, StuntDouble* sd1,  StuntDouble* sd2) {
    
    // Vectors v1x, v1y, and v1z are the body-fixed axes on the
    // molecule in frame 1 in the laboratory frame.

    // Vectors v2x, v2y, and v2z are the body-fixed axes on the
    // molecule in frame 2 in the laboratory frame.

    // Vectors u1 & u2 are the first OH bond vector in frames 1 & 2
    // respectively.  Here we assume SPC/E geometry.

    // Vectors w1 & w2 are the second OH bond vector in frames 1 & 2
    // respectively.  Here we assume SPC/E geometry.

    Vector3d v1x = sd1->getA(frame1).getRow(0);
    Vector3d v2x = sd2->getA(frame2).getRow(0);

    Vector3d v1y = sd1->getA(frame1).getRow(1);
    Vector3d v2y = sd2->getA(frame2).getRow(1);

    Vector3d v1z = sd1->getA(frame1).getRow(2);
    Vector3d v2z = sd2->getA(frame2).getRow(2);

    Vector3d u1 = 0.81649 * v1y + 0.57736 * v1z;
    Vector3d u2 = 0.81649 * v2y + 0.57736 * v2z;

    Vector3d w1 = -0.81649 * v1y + 0.57736 * v1z;
    Vector3d w2 = -0.81649 * v2y + 0.57736 * v2z;

    
    RealType vprod = legendre_.evaluate(dot(v1z, v2z)/(v1z.length()*v2z.length()));
    RealType uprod = legendre_.evaluate(dot(u1, u2)/(u1.length()*u2.length()));
    RealType wprod = legendre_.evaluate(dot(w1, w2)/(w1.length()*w2.length()));

    return Vector3d(vprod, uprod, wprod);

  }


  void COHZ::validateSelection(const SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;    
    for (sd = seleMan1_.beginSelected(i); sd != NULL; sd = seleMan1_.nextSelected(i)) {
      if (!sd->isDirectionalAtom()) {
	sprintf(painCave.errMsg,
                "LegendreCorrFunc::validateSelection Error: selected atoms are not Directional\n");
	painCave.isFatal = 1;
	simError();        
      }
    }
    
  }

  void COHZ::writeCorrelate() {
    std::ofstream ofs1(getOutputFileName1().c_str());
    std::ofstream ofs2(getOutputFileName2().c_str());

    if (ofs1.is_open()) {

      ofs1 << "#" << getCorrFuncType() << "\n";
      ofs1 << "#time\tPn(costheta_z)\n";

      for (int i = 0; i < nTimeBins_; ++i) {

        ofs1 << time_[i];

        for (int j = 0; j < nZBins_; ++j) {          
          ofs1 << "\t" << 0.5*(histogram_[i][j](1) +  histogram_[i][j](2));
        }
        ofs1 << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "cOHz::writeCorrelate Error: fail to open %s\n", getOutputFileName1().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs1.close();    

    if (ofs2.is_open()) {

      ofs2 << "#" << getCorrFuncType() << "\n";
      ofs2 << "#time\tPn(costheta_z)\n";

      for (int i = 0; i < nTimeBins_; ++i) {

        ofs2 << time_[i];

        for (int j = 0; j < nZBins_; ++j) {          
          ofs2 << "\t" << histogram_[i][j](0);
        }
        ofs2 << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "cOHz::writeCorrelate Error: fail to open %s\n", getOutputFileName2().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs2.close();    
  }
}
