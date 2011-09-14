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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#include "applications/dynamicProps/LegendreCorrFunc.hpp"
#include "math/LegendrePolynomial.hpp"
#include "utils/simError.h"

namespace OpenMD {
  LegendreCorrFunc::LegendreCorrFunc(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, int order, long long int memSize)
    : ParticleTimeCorrFunc(info, filename, sele1, sele2, DataStorage::dslAmat | DataStorage::dslElectroFrame, memSize){

      setCorrFuncType("Legendre Correlation Function");
      setOutputName(getPrefix(dumpFilename_) + ".lcorr");
      histogram_.resize(nTimeBins_); 
      count_.resize(nTimeBins_);
      nSelected_ =   seleMan1_.getSelectionCount();  
      std::cout << "sele1 = " << sele1 << "\n";
      std::cout << "sele2 = " << sele2 << "\n";
      std::cout << "nSelected = " << nSelected_ << "\n";
      assert(  nSelected_ == seleMan2_.getSelectionCount());

      LegendrePolynomial polynomial(order);
      legendre_ = polynomial.getLegendrePolynomial(order);
    }

  void LegendreCorrFunc::correlateFrames(int frame1, int frame2) {
    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    RealType time1 = snapshot1->getTime();
    RealType time2 = snapshot2->getTime();

    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);
    count_[timeBin] += nSelected_;    

    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd2;

    for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan2_.beginSelected(j);
         sd1 != NULL && sd2 != NULL;
         sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j)) {

      Vector3d corrVals = calcCorrVals(frame1, frame2, sd1, sd2);
      histogram_[timeBin] += corrVals; 
    }
    
  }

  void LegendreCorrFunc::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= count_[i];
      }
    }
  }

  void LegendreCorrFunc::preCorrelate() {
    // Fill the histogram with empty Vector3d:
    std::fill(histogram_.begin(), histogram_.end(), Vector3d(0.0));
    // count array set to zero
    std::fill(count_.begin(), count_.end(), 0);
  }



  Vector3d LegendreCorrFunc::calcCorrVals(int frame1, int frame2, StuntDouble* sd1,  StuntDouble* sd2) {
    
    // The lab frame vector corresponding to the body-fixed 
    // z-axis is simply the second column of A.transpose()
    // or, identically, the second row of A itself.
    // Similar identites give the 0th and 1st rows of A for
    // the lab vectors associated with body-fixed x- and y- axes.

    Vector3d v1x = sd1->getA(frame1).getRow(0);
    Vector3d v2x = sd2->getA(frame2).getRow(0);
    Vector3d v1y = sd1->getA(frame1).getRow(1);
    Vector3d v2y = sd2->getA(frame2).getRow(1);
    Vector3d v1z = sd1->getA(frame1).getRow(2);
    Vector3d v2z = sd2->getA(frame2).getRow(2);

    RealType uxprod = legendre_.evaluate(dot(v1x, v2x)/(v1x.length()*v2x.length()));
    RealType uyprod = legendre_.evaluate(dot(v1y, v2y)/(v1y.length()*v2y.length()));
    RealType uzprod = legendre_.evaluate(dot(v1z, v2z)/(v1z.length()*v2z.length()));

    return Vector3d(uxprod, uyprod, uzprod);

  }


  void LegendreCorrFunc::validateSelection(const SelectionManager& seleMan) {
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

  void LegendreCorrFunc::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#time\tPn(costheta_x)\tPn(costheta_y)\tPn(costheta_z)\n";

      for (int i = 0; i < nTimeBins_; ++i) {
        ofs << time_[i] << "\t" << 
          histogram_[i](0) << "\t" <<
          histogram_[i](1) << "\t" <<
          histogram_[i](2) << "\t" << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "LegendreCorrFunc::writeCorrelate Error: fail to open %s\n", getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs.close();    
  }
}
