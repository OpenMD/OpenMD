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
#include "applications/staticProps/TwoDGofR.hpp"
#include "utils/simError.h"

namespace OpenMD {

  TwoDGofR::TwoDGofR(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, RealType len, RealType dz, int nrbins)
    : RadialDistrFunc(info, filename, sele1, sele2), len_(len), nRBins_(nrbins){

      deltaR_ = len_ /nRBins_;

      deltaZ_ = dz;
    
      histogram_.resize(nRBins_);
      avgTwoDGofR_.resize(nRBins_);

      setOutputName(getPrefix(filename) + ".TwoDGofR");
    }


  void TwoDGofR::preProcess() {
    std::fill(avgTwoDGofR_.begin(), avgTwoDGofR_.end(), 0.0);    
  }

  void TwoDGofR::initalizeHistogram() {
    std::fill(histogram_.begin(), histogram_.end(), 0);
  }


  void TwoDGofR::processHistogram() {

    int nPairs = getNPairs();

    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();

    RealType volume = hmat(0,0) * hmat(1,1) * deltaZ_;

    RealType pairDensity = nPairs /volume * 2.0;
    RealType pairConstant = (NumericConstant::PI * pairDensity);

    for(unsigned int i = 0 ; i < histogram_.size(); ++i){

      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice = deltaZ_ * (( rUpper*rUpper ) - ( rLower*rLower ));
      RealType nIdeal = volSlice * pairConstant;
      
      avgTwoDGofR_[i] += histogram_[i] / nIdeal;    
    }

  }

  void TwoDGofR::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    if (sd1 == sd2) {
      return;
    }
    
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = sqrt(r12.x()*r12.x() + r12.y()*r12.y());

    if (distance < len_) {
      int whichBin = distance / deltaR_;
      histogram_[whichBin] += 2;
    }
  }


  void TwoDGofR::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#2D radial distribution function\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      rdfStream << "#r\tcorrValue\n";
      for (unsigned int i = 0; i < avgTwoDGofR_.size(); ++i) {
	RealType r = deltaR_ * (i + 0.5);
	rdfStream << r << "\t" << avgTwoDGofR_[i]/nProcessed_ << "\n";
      }
        
    } else {

      sprintf(painCave.errMsg, "TwoDGofR: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();
  }

}

