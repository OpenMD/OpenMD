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
#include "applications/staticProps/GofZ.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofZ::GofZ(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, RealType len, int nrbins)
    : RadialDistrFunc(info, filename, sele1, sele2), len_(len), nRBins_(nrbins){

      deltaZ_ = len_ /nRBins_;
      rC_ = len_ / 2.0;
    
      histogram_.resize(nRBins_);
      avgGofz_.resize(nRBins_);

      setOutputName(getPrefix(filename) + ".gofz");
    }


  void GofZ::preProcess() {
    std::fill(avgGofz_.begin(), avgGofz_.end(), 0.0);    
  }

  void GofZ::initalizeHistogram() {
    std::fill(histogram_.begin(), histogram_.end(), 0);
  }


  void GofZ::processHistogram() {

    int nPairs = getNPairs();
    RealType volume = info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = nPairs /volume * 2.0;
    RealType pairConstant = ( 2.0 * NumericConstant::PI * pairDensity );

    for(unsigned int i = 0 ; i < histogram_.size(); ++i){

      RealType zLower = i * deltaZ_;
      RealType zUpper = zLower + deltaZ_;
      RealType cylSlice = (zUpper - zLower)*rC_*rC_;
      RealType nIdeal = cylSlice * pairConstant;

      avgGofz_[i] += histogram_[i] / nIdeal;    
    }

  }

  void GofZ::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    RealType thisZ;

    if (sd1 == sd2) {
      return;
    }
    
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();
    RealType z2 = r12.z()*r12.z();
    RealType xydist = sqrt(distance*distance - z2);

    if (xydist < rC_) {
      thisZ = abs(r12.z());
      int whichBin = thisZ / deltaZ_;
      histogram_[whichBin] += 2;
    }
  }


  void GofZ::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#z-separation distribution function\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      rdfStream << "#r\tcorrValue\n";
      for (unsigned int i = 0; i < avgGofz_.size(); ++i) {
	RealType z = deltaZ_ * (i + 0.5);
	rdfStream << z << "\t" << avgGofz_[i]/nProcessed_ << "\n";
      }
        
    } else {

      sprintf(painCave.errMsg, "GofZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();
  }

}

