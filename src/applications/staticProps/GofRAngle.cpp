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
#include "applications/staticProps/GofRAngle.hpp"
#include "primitives/Atom.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofRAngle::GofRAngle(SimInfo* info, const std::string& filename, const std::string& sele1, 
		       const std::string& sele2, RealType len, int nrbins, int nangleBins)
    : RadialDistrFunc(info, filename, sele1, sele2), len_(len), nRBins_(nrbins), nAngleBins_(nangleBins){

      deltaR_ = len_ /(double) nRBins_;
      deltaCosAngle_ = 2.0 / (double)nAngleBins_;    
      histogram_.resize(nRBins_);
      avgGofr_.resize(nRBins_);
      for (int i = 0 ; i < nRBins_; ++i) {
        histogram_[i].resize(nAngleBins_);
        avgGofr_[i].resize(nAngleBins_);
      }
    }


  void GofRAngle::preProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      std::fill(avgGofr_[i].begin(), avgGofr_[i].end(), 0);
    }
  }

  void GofRAngle::initalizeHistogram() {
    npairs_ = 0;
    for (unsigned int i = 0; i < histogram_.size(); ++i){
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0);
    }
  }

  void GofRAngle::processHistogram() {
    int nPairs = getNPairs();
    RealType volume = info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = nPairs /volume;
    RealType pairConstant = ( 4.0 * NumericConstant::PI * pairDensity ) / 3.0;

    for(unsigned int i = 0 ; i < histogram_.size(); ++i){

      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice = ( rUpper * rUpper * rUpper ) - ( rLower * rLower * rLower );
      RealType nIdeal = volSlice * pairConstant;

      for (unsigned int j = 0; j < histogram_[i].size(); ++j){
	avgGofr_[i][j] += histogram_[i][j] / nIdeal;    
      }
    }

  }

  void GofRAngle::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    if (sd1 == sd2) {
      return;
    }
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();
    int whichRBin = distance / deltaR_;

    if (distance <= len_) {

      RealType cosAngle = evaluateAngle(sd1, sd2);
      RealType halfBin = (nAngleBins_ - 1) * 0.5;
      int whichThetaBin = halfBin * (cosAngle + 1.0);
      ++histogram_[whichRBin][whichThetaBin];
        
      ++npairs_;
    }
  }

  void GofRAngle::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#radial distribution function\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      rdfStream << "#nRBins = " << nRBins_ << "\t maxLen = " << len_ << "deltaR = " << deltaR_ <<"\n";
      rdfStream << "#nAngleBins =" << nAngleBins_ << "deltaCosAngle = " << deltaCosAngle_ << "\n";
      for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
	RealType r = deltaR_ * (i + 0.5);

	for(unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
	  RealType cosAngle = -1.0 + (j + 0.5)*deltaCosAngle_;
	  rdfStream << avgGofr_[i][j]/nProcessed_ << "\t";
	}

	rdfStream << "\n";
      }
        
    } else {
      sprintf(painCave.errMsg, "GofRAngle: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();
  }

  RealType GofRTheta::evaluateAngle(StuntDouble* sd1, StuntDouble* sd2) {
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
  
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    r12.normalize();

    AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();
    MultipoleAdapter ma1 = MultipoleAdapter(atype1);
    Vector3d vec;
    if (ma1.isDipole() )
      vec = sd1->getDipole();
    else 
      vec = sd1->getA().transpose() * V3Z;
    vec.normalize();    

    return dot(r12, vec);
  }

  RealType GofROmega::evaluateAngle(StuntDouble* sd1, StuntDouble* sd2) {

    AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();
    AtomType* atype2 = static_cast<Atom*>(sd2)->getAtomType();

    MultipoleAdapter ma1 = MultipoleAdapter(atype1);
    MultipoleAdapter ma2 = MultipoleAdapter(atype2);

    Vector3d v1, v2;

    if (ma1.isDipole() )
      v1 = sd1->getDipole();
    else 
      v1 = sd1->getA().transpose() * V3Z;

    if (ma2.isDipole() )
      v2 = sd2->getDipole();
    else 
      v2 = sd2->getA().transpose() * V3Z;

    v1.normalize();
    v2.normalize();
    return dot(v1, v2);
  }


}


