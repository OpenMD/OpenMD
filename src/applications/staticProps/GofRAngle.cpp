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

#include <algorithm>
#include <fstream>
#include "applications/staticProps/GofRAngle.hpp"
#include "primitives/Atom.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  
  GofRAngle::GofRAngle(SimInfo* info, const std::string& filename, 
                       const std::string& sele1, 
		       const std::string& sele2, 
                       RealType len, int nrbins, int nangleBins)
    : RadialDistrFunc(info, filename, sele1, sele2), nAngleBins_(nangleBins), 
      len_(len), nRBins_(nrbins), 
      doSele3_(false), seleMan3_(info), evaluator3_(info) {
    
    deltaR_ = len_ /(double) nRBins_;
    deltaCosAngle_ = 2.0 / (double)nAngleBins_;    
    histogram_.resize(nRBins_);
    avgGofr_.resize(nRBins_);
    for (int i = 0 ; i < nRBins_; ++i) {
      histogram_[i].resize(nAngleBins_);
      avgGofr_[i].resize(nAngleBins_);
    }
  }
  
  GofRAngle::GofRAngle(SimInfo* info, const std::string& filename, 
                       const std::string& sele1, 
		       const std::string& sele2, 
                       const std::string& sele3,
                       RealType len, int nrbins, int nangleBins)
    : RadialDistrFunc(info, filename, sele1, sele2), nAngleBins_(nangleBins),
      len_(len), nRBins_(nrbins), doSele3_(true), selectionScript3_(sele3),
      seleMan3_(info), evaluator3_(info) {

    deltaR_ = len_ /(double) nRBins_;
    deltaCosAngle_ = 2.0 / (double)nAngleBins_;    
    histogram_.resize(nRBins_);
    avgGofr_.resize(nRBins_);
    for (int i = 0 ; i < nRBins_; ++i) {
      histogram_[i].resize(nAngleBins_);
      avgGofr_[i].resize(nAngleBins_);
    }

    evaluator3_.loadScriptString(sele3);      
    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }

  }
  
  void GofRAngle::processNonOverlapping( SelectionManager& sman1, 
                                         SelectionManager& sman2) {
    StuntDouble* sd1;
    StuntDouble* sd2;
    StuntDouble* sd3;
    int i;    
    int j;
    int k;
    
    // This is the same as a non-overlapping pairwise loop structure:
    // for (int i = 0;  i < ni ; ++i ) {
    //   for (int j = 0; j < nj; ++j) {} 
    // }

    if (doSele3_) {
      if  (evaluator3_.isDynamic()) {
        seleMan3_.setSelectionSet(evaluator3_.evaluate());
      }
      if (sman1.getSelectionCount() != seleMan3_.getSelectionCount() ) {
        RadialDistrFunc::processNonOverlapping( sman1, sman2 );
      }

      for (sd1 = sman1.beginSelected(i), sd3 = seleMan3_.beginSelected(k); 
           sd1 != NULL && sd3 != NULL; 
           sd1 = sman1.nextSelected(i), sd3 = seleMan3_.nextSelected(k)) {
        for (sd2 = sman2.beginSelected(j); sd2 != NULL; 
             sd2 = sman2.nextSelected(j)) {
          collectHistogram(sd1, sd2, sd3);
        }
      }
    } else {
      RadialDistrFunc::processNonOverlapping( sman1, sman2 );
    }
  }

  void GofRAngle::processOverlapping( SelectionManager& sman) {
    StuntDouble* sd1;
    StuntDouble* sd2;
    StuntDouble* sd3;
    int i;    
    int j;
    int k;

    // This is the same as a pairwise loop structure:
    // for (int i = 0;  i < n-1 ; ++i ) {
    //   for (int j = i + 1; j < n; ++j) {} 
    // }
    
    if (doSele3_) {
      if  (evaluator3_.isDynamic()) {
        seleMan3_.setSelectionSet(evaluator3_.evaluate());
      }
      if (sman.getSelectionCount() != seleMan3_.getSelectionCount() ) {
        RadialDistrFunc::processOverlapping( sman);
      }
      for (sd1 = sman.beginSelected(i), sd3 = seleMan3_.beginSelected(k); 
           sd1 != NULL && sd3 != NULL; 
           sd1 = sman.nextSelected(i), sd3 = seleMan3_.nextSelected(k)) {
        for (j  = i, sd2 = sman.nextSelected(j); sd2 != NULL; 
             sd2 = sman.nextSelected(j)) {
          collectHistogram(sd1, sd2, sd3);
        }            
      }
    } else {
      RadialDistrFunc::processOverlapping( sman);
    }    
  }
  

  void GofRAngle::preProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      std::fill(avgGofr_[i].begin(), avgGofr_[i].end(), 0);
    }
  }

  void GofRAngle::initializeHistogram() {
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
      RealType volSlice = ( rUpper * rUpper * rUpper ) - 
        ( rLower * rLower * rLower );
      RealType nIdeal = volSlice * pairConstant;

      for (unsigned int j = 0; j < histogram_[i].size(); ++j){
	avgGofr_[i][j] += histogram_[i][j] / nIdeal;    
      }
    }

  }
  
  void GofRTheta::processHistogram() {
    int nPairs = getNPairs();
    RealType volume = info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = nPairs /volume;
    RealType pairConstant = ( 4.0 * NumericConstant::PI * pairDensity ) /
      (3.0 * (double)nAngleBins_ );

    for(unsigned int i = 0 ; i < histogram_.size(); ++i){

      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice = ( rUpper * rUpper * rUpper ) - 
        ( rLower * rLower * rLower );
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
    int whichRBin = int(distance / deltaR_);

    if (distance <= len_) {

      RealType cosAngle = evaluateAngle(sd1, sd2);
      RealType halfBin = (nAngleBins_ - 1) * 0.5;
      int whichThetaBin = int(halfBin * (cosAngle + 1.0));
      ++histogram_[whichRBin][whichThetaBin];
        
      ++npairs_;
    }
  }

  void GofRAngle::collectHistogram(StuntDouble* sd1, StuntDouble* sd2, 
                                   StuntDouble* sd3) {

    if (sd1 == sd2) {
      return;
    }

    Vector3d p1 = sd1->getPos();
    Vector3d p3 = sd3->getPos();

    Vector3d c = 0.5 * (p1 + p3);
    Vector3d r13 = p3 - p1;

    Vector3d r12 = sd2->getPos() - c;
  
    if (usePeriodicBoundaryConditions_) {
      currentSnapshot_->wrapVector(r12);
      currentSnapshot_->wrapVector(r13);
    }
    
    RealType distance = r12.length();
    int whichRBin = int(distance / deltaR_);

    if (distance <= len_) {

      RealType cosAngle = evaluateAngle(sd1, sd2, sd3);
      RealType halfBin = (nAngleBins_ - 1) * 0.5;
      int whichThetaBin = int(halfBin * (cosAngle + 1.0));
      ++histogram_[whichRBin][whichThetaBin];
        
      ++npairs_;
    }
  }

  void GofRAngle::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#radial distribution function\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      rdfStream <<  "selection2: (" << selectionScript2_ << ")";
      if (doSele3_) {
        rdfStream << "\tselection3: (" << selectionScript3_ << ")\n";
      } else {
        rdfStream << "\n";
      }
      rdfStream << "#nRBins = " << nRBins_ << "\tmaxLen = " 
                << len_ << "\tdeltaR = " << deltaR_ <<"\n";
      rdfStream << "#nAngleBins =" << nAngleBins_ << "\tdeltaCosAngle = " 
                << deltaCosAngle_ << "\n";
      for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
	// RealType r = deltaR_ * (i + 0.5);

	for(unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
	  // RealType cosAngle = -1.0 + (j + 0.5)*deltaCosAngle_;
	  rdfStream << avgGofr_[i][j]/nProcessed_ << "\t";
	}

	rdfStream << "\n";
      }
        
    } else {
      sprintf(painCave.errMsg, "GofRAngle: unable to open %s\n", 
              outputFilename_.c_str());
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

    Vector3d vec;    
    
    if (!sd1->isDirectional()) {
      sprintf(painCave.errMsg, 
              "GofRTheta: attempted to use a non-directional object: %s\n", 
              sd1->getType().c_str());
      painCave.isFatal = 1;
      simError();  
    }

    if (sd1->isAtom()) {
      AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();
      MultipoleAdapter ma1 = MultipoleAdapter(atype1);
      
      if (ma1.isDipole() )
        vec = sd1->getDipole();
      else 
        vec = sd1->getA().transpose() * V3Z;
    } else {
      vec = sd1->getA().transpose() * V3Z;
    }

    vec.normalize();    
      
    return dot(r12, vec);
  }

  RealType GofRTheta::evaluateAngle(StuntDouble* sd1, StuntDouble* sd2, 
                                    StuntDouble* sd3) {
    Vector3d p1 = sd1->getPos();
    Vector3d p3 = sd3->getPos();

    Vector3d c = 0.5 * (p1 + p3);
    Vector3d r13 = p3 - p1;

    Vector3d r12 = sd2->getPos() - c;
  
    if (usePeriodicBoundaryConditions_) {
      currentSnapshot_->wrapVector(r12);
      currentSnapshot_->wrapVector(r13);
    }

    r12.normalize();
    r13.normalize();
          
    return dot(r12, r13);
  }

  RealType GofROmega::evaluateAngle(StuntDouble* sd1, StuntDouble* sd2) {
    Vector3d v1, v2;

    if (!sd1->isDirectional()) {
      sprintf(painCave.errMsg, 
              "GofROmega: attempted to use a non-directional object: %s\n", 
              sd1->getType().c_str());
      painCave.isFatal = 1;
      simError();  
    }

    if (sd1->isAtom()){
      AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();
      MultipoleAdapter ma1 = MultipoleAdapter(atype1);
      if (ma1.isDipole() )
        v1 = sd1->getDipole();
      else 
        v1 = sd1->getA().transpose() * V3Z;
    } else {
      v1 = sd1->getA().transpose() * V3Z;
    }
    
    if (!sd2->isDirectional()) {
      sprintf(painCave.errMsg, 
              "GofROmega attempted to use a non-directional object: %s\n", 
              sd2->getType().c_str());
      painCave.isFatal = 1;
      simError();  
    }

    if (sd2->isAtom()) {
      AtomType* atype2 = static_cast<Atom*>(sd2)->getAtomType();
      MultipoleAdapter ma2 = MultipoleAdapter(atype2);
           
      if (ma2.isDipole() )
        v2 = sd2->getDipole();
      else 
        v2 = sd2->getA().transpose() * V3Z;
    } else {
      v2 = sd2->getA().transpose() * V3Z;
    }
      
    v1.normalize();
    v2.normalize();
    return dot(v1, v2);
  }

  RealType GofROmega::evaluateAngle(StuntDouble* sd1, StuntDouble* sd2, 
                                    StuntDouble* sd3) {

    Vector3d v1;
    Vector3d v2;
    
    v1 = sd3->getPos() - sd1->getPos();
    if (usePeriodicBoundaryConditions_) 
      currentSnapshot_->wrapVector(v1);
    
    if (!sd2->isDirectional()) {
      sprintf(painCave.errMsg, 
              "GofROmega: attempted to use a non-directional object: %s\n", 
              sd2->getType().c_str());
      painCave.isFatal = 1;
      simError();  
    }

    if (sd2->isAtom()) {
      AtomType* atype2 = static_cast<Atom*>(sd2)->getAtomType();
      MultipoleAdapter ma2 = MultipoleAdapter(atype2);
           
      if (ma2.isDipole() )
        v2 = sd2->getDipole();
      else 
        v2 = sd2->getA().transpose() * V3Z;
    } else {
      v2 = sd2->getA().transpose() * V3Z;
    }
      
    v1.normalize();
    v2.normalize();
    return dot(v1, v2);
  }
}


