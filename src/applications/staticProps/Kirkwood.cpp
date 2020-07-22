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
#include <sstream>
#include "applications/staticProps/Kirkwood.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include "types/MultipoleAdapter.hpp"

namespace OpenMD {

  Kirkwood::Kirkwood(SimInfo* info, const std::string& filename,
                     const std::string& sele1, const std::string& sele2,
                     RealType len, int nrbins)
    : RadialDistrFunc(info, filename, sele1, sele2, nrbins), len_(len) {
    
    setAnalysisType("Distance-dependent Kirkwood G-factor");
    setOutputName(getPrefix(filename) + ".kirkwood");
    
    deltaR_ = len_ /nBins_;
    
    histogram_.resize(nBins_);
    avgKirkwood_.resize(nBins_);
    std::stringstream params;
    params << " len = " << len_
           << ", nrbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString( paramString );
  }
  

  void Kirkwood::preProcess() {
    std::fill(avgKirkwood_.begin(), avgKirkwood_.end(), 0.0);    
  }

  void Kirkwood::initializeHistogram() {
    std::fill(histogram_.begin(), histogram_.end(), 0);
  }


  void Kirkwood::processHistogram() {
    int nSelected1 = seleMan1_.getSelectionCount();
    for(unsigned int i = 0 ; i < histogram_.size(); ++i){
      avgKirkwood_[i] += histogram_[i] / nSelected1;    
    }
  }

  void Kirkwood::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    if (sd1 == sd2) {
      return;
    }
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();

    AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();          
    AtomType* atype2 = static_cast<Atom*>(sd2)->getAtomType();          

    MultipoleAdapter ma1 = MultipoleAdapter(atype1);
    MultipoleAdapter ma2 = MultipoleAdapter(atype2);

    Vector3d d1(0.0);
    Vector3d d2(0.0);
    RealType dotProduct(0.0);

    if (ma1.isDipole()) {
      d1 = sd1->getDipole();
      d1.normalize();
      if (ma2.isDipole()) {
        d2 = sd2->getDipole();
        d2.normalize();
        dotProduct = dot(d1,d2);
      }
    }
    
    if (distance < len_) {
      int whichBin = int(distance / deltaR_);
      // each dipole pair contributes to all of the radii that contain it.
      for (unsigned int i = whichBin; i < nBins_; i++) {
        histogram_[i] += dotProduct;
      }
    }
  }


  void Kirkwood::writeRdf() {
    std::ofstream ofs(outputFilename_.c_str());
    if (ofs.is_open()) {
      Revision r;
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_ ;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#r\tcorrValue\n";
      for (unsigned int i = 0; i < avgKirkwood_.size(); ++i) {
	RealType r = deltaR_ * (i + 0.5);
	ofs << r << "\t" << avgKirkwood_[i]/nProcessed_ << "\n";
      }
    } else {
      
      sprintf(painCave.errMsg, "Kirkwood: unable to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    ofs.close();
  }

  KirkwoodQuadrupoles::KirkwoodQuadrupoles(SimInfo* info,
                                           const std::string& filename,
                                           const std::string& sele1,
                                           const std::string& sele2,
                                           RealType len, int nrbins)
    : Kirkwood(info, filename, sele1, sele2, len, nrbins) {
    
    setAnalysisType("Distance-dependent Kirkwood G-factor for quadrupoles");
    setOutputName(getPrefix(filename) + ".kirkwoodQ");
  }
  
  void KirkwoodQuadrupoles::collectHistogram(StuntDouble* sd1,
                                             StuntDouble* sd2) {
    
    if (sd1 == sd2) {
      return;
    }
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    if (usePeriodicBoundaryConditions_)
      currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();

    AtomType* atype1 = static_cast<Atom*>(sd1)->getAtomType();          
    AtomType* atype2 = static_cast<Atom*>(sd2)->getAtomType();          

    MultipoleAdapter ma1 = MultipoleAdapter(atype1);
    MultipoleAdapter ma2 = MultipoleAdapter(atype2);

    Mat3x3d Q1(0.0);
    Mat3x3d Q2(0.0);
    RealType trQ1(0.0);
    RealType trQ2(0.0);
    RealType Q1dQ1(0.0);
    RealType Q2dQ2(0.0);
    
    RealType quadrupoleProduct(0.0);

    // Similar to the dipole case, but the effective quadrupole moment
    // is defined (in our electrostatics work) as:
    //    sqrt (3 Q:Q - Tr(Q)^2 )
    // so normalization is a bit different.  Here : denotes a
    // contraction (double dot product) of the quadrupole tensor.
    
    if (ma1.isQuadrupole()) {
      Q1 = sd1->getQuadrupole();      
      trQ1 = Q1.trace();
      Q1dQ1 = doubleDot(Q1, Q1);
      Q1 /= sqrt(3.0 * Q1dQ1 - trQ1*trQ1);
      // recompute the trace after normalizing:
      trQ1 = Q1.trace();
       
      if (ma2.isQuadrupole()) {
        Q2 = sd2->getQuadrupole();      
        trQ2 = Q2.trace();
        Q2dQ2 = doubleDot(Q2, Q2);
        Q2 /= sqrt(3.0 * Q2dQ2 - trQ2*trQ2);
        // recompute the trace after normalizing:
        trQ2 = Q2.trace();
        
        quadrupoleProduct = 3.0 * doubleDot(Q1, Q2) - trQ1 * trQ2; 
      }
    }
    
    if (distance < len_) {
      int whichBin = int(distance / deltaR_);
      // each dipole pair contributes to all of the radii that contain it.
      for (unsigned int i = whichBin; i < nBins_; i++) {
        histogram_[i] += quadrupoleProduct;
      }
    }
  }
}

