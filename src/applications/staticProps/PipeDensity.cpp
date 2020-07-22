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

/* Calculates PipeDensity, rho(axis2,axis3) in the box */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/PipeDensity.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {
  
  PipeDensity::PipeDensity(SimInfo* info, const std::string& filename, 
                           const std::string& sele, int nbins, int nbins2,
                           int axis)
    : StaticAnalyser(info, filename, nbins2), selectionScript_(sele), 
      evaluator_(info), seleMan_(info), nBins2_(nbins), axis_(axis) {
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }       
    
    // fixed number of bins

    sliceSDLists_.resize(nBins2_);
    density_.resize(nBins2_);
    for (unsigned int i = 0 ; i < nBins2_; ++i) {
      sliceSDLists_[i].resize(nBins_);
      density_[i].resize(nBins_);
    }

    // Compute complementary axes to the privileged axis
    axis1_ = (axis_ + 1) % 3;
    axis2_ = (axis_ + 2) % 3;
    
    // Set the axis labels for the non-privileged axes
    switch(axis_) {
    case 0:
      axisLabel1_ = "y";
      axisLabel2_ = "z";
      break;
    case 1:
      axisLabel1_ = "z";
      axisLabel2_ = "x";            
      break;
    case 2:
    default:
      axisLabel1_ = "x";
      axisLabel2_ = "y";
      break;
    }

    setOutputName(getPrefix(filename) + ".PipeDensity");
  }

  void PipeDensity::process() {
    StuntDouble* sd;
    int ii;

    bool usePeriodicBoundaryConditions_ = 
      info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (unsigned int i = 0; i < nBins2_; i++) {
        for (unsigned int j = 0; j < nBins_; j++) {          
          sliceSDLists_[i][j].clear();
        }
      }

      RealType sliceVolume = currentSnapshot_->getVolume() /(nBins2_ * nBins_);
      Mat3x3d hmat = currentSnapshot_->getHmat();

      RealType halfBox1_ = hmat(axis1_,axis1_) / 2.0;      
      RealType halfBox2_ = hmat(axis2_,axis2_) / 2.0;      

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      //wrap the stuntdoubles into a cell      
      for (sd = seleMan_.beginSelected(ii); sd != NULL; 
	   sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_)
          currentSnapshot_->wrapVector(pos);
        sd->setPos(pos);
      }
      
      //determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(ii); sd != NULL; 
	   sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        // shift molecules by half a box to have bins start at 0
        int binNo1 = int(nBins2_ * (halfBox1_ + pos[axis1_]) /
                         hmat(axis1_,axis1_));
        int binNo2 = int(nBins_  * (halfBox2_ + pos[axis2_]) /
                         hmat(axis2_,axis2_));
        sliceSDLists_[binNo1][binNo2].push_back(sd);
      }

      //loop over the slices to calculate the densities
      for (unsigned int i = 0; i < nBins2_; i++) {
        for (unsigned int j = 0; j < nBins_; j++) {

          RealType totalMass = 0;
          for (unsigned int k = 0; k < sliceSDLists_[i][j].size(); ++k) {
            totalMass += sliceSDLists_[i][j][k]->getMass();
          }
          density_[i][j] += totalMass/sliceVolume;
        }
      }
    }
      
    writeDensity();

  }
    
  void PipeDensity::writeDensity() {

    std::vector<RealType>::iterator j;
    std::ofstream rdfStream(outputFilename_.c_str());

    if (rdfStream.is_open()) {
      rdfStream << "#PipeDensity\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#density (" << axisLabel1_ << "," << axisLabel2_ << ")\n";
      for (unsigned int i = 0; i < density_.size(); ++i) {
        for (unsigned int j = 0; j < density_[i].size(); ++j) {          
          rdfStream << Constants::densityConvert * density_[i][j] / nProcessed_;
          rdfStream << "\t";
        }
        rdfStream << "\n";
      }
      
    } else {
      
      sprintf(painCave.errMsg, "PipeDensity: unable to open %s\n", 
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    rdfStream.close();
  }
  
}

