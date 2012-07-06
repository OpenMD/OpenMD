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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by Charles F. Vardeman II on 11/26/05.
 *  @author  Charles F. Vardeman II 
 *  @version $Id$
 *
 */

/* Calculates Rho(Z) for density profile of liquid slab. */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/RhoZ.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
namespace OpenMD {
  
  RhoZ::RhoZ(SimInfo* info, const std::string& filename, const std::string& sele, int nzbins)
    : StaticAnalyser(info, filename), selectionScript_(sele),  evaluator_(info), seleMan_(info), nZBins_(nzbins){

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }       
    
    // fixed number of bins

    sliceSDLists_.resize(nZBins_);
    density_.resize(nZBins_);
    
    setOutputName(getPrefix(filename) + ".RhoZ");
  }

  void RhoZ::process() {
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }
      }

      int i;    
      for (i=0; i < nZBins_; i++) {
        sliceSDLists_[i].clear();
      }

      RealType sliceVolume = currentSnapshot_->getVolume() /nZBins_;
      Mat3x3d hmat = currentSnapshot_->getHmat();
      zBox_.push_back(hmat(2,2));
      
      RealType halfBoxZ_ = hmat(2,2) / 2.0;      

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      //wrap the stuntdoubles into a cell      
      for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
        Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_)
          currentSnapshot_->wrapVector(pos);
        sd->setPos(pos);
      }
      
      //determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
        Vector3d pos = sd->getPos();
        // shift molecules by half a box to have bins start at 0
        int binNo = int(nZBins_ * (halfBoxZ_ + pos.z()) / hmat(2,2));
        sliceSDLists_[binNo].push_back(sd);
      }

      //loop over the slices to calculate the densities
      for (i = 0; i < nZBins_; i++) {
        RealType totalMass = 0;
        for (unsigned int k = 0; k < sliceSDLists_[i].size(); ++k) {
          totalMass += sliceSDLists_[i][k]->getMass();
        }
        density_[i] += totalMass/sliceVolume;
      }
    }
    
    writeDensity();

  }
  
  
  
  void RhoZ::writeDensity() {

    // compute average box length:
    std::vector<RealType>::iterator j;
    RealType zSum = 0.0;
    for (j = zBox_.begin(); j != zBox_.end(); ++j) {
      zSum += *j;       
    }
    RealType zAve = zSum / zBox_.size();

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#RhoZ\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#z\tdensity\n";
      for (unsigned int i = 0; i < density_.size(); ++i) {
        RealType z = zAve * (i+0.5)/density_.size();
        rdfStream << z << "\t" << 1.660535*density_[i]/nProcessed_ << "\n";
      }
      
    } else {
      
      sprintf(painCave.errMsg, "RhoZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    rdfStream.close();
  }
  
}

