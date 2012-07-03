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
 */

/* Calculates Rho(theta) */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/pAngle.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {
  
  pAngle::pAngle(SimInfo* info, const std::string& filename, 
                 const std::string& sele, int nthetabins)
    : StaticAnalyser(info, filename), selectionScript_(sele),  
      evaluator_(info), seleMan_(info), nThetaBins_(nthetabins){
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }            
    
    count_.resize(nThetaBins_);
    histogram_.resize(nThetaBins_);
    
    setOutputName(getPrefix(filename) + ".pAngle");
  }
  
  void pAngle::process() {
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    int i;

    Thermo thermo(info_);
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);
    
    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }
      }
      
      Vector3d CenterOfMass = thermo.getCom();      

      if  (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      int runningTot = 0;
      for (sd = seleMan_.beginSelected(i); sd != NULL; 
           sd = seleMan_.nextSelected(i)) {
        
        Vector3d pos = sd->getPos();
        
        Vector3d r1 = CenterOfMass - pos;
        // only do this if the stunt double actually has a vector associated
        // with it
        if (sd->isDirectional()) {
          Vector3d dipole = sd->getA().getColumn(2);
          RealType distance = r1.length();
          
          dipole.normalize();
          r1.normalize();
          RealType cosangle = dot(r1, dipole);
                          
          int binNo = int(nThetaBins_ * (1.0 + cosangle) / 2.0);
          count_[binNo]++;
        }
        
      }
    }
    processHistogram();
    writeProbs();
    
  }
  
  void pAngle::processHistogram() {
    
    int atot = 0;
    for(int i = 0; i < count_.size(); ++i) 
      atot += count_[i];
    
    for(int i = 0; i < count_.size(); ++i) {
      histogram_[i] = double(count_[i] / double(atot));
    }    
  }
  
  
  void pAngle::writeProbs() {

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#pAngle\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#cos(theta)\tp(cos(theta))\n";
      RealType dct = 2.0 / histogram_.size();
      for (int i = 0; i < histogram_.size(); ++i) {
        RealType ct = -1.0 + (2.0 * i + 1) / (histogram_.size());
        rdfStream << ct << "\t" << histogram_[i]/dct << "\n";
      }
      
    } else {
      
      sprintf(painCave.errMsg, "pAngle: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    rdfStream.close();
  }
  
}

