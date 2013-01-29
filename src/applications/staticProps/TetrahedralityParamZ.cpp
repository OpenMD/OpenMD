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
 * [6]  Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 */
 
#include "applications/staticProps/TetrahedralityParamZ.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
namespace OpenMD {
  TetrahedralityParamZ::TetrahedralityParamZ(SimInfo* info,  
                                             const std::string& filename, 
                                             const std::string& sele,
                                             double rCut, int nzbins) 
    : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), 
      seleMan_(info), nZBins_(nzbins) {
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // Set up cutoff radius:    
    rCut_ = rCut;

    // fixed number of bins
    sliceQ_.resize(nZBins_);
    sliceCount_.resize(nZBins_);    
    std::fill(sliceQ_.begin(), sliceQ_.end(), 0.0);
    std::fill(sliceCount_.begin(), sliceCount_.end(), 0);
    
    setOutputName(getPrefix(filename) + ".Qz");
  }
  
  TetrahedralityParamZ::~TetrahedralityParamZ() {
    sliceQ_.clear();
    sliceCount_.clear();
    zBox_.clear();
  }
    
  void TetrahedralityParamZ::process() {
    Molecule* mol;
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sdi;
    StuntDouble* sdj;
    RigidBody* rb;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule::RigidBodyIterator rbIter;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType,StuntDouble*> > myNeighbors;
    int isd;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      Mat3x3d hmat = currentSnapshot_->getHmat();
      zBox_.push_back(hmat(2,2));
      
      RealType halfBoxZ_ = hmat(2,2) / 2.0;      

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      // update the positions of atoms which belong to the rigidbodies
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL;
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }
      }
      
      // outer loop is over the selected StuntDoubles:
      for (sd = seleMan_.beginSelected(isd); sd != NULL;
           sd = seleMan_.nextSelected(isd)) {
        
        myIndex = sd->getGlobalIndex();

        Qk = 1.0;	  
        myNeighbors.clear();
        
        // inner loop is over all StuntDoubles in the system:
        
        for (mol = info_->beginMolecule(mi); mol != NULL; 
             mol = info_->nextMolecule(mi)) {

          for (sd2 = mol->beginIntegrableObject(ioi); sd2 != NULL; 
               sd2 = mol->nextIntegrableObject(ioi)) {
            
            if (sd2->getGlobalIndex() != myIndex) {
              
              vec = sd->getPos() - sd2->getPos();       
              
              if (usePeriodicBoundaryConditions_) 
                currentSnapshot_->wrapVector(vec);
              
              r = vec.length();             

              // Check to see if neighbor is in bond cutoff 
              
              if (r < rCut_) {                
                myNeighbors.push_back(std::make_pair(r,sd2));
              }
            }
          }
        }
        
        // Sort the vector using predicate and std::sort
        std::sort(myNeighbors.begin(), myNeighbors.end());
        
        // Use only the 4 closest neighbors to do the rest of the work:
        
        int nbors =  myNeighbors.size()> 4 ? 4 : myNeighbors.size();
        int nang = int (0.5 * (nbors * (nbors - 1)));
        
        rk = sd->getPos();

        for (int i = 0; i < nbors-1; i++) {       

          sdi = myNeighbors[i].second;
          ri = sdi->getPos();
          rik = rk - ri;
          if (usePeriodicBoundaryConditions_) 
            currentSnapshot_->wrapVector(rik);
          
          rik.normalize();

          for (int j = i+1; j < nbors; j++) {       

            sdj = myNeighbors[j].second;
            rj = sdj->getPos();
            rkj = rk - rj;
            if (usePeriodicBoundaryConditions_) 
              currentSnapshot_->wrapVector(rkj);
            rkj.normalize();
            
            cospsi = dot(rik,rkj);           

            // Calculates scaled Qk for each molecule using calculated
            // angles from 4 or fewer nearest neighbors.
            Qk -=  (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);            
          }
        }

        if (nang > 0) {
          if (usePeriodicBoundaryConditions_)
            currentSnapshot_->wrapVector(rk);
          
          int binNo = int(nZBins_ * (halfBoxZ_ + rk.z()) / hmat(2,2));
          sliceQ_[binNo] += Qk;
          sliceCount_[binNo] += 1;
        }  
      }
    }
    writeQz();
  }
  
  void TetrahedralityParamZ::writeQz() {

    // compute average box length:
    
    RealType zSum = 0.0;
    for (std::vector<RealType>::iterator j = zBox_.begin(); 
         j != zBox_.end(); ++j) {
      zSum += *j;       
    }
    RealType zAve = zSum / zBox_.size();

    std::ofstream qZstream(outputFilename_.c_str());
    if (qZstream.is_open()) {
      qZstream << "#Tetrahedrality Parameters (z)\n";
      qZstream << "#nFrames:\t" << zBox_.size() << "\n";
      qZstream << "#selection: (" << selectionScript_ << ")\n";
      qZstream << "#z\tQk\n";
      for (unsigned int i = 0; i < sliceQ_.size(); ++i) {
        RealType z = zAve * (i+0.5) / sliceQ_.size();
        qZstream << z << "\t" << sliceQ_[i] / sliceCount_[i] << "\n";
      }
      
    } else {      
      sprintf(painCave.errMsg, "TetrahedralityParamZ: unable to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }    
    qZstream.close();
  }
}



