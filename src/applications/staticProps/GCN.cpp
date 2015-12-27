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
#include "applications/staticProps/GCN.hpp"
#include "io/DumpReader.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GCN::GCN(SimInfo* info, const std::string& filename,
           const std::string& sele1, const std::string& sele2,
           RealType rCut, int bins):
    StaticAnalyser(info, filename),
    sele1_(sele1), seleMan1_(info), evaluator1_(info),
    sele2_(sele2), seleMan2_(info), evaluator2_(info) {
    
    setOutputName(getPrefix(filename) + ".gcn");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      selectionCount1_ = seleMan1_.getSelectionCount();
    }
    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
      selectionCount2_ = seleMan2_.getSelectionCount();
    }

    bins_ = bins; //Treat bins as division, bins per single integer

    setOutputName(getPrefix(filename) + ".gcn");

    rCut_ = rCut;
    nnMax_ = 12;
    // Overestimate the histogram size (still protected in process
    // function below):
    hBins_ = bins_ * (nnMax_ * 1.5);
  }

  GCN::~GCN() {
    histogram_.clear();
  }

  void GCN::process() {
    SelectionManager common(info_);
    
    std::vector<std::vector<int> > listNN;
    std::vector<int> globalToLocal;

    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd1;
    StuntDouble* sd2;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    
    int iterator1;
    int iterator2;
    unsigned int mapIndex1(0);
    unsigned int mapIndex2(0);
    unsigned int tempIndex(0);
    unsigned int binIndex(0);
    RealType gcn(0.0);
    RealType binValue(0.0);
    Vector3d pos1;
    Vector3d pos2;
    Vector3d diff;
    RealType distance;

    std::ofstream os(getOutputFileName().c_str());
    os << "#Generalized Coordinate Number\n";
    
    //First have to calculate lists of nearest neighbors (listNN_):
    
    for(int istep = 0; istep < nFrames; istep += step_){
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      for(mol = info_->beginMolecule(mi); mol != NULL;
          mol = info_->nextMolecule(mi)){
	for(rb = mol->beginRigidBody(rbIter); rb != NULL;
            rb = mol->nextRigidBody(rbIter)){
	  rb->updateAtoms();
	}
      }

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
	selectionCount1_ = seleMan1_.getSelectionCount();
      }
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
	selectionCount2_ = seleMan2_.getSelectionCount();
      }

      // We need a common selection set:
      common = seleMan1_ | seleMan2_;
      int commonCount = common.getSelectionCount();

      globalToLocal.clear();
      globalToLocal.resize(info_->getNGlobalAtoms(),-1);
      for (unsigned int i = 0; i < listNN.size(); i++)         
        listNN.at(i).clear();
      listNN.clear();
      listNN.resize(commonCount);
      histogram_.clear();
      histogram_.resize(hBins_, 0.0);

      mapIndex1 = 0;
      for(sd1 = common.beginSelected(iterator1); sd1 != NULL;
          sd1 = common.nextSelected(iterator1)) {
        
	globalToLocal.at(sd1->getGlobalIndex()) = mapIndex1;

        pos1 = sd1->getPos();

	mapIndex2 = 0;
 	for(sd2 = common.beginSelected(iterator2); sd2 != NULL;
            sd2 = common.nextSelected(iterator2)) {
          
	  if (mapIndex1 < mapIndex2) {
            pos2 = sd2->getPos();
            diff = pos2 - pos1;
	    if (usePeriodicBoundaryConditions_) {
	      currentSnapshot_->wrapVector(diff);
	    }
            distance = diff.length();
	    if (distance < rCut_) {
              listNN.at(mapIndex1).push_back(mapIndex2);
	      listNN.at(mapIndex2).push_back(mapIndex1);
	    }
	  }
	  mapIndex2++;
	}
	mapIndex1++;
      }
      
      // Fill up the histogram with gcn values
      for(sd1 = seleMan1_.beginSelected(iterator1); sd1 != NULL;
          sd1 = seleMan1_.nextSelected(iterator1)){
            
	mapIndex1 = globalToLocal.at(sd1->getGlobalIndex());
	gcn = 0.0;
	for(unsigned int i = 0; i < listNN.at(mapIndex1).size(); i++){
          // tempIndex is the index of one of i's nearest neighbors
	  tempIndex = listNN.at(mapIndex1).at(i);
	  gcn += listNN.at(tempIndex).size();
	}

        gcn = gcn / nnMax_;
        binIndex = int(gcn*bins_);
        if (binIndex < histogram_.size()) {
          histogram_.at(binIndex) += 1;
        } else {
          cerr << "In frame " <<  istep <<  ", object "
               << sd1->getGlobalIndex() << " has GCN value = " << gcn << "\n";
        }
      }
      
      os << "#Selection Count: " << selectionCount1_ << "\n";
      os << "#Frame " << istep << "\n";
      for(unsigned int n = 0; n < histogram_.size(); n++){
	binValue = n/(1.0*bins_);
	os << binValue << "\t"
           << (histogram_.at(n)/(selectionCount1_*1.0))
           << "\n";
      } 
    }   
    os.close();
  }
}
