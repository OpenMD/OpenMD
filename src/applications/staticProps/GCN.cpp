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

  GCN::GCN(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, int bins):
    StaticAnalyser(info, filename),
    sele1_(sele1), seleMan1_(info), evaluator1_(info),
    sele2_(sele2), seleMan2_(info), evaluator2_(info) {

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


      filename_ = filename;
      bins_ = bins; //Treat bins as division, bins per single integer

      setOutputName(getPrefix(filename) + ".gcn");

      nnMax_ = 12.0;
      solShell_ = 3.3;
      cout << "Bins per integer: " << bins_ << "\n";
      hBins_ = bins_*(nnMax_+2); //temporary so 0.0 to 12.0 by 0.1
      cout << "Bins: " << hBins_ << "\n";
    }


  GCN::~GCN(){


  }


  void GCN::process() {
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd1;
    StuntDouble* sd2;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    std::ofstream gcnStream;
    setOutputName(getPrefix(filename_) + ".gcn");
    gcnStream.open(outputFilename_.c_str());
    
    DumpReader reader(info_, dumpFilename_);
    
    int nFrames = reader.getNFrames();
    int iterator1;
    int iterator2;
    int mapIndex1 = 0;
    int mapIndex2 = 0;
    int tempIndex = 0;
    int binIndex = 0;
    RealType gcn = 0.0;
    RealType binValue = 0.0;

    gcnStream << "#Generalized Coordinate Number\n";
    //First have to calculate nearestNeighbors and listNN
    
    for(int istep = 0; istep < nFrames; istep += step_){
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      for(mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)){
	for(rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)){
	  rb->updateAtoms();
	}
      }

      if  (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
	selectionCount1_ = seleMan1_.getSelectionCount();
      }
      if  (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
	selectionCount2_ = seleMan2_.getSelectionCount();
      }

      globalToLocal_.clear();
      globalToLocal_.resize(info_->getNGlobalAtoms(),-1);
      listNN_.clear();
      listNN_.resize(selectionCount2_);
      nearestNeighbors_.clear();
      nearestNeighbors_.resize(selectionCount2_, 0);
      histogram_.clear();
      histogram_.resize(hBins_, 0.0);
      
      mapIndex1 = 0;
      for(sd1 = seleMan2_.beginSelected(iterator1); sd1 != NULL; sd1 = seleMan2_.nextSelected(iterator1)){
	globalToLocal_[sd1->getGlobalIndex()] = mapIndex1;

	mapIndex2 = 0;
 	for(sd2 = seleMan2_.beginSelected(iterator2); sd2 != NULL; sd2 = seleMan2_.nextSelected(iterator2)){
	  if(mapIndex1 >= mapIndex2){
	  
	  }else{
	    Vector3d pos1 = sd1->getPos();
	    Vector3d pos2 = sd2->getPos();
	    Vector3d diff = pos2-pos1;
	    if(usePeriodicBoundaryConditions_){
	      currentSnapshot_->wrapVector(diff);
	    }
	    double distance = diff.length();
	    if(distance < solShell_){
	      listNN_[mapIndex2].push_back(mapIndex1);
	      listNN_[mapIndex1].push_back(mapIndex2);
	      nearestNeighbors_[mapIndex1] += 1;
	      nearestNeighbors_[mapIndex2] += 1;
	    }
	  }
	  mapIndex2++;
	}
	mapIndex1++;
      }	
      // Fill up the histogram with gcn values
      for(sd1 = seleMan1_.beginSelected(iterator1); sd1 != NULL; sd1 = seleMan1_.nextSelected(iterator1)){
	mapIndex1 = globalToLocal_[sd1->getGlobalIndex()];
	if(mapIndex1 == -1){
	  cerr << "mapindex1: -1\n"
	}
	gcn = 0.0;
	for(int i = 0; i < listNN_[mapIndex1].size(); i++){
	  tempIndex = listNN_[mapIndex1][i];
	  gcn += nearestNeighbors_[tempIndex];
	}
	gcn = gcn / nnMax_;
	binIndex = int(gcn*bins_);
	histogram_[binIndex] += 1;

      }
      gcnStream << "#Selection Count: " << selectionCount1_ << "\n";
      gcnStream << "#Frame " << istep << "\n";
      for(int n = 0; n < histogram_.size(); n++){
	binValue = n/(1.0*bins_);
	gcnStream << binValue << "  " << (histogram_[n]/(selectionCount1_*1.0)) << "\n";
      } 
    }   
    gcnStream.close();
  }
}
