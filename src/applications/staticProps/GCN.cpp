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

  GCN::GCN(SimInfo* info, const std::string& filename, const std::string& sele, int bins, int zRegions)
    : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan1_(info){


      evaluator_.loadScriptString(sele);
      if(!evaluator_.isDynamic()){
	seleMan1_.setSelectionSet(evaluator_.evaluate());
      }

      selectionCount_ = seleMan1_.getSelectionCount();

      filename_ = filename;
      regions_ = zRegions;
      bins_ = bins; //Treat bins as division, bins per single integer

      setOutputName(getPrefix(filename) + ".gcn");

      nnMax_ = 12.0;
      solShell_ = 3.3;
      cout << "Bins per integer: " << bins_ << "\n";
      bins_ = bins_*nnMax_; //temporary so 0.0 to 12.0 by 0.1
      cout << "Bins: " << bins_ << "\n";
      cout << "ZRegions: " << regions_ << "\n";
    }


  GCN::~GCN(){
    for(int i = 0; i < nearestNeighbors_.size(); i++){
      nearestNeighbors_[i].clear();
    }
    nearestNeighbors_.clear();
    for(int i = 0; i < listNN_.size(); i++){
      for(int j = 0; j < listNN_[i].size(); j++){
	listNN_[i][j].clear();
      }
      listNN_[i].clear();
    }
    listNN_.clear();

    for(int i = 0; i < histogram_.size(); i++){
      for(int j = 0; j < histogram_[i].size(); j++){
	histogram_[i][j].clear();
      }
      histogram_[i].clear();
    } 
    histogram_.clear();

    for(int i = 0; i < indexMapping_.size(); i++){
      indexMapping_[i].clear();
    }
    indexMapping_.clear();
  }


  void GCN::process() {
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd1;
    StuntDouble* sd2;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    std::vector<RealType> minZ;
    std::vector<RealType> maxZ;
    std::vector< std::vector< RealType > > zPos;

    zPos.resize(nFrames);
    minZ.resize(nFrames);
    maxZ.resize(nFrames);
    histogram_.resize(nFrames);
    nearestNeighbors_.resize(nFrames);
    listNN_.resize(nFrames);
    indexMapping_.resize(nFrames);

    for(int i = 0; i < nFrames; i++){
      nearestNeighbors_[i].resize(selectionCount_);
      listNN_[i].resize(selectionCount_);
      zPos[i].resize(selectionCount_);
      indexMapping_[i].resize(selectionCount_);
    }
  
    for(int i = 0; i < nFrames; i++){
      histogram_[i].resize(regions_);
      for(int j = 0; j < regions_; j++){
	//histogram_[i].resize(selectionCount_);
	histogram_[i][j].resize(bins_);
      }
    }

    int iterator1;
    int iterator2;
    int gIndex = 0;
    int mapIndex1 = 0;
    int mapIndex2 = 0;

    /* Arrays are sized correctly
     * Histogram[frames][regions][bins]
     * nearestNeighbors[frames][selectionCount] ([0] = 9, [1] = 12, [globalIndex] = nn)
     * listNN[frames][selectionCount][0-13] ([0] ?= 199, [1] = 252, [2] = globalIndex of nn
     *
     * method, loop over listNN and use listNN_[frames][i][0] as index into nearestNeigbhors_[frames][listNN_[frames][i][j]] when calculating gcn
     */


    //First have to calculate nearestNeighbors and listNN
    for(int istep = 0; istep < nFrames; istep += step_){
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for(mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)){
	for(rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)){
	  rb->updateAtoms();
	}
      }
      minZ[istep] = 100.00; //Bad Joe, bad
      maxZ[istep] = -100.00;
      gIndex = 0;
      mapIndex1 = 0;
      //Double looping, TODO get sd2 to start later, will need to double up on nearestNeighbors and stuff
      for(sd1 = seleMan1_.beginSelected(iterator1); sd1 != NULL; sd1 = seleMan1_.nextSelected(iterator1)){
	gIndex = sd1->getGlobalIndex();
	indexMapping_[istep][mapIndex1] = gIndex;//So the first entry in indexMapping may actually refer to Platinum 337, and that is okay
	zPos[istep][mapIndex1] = sd1->getPos()[2];
	if(zPos[istep][mapIndex1] < minZ[istep]){
	  minZ[istep] = zPos[istep][mapIndex1];
	}
	if(zPos[istep][mapIndex2] > maxZ[istep]){
	  maxZ[istep] = zPos[istep][mapIndex1];
	}
	mapIndex2 = 0;
	for(sd2 = seleMan1_.beginSelected(iterator2); sd2 != NULL; sd2 = seleMan1_.nextSelected(iterator2)){
	  //if(sd1->getGlobalIndex() == sd2->getGlobalIndex()){
	  if(mapIndex2 > selectionCount_){
	    cout << "mapIndex2: " << mapIndex2 << "\n";
	  }
	  if(mapIndex1 == mapIndex2){
	    //pass, same stuntDouble
	  }else{
	    Vector3d pos1 = sd1->getPos();
	    Vector3d pos2 = sd2->getPos();
	    Vector3d diff = pos2-pos1;
	    if(usePeriodicBoundaryConditions_){
	      currentSnapshot_->wrapVector(diff);
	    }
	    double distance = diff.length();
	    if(distance < solShell_){
	      listNN_[istep][mapIndex1].push_back(mapIndex2);
	      nearestNeighbors_[istep][mapIndex1] += 1;
	    }
	  }
	  mapIndex2++;
	}
	//cout << "istep: " << istep << "\tmapIndex1: " << mapIndex1 << "\n";
	//for(int k = 0; k < listNN_[istep][mapIndex1].size(); k++){
	//  cout << "k: " << listNN_[istep][mapIndex1][k] << "\n";
	//}
	mapIndex1++;
      }
    }
    cout << "Filled nearestNeighbors_ and listNN_\n\n";
    
    for(int test = 0; test < 10; test++){
      cout << "\tGlobal Index: " << indexMapping_[0][test] << "\t" << test << "\t" << nearestNeighbors_[0][test] << "\n";
    }    
    cout << "^ First 10 indices\n";
    cout << "selection Count: " << selectionCount_ << "\n";

    double regionShift = 0.0;
    int binIndex = 0; 
    RealType gcn = 0.0;
    int regionIndex = 0;
    int tempIndex = 0;
    // Need to check, but listNN_ and nearestNeighbors_ should be properly filled
    for(int istep = 0; istep < nFrames; istep += step_){
      regionShift = (maxZ[istep]-minZ[istep])/regions_;
      cout << "Frame: " << istep << "\n";
      for(int i = 0; i < selectionCount_; i++){
	//cout << "i: " << i << "\n";
	gcn = 0.0;
	for(int j = 0; j < listNN_[istep][i].size(); j++){
	  //cout << "after for listNN\n";
	  tempIndex = listNN_[istep][i][j];
	  if(tempIndex > selectionCount_){
	    for(int k = 0; k < listNN_[istep][i].size(); k++){
	      cout << "\t" << listNN_[istep][i][k] << "\n";
	    }
	    cout << "istep: " << istep << "\tatom: " << i << "\tsize: " << listNN_[istep][i].size() << "\n";
	    cout << "j: " << j << "\n";
	    cout << "tempIndex: " << tempIndex << "\n";
	  }
	  //cout << "tempIndex: " << tempIndex << "\n";
	  gcn += nearestNeighbors_[istep][tempIndex];
	  //cout << "after gcn addition\n\n";
	}
	gcn = gcn / nnMax_;
	regionIndex = int((zPos[istep][i] - minZ[istep])/regionShift);
	if(regionIndex >= regions_){
	  regionIndex = regions_ - 1;
	  cout << (zPos[istep][i] - minZ[istep])/regionShift << "\n";
	}
	binIndex = int(gcn*10);
	histogram_[istep][regionIndex][binIndex] += 1;
      }
    }
  }
}
