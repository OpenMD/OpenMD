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
#include "applications/staticProps/CoordinationNumber.hpp"
#include "io/DumpReader.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"

namespace OpenMD {

  CoordinationNumber::CoordinationNumber(SimInfo* info,
                                         const std::string& filename,
                                         const std::string& sele1,
                                         const std::string& sele2,
                                         RealType rCut, int bins):
    StaticAnalyser(info, filename, bins), rCut_(rCut), bins_(bins),
    sele1_(sele1), seleMan1_(info), evaluator1_(info),
    sele2_(sele2), seleMan2_(info), evaluator2_(info) {

    setAnalysisType("Coordination Number Distribution");
    setOutputName(getPrefix(filename) + ".cn");
    
    nnMax_ = 12;
    RealType binMax_ = nnMax_ * 1.5;
    delta_ = binMax_ / bins_;
    
    std::stringstream params;
    params << " rcut = " << rCut_
           << ", nbins = " << bins_
           << ", nnMax = " << nnMax_;
    const std::string paramString = params.str();
    setParameterString( paramString );

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
  }

  CoordinationNumber::~CoordinationNumber() {
    histogram_.clear();
  }

  void CoordinationNumber::process() {
    SelectionManager common(info_);
    
    std::vector<std::vector<int> > listNN;
    std::vector<int> globalToLocal;

    StuntDouble* sd1;
    StuntDouble* sd2;

    Snapshot* currentSnapshot_;
    bool usePeriodicBoundaryConditions_ =
      info_->getSimParams()->getUsePeriodicBoundaryConditions();

    int iterator1;
    int iterator2;
    unsigned int mapIndex1(0);
    unsigned int mapIndex2(0);
    unsigned int whichBin(0);
    RealType cn(0.0);
    Vector3d pos1;
    Vector3d pos2;
    Vector3d diff;
    RealType distance;

    histogram_.clear();
    histogram_.resize(bins_, 0.0);

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    count_ = 0;

    //First have to calculate lists of nearest neighbors (listNN_):
    
    for(int istep = 0; istep < nFrames; istep += step_){
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
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
      globalToLocal.resize(info_->getNGlobalAtoms() +
                           info_->getNGlobalRigidBodies(), -1);
      for (unsigned int i = 0; i < listNN.size(); i++)         
        listNN.at(i).clear();
      listNN.clear();
      listNN.resize(commonCount);
      
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
	    if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(diff);
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
      
      // Fill up the histogram with cn values

      for(sd1 = seleMan1_.beginSelected(iterator1); sd1 != NULL;
          sd1 = seleMan1_.nextSelected(iterator1)){
            
	mapIndex1 = globalToLocal.at(sd1->getGlobalIndex());

        cn = computeCoordination(mapIndex1, listNN);        
        whichBin = int(cn / delta_);
        
        if (whichBin < histogram_.size()) {
          histogram_[whichBin] += 1;
        } else {
          sprintf(painCave.errMsg, "Coordination Number: Error: "
                  "In frame, %d, object %d has CN %lf outside of range max.\n",
                  istep, sd1->getGlobalIndex(), cn );
          painCave.isFatal = 1;
          simError();  
        }
      }
      count_ += selectionCount1_;
    }

    for(unsigned int n = 0; n < histogram_.size(); n++){
      if (count_ > 0) 
        histogram_[n] /= RealType(count_);
      else
        histogram_[n] = 0.0;               
    } 
  
    writeOutput();
  }

  RealType CoordinationNumber::computeCoordination(int a,
                                                   vector<vector<int> > nl) {
    return RealType(nl.at(a).size());
  }  
    
  void CoordinationNumber::writeOutput() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    
    if (ofs.is_open()) {
      
      Revision r;
      RealType binValue(0.0);
      
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << sele1_ ;
      ofs << "\"\tselection script2: \"" << sele2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      
      for(unsigned int n = 0; n < histogram_.size(); n++){
        binValue = n * delta_;
        ofs << binValue << "\t"
            << histogram_[n] / delta_
            << "\n";
      } 
    }   
    ofs.close();
  }

  SCN::SCN(SimInfo* info,
           const std::string& filename,
           const std::string& sele1,
           const std::string& sele2,
           RealType rCut, int bins):
    CoordinationNumber(info, filename, sele1, sele2, rCut, bins) {
    setOutputName(getPrefix(filename) + ".scn");
    setAnalysisType("Secondary Coordination Number");
  }
  
  SCN::~SCN() {
  }
 
  RealType SCN::computeCoordination(int a, vector<vector<int> > nl) {
    RealType scn = 0.0;
    int b;
    
    if (nl.at(a).size() != 0) {
      for(unsigned int i = 0; i < nl.at(a).size(); i++){
        // b is the index of one of i's nearest neighbors
        b = nl.at(a).at(i);
        scn += nl.at(b).size();
      }
      return scn / RealType(nl.at(a).size());
    } else {
      return 0.0;
    }
    
  }

  GCN::GCN(SimInfo* info,
           const std::string& filename,
           const std::string& sele1,
           const std::string& sele2,
           RealType rCut, int bins):
    CoordinationNumber(info, filename, sele1, sele2, rCut, bins) {
    setOutputName(getPrefix(filename) + ".gcn");
    setAnalysisType("Generalized Coordination Number");
  }

  GCN::~GCN() {
  }

  RealType GCN::computeCoordination(int a, vector<vector<int> > nl) {
    RealType gcn = 0.0;
    int b;
    for(unsigned int i = 0; i < nl.at(a).size(); i++){
      // b is the index of one of i's nearest neighbors
      b = nl.at(a).at(i);
      gcn += nl.at(b).size();
    }
    return gcn / RealType(nnMax_);
  }
}
