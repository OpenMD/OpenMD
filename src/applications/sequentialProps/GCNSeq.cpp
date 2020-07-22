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
#include "applications/sequentialProps/GCNSeq.hpp"
#include "io/DumpReader.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"

namespace OpenMD {

  GCNSeq::GCNSeq(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 RealType rCut, int bins):
    SequentialAnalyzer(info, filename, sele1, sele2), rCut_(rCut), bins_(bins) {

    setSequenceType("Generalized Coordination Number Distribution");
    setOutputName(getPrefix(filename) + ".gcnSeq");

    nnMax_ = 12;
    RealType binMax_ = nnMax_ * 1.5;
    delta_ = binMax_ / bins_;
    usePBC_ = info->getSimParams()->getUsePeriodicBoundaryConditions();

    std::stringstream params;
    params << " rcut = " << rCut_
           << ", nbins = " << bins_
           << ", max neighbors = " << nnMax_;
    const std::string paramString = params.str();
    setParameterString( paramString );
  }

  GCNSeq::~GCNSeq() {
    histogram_.clear();
  }

  void GCNSeq::doFrame(int istep) {
    SelectionManager common(info_);
    
    std::vector<std::vector<int> > listNN;
    std::vector<int> globalToLocal;

    StuntDouble* sd1;
    StuntDouble* sd2;
        
    int iterator1;
    int iterator2;
    unsigned int mapIndex1(0);
    unsigned int mapIndex2(0);
    unsigned int tempIndex(0);
    unsigned int whichBin(0);
    RealType gcn(0.0);
    Vector3d pos1;
    Vector3d pos2;
    Vector3d diff;
    RealType distance;
    
    //First have to calculate lists of nearest neighbors (listNN_):
             
    selectionCount1_ = seleMan1_.getSelectionCount();
    selectionCount2_ = seleMan2_.getSelectionCount();
    
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
    std::vector<RealType> histo;
    histo.resize(bins_, 0.0);
    
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
          if (usePBC_) currentSnapshot_->wrapVector(diff);
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
      whichBin = int(gcn / delta_);
      if (whichBin < histo.size()) {
        histo[whichBin] += 1;
      } else {
        cerr << "In frame " <<  istep <<  ", object "
             << sd1->getGlobalIndex() << " has GCN value = " << gcn << "\n";
      }
    }

    for(unsigned int n = 0; n < histo.size(); n++){
      if (selectionCount1_ > 0)
        histo[n] /= RealType(selectionCount1_);
      else
        histo[n] = 0.0;
    }

    count_.push_back(selectionCount1_);
    histogram_.push_back(histo);
  }

  void GCNSeq::writeSequence() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    
    if (ofs.is_open()) {
      
      Revision r;
      RealType binValue(0.0);
      
      ofs << "# " << getSequenceType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_ ;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      
      ofs << "#time\tvalue\n";

      for (unsigned int i = 0; i < times_.size(); ++i) {
        ofs << "#Frame " << i << "\n";
        ofs << "#Selection 1 Count: " << count_[i] << "\n";

        for(unsigned int n = 0; n < histogram_[i].size(); n++){
          binValue = n * delta_;
          ofs << binValue << "\t"
              << histogram_[i][n]
              << "\n";
        }
      }
    } else {
      sprintf(painCave.errMsg,
              "GCN::writeSequence Error: failed to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();        
    }
    
    ofs.close();    
  }
}
