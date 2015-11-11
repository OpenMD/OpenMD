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

#include "applications/dynamicProps/RCorrFunc.hpp"

namespace OpenMD {
  RCorrFunc::RCorrFunc(SimInfo* info, const std::string& filename, 
                       const std::string& sele1, const std::string& sele2)
    : AutoCorrFunc(info, filename, sele1, sele2, 
                   DataStorage::dslPosition | DataStorage::dslAmat){
    
    setCorrFuncType("Mean Square Displacement");
    setOutputName(getPrefix(dumpFilename_) + ".rcorr");

    positions_.resize(nFrames_);
  }

  RCorrFuncZ::RCorrFuncZ(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2,
                         int nZbins)
    : AutoCorrFunc(info, filename, sele1, sele2, 
                   DataStorage::dslPosition | DataStorage::dslAmat){

    setCorrFuncType("Mean Square Displacement binned by Z");
    setOutputName(getPrefix(dumpFilename_) + ".rcorrZ");
    positions_.resize(nFrames_);
    zBins_.resize(nFrames_);
    nZBins_ = nZbins;
 
    histograms_.resize(nTimeBins_);
    counts_.resize(nTimeBins_);
    for (int i = 0; i < nTimeBins_; i++) {
      histograms_[i].resize(nZBins_);
      counts_[i].resize(nZBins_);
      std::fill(histograms_[i].begin(), histograms_[i].end(), 0.0);
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
  }
  
  RCorrFuncR::RCorrFuncR(SimInfo* info, const std::string& filename, 
                         const std::string& sele1, const std::string& sele2)
    : AutoCorrFunc(info, filename, sele1, sele2, 
                   DataStorage::dslPosition | DataStorage::dslAmat){
    
    // Turn on COM calculation in reader:
    bool ncp = true;
    reader_->setNeedCOMprops(ncp);    
    setCorrFuncType("MSD (radial projection)");
    setOutputName(getPrefix(dumpFilename_) + ".r_rcorr");
    positions_.resize(nFrames_);
  }

  int RCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    positions_[frame].push_back( sd->getPos() );
    return positions_[frame].size() - 1;
  }
  RealType RCorrFunc::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    Vector3d diff = positions_[frame2][id2] - positions_[frame1][id1];
    return diff.lengthSquare();
  }


  void RCorrFuncZ::computeFrame(int istep) {
    hmat_ = currentSnapshot_->getHmat();
    halfBoxZ_ = hmat_(2,2) / 2.0;      

    StuntDouble* sd;

    int isd1, isd2;
    unsigned int index;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    
    if (uniqueSelections_ && evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }      
    
    for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
         sd = seleMan1_.nextSelected(isd1)) {

      index = computeProperty1(istep, sd);        
      if (index == sele1ToIndex_[istep].size()) {
        sele1ToIndex_[istep].push_back(sd->getGlobalIndex());
      } else {
        sele1ToIndex_[istep].resize(index+1);
        sele1ToIndex_[istep][index] = sd->getGlobalIndex();
      }
    }
   
    if (uniqueSelections_) { 
      for (sd = seleMan2_.beginSelected(isd2); sd != NULL;
           sd = seleMan2_.nextSelected(isd2)) {        

        index = computeProperty1(istep, sd);

        if (index == sele2ToIndex_[istep].size()) {
          sele2ToIndex_[istep].push_back(sd->getGlobalIndex());
        } else {
          sele2ToIndex_[istep].resize(index+1);
          sele2ToIndex_[istep][index] = sd->getGlobalIndex();
        }
      }
    }
  }

  int RCorrFuncZ::computeProperty1(int frame, StuntDouble* sd) {
    Vector3d pos = sd->getPos();
    // we need the raw (not wrapped) positions for RMSD:
    positions_[frame].push_back( sd->getPos() );

    if (info_->getSimParams()->getUsePeriodicBoundaryConditions()) {
      currentSnapshot_->wrapVector(pos);
    }    
    int zBin = int(nZBins_ * (halfBoxZ_ + pos.z()) / hmat_(2,2));
    zBins_[frame].push_back(zBin);
    
    return positions_[frame].size() - 1;
  }

  void RCorrFuncZ::correlateFrames(int frame1, int frame2, int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;
    
    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

    RealType corrVal(0.0);

    s1 = sele1ToIndex_[frame1];

    if (uniqueSelections_) 
       s2 = sele2ToIndex_[frame2];
    else
       s2 = sele1ToIndex_[frame2];

    for (i1 = s1.begin(), i2 = s2.begin();
         i1 != s1.end() && i2 != s2.end(); ++i1, ++i2){
      
      // If the selections are dynamic, they might not have the
      // same objects in both frames, so we need to roll either of
      // the selections until we have the same object to
      // correlate.

      while ( i1 != s1.end() && *i1 < *i2 ) {
        ++i1;
      }
      
      while ( i2 != s2.end() && *i2 < *i1 ) {
        ++i2;
      }
          
      if ( i1 == s1.end() || i2 == s2.end() ) break;

      corrVal = calcCorrVal(frame1, frame2, i1 - s1.begin(), i2 - s2.begin(),
                            timeBin);
    }
  }
  
  RealType RCorrFuncZ::calcCorrVal(int frame1, int frame2, int id1, int id2,
                                   int timeBin) {

    int zBin1 = zBins_[frame1][id1];
    int zBin2 = zBins_[frame2][id2];
    
    if (zBin1 == zBin2) {
      Vector3d diff = positions_[frame2][id2] - positions_[frame1][id1];
      histograms_[timeBin][zBin1] += diff.lengthSquare();
      counts_[timeBin][zBin1]++;
    }
    return 0.0;
  }
  
  void RCorrFuncZ::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
      for (int j = 0; j < nZBins_; ++j) {        
        if (counts_[i][j] > 0) {
          histograms_[i][j] /= counts_[i][j];
        } else {
          histograms_[i][j] = 0;
        }
      }
    }
  }
  void RCorrFuncZ::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {

      ofs << "#" << getCorrFuncType() << "\n";
      ofs << "#selection script1: \"" << selectionScript1_ ;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs << "#extra information: " << extra_ << "\n";
      ofs << "#time\tcorrVal\n";

      for (int i = 0; i < nTimeBins_; ++i) {

        ofs << times_[i] - times_[0];

        for (int j = 0; j < nZBins_; ++j) {          
          ofs << "\t" << histograms_[i][j];
        }
        ofs << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "RCorrFuncZ::writeCorrelate Error: fail to open %s\n",
              getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs.close();    
  }


  int RCorrFuncR::computeProperty1(int frame, StuntDouble* sd) {
    // get the radial vector from the frame's center of mass:
    Vector3d coord_t = sd->getPos() - sd->getCOM();
    
    positions_[frame].push_back( coord_t.length() );
    return positions_[frame].size() - 1;
  }
  RealType RCorrFuncR::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType dr;
    dr  = positions_[frame2][id2] - positions_[frame1][id1];
    return dr * dr;
  }
}

