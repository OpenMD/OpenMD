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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). 
 */

#include "analysis/NonSpatialStatistics.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#ifdef _MSC_VER
#define isnan(x) _isnan((x))
#define isinf(x) (!_finite(x) && !_isnan(x))
#endif

namespace OpenMD {
  
  NonSpatialStatistics::NonSpatialStatistics(SimInfo* info, 
                                       const string& sele, int nbins)
    : StaticAnalyser(info, nbins), selectionScript_(sele),
      evaluator_(info), seleMan_(info) {
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    string prefixFileName = info->getPrefixFileName();
    setOutputName(prefixFileName + ".nspst");
  }

  NonSpatialStatistics::NonSpatialStatistics(SimInfo* info, 
					     const string& sele1,
					     const string& sele2,
					     int nbins)
    : StaticAnalyser(info, nbins), selectionScript_(sele1),
      evaluator_(info), seleMan_(info) {

    SelectionEvaluator evaluator2(info);
    SelectionManager seleMan2(info);
 
    evaluator_.loadScriptString(sele1);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    evaluator2.loadScriptString(sele2);
    if (!evaluator2.isDynamic()) {
      seleMan2.setSelectionSet(evaluator2.evaluate());
    }

    string prefixFileName = info->getPrefixFileName();
    setOutputName(prefixFileName + ".nspst");
  }

  NonSpatialStatistics::NonSpatialStatistics(SimInfo* info, 
					     const string& sele1,
					     const string& sele2,
					     const string& sele3,
					     int nbins)
    : StaticAnalyser(info, nbins), selectionScript_(sele1),
      evaluator_(info), seleMan_(info) {

    SelectionEvaluator evaluator2(info);
    SelectionManager seleMan2(info);
    
    SelectionEvaluator evaluator3(info);
    SelectionManager seleMan3(info);
    
    evaluator_.loadScriptString(sele1);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    evaluator2.loadScriptString(sele2);
    if (!evaluator2.isDynamic()) {
      seleMan2.setSelectionSet(evaluator2.evaluate());
    }

    evaluator3.loadScriptString(sele3);
    if (!evaluator3.isDynamic()) {
      seleMan3.setSelectionSet(evaluator3.evaluate());
    }

    string prefixFileName = info->getPrefixFileName();
    setOutputName(prefixFileName + ".nspst");
  }

  NonSpatialStatistics::~NonSpatialStatistics() {
    // vector<OutputData*>::iterator i;
    // OutputData* outputData;
    
    // for(outputData = beginOutputData(i); outputData; 
    //     outputData = nextOutputData(i)) {
    //   delete outputData;
    // }
    // data_.clear();

    // delete counts_;
  }

  


  void NonSpatialStatistics::processDump() {
    string dumpFileName_ = info_->getDumpFileName();
    DumpReader reader(info_, dumpFileName_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      processFrame(istep);
    }
    processHistogram();
    writeOutput();
  }


  void NonSpatialStatistics::processFrame(int istep) {
    StuntDouble* sd;
    int i;
        
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // loop over the selected atoms:
    
    for (sd = seleMan_.beginSelected(i); sd != NULL; 
         sd = seleMan_.nextSelected(i)) {
      
      // figure out where that object is:
      
      Vector3d pos = sd->getPos();
      
      int bin = getBin(pos);
    
      
      // forward the work of statistics on to the subclass:
      
      processStuntDouble( sd, bin );

      dynamic_cast<Accumulator *>(counts_->accumulator[bin])->add(1);
    }
  }


  int NonSpatialStatistics::getBin(Vector3d pos) {
    unsigned int zero = 0;
    return zero;  
  }
  
  
  
}

