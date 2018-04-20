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
#ifndef ANALYSIS_STATICANALYSER_HPP
#define ANALYSIS_STATICANALYSER_HPP

#include <string>
#include "brains/SimInfo.hpp"
#include "brains/Snapshot.hpp"
#include "utils/Accumulator.hpp"

namespace OpenMD {
  enum OutputDataType {
    odtReal,
    odtVector3,
    odtArray2d,
    odtUnknownDataType
  };
  
  enum OutputDataHandling {
    odhAverage,
    odhTotal,
    odhLastValue,
    odhMax,
    odhUnknownDataHandling
  };
  
  struct OutputData {
    string title;
    string units;
    OutputDataType dataType;
    OutputDataHandling dataHandling;
    vector<BaseAccumulator*> accumulator;
    vector<vector<BaseAccumulator*> > accumulatorArray2d;
  };

  class StaticAnalyser{
  public:
    StaticAnalyser(SimInfo* info, unsigned int nbins);
    
    virtual ~StaticAnalyser() {}
    virtual void processFrame(int frame)=0;
    virtual void processDump()=0;

    void setOutputName(const std::string& filename) {
      outputFilename_ = filename;
    }
    
    const std::string& getOutputFileName() const {
      return outputFilename_;
    }
        
    void setStep(int step) {
      assert(step > 0);
      step_ =step;    
    }

    int getStep() { return step_;}
    
    const std::string& getAnalysisType() const {
      return analysisType_;
    }

    void setAnalysisType(const std::string& type) {
      analysisType_ = type;
    }
    
    void setParameterString(const std::string& params) {
      paramString_ = params;
    }

  protected:
    virtual void writeOutput();
    virtual void processHistogram();
    virtual void writeData(ostream& os, OutputData* dat, unsigned int bin);
    virtual void writeErrorBars(ostream& os, OutputData* dat, unsigned int bin);
    OutputData* beginOutputData(vector<OutputData*>::iterator& i);
    OutputData* nextOutputData(vector<OutputData*>::iterator& i);

    SimInfo* info_;
    std::string outputFilename_;
    int step_;
    std::string analysisType_;
    std::string paramString_;
    
    unsigned int nBins_;
    OutputData* counts_;
    vector<OutputData*> data_;

  };
}
#endif
