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
#ifndef APPLICATIONS_DYNAMICPROPS_MULTIPASSCORRFUNC_HPP
#define APPLICATIONS_DYNAMICPROPS_MULTIPASSCORRFUNC_HPP

#include <string>
#include <vector>

#include "applications/dynamicProps/DynamicProperty.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  //! Computes a correlation function by scanning a trajectory once to precompute quantities to be correlated

  class MultipassCorrFunc : public DynamicProperty {
  public:
    MultipassCorrFunc(SimInfo* info, const std::string& filename, 
                      const std::string& sele1, const std::string& sele2,
                      int storageLayout);
    
    virtual ~MultipassCorrFunc(){ }    
    virtual void doCorrelate();

    const std::string& getCorrFuncType() const {
      return corrFuncType_;
    }

    void setCorrFuncType(const std::string& type) {
      corrFuncType_ = type;
    }

    void setParameterString(const std::string& params) {
      paramString_ = params;
    }
            
  protected:
        
    virtual void preCorrelate();
    virtual void correlation();
    virtual void postCorrelate();
    virtual void computeFrame(int frame);
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) = 0;
    virtual void correlateFrames(int frame1, int frame2, int timeBin);
    virtual RealType calcCorrVal(int frame1, int frame2, int id1, int id2) = 0;

    int storageLayout_;

    RealType deltaTime_;
    int nTimeBins_;
    int nFrames_;
    std::vector<RealType> histogram_;
    std::vector<int> count_;
    std::vector<RealType> times_;
    bool uniqueSelections_;
        
    SimInfo* info_;
    DumpReader* reader_;
    std::string dumpFilename_;        
    SelectionManager seleMan1_;
    SelectionManager seleMan2_;

    virtual void writeCorrelate();
    virtual void validateSelection(const SelectionManager& seleMan) {}

    Snapshot* currentSnapshot_;
    
    std::string selectionScript1_;
    std::string selectionScript2_;
        
    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;

    bool autoCorrFunc_;

    std::string corrFuncType_;
    std::string paramString_;

    std::vector<std::vector<int> > sele1ToIndex_;
    std::vector<std::vector<int> > sele2ToIndex_;

  };


  class AutoCorrFunc : public MultipassCorrFunc {
  public:
    AutoCorrFunc(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 int storageLayout) :
      MultipassCorrFunc(info, filename, sele1, sele2, storageLayout) {
      autoCorrFunc_ = true;
    }
        
  protected:
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) { return -1; }
  };
  
  class CrossCorrFunc : public MultipassCorrFunc {
  public:
    CrossCorrFunc(SimInfo* info, const std::string& filename,
                  const std::string& sele1, const std::string& sele2,
                  int storageLayout) :
      MultipassCorrFunc(info, filename, sele1, sele2, storageLayout) {
      autoCorrFunc_ = false;
    }
        
  protected:
    virtual int computeProperty1(int frame, StuntDouble* sd) = 0;
    virtual int computeProperty2(int frame, StuntDouble* sd) = 0;
  };
}
#endif
