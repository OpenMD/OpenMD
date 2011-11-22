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
 */
#ifndef APPLICATIONS_DYNAMICPROPS_TIMECORRFUNC_HPP
#define APPLICATIONS_DYNAMICPROPS_TIMECORRFUNC_HPP

#include <string>
#include <vector>

#include "brains/SimInfo.hpp"
#include "brains/BlockSnapshotManager.hpp"

#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  /**
   * @class TimeCorrFunc TimeCorrFunc.hpp "applications/dynamicProps/TimeCorrFunc"
   * @brief Base class for Correlation function
   */
 
  class TimeCorrFunc {
  public:
    TimeCorrFunc(SimInfo* info, const std::string& filename, 
		 const std::string& sele1, const std::string& sele2,
                 int storageLayout, long long int memSize);
        
    void doCorrelate();


    void setOutputName(const std::string& filename) {
      outputFilename_ = filename;
    }

    const std::string& getOutputFileName() const {
      return outputFilename_;
    }


    const std::string& getCorrFuncType() const {
      return corrFuncType_;
    }

    void setCorrFuncType(const std::string& type) {
      corrFuncType_ = type;
    }

    void setExtraInfo(const std::string& extra) {
      extra_ = extra;
    }
            
  protected:
        
    virtual void preCorrelate();        
    virtual void postCorrelate();
    virtual void updateFrame(int frame);

    RealType deltaTime_;
    int nTimeBins_;
    std::vector<RealType> histogram_;
    std::vector<int> count_;
    std::vector<RealType> time_;
        
    SimInfo* info_;
    int storageLayout_;
    long long int memSize_;
    std::string dumpFilename_;        
    SelectionManager seleMan1_;
    SelectionManager seleMan2_;          

    BlockSnapshotManager* bsMan_;       
        
  private:

    void correlateBlocks(int block1, int block2);
    virtual void correlateFrames(int frame1, int frame2) = 0;       
        
    virtual void writeCorrelate();

    virtual void validateSelection(const SelectionManager& seleMan) {}


    std::string selectionScript1_;
    std::string selectionScript2_;
        
    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;
 

    std::string outputFilename_;

    std::string corrFuncType_;
    std::string extra_;
  };

}
#endif
