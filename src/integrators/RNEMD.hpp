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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
/**
 * @file RNEMD.hpp
 * @author gezelter
 * @date 03/13/2009
 * @time 15:56pm
 * @version 1.0
 */

#ifndef INTEGRATORS_RNEMD_HPP
#define INTEGRATORS_RNEMD_HPP
#include "brains/SimInfo.hpp"
#include "math/RandNumGen.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include <iostream>

namespace OpenMD {

  /**
   * @class RNEMD RNEMD.hpp "integrators/RNEMD.hpp"
   * @todo document
   */
  class RNEMD {
  public:
    RNEMD(SimInfo* info);
    virtual ~RNEMD();
    
    void doRNEMD();
    void doSwap();
    void doScale();
    void collectData();
    void getStarted();
    void getStatus();
    void set_RNEMD_exchange_time(RealType exchangeTime) {
      exchangeTime_ = exchangeTime;
    }
    void set_RNEMD_nBins(int nbins) { nBins_ = nbins; }
    void set_RNEMD_logWidth(int logWidth) { rnemdLogWidth_ = logWidth; }
    void set_RNEMD_exchange_total(RealType et) { exchangeSum_ = et; }
    void set_RNEMD_target_flux(RealType targetFlux) {targetFlux_ = targetFlux;}
    RealType get_RNEMD_exchange_total() { return exchangeSum_; }

  private:

    enum RNEMDTypeEnum {
      rnemdKineticSwap,
      rnemdKineticScale,
      rnemdPxScale,
      rnemdPyScale,
      rnemdPzScale,
      rnemdPx,
      rnemdPy,
      rnemdPz,
      rnemdUnknown
    };
    
    SimInfo* info_;
    RandNumGen* randNumGen_;
    std::map<std::string, RNEMDTypeEnum> stringToEnumMap_;
    RNEMDTypeEnum rnemdType_;
    std::string rnemdObjectSelection_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;
    bool usePeriodicBoundaryConditions_;
    bool output3DTemp_;
    int nBins_;
    int midBin_;
    int rnemdLogWidth_;
    RealType zShift_;
    RealType exchangeTime_;
    RealType targetFlux_;
    RealType exchangeSum_;
    int failTrialCount_;
    int failRootCount_;
    std::ofstream rnemdLog_;
    // keeps track of what's being averaged
    std::vector<RealType> valueHist_;
    std::vector<int> valueCount_, xyzTempCount_;
    // keeps track of the number of degrees of freedom being averaged
    std::vector<RealType> xTempHist_, yTempHist_, zTempHist_;
    std::ofstream xTempLog_, yTempLog_, zTempLog_;
  };

}
#endif //INTEGRATORS_RNEMD_HPP
