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

using namespace std;
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
    void doShiftScale();
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
    void set_RNEMD_target_JzKE(RealType targetJzKE) {targetJzKE_ = targetJzKE;}
    void set_RNEMD_target_jzpx(RealType targetJzpx) {targetJzpx_ = targetJzpx;}
    void set_RNEMD_target_jzpy(RealType targetJzpy) {targetJzpy_ = targetJzpy;}
    void set_RNEMD_target_jzpz(RealType targetJzpz) {targetJzpz_ = targetJzpz;}
    RealType get_RNEMD_exchange_total() { return exchangeSum_; }

  private:

    enum RNEMDTypeEnum {
      rnemdKineticSwap,
      rnemdKineticScale,
      rnemdKineticScaleVAM,
      rnemdKineticScaleAM,
      rnemdPxScale,
      rnemdPyScale,
      rnemdPzScale,
      rnemdPx,
      rnemdPy,
      rnemdPz,
      rnemdShiftScaleV,
      rnemdShiftScaleVAM,
      rnemdUnknown
    };
    
    SimInfo* info_;
    RandNumGen* randNumGen_;
    map<string, RNEMDTypeEnum> stringToEnumMap_;
    RNEMDTypeEnum rnemdType_;
    string rnemdObjectSelection_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;
    bool usePeriodicBoundaryConditions_;
    bool outputTemp_;
    bool outputVx_;
    bool outputVy_;
    bool output3DTemp_;
    bool outputRotTemp_;
    int nBins_; /**< The number of bins to divide the simulation box into.  */
    /*!
      The middle bin for the RNEMD method. midBin_ = nBins_/2;
      Depending on the setting of the flux, this box should contain the minimum energy (temperature)
      within the simulation. 
    */
    int midBin_;
    int rnemdLogWidth_; /**< Number of elements to print out in logs */
    RealType zShift_;
    RealType exchangeTime_;
    RealType targetFlux_;
    RealType targetJzKE_;
    RealType targetJzpx_;
    RealType targetJzpy_;
    RealType targetJzpz_;
    Vector3d jzp_, njzp_;
    RealType exchangeSum_;
    int failTrialCount_;
    int failRootCount_;
    ofstream tempLog_, vxzLog_, vyzLog_;
    ofstream xTempLog_, yTempLog_, zTempLog_, rotTempLog_;
    // keeps track of what's being averaged
    vector<RealType> tempHist_, pxzHist_, pyzHist_, mHist_;
    vector<RealType> xTempHist_, yTempHist_, zTempHist_, rotTempHist_;
    // keeps track of the number of degrees of freedom being averaged
    //vector<int> pxzCount_, pyzCount_;
    vector<int> tempCount_, xyzTempCount_, rotTempCount_;
  };

}
#endif //INTEGRATORS_RNEMD_HPP
