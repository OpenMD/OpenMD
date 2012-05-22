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
 
/**
 * @file NPT.hpp
 * @author tlin
 * @date 11/19/2004
 * @time 13:25am
 * @version 1.0
 */

#ifndef INTEGRATORS_NPT_HPP
#define INTEGRATORS_NPT_HPP

#include "integrators/VelocityVerletIntegrator.hpp"

namespace OpenMD {
  class NPT : public VelocityVerletIntegrator {
  public:

    NPT(SimInfo * info);
    virtual ~NPT();


    int getMaxIterationNumber() {
      return maxIterNum_;
    }

    void setMaxIterationNumber(int maxIter) {
      maxIterNum_ = maxIter;
    }
    RealType getTauThermostat() {
      return tauThermostat;
    }

    void setTauThermostat(RealType tt) {
      tauThermostat = tt;
    }

    RealType getTauBarostat() {
      return tauBarostat;
    }
    void setTauBarostat(RealType tb) {
      tauBarostat = tb;
    }

    RealType getTargetTemp() {
      return targetTemp;
    }
            
    void setTargetTemp(RealType tt) {
      targetTemp = tt;
    }

    RealType getTargetPressure() {
      return targetTemp;
    }
            
    void setTargetPressure(RealType tp) {
      targetPressure = tp;
    }

    RealType getChiTolerance() {
      return chiTolerance;
    }
            
    void setChiTolerance(RealType tol) {
      chiTolerance = tol;
    }

    RealType getEtaTolerance() {
      return etaTolerance;
    }
            
    void setEtaTolerance(RealType tol) {
      etaTolerance = tol;
    }

  protected:

    virtual void integrateStep() {
      needStress= true;
      VelocityVerletIntegrator::integrateStep();
    }

    virtual void doUpdateSizes();

    virtual void resetIntegrator();

    virtual void resetEta();
    
    RealType NkBT;
    RealType fkBT;

    RealType tt2;            
    RealType tb2;
            
    RealType instaTemp;
    RealType instaPress;
    RealType instaVol;


    // targetTemp, targetPressure, and tauBarostat must be set.
    // One of qmass or tauThermostat must be set;

    RealType targetTemp;
    RealType targetPressure;
    RealType tauThermostat;
    RealType tauBarostat;

    std::vector<Vector3d> oldPos;
    std::vector<Vector3d> oldVel;
    std::vector<Vector3d> oldJi;

    RealType etaTolerance;
       
    RealType chi;
    RealType integralOfChidt;
    Mat3x3d press;
            
  private:

    virtual void moveA();
    virtual void moveB();
            
    virtual void getVelScaleA(Vector3d& sc, const Vector3d& vel) = 0;
            
    virtual void getVelScaleB(Vector3d& sc, int index) = 0;
            
    virtual void getPosScale(const Vector3d& pos, const Vector3d& COM, 
			     int index, Vector3d& sc) = 0;

    virtual void calcVelScale() = 0;

    virtual bool etaConverged() = 0;

    virtual void evolveEtaA() = 0;

    virtual void evolveEtaB() = 0;

    virtual void scaleSimBox() = 0;
            
    virtual RealType calcConservedQuantity() = 0;      

    virtual void loadEta() = 0;
    virtual void saveEta() = 0;
            
    int maxIterNum_;

    RealType chiTolerance;    
  };

}      //end namespace OpenMD

#endif //INTEGRATORS_NPT_HPP
