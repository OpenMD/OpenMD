/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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

namespace oopse {
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
    double getTauThermostat() {
      return tauThermostat;
    }

    void setTauThermostat(double tt) {
      tauThermostat = tt;
    }

    double getTauBarostat() {
      return tauBarostat;
    }
    void setTauBarostat(double tb) {
      tauBarostat = tb;
    }

    double getTargetTemp() {
      return targetTemp;
    }
            
    void setTargetTemp(double tt) {
      targetTemp = tt;
    }

    double getTargetPressure() {
      return targetTemp;
    }
            
    void setTargetPressure(double tp) {
      targetPressure = tp;
    }

    double getChiTolerance() {
      return chiTolerance;
    }
            
    void setChiTolerance(double tol) {
      chiTolerance = tol;
    }

    double getEtaTolerance() {
      return etaTolerance;
    }
            
    void setEtaTolerance(double tol) {
      etaTolerance = tol;
    }

  protected:

    virtual void integrateStep() {
      needStress= true;
      VelocityVerletIntegrator::integrateStep();
    }

    virtual void doUpdate();

    virtual void resetIntegrator();

    virtual void resetEta();
    
    double NkBT;
    double fkBT;

    double tt2;            
    double tb2;
            
    double instaTemp;
    double instaPress;
    double instaVol;


    // targetTemp, targetPressure, and tauBarostat must be set.
    // One of qmass or tauThermostat must be set;

    double targetTemp;
    double targetPressure;
    double tauThermostat;
    double tauBarostat;

    std::vector<Vector3d> oldPos;
    std::vector<Vector3d> oldVel;
    std::vector<Vector3d> oldJi;

    double etaTolerance;
       
    double chi;
    double integralOfChidt;
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
            
    virtual double calcConservedQuantity() = 0;      

    virtual void loadEta() = 0;
    virtual void saveEta() = 0;
            
    int maxIterNum_;

    double chiTolerance;    
  };

}      //end namespace oopse

#endif //INTEGRATORS_NPT_HPP
