/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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

/**
 * @file NPT.hpp
 * @author tlin
 * @date 11/19/2004
 * @version 1.0
 */

#ifndef INTEGRATORS_NPT_HPP
#define INTEGRATORS_NPT_HPP

#include "integrators/VelocityVerletIntegrator.hpp"

namespace OpenMD {
  class NPT : public VelocityVerletIntegrator {
  public:
    NPT(SimInfo* info);

    int getMaxIterationNumber() { return maxIterNum_; }

    void setMaxIterationNumber(int maxIter) { maxIterNum_ = maxIter; }
    RealType getTauThermostat() { return tauThermostat; }

    void setTauThermostat(RealType tt) { tauThermostat = tt; }

    RealType getTauBarostat() { return tauBarostat; }
    void setTauBarostat(RealType tb) { tauBarostat = tb; }

    RealType getTargetTemp() { return targetTemp; }

    void setTargetTemp(RealType tt) { targetTemp = tt; }

    RealType getTargetPressure() { return targetTemp; }

    void setTargetPressure(RealType tp) { targetPressure = tp; }

    RealType getChiTolerance() { return chiTolerance; }

    void setChiTolerance(RealType tol) { chiTolerance = tol; }

    RealType getEtaTolerance() { return etaTolerance; }

    void setEtaTolerance(RealType tol) { etaTolerance = tol; }

  protected:
    virtual void step() {
      needVirial = true;
      VelocityVerletIntegrator::step();
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

    pair<RealType, RealType> thermostat;
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
    RealType chiTolerance;

  protected:
    int maxIterNum_;
  };
}  // namespace OpenMD

#endif  // INTEGRATORS_NPT_HPP
