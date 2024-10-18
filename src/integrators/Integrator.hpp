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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#ifndef INTEGRATORS_INTEGRATOR_HPP
#define INTEGRATORS_INTEGRATOR_HPP

#include <memory>

#include "brains/ForceManager.hpp"
#include "brains/Stats.hpp"
#include "brains/Velocitizer.hpp"
#include "constraints/Rattle.hpp"
#include "flucq/FluctuatingChargePropagator.hpp"
#include "integrators/DLM.hpp"
#include "integrators/RotationAlgorithm.hpp"
#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"
#include "rnemd/RNEMD.hpp"
#include "utils/ProgressBar.hpp"

namespace OpenMD {
  /**
   * @brief Declaration of the Integrator base class, which
   * all other integrators inherit from.
   *
   * It provides an abstract integrate() function which will be called
   * by the host application to integrate the system forward in time.
   *
   * Several convenience functions are also provided that are commonly
   * needed by subclasses.
   */
  class Integrator {
  public:
    /**
     * @brief Default Destructor
     */
    virtual ~Integrator();

    void integrate();
    void updateSizes();
    void setVelocitizer(std::unique_ptr<Velocitizer> velocitizer);
    void setFluctuatingChargePropagator(FluctuatingChargePropagator* prop);
    void setRotationAlgorithm(RotationAlgorithm* algo);
    void setRNEMD(std::unique_ptr<RNEMD::RNEMD> rnemd);

  protected:
    Integrator(SimInfo* info);

    virtual void initialize();
    virtual void preStep();
    /** @brief Computes an integration step from t to t+dt
     *
     * This function must be implemented by any subclasses, and computes a
     * single integration step from the current time (t) to (t+dt).
     */
    virtual void step() = 0;
    virtual void calcForce();
    virtual void postStep();
    virtual void finalize();
    virtual void resetIntegrator() {}
    virtual void doUpdateSizes() {}
    void saveConservedQuantity();

    RealType dt, dt2;
    RealType runTime;
    RealType sampleTime;
    RealType statusTime;
    RealType thermalTime;
    RealType resetTime;
    RealType RNEMD_exchangeTime;
    RealType currSample;
    RealType currStatus;
    RealType currThermal;
    RealType currReset;
    RealType currRNEMD;

    SimInfo* info_ {nullptr};
    Globals* simParams {nullptr};
    ForceManager* forceMan_ {nullptr};
    RotationAlgorithm* rotAlgo_ {nullptr};
    FluctuatingChargePropagator* flucQ_ {nullptr};
    Rattle* rattle_ {nullptr};
    std::unique_ptr<Velocitizer> velocitizer_ {nullptr};
    std::unique_ptr<RNEMD::RNEMD> rnemd_ {nullptr};

    bool needPotential {false};
    bool needVirial {false};
    bool needReset {false};
    bool needVelocityScaling {false};
    bool useRNEMD {false};

    RealType targetScalingTemp;

    Stats* stats {nullptr};
    DumpWriter* dumpWriter {nullptr};
    StatWriter* statWriter {nullptr};
    Thermo thermo;

    Snapshot* snap {nullptr};
    ProgressBarPtr progressBar {nullptr};

    const RealType dtEps = 1.0e-4;

  private:
    virtual RealType calcConservedQuantity() = 0;
    virtual DumpWriter* createDumpWriter();
    virtual StatWriter* createStatWriter();
  };
}  // namespace OpenMD

#endif  // INTEGRATORS_INTEGRATOR_HPP
