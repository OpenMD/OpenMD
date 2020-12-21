/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */
 
#ifndef INTEGRATORS_INTEGRATOR_HPP
#define INTEGRATORS_INTEGRATOR_HPP

#include <memory>

#include "brains/ForceManager.hpp"
#include "restraints/ThermoIntegrationForceManager.hpp"
#include "brains/Stats.hpp"
#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"
#include "integrators/RotationAlgorithm.hpp"
#include "flucq/FluctuatingChargePropagator.hpp"
#include "brains/Velocitizer.hpp"
#include "rnemd/RNEMD.hpp"
#include "constraints/Rattle.hpp"
#include "integrators/DLM.hpp"
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
    void setForceManager(ForceManager* forceMan);
    void setVelocitizer(VelocitizerPtr velocitizer);
    void setFluctuatingChargePropagator(FluctuatingChargePropagator* prop);
    void setRotationAlgorithm(RotationAlgorithm* algo);
    void setRNEMD(RNEMD* rnemd);

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
    VelocitizerPtr velocitizer_ {nullptr};
    RNEMD* rnemd_ {nullptr};

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
    
  private:        
    virtual RealType calcConservedQuantity() = 0;
    virtual DumpWriter* createDumpWriter();
    virtual StatWriter* createStatWriter();

    ProgressBarPtr progressBar {nullptr};
  };
}

#endif // INTEGRATORS_INTEGRATOR_HPP
