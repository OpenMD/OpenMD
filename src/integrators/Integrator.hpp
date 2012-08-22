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
 * @file Integrator.hpp
 * @author tlin
 * @date 11/08/2004
 * @time 13:25am
 * @version 1.0
 */

#ifndef INTEGRATORS_INTEGRATOR_HPP
#define INTEGRATORS_INTEGRATOR_HPP

#include "brains/ForceManager.hpp"
#include "restraints/ThermoIntegrationForceManager.hpp"
#include "brains/Stats.hpp"
#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"
#include "integrators/RotationAlgorithm.hpp"
#include "flucq/FluctuatingChargePropagator.hpp"
#include "integrators/Velocitizer.hpp"
#include "rnemd/RNEMD.hpp"
#include "constraints/Rattle.hpp"

namespace OpenMD {


  /**
   * @class Integrator Integrator.hpp "integrators/Integrator.hpp"
   * @brief Base class of Integrator
   * @todo document
   */
  class Integrator {
  public:

    virtual ~Integrator();

    //avoid public virtual function        
    void integrate() {
      doIntegrate();
    }

    void updateSizes() {
      doUpdateSizes();
      flucQ_->updateSizes();
    }

    void setForceManager(ForceManager* forceMan) {

      if (forceMan_ != forceMan && forceMan_  != NULL) {
	delete forceMan_;
      }
      forceMan_ = forceMan;
      // forward this on:
      if (flucQ_ != NULL) {
        flucQ_->setForceManager(forceMan_);
      }
    }

    void setVelocitizer(Velocitizer* velocitizer) {
      if (velocitizer_ != velocitizer && velocitizer_ != NULL) {
	delete velocitizer_;
      }
      velocitizer_ = velocitizer;
    }

    void setFluctuatingChargePropagator(FluctuatingChargePropagator* prop) {
      if (prop != flucQ_ && flucQ_ != NULL){            
        delete flucQ_;
      }            
      flucQ_ = prop;
      if (forceMan_ != NULL) {
        flucQ_->setForceManager(forceMan_);
      }
    }

    void setRotationAlgorithm(RotationAlgorithm* algo) {
      if (algo != rotAlgo_ && rotAlgo_ != NULL){            
        delete rotAlgo_;
      }
            
      rotAlgo_ = algo;
    }

    void setRNEMD(RNEMD* rnemd) {
      if (rnemd_ != rnemd && rnemd_  != NULL) {
	delete rnemd_;
      }
      rnemd_ = rnemd;
    }
    
  protected:

    Integrator(SimInfo* info);

    virtual void doIntegrate() = 0;

    virtual void doUpdateSizes() {}
        
    void saveConservedQuantity() {
      snap->setConservedQuantity( calcConservedQuantity() );
    }
        
    SimInfo* info_;
    Globals* simParams;
    ForceManager* forceMan_;
    RotationAlgorithm* rotAlgo_;
    FluctuatingChargePropagator* flucQ_;
    Rattle* rattle_;
    Velocitizer* velocitizer_;
    RNEMD* rnemd_;

    bool needPotential;
    bool needStress;
    bool needReset;    
    bool needVelocityScaling;
    RealType targetScalingTemp;

    bool useRNEMD;    
    
    Stats* stats;
    DumpWriter* dumpWriter;
    StatWriter* statWriter;
    Thermo thermo;

    RealType runTime;
    RealType sampleTime;
    RealType statusTime;
    RealType thermalTime;
    RealType resetTime;
    RealType RNEMD_exchangeTime;
    RealType dt;

    Snapshot* snap; // During the integration, the address of snap Will not change

        
  private:
        
    virtual RealType calcConservedQuantity() = 0;
        
    virtual DumpWriter* createDumpWriter() = 0;

    virtual StatWriter* createStatWriter() = 0;
  };

    
}
#endif //INTEGRATORS_INTEGRATOR_HPP
