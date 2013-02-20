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
 
/**
 * @file VelocityVerletIntegrator.hpp
 * @author tlin
 * @date 11/08/2004
 * @version 1.0
 */

#ifndef INTEGRATORS_VELOCITYVERLETINTEGRATOR_HPP
#define INTEGRATORS_VELOCITYVERLETINTEGRATOR_HPP

#include "integrators/Integrator.hpp"
#include "integrators/RotationAlgorithm.hpp"
#include "flucq/FluctuatingChargePropagator.hpp"
#include "constraints/Rattle.hpp"
#include "utils/ProgressBar.hpp"

namespace OpenMD {

  /**
   * @class VelocityVerletIntegrator VelocityVerletIntegrator.hpp "integrators/VelocityVerletIntegrator.hpp"
   * @brief  Velocity-Verlet Family Integrator
   * Template pattern is used in VelocityVerletIntegrator class. 
   */
  class VelocityVerletIntegrator : public Integrator {
  public:
    virtual ~VelocityVerletIntegrator();
        
  protected:

    VelocityVerletIntegrator(SimInfo* info);
    virtual void doIntegrate();
    virtual void initialize();
    virtual void preStep();
    virtual void integrateStep();        
    virtual void postStep();
    virtual void finalize();
    virtual void resetIntegrator() {}
    
    RealType dt2;
    RealType currSample;
    RealType currStatus;
    RealType currThermal;
    RealType currReset;
    RealType currRNEMD;
        
  private:
        
    virtual void calcForce();    
    virtual void moveA() = 0;
    virtual void moveB() = 0;
    virtual RealType calcConservedQuantity() = 0;
    virtual DumpWriter* createDumpWriter();
    virtual StatWriter* createStatWriter();

    ProgressBar* progressBar;

  };

} //end namespace OpenMD
#endif //INTEGRATORS_VELOCITYVERLETINTEGRATOR_HPP
