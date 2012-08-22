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
 
#ifndef INTEGRATOR_LDFORCEMANAGER_HPP
#define INTEGRATOR_LDFORCEMANAGER_HPP

#include "brains/ForceManager.hpp"
#include "primitives/Molecule.hpp"
#include "math/SeqRandNumGen.hpp"
#include "hydrodynamics/Shape.hpp"
#include "integrators/Velocitizer.hpp"

namespace OpenMD {
   
  struct SDShape{
    StuntDouble* sd;
    Shape* shape;
  };
    
  /**
   * @class LDForceManager
   * Force manager for Lagevin Dynamics applying friction and random 
   * forces as well as torques.
   */
  class LDForceManager : public ForceManager{
    
  public:
    LDForceManager(SimInfo * info);
    
    int getMaxIterationNumber() {
      return maxIterNum_;
    }
        
    void setMaxIterationNumber(int maxIter) {
      maxIterNum_ = maxIter;
    }

    RealType getForceTolerance() {
      return forceTolerance_;
    }

    void setForceTolerance(RealType tol) {
      forceTolerance_ = tol;
    }

    RealType getDt2() {
      return dt2_;
    }

    void setDt2(RealType dt2) {
      dt2_ = dt2;
    }


  protected:
    virtual void postCalculation();
    
  private:
    std::map<std::string, HydroProp*> parseFrictionFile(const std::string& filename);    
    void genRandomForceAndTorque(Vector3d& force, Vector3d& torque, unsigned int index, RealType variance);
    std::vector<HydroProp*> hydroProps_;
    SeqRandNumGen randNumGen_;    
    RealType variance_;
    RealType langevinBufferRadius_;
    RealType frozenBufferRadius_;
    bool sphericalBoundaryConditions_;
    Globals* simParams;
    Velocitizer* veloMunge;
    // convergence parameters:
    int maxIterNum_;
    RealType forceTolerance_;
    RealType dt2_;
  };
  
} //end namespace OpenMD
#endif //BRAINS_FORCEMANAGER_HPP

