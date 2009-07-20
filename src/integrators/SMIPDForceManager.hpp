/*
 * Copyright (c) 2008 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef INTEGRATOR_SMIPDFORCEMANAGER_HPP
#define INTEGRATOR_SMIPDFORCEMANAGER_HPP

#include "brains/ForceManager.hpp"
#include "brains/Thermo.hpp"
#include "primitives/Molecule.hpp"
#include "integrators/Velocitizer.hpp"
#include "math/Hull.hpp"
#include "math/Triangle.hpp"
#include "math/SeqRandNumGen.hpp"

namespace oopse {
   
  /**
   * @class SMIDForceManager
   * Force manager for Surface Mesh Implicit Pressure Dynamics
   * applying friction and random forces as well as torques. 
   * Stochasitc force is determined by area of surface triangles 
   * on the convex hull. See: Kohanoff et al. CHEMPHYSCHEM 2005, 6, 1848-1852.
   */
  class SMIPDForceManager : public ForceManager{
    
  public:
    SMIPDForceManager(SimInfo * info);
    
  protected:
    virtual void postCalculation(bool needStress);
    
  private:
    std::vector<RealType> genTriangleForces(int nTriangles, RealType variance);

    Globals* simParams;
    Thermo* thermo;
    SeqRandNumGen randNumGen_;    
    Velocitizer* veloMunge;

    RealType dt_;
    RealType targetTemp_;
    RealType targetPressure_; 
    RealType thermalConductivity_;
    RealType thermalLength_;

    RealType variance_;

    Hull* surfaceMesh_;
    std::vector<StuntDouble*> localSites_;
  };
  
} //end namespace oopse
#endif //SIMPD_FORCEMANAGER

