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
 * @file ForceManager.hpp
 * @author tlin
 * @date 11/09/2004
 * @time 10:36am
 * @version 1.0
 */

#ifndef BRAINS_FORCEMANAGER_HPP
#define BRAINS_FORCEMANAGER_HPP

#include "brains/SimInfo.hpp"
#include "primitives/Molecule.hpp"
namespace oopse {
  /**
   * @class ForceManager ForceManager.hpp "brains/ForceManager.hpp"
   * ForceManager is responsible for calculating the short range
   * interactions (C++) and long range interactions (Fortran). If the
   * Fortran side is not set up before the force calculation, call
   * SimInfo's update function to settle it down.
   *
   * @note the reason we delay fortran side's setup is that some
   * applications (Dump2XYZ etc.) may not need force calculation, so why
   * bother?
   */
  class ForceManager {

  public:
    ForceManager(SimInfo * info) : info_(info) {}
        
    virtual ~ForceManager() {}

    // public virtual functions should be avoided
    /**@todo needs refactoring */
    virtual void calcForces(bool needPotential, bool needStress);

    virtual void init() {}
  protected:

    virtual void preCalculation();
        
    virtual void calcShortRangeInteraction();

    virtual void calcLongRangeInteraction(bool needPotential, bool needStress);

    virtual void postCalculation(bool needStress);
 
    SimInfo * info_;        

    std::map<Bend*, BendDataSet> bendDataSets;
    std::map<Torsion*, TorsionDataSet> torsionDataSets;
    Mat3x3d tau;
    
  };

} //end namespace oopse
#endif //BRAINS_FORCEMANAGER_HPP
