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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
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
#include "nonbonded/InteractionManager.hpp"
#include "parallel/Decomposition.hpp"

#define PREPAIR_LOOP 0
#define PAIR_LOOP 1

using namespace std;
namespace OpenMD {
  /**
   * @class ForceManager ForceManager.hpp "brains/ForceManager.hpp"
   * ForceManager is responsible for calculating the short range
   * interactions and long range interactions.
   *
   * @note the reason we delay some of the setup is that some
   * applications (Dump2XYZ etc.) may not need force calculation, so
   * why bother?
   */
  class ForceManager {

  public:
    
    ForceManager(SimInfo * info);                          
    virtual ~ForceManager() {}
    virtual void calcForces();
    virtual void init() {};

  protected:

    virtual void preCalculation();        
    virtual void shortRangeInteractions();
    virtual void longRangeInteractions();
    virtual void postCalculation();
 
    SimInfo * info_;        
    map<Bend*, BendDataSet> bendDataSets;
    map<Torsion*, TorsionDataSet> torsionDataSets;
    map<Inversion*, InversionDataSet> inversionDataSets;
    Mat3x3d tau;

    vector<pair<int, int> > neighborList_;

    InteractionManager* interactionMan_;
    Decomposition* decomp_;
    SwitchingFunction* swfun_;
    vector<pair<int, int> > neighborList;
    map< pair<int, int>, pair<RealType, RealType> > groupCutoffMap;
    
  };
  
} 
#endif //BRAINS_FORCEMANAGER_HPP
