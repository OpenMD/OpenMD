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
 * @file ForceManager.hpp
 * @author tlin
 * @date 11/09/2004
 * @version 1.0
 */

#ifndef BRAINS_FORCEMANAGER_HPP
#define BRAINS_FORCEMANAGER_HPP

#include "brains/SimInfo.hpp"
#include "primitives/Molecule.hpp"
#include "nonbonded/Cutoffs.hpp"
#include "nonbonded/SwitchingFunction.hpp"
#include "nonbonded/InteractionManager.hpp"
#include "perturbations/Perturbation.hpp"
#include "parallel/ForceDecomposition.hpp"

#define PREPAIR_LOOP 0
#define PAIR_LOOP 1

using namespace std;
namespace OpenMD {
  /**
   * @class ForceManager ForceManager.hpp "brains/ForceManager.hpp"
   * ForceManager is responsible for calculating both the short range
   * (bonded) interactions and long range (non-bonded) interactions.
   *
   * @note the reason we delay some of the setup is that
   * initialization must wait until after the force field has been
   * parsed so that the atom types are known.
   */
  class ForceManager {

  public:
    ForceManager(SimInfo * info);                          
    virtual ~ForceManager() {}
    virtual void calcForces();
    void initialize();

  protected: 
    bool initialized_; 
    bool doParticlePot_;
    bool doElectricField_;
    bool doHeatFlux_;
    bool doLongRangeCorrections_;

    virtual void setupCutoffs();
    virtual void preCalculation();        
    virtual void shortRangeInteractions();
    virtual void longRangeInteractions();
    virtual void postCalculation();

    SimInfo* info_;        
    ForceField* forceField_;
    InteractionManager* interactionMan_;
    ForceDecomposition* fDecomp_;
    SwitchingFunction* switcher_;

    SwitchingFunctionType sft_;/**< Type of switching function in use */
    RealType rCut_;            /**< cutoff radius for non-bonded interactions */
    RealType rSwitch_;         /**< inner radius of switching function */
    CutoffMethod cutoffMethod_;/**< Cutoff Method for most non-bonded interactions */
    CutoffPolicy cutoffPolicy_;/**< Cutoff Policy for non-bonded interactions */

    set<AtomType*> atomTypes_;
    vector<pair<AtomType*, AtomType*> > interactions_;
    map<Bend*, BendDataSet> bendDataSets;
    map<Torsion*, TorsionDataSet> torsionDataSets;
    map<Inversion*, InversionDataSet> inversionDataSets;
    vector<pair<int, int> > neighborList;

    vector<RealType> vdwScale_;
    vector<RealType> electrostaticScale_;

    Mat3x3d stressTensor;

    vector<Perturbation*> perturbations_;
  };
} 
#endif //BRAINS_FORCEMANAGER_HPP
