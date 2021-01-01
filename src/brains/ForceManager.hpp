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
#include "brains/Thermo.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

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
    virtual ~ForceManager();
    virtual void calcForces();
    virtual void calcSelectedForces(Molecule* mol1, Molecule* mol2);
    void setDoElectricField(bool def) { doElectricField_ = def; }
    void initialize();

  protected: 
    bool initialized_ {false}; 
    bool doParticlePot_ {false};
    bool doElectricField_ {false};
    bool doSitePotential_ {false};
    bool doHeatFlux_ {false};
    bool doLongRangeCorrections_ {false};
    bool usePeriodicBoundaryConditions_ {false};
    bool useSurfaceTerm_ {false};
    bool useSlabGeometry_ {false};
    int axis_ {false};

    virtual void setupCutoffs();
    virtual void preCalculation();        
    virtual void shortRangeInteractions();
    virtual void longRangeInteractions();
    virtual void postCalculation();

    virtual void selectedPreCalculation(Molecule* mol1, Molecule* mol2);        
    virtual void selectedShortRangeInteractions(Molecule* mol1, Molecule* mol2);
    virtual void selectedLongRangeInteractions(Molecule* mol1, Molecule* mol2);
    virtual void selectedPostCalculation(Molecule* mol1, Molecule* mol2);
    
    SimInfo* info_ {nullptr};        
    ForceField* forceField_ {nullptr};
    InteractionManager* interactionMan_ {nullptr};
    ForceDecomposition* fDecomp_ {nullptr};
    SwitchingFunction* switcher_ {nullptr};
    Thermo* thermo {nullptr};

    SwitchingFunctionType sft_;/**< Type of switching function in use */
    RealType rCut_;            /**< cutoff radius for non-bonded interactions */
    RealType rCutSq_;
    RealType rSwitch_;         /**< inner radius of switching function */
    CutoffMethod cutoffMethod_;/**< Cutoff Method for most non-bonded interactions */

    set<AtomType*> atomTypes_;
    vector<pair<AtomType*, AtomType*> > interactions_;
    map<Bend*, BendDataSet> bendDataSets;
    map<Torsion*, TorsionDataSet> torsionDataSets;
    map<Inversion*, InversionDataSet> inversionDataSets;
    //vector<pair<int, int> > neighborList_;
    vector<int> neighborList_;
    vector<int> point_;

    vector<RealType> vdwScale_;
    vector<RealType> electrostaticScale_;

    Mat3x3d virialTensor;

    vector<Perturbation*> perturbations_;

    bool doPotentialSelection_ {false};
    string selectionScript_;
    SelectionManager seleMan_;
    SelectionEvaluator evaluator_;
    // And all of the variables and structures for long range interactions:
    
    InteractionData idat;
    SelfData sdat;

  };
} 
#endif //BRAINS_FORCEMANAGER_HPP
