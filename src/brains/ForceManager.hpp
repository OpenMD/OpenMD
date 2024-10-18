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

/**
 * @file ForceManager.hpp
 * @author tlin
 * @date 11/09/2004
 * @version 1.0
 */

#ifndef BRAINS_FORCEMANAGER_HPP
#define BRAINS_FORCEMANAGER_HPP

#include "brains/ForceModifier.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "nonbonded/Cutoffs.hpp"
#include "nonbonded/InteractionManager.hpp"
#include "nonbonded/SwitchingFunction.hpp"
#include "parallel/ForceDecomposition.hpp"
#include "primitives/Molecule.hpp"
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
    ForceManager(SimInfo* info);
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

    SwitchingFunctionType sft_; /**< Type of switching function in use */
    RealType rCut_; /**< cutoff radius for non-bonded interactions */
    RealType rCutSq_;
    RealType rSwitch_; /**< inner radius of switching function */
    CutoffMethod
        cutoffMethod_; /**< Cutoff Method for most non-bonded interactions */

    AtomTypeSet atomTypes_;
    std::vector<pair<AtomType*, AtomType*>> interactions_;
    std::map<Bend*, BendDataSet> bendDataSets;
    std::map<Torsion*, TorsionDataSet> torsionDataSets;
    std::map<Inversion*, InversionDataSet> inversionDataSets;
    std::vector<int> neighborList_;
    std::vector<int> point_;
    std::vector<Vector3d> savedPositions_;

    std::vector<RealType> vdwScale_;
    std::vector<RealType> electrostaticScale_;

    Mat3x3d virialTensor;

    std::vector<ForceModifier*> forceModifiers_;

    bool doPotentialSelection_ {false};
    std::string selectionScript_;
    SelectionManager seleMan_;
    SelectionEvaluator evaluator_;

    // And all of the variables and structures for long range interactions:
    InteractionData idat;
    SelfData sdat;
  };
}  // namespace OpenMD

#endif  // BRAINS_FORCEMANAGER_HPP
