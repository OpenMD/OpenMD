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
 
#ifndef PARALLEL_FORCEDECOMPOSITION_HPP
#define PARALLEL_FORCEDECOMPOSITION_HPP

#include "brains/SimInfo.hpp"
#include "brains/SnapshotManager.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "nonbonded/InteractionManager.hpp"
#include "utils/Tuple.hpp"

using namespace std;
namespace OpenMD { 

  /**
   * @class ForceDecomposition 
   *
   * ForceDecomposition is an interface for passing out and collecting
   * information from many processors at various stages of the main
   * non-bonded ForceLoop.
   *
   * The pairwise force calculation has an outer-running loop (the "I"
   * loop) and an inner-running loop (the "J" loop).  In parallel
   * decompositions, these loop over different groups of atoms on
   * different processors.  Between each set of computations on the
   * local processor, data must be exchanged among the processors.
   * This can happen at different times in the calculation:
   *
   *  distributeInitialData      (parallel communication - one time only)
   *  distributeData             (parallel communication - every ForceLoop)
   *
   *  loop iLoop over nLoops     (nLoops may be 1, 2, or until self consistent)
   *  |  loop over i
   *  |  | loop over j
   *  |  | | localComputation
   *  |  |  end
   *  |  end
   *  |  if (nLoops > 1):
   *  |  |   collectIntermediateData    (parallel communication)
   *  |  |   distributeIntermediateData (parallel communication)
   *  |  endif
   *  end
   * collectData                        (parallel communication)
   * loop over i
   * | localComputation
   * end
   * collectSelfData                    (parallel communication)
   *
   * ForceDecomposition provides the interface for ForceLoop to do the
   * communication steps and to iterate using the correct set of atoms
   * and cutoff groups.
   */
  class ForceDecomposition {
  public:

    ForceDecomposition(SimInfo* info, InteractionManager* iMan);
    virtual ~ForceDecomposition() {}

    virtual void setSnapshot(Snapshot* snap) {snap_ = snap;}
    
    virtual void distributeInitialData() = 0;
    virtual void distributeData() = 0;
    virtual void zeroWorkArrays() = 0;
    virtual void collectIntermediateData() = 0;
    virtual void distributeIntermediateData() = 0;
    virtual void collectData() = 0;
    virtual void collectSelfData() = 0;
    virtual potVec getSelfPotential() { return selfPot; }
    virtual potVec getPairwisePotential() { return pairwisePot; }
    virtual potVec getExcludedPotential() { return excludedPot; }
    virtual potVec getSelectedPotential() { return selectedPot; }
    virtual potVec getExcludedSelfPotential() { return excludedSelfPot; }
    virtual potVec getSelectedSelfPotential() { return selectedSelfPot; }

    // neighbor list routines
    virtual bool checkNeighborList();
    virtual void buildNeighborList(vector<int>& neighborList, vector<int>& point) = 0;

    void setCutoffRadius(RealType rCut);
    
    // group bookkeeping
    virtual Vector3d& getGroupVelocityColumn(int atom2) = 0;

    // Group->atom bookkeeping
    virtual vector<int>& getAtomsInGroupRow(int cg1) = 0; 
    virtual vector<int>& getAtomsInGroupColumn(int cg2) = 0;

    virtual Vector3d getAtomToGroupVectorRow(int atom1, int cg1) = 0;
    virtual Vector3d getAtomToGroupVectorColumn(int atom2, int cg2) = 0;
    virtual RealType& getMassFactorRow(int atom1) = 0;
    virtual RealType& getMassFactorColumn(int atom2) = 0;

    // spatial data
    virtual Vector3d getIntergroupVector(int cg1, int cg2) = 0;
    virtual Vector3d getInteratomicVector(int atom1, int atom2) = 0;
       
    // atom bookkeeping
    virtual int& getNAtomsInRow() = 0;
    virtual vector<int>& getExcludesForAtom(int atom1) = 0;
    virtual bool skipAtomPair(int atom1, int atom2, int cg1, int cg2) = 0;
    virtual bool excludeAtomPair(int atom1, int atom2) = 0;
    virtual int getGlobalIDRow(int atom1) = 0;
    virtual int getGlobalIDCol(int atom2) = 0;
    virtual int getGlobalID(int atom1) = 0;
    
    virtual int getTopologicalDistance(int atom1, int atom2) = 0;
    virtual void addForceToAtomRow(int atom1, Vector3d fg) = 0;
    virtual void addForceToAtomColumn(int atom2, Vector3d fg) = 0;
    virtual Vector3d& getAtomVelocityColumn(int atom2) = 0;

    // filling & unpacking interaction data from Snapshot or Row/Col data
    virtual void fillInteractionData(InteractionData &idat, int atom1,
				     int atom2, bool newAtom1 = true) = 0;
    virtual void unpackPrePairData(InteractionData &idat, int atom1,
				   int atom2) = 0;
    virtual void unpackInteractionData(InteractionData &idat, int atom1,
				       int atom2) = 0;

    // filling & unpacking self data from Snapshot data
    virtual void fillSelfData(SelfData &sdat, int atom);
    virtual void unpackSelfData(SelfData &sdat, int atom);

    virtual void fillPreForceData(SelfData &sdat, int atom);
    virtual void unpackPreForceData(SelfData &sdat, int atom);
    
    virtual void addToHeatFlux(Vector3d hf);
    virtual void setHeatFlux(Vector3d hf);
    
  protected:
    SimInfo* info_ {nullptr};   
    SnapshotManager* sman_;    
    Snapshot* snap_;
    ForceField* ff_;
    InteractionManager* interactionMan_;

    int storageLayout_;
    bool needVelocities_;
    bool usePeriodicBoundaryConditions_;
    RealType skinThickness_;   /**< Verlet neighbor list skin thickness */
    RealType rCut_;
    RealType rList_;
    RealType rListSq_;

    vector<int> idents;
    vector<int> regions;
    potVec pairwisePot;
    potVec selfPot;
    potVec excludedPot;
    potVec excludedSelfPot;
    potVec selectedPot;
    potVec selectedSelfPot;

    /** 
     * The topological distance between two atomic sites is handled
     * via two vector structures for speed.  These structures agnostic
     * regarding the parallel decomposition.  The index for
     * toposForAtom could be local or row, while the values could be
     * local or column.  It will be up to the specific decomposition
     * method to fill these.
     */
    vector<vector<int> > toposForAtom; 
    vector<vector<int> > topoDist;                                       
    vector<vector<int> > excludesForAtom;
    vector<vector<int> > groupList_;
    vector<RealType> massFactors;
    vector<AtomType*> atypesLocal;

    vector<Vector3i> cellOffsets_;
    Vector3i nCells_;
    vector<vector<int> > cellList_;
    vector<Vector3d> saved_CG_positions_;
  };    
}
#endif
