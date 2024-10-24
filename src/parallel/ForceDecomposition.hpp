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

#ifndef PARALLEL_FORCEDECOMPOSITION_HPP
#define PARALLEL_FORCEDECOMPOSITION_HPP

#include "brains/SimInfo.hpp"
#include "brains/SnapshotManager.hpp"
#include "nonbonded/InteractionManager.hpp"
#include "nonbonded/NonBondedInteraction.hpp"

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

    virtual void setSnapshot(Snapshot* snap) { snap_ = snap; }

    virtual void distributeInitialData()      = 0;
    virtual void distributeData()             = 0;
    virtual void zeroWorkArrays()             = 0;
    virtual void collectIntermediateData()    = 0;
    virtual void distributeIntermediateData() = 0;
    virtual void collectData()                = 0;
    virtual void collectSelfData()            = 0;
    virtual potVec getSelfPotential() { return selfPot; }
    virtual potVec getPairwisePotential() { return pairwisePot; }
    virtual potVec getExcludedPotential() { return excludedPot; }
    virtual potVec getSelectedPotential() { return selectedPot; }
    virtual potVec getExcludedSelfPotential() { return excludedSelfPot; }
    virtual potVec getSelectedSelfPotential() { return selectedSelfPot; }

    // neighbor list routines
    virtual bool checkNeighborList(vector<Vector3d> savedPositions);
    virtual void buildNeighborList(vector<int>& neighborList,
                                   vector<int>& point,
                                   vector<Vector3d>& savedPositions) = 0;

    void setCutoffRadius(RealType rCut);

    // group bookkeeping
    virtual Vector3d& getGroupVelocityColumn(int atom2) = 0;

    // Group->atom bookkeeping
    virtual vector<int>& getAtomsInGroupRow(int cg1)    = 0;
    virtual vector<int>& getAtomsInGroupColumn(int cg2) = 0;

    virtual Vector3d getAtomToGroupVectorRow(int atom1, int cg1)    = 0;
    virtual Vector3d getAtomToGroupVectorColumn(int atom2, int cg2) = 0;
    virtual RealType& getMassFactorRow(int atom1)                   = 0;
    virtual RealType& getMassFactorColumn(int atom2)                = 0;

    // spatial data
    virtual Vector3d getIntergroupVector(int cg1, int cg2)      = 0;
    virtual Vector3d getInteratomicVector(int atom1, int atom2) = 0;

    // atom bookkeeping
    virtual int& getNAtomsInRow()                                     = 0;
    virtual vector<int>& getExcludesForAtom(int atom1)                = 0;
    virtual bool skipAtomPair(int atom1, int atom2, int cg1, int cg2) = 0;
    virtual bool excludeAtomPair(int atom1, int atom2)                = 0;
    virtual int getGlobalIDRow(int atom1)                             = 0;
    virtual int getGlobalIDCol(int atom2)                             = 0;
    virtual int getGlobalID(int atom1)                                = 0;

    virtual int getTopologicalDistance(int atom1, int atom2)  = 0;
    virtual void addForceToAtomRow(int atom1, Vector3d fg)    = 0;
    virtual void addForceToAtomColumn(int atom2, Vector3d fg) = 0;
    virtual Vector3d& getAtomVelocityColumn(int atom2)        = 0;

    // filling & unpacking interaction data from Snapshot or Row/Col data
    virtual void fillInteractionData(InteractionData& idat, int atom1,
                                     int atom2, bool newAtom1 = true) = 0;
    virtual void unpackPrePairData(InteractionData& idat, int atom1,
                                   int atom2)                         = 0;
    virtual void unpackInteractionData(InteractionData& idat, int atom1,
                                       int atom2)                     = 0;

    // filling & unpacking self data from Snapshot data
    virtual void fillSelfData(SelfData& sdat, int atom);
    virtual void unpackSelfData(SelfData& sdat, int atom);

    virtual void fillPreForceData(SelfData& sdat, int atom);
    virtual void unpackPreForceData(SelfData& sdat, int atom);

    virtual void addToHeatFlux(Vector3d hf);
    virtual void setHeatFlux(Vector3d hf);

  protected:
    SimInfo* info_ {nullptr};
    SnapshotManager* sman_;
    Snapshot* snap_;
    ForceField* ff_;
    InteractionManager* interactionMan_;

    int atomStorageLayout_;
    int rigidBodyStorageLayout_;
    int cutoffGroupStorageLayout_;
    bool needVelocities_;
    bool usePeriodicBoundaryConditions_;
    RealType skinThickness_; /**< Verlet neighbor list skin thickness */
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
    vector<vector<int>> toposForAtom;
    vector<vector<int>> topoDist;
    vector<vector<int>> excludesForAtom;
    vector<vector<int>> groupList_;
    vector<RealType> massFactors;
    vector<AtomType*> atypesLocal;

    vector<Vector3i> cellOffsets_;
    Vector3i nCells_;
    vector<vector<int>> cellList_;
  };
}  // namespace OpenMD

#endif
