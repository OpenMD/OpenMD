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
 
#ifndef PARALLEL_FORCEMATRIXDECOMPOSITION_HPP
#define PARALLEL_FORCEMATRIXDECOMPOSITION_HPP

#include "parallel/ForceDecomposition.hpp"
#include "math/SquareMatrix3.hpp"
#include "brains/Snapshot.hpp"

#ifdef IS_MPI
#include "parallel/Communicator.hpp"
#endif

using namespace std;
namespace OpenMD {
  
  class ForceMatrixDecomposition : public ForceDecomposition {
  public:
    ForceMatrixDecomposition(SimInfo* info, InteractionManager* iMan);
    ~ForceMatrixDecomposition();

    void distributeInitialData();
    void zeroWorkArrays();
    void distributeData();
    void collectIntermediateData();
    void distributeIntermediateData();
    void collectSelfData();
    void collectData();

    // neighbor list routines
    void buildNeighborList(vector<int>& neighborList, vector<int>& point);

    // group bookkeeping
    Vector3d& getGroupVelocityColumn(int cg2);

    // Group->atom bookkeeping
    vector<int>& getAtomsInGroupRow(int cg1);
    vector<int>& getAtomsInGroupColumn(int cg2);
    Vector3d getAtomToGroupVectorRow(int atom1, int cg1);
    Vector3d getAtomToGroupVectorColumn(int atom2, int cg2);
    RealType& getMassFactorRow(int atom1);
    RealType& getMassFactorColumn(int atom2);

    // spatial data
    Vector3d getIntergroupVector(int cg1, int cg2);
    Vector3d getInteratomicVector(int atom1, int atom2);
       
    // atom bookkeeping
    int& getNAtomsInRow();
    int getTopologicalDistance(int atom1, int atom2);
    vector<int>& getExcludesForAtom(int atom1); 
    bool skipAtomPair(int atom1, int atom2, int cg1, int cg2);
    bool excludeAtomPair(int atom1, int atom2);
    int getGlobalIDRow(int atom1);
    int getGlobalIDCol(int atom1);
    int getGlobalID(int atom1);
    void addForceToAtomRow(int atom1, Vector3d fg);
    void addForceToAtomColumn(int atom2, Vector3d fg);
    Vector3d& getAtomVelocityColumn(int atom2);

    // filling interaction blocks with pointers
    void fillInteractionData(InteractionData &idat, int atom1, int atom2, bool newAtom1 = true);
    void unpackInteractionData(InteractionData &idat, int atom1, int atom2);
    void unpackPrePairData(InteractionData &idat, int atom1, int atom2);

  private:     
    int nLocal_;
    int nGroups_;
    vector<int> AtomLocalToGlobal;
    vector<int> cgLocalToGlobal;
    vector<RealType> groupCutoff;
    vector<int> groupToGtype;

#ifdef IS_MPI    
    DataStorage atomRowData;
    DataStorage atomColData;
    DataStorage cgRowData;
    DataStorage cgColData;

    int nAtomsInRow_; 
    int nAtomsInCol_;
    int nGroupsInRow_;
    int nGroupsInCol_;

    Communicator<Row> rowComm;
    Communicator<Column> colComm;

    Plan<int>* AtomPlanIntRow {nullptr};
    Plan<RealType>* AtomPlanRealRow {nullptr}; 
    Plan<Vector3d>* AtomPlanVectorRow {nullptr}; 
    Plan<Mat3x3d>*  AtomPlanMatrixRow {nullptr}; 
    Plan<potVec>* AtomPlanPotRow {nullptr};

    Plan<int>* AtomPlanIntColumn {nullptr};
    Plan<RealType>* AtomPlanRealColumn {nullptr}; 
    Plan<Vector3d>* AtomPlanVectorColumn {nullptr}; 
    Plan<Mat3x3d>*  AtomPlanMatrixColumn {nullptr};
    Plan<potVec>* AtomPlanPotColumn {nullptr}; 

    Plan<int>* cgPlanIntRow {nullptr};
    Plan<Vector3d>* cgPlanVectorRow {nullptr}; 
    Plan<int>* cgPlanIntColumn {nullptr};
    Plan<Vector3d>* cgPlanVectorColumn {nullptr}; 

    // work arrays for assembling potential energy
    vector<potVec> pot_row;
    vector<potVec> pot_col;

    vector<potVec> expot_row;
    vector<potVec> expot_col;

    vector<potVec> selepot_row;
    vector<potVec> selepot_col;

    vector<int> identsRow;
    vector<int> identsCol;

    vector<int> regionsRow;
    vector<int> regionsCol;

    vector<AtomType*> atypesRow;
    vector<AtomType*> atypesCol;

    vector<int> AtomRowToGlobal;
    vector<int> AtomColToGlobal;

public:
    vector<int> cgRowToGlobal;
    vector<int> cgColToGlobal;

private:
    vector<vector<int> > cellListRow_;
    vector<vector<int> > cellListCol_;

    vector<vector<int> > groupListRow_;
    vector<vector<int> > groupListCol_;

    vector<RealType> massFactorsRow;
    vector<RealType> massFactorsCol;

    vector<int> regionRow;
    vector<int> regionCol;
#endif

  };

}
#endif

