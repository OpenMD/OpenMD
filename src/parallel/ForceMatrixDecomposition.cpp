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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#include "parallel/ForceMatrixDecomposition.hpp"
#include "math/SquareMatrix3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "brains/SnapshotManager.hpp"
#include "brains/PairList.hpp"

using namespace std;
namespace OpenMD {

  ForceMatrixDecomposition::ForceMatrixDecomposition(SimInfo* info, InteractionManager* iMan) : ForceDecomposition(info, iMan) {

    // In a parallel computation, row and colum scans must visit all
    // surrounding cells (not just the 14 upper triangular blocks that
    // are used when the processor can see all pairs)
#ifdef IS_MPI
    cellOffsets_.clear();
    cellOffsets_.push_back( Vector3i(-1,-1,-1) );
    cellOffsets_.push_back( Vector3i( 0,-1,-1) );
    cellOffsets_.push_back( Vector3i( 1,-1,-1) );                          
    cellOffsets_.push_back( Vector3i(-1, 0,-1) );
    cellOffsets_.push_back( Vector3i( 0, 0,-1) );
    cellOffsets_.push_back( Vector3i( 1, 0,-1) );
    cellOffsets_.push_back( Vector3i(-1, 1,-1) );
    cellOffsets_.push_back( Vector3i( 0, 1,-1) );      
    cellOffsets_.push_back( Vector3i( 1, 1,-1) );
    cellOffsets_.push_back( Vector3i(-1,-1, 0) );
    cellOffsets_.push_back( Vector3i( 0,-1, 0) );
    cellOffsets_.push_back( Vector3i( 1,-1, 0) );
    cellOffsets_.push_back( Vector3i(-1, 0, 0) );       
    cellOffsets_.push_back( Vector3i( 0, 0, 0) );
    cellOffsets_.push_back( Vector3i( 1, 0, 0) );
    cellOffsets_.push_back( Vector3i(-1, 1, 0) );
    cellOffsets_.push_back( Vector3i( 0, 1, 0) );
    cellOffsets_.push_back( Vector3i( 1, 1, 0) );
    cellOffsets_.push_back( Vector3i(-1,-1, 1) );
    cellOffsets_.push_back( Vector3i( 0,-1, 1) );
    cellOffsets_.push_back( Vector3i( 1,-1, 1) );
    cellOffsets_.push_back( Vector3i(-1, 0, 1) );
    cellOffsets_.push_back( Vector3i( 0, 0, 1) );
    cellOffsets_.push_back( Vector3i( 1, 0, 1) );
    cellOffsets_.push_back( Vector3i(-1, 1, 1) );
    cellOffsets_.push_back( Vector3i( 0, 1, 1) );
    cellOffsets_.push_back( Vector3i( 1, 1, 1) );
#endif    
  }


  /**
   * distributeInitialData is essentially a copy of the older fortran 
   * SimulationSetup
   */
  void ForceMatrixDecomposition::distributeInitialData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
    ff_ = info_->getForceField();
    nLocal_ = snap_->getNumberOfAtoms();
   
    nGroups_ = info_->getNLocalCutoffGroups();
    // gather the information for atomtype IDs (atids):
    idents = info_->getIdentArray();
    regions = info_->getRegions();
    AtomLocalToGlobal = info_->getGlobalAtomIndices();
    cgLocalToGlobal = info_->getGlobalGroupIndices();
    vector<int> globalGroupMembership = info_->getGlobalGroupMembership();

    massFactors = info_->getMassFactors();

    PairList* excludes = info_->getExcludedInteractions();
    PairList* oneTwo = info_->getOneTwoInteractions();
    PairList* oneThree = info_->getOneThreeInteractions();
    PairList* oneFour = info_->getOneFourInteractions();
    
    if (needVelocities_) 
      snap_->cgData.setStorageLayout(DataStorage::dslPosition | 
                                     DataStorage::dslVelocity);
    else 
      snap_->cgData.setStorageLayout(DataStorage::dslPosition);
    
#ifdef IS_MPI
 
    MPI_Comm row = rowComm.getComm();
    MPI_Comm col = colComm.getComm();

    AtomPlanIntRow = new Plan<int>(row, nLocal_);
    AtomPlanRealRow = new Plan<RealType>(row, nLocal_);
    AtomPlanVectorRow = new Plan<Vector3d>(row, nLocal_);
    AtomPlanMatrixRow = new Plan<Mat3x3d>(row, nLocal_);
    AtomPlanPotRow = new Plan<potVec>(row, nLocal_);

    AtomPlanIntColumn = new Plan<int>(col, nLocal_);
    AtomPlanRealColumn = new Plan<RealType>(col, nLocal_);
    AtomPlanVectorColumn = new Plan<Vector3d>(col, nLocal_);
    AtomPlanMatrixColumn = new Plan<Mat3x3d>(col, nLocal_);
    AtomPlanPotColumn = new Plan<potVec>(col, nLocal_);

    cgPlanIntRow = new Plan<int>(row, nGroups_);
    cgPlanVectorRow = new Plan<Vector3d>(row, nGroups_);
    cgPlanIntColumn = new Plan<int>(col, nGroups_);
    cgPlanVectorColumn = new Plan<Vector3d>(col, nGroups_);

    nAtomsInRow_ = AtomPlanIntRow->getSize();
    nAtomsInCol_ = AtomPlanIntColumn->getSize();
    nGroupsInRow_ = cgPlanIntRow->getSize();
    nGroupsInCol_ = cgPlanIntColumn->getSize();

    // Modify the data storage objects with the correct layouts and sizes:
    atomRowData.resize(nAtomsInRow_);
    atomRowData.setStorageLayout(storageLayout_);
    atomColData.resize(nAtomsInCol_);
    atomColData.setStorageLayout(storageLayout_);
    cgRowData.resize(nGroupsInRow_);
    cgRowData.setStorageLayout(DataStorage::dslPosition);
    cgColData.resize(nGroupsInCol_);
    if (needVelocities_)
      // we only need column velocities if we need them.
      cgColData.setStorageLayout(DataStorage::dslPosition |
                                 DataStorage::dslVelocity);
    else     
      cgColData.setStorageLayout(DataStorage::dslPosition);
      
    identsRow.resize(nAtomsInRow_);
    identsCol.resize(nAtomsInCol_);
    
    AtomPlanIntRow->gather(idents, identsRow);
    AtomPlanIntColumn->gather(idents, identsCol);

    regionsRow.resize(nAtomsInRow_);
    regionsCol.resize(nAtomsInCol_);
    
    AtomPlanIntRow->gather(regions, regionsRow);
    AtomPlanIntColumn->gather(regions, regionsCol);
    
    // allocate memory for the parallel objects
    atypesRow.resize(nAtomsInRow_);
    atypesCol.resize(nAtomsInCol_);

    for (int i = 0; i < nAtomsInRow_; i++) 
      atypesRow[i] = ff_->getAtomType(identsRow[i]);
    for (int i = 0; i < nAtomsInCol_; i++) 
      atypesCol[i] = ff_->getAtomType(identsCol[i]);         

    pot_row.resize(nAtomsInRow_);
    pot_col.resize(nAtomsInCol_);

    expot_row.resize(nAtomsInRow_);
    expot_col.resize(nAtomsInCol_);

    AtomRowToGlobal.resize(nAtomsInRow_);
    AtomColToGlobal.resize(nAtomsInCol_);
    AtomPlanIntRow->gather(AtomLocalToGlobal, AtomRowToGlobal);
    AtomPlanIntColumn->gather(AtomLocalToGlobal, AtomColToGlobal);

    cgRowToGlobal.resize(nGroupsInRow_);
    cgColToGlobal.resize(nGroupsInCol_);
    cgPlanIntRow->gather(cgLocalToGlobal, cgRowToGlobal);
    cgPlanIntColumn->gather(cgLocalToGlobal, cgColToGlobal);

    massFactorsRow.resize(nAtomsInRow_);
    massFactorsCol.resize(nAtomsInCol_);
    AtomPlanRealRow->gather(massFactors, massFactorsRow);
    AtomPlanRealColumn->gather(massFactors, massFactorsCol);

    groupListRow_.clear();
    groupListRow_.resize(nGroupsInRow_);
    for (int i = 0; i < nGroupsInRow_; i++) {
      int gid = cgRowToGlobal[i];
      for (int j = 0; j < nAtomsInRow_; j++) {
        int aid = AtomRowToGlobal[j];
        if (globalGroupMembership[aid] == gid)
          groupListRow_[i].push_back(j);
      }      
    }

    groupListCol_.clear();
    groupListCol_.resize(nGroupsInCol_);
    for (int i = 0; i < nGroupsInCol_; i++) {
      int gid = cgColToGlobal[i];
      for (int j = 0; j < nAtomsInCol_; j++) {
        int aid = AtomColToGlobal[j];
        if (globalGroupMembership[aid] == gid)
          groupListCol_[i].push_back(j);
      }      
    }

    excludesForAtom.clear();
    excludesForAtom.resize(nAtomsInRow_);
    toposForAtom.clear();
    toposForAtom.resize(nAtomsInRow_);
    topoDist.clear();
    topoDist.resize(nAtomsInRow_);
    for (int i = 0; i < nAtomsInRow_; i++) {
      int iglob = AtomRowToGlobal[i];

      for (int j = 0; j < nAtomsInCol_; j++) {
        int jglob = AtomColToGlobal[j];

        if (excludes->hasPair(iglob, jglob)) 
          excludesForAtom[i].push_back(j);       
        
        if (oneTwo->hasPair(iglob, jglob)) {
          toposForAtom[i].push_back(j);
          topoDist[i].push_back(1);
        } else {
          if (oneThree->hasPair(iglob, jglob)) {
            toposForAtom[i].push_back(j);
            topoDist[i].push_back(2);
          } else {
            if (oneFour->hasPair(iglob, jglob)) {
              toposForAtom[i].push_back(j);
              topoDist[i].push_back(3);
            }
          }
        }
      }      
    }

#else
    excludesForAtom.clear();
    excludesForAtom.resize(nLocal_);
    toposForAtom.clear();
    toposForAtom.resize(nLocal_);
    topoDist.clear();
    topoDist.resize(nLocal_);

    for (int i = 0; i < nLocal_; i++) {
      int iglob = AtomLocalToGlobal[i];

      for (int j = 0; j < nLocal_; j++) {
        int jglob = AtomLocalToGlobal[j];

        if (excludes->hasPair(iglob, jglob)) 
          excludesForAtom[i].push_back(j);              
        
        if (oneTwo->hasPair(iglob, jglob)) {
          toposForAtom[i].push_back(j);
          topoDist[i].push_back(1);
        } else {
          if (oneThree->hasPair(iglob, jglob)) {
            toposForAtom[i].push_back(j);
            topoDist[i].push_back(2);
          } else {
            if (oneFour->hasPair(iglob, jglob)) {
              toposForAtom[i].push_back(j);
              topoDist[i].push_back(3);
            }
          }
        }
      }      
    }
#endif

    // allocate memory for the parallel objects
    atypesLocal.resize(nLocal_);

    for (int i = 0; i < nLocal_; i++) 
      atypesLocal[i] = ff_->getAtomType(idents[i]);

    groupList_.clear();
    groupList_.resize(nGroups_);
    for (int i = 0; i < nGroups_; i++) {
      int gid = cgLocalToGlobal[i];
      for (int j = 0; j < nLocal_; j++) {
        int aid = AtomLocalToGlobal[j];
        if (globalGroupMembership[aid] == gid) {
          groupList_[i].push_back(j);
        }
      }      
    }


    createGtypeCutoffMap();

  }
   
  void ForceMatrixDecomposition::createGtypeCutoffMap() {
    
    GrCut.clear();
    GrCutSq.clear();
    GrlistSq.clear();

    RealType tol = 1e-6;
    largestRcut_ = 0.0;
    int atid;
    set<AtomType*> atypes = info_->getSimulatedAtomTypes();
    
    map<int, RealType> atypeCutoff;
      
    for (set<AtomType*>::iterator at = atypes.begin(); 
         at != atypes.end(); ++at){
      atid = (*at)->getIdent();
      if (userChoseCutoff_) 
        atypeCutoff[atid] = userCutoff_;
      else
        atypeCutoff[atid] = interactionMan_->getSuggestedCutoffRadius(*at);
    }
    
    vector<RealType> gTypeCutoffs;
    // first we do a single loop over the cutoff groups to find the
    // largest cutoff for any atypes present in this group.
#ifdef IS_MPI
    vector<RealType> groupCutoffRow(nGroupsInRow_, 0.0);
    groupRowToGtype.resize(nGroupsInRow_);
    for (int cg1 = 0; cg1 < nGroupsInRow_; cg1++) {
      vector<int> atomListRow = getAtomsInGroupRow(cg1);
      for (vector<int>::iterator ia = atomListRow.begin(); 
           ia != atomListRow.end(); ++ia) {            
        int atom1 = (*ia);
        atid = identsRow[atom1];
        if (atypeCutoff[atid] > groupCutoffRow[cg1]) {
          groupCutoffRow[cg1] = atypeCutoff[atid];
        }
      }

      bool gTypeFound = false;
      for (int gt = 0; gt < gTypeCutoffs.size(); gt++) {
        if (abs(groupCutoffRow[cg1] - gTypeCutoffs[gt]) < tol) {
          groupRowToGtype[cg1] = gt;
          gTypeFound = true;
        } 
      }
      if (!gTypeFound) {
        gTypeCutoffs.push_back( groupCutoffRow[cg1] );
        groupRowToGtype[cg1] = gTypeCutoffs.size() - 1;
      }
      
    }
    vector<RealType> groupCutoffCol(nGroupsInCol_, 0.0);
    groupColToGtype.resize(nGroupsInCol_);
    for (int cg2 = 0; cg2 < nGroupsInCol_; cg2++) {
      vector<int> atomListCol = getAtomsInGroupColumn(cg2);
      for (vector<int>::iterator jb = atomListCol.begin(); 
           jb != atomListCol.end(); ++jb) {            
        int atom2 = (*jb);
        atid = identsCol[atom2];
        if (atypeCutoff[atid] > groupCutoffCol[cg2]) {
          groupCutoffCol[cg2] = atypeCutoff[atid];
        }
      }
      bool gTypeFound = false;
      for (int gt = 0; gt < gTypeCutoffs.size(); gt++) {
        if (abs(groupCutoffCol[cg2] - gTypeCutoffs[gt]) < tol) {
          groupColToGtype[cg2] = gt;
          gTypeFound = true;
        } 
      }
      if (!gTypeFound) {
        gTypeCutoffs.push_back( groupCutoffCol[cg2] );
        groupColToGtype[cg2] = gTypeCutoffs.size() - 1;
      }
    }
#else

    vector<RealType> groupCutoff(nGroups_, 0.0);
    groupToGtype.resize(nGroups_);
    for (int cg1 = 0; cg1 < nGroups_; cg1++) {
      groupCutoff[cg1] = 0.0;
      vector<int> atomList = getAtomsInGroupRow(cg1);
      for (vector<int>::iterator ia = atomList.begin(); 
           ia != atomList.end(); ++ia) {            
        int atom1 = (*ia);
        atid = idents[atom1];
        if (atypeCutoff[atid] > groupCutoff[cg1]) 
          groupCutoff[cg1] = atypeCutoff[atid];
      }
      
      bool gTypeFound = false;
      for (unsigned int gt = 0; gt < gTypeCutoffs.size(); gt++) {
        if (abs(groupCutoff[cg1] - gTypeCutoffs[gt]) < tol) {
          groupToGtype[cg1] = gt;
          gTypeFound = true;
        } 
      }
      if (!gTypeFound) {      
        gTypeCutoffs.push_back( groupCutoff[cg1] );
        groupToGtype[cg1] = gTypeCutoffs.size() - 1;
      }      
    }
#endif

    // Now we find the maximum group cutoff value present in the simulation

    RealType groupMax = *max_element(gTypeCutoffs.begin(), 
                                     gTypeCutoffs.end());

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &groupMax, 1, MPI_REALTYPE, 
                  MPI_MAX, MPI_COMM_WORLD);
#endif
    
    RealType tradRcut = groupMax;

    GrCut.resize( gTypeCutoffs.size() );
    GrCutSq.resize( gTypeCutoffs.size() );
    GrlistSq.resize( gTypeCutoffs.size() );


    for (unsigned int i = 0; i < gTypeCutoffs.size();  i++) {
      GrCut[i].resize( gTypeCutoffs.size() , 0.0);
      GrCutSq[i].resize( gTypeCutoffs.size(), 0.0 );
      GrlistSq[i].resize( gTypeCutoffs.size(), 0.0 );

      for (unsigned int j = 0; j < gTypeCutoffs.size();  j++) {       
        RealType thisRcut;
        switch(cutoffPolicy_) {
        case TRADITIONAL:
          thisRcut = tradRcut;
          break;
        case MIX:
          thisRcut = 0.5 * (gTypeCutoffs[i] + gTypeCutoffs[j]);
          break;
        case MAX:
          thisRcut = max(gTypeCutoffs[i], gTypeCutoffs[j]);
          break;
        default:
          sprintf(painCave.errMsg,
                  "ForceMatrixDecomposition::createGtypeCutoffMap " 
                  "hit an unknown cutoff policy!\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();
          break;
        }

        GrCut[i][j] = thisRcut;
        if (thisRcut > largestRcut_) largestRcut_ = thisRcut;
        GrCutSq[i][j] = thisRcut * thisRcut;
        GrlistSq[i][j] = pow(thisRcut + skinThickness_, 2);

        // pair<int,int> key = make_pair(i,j);
        // gTypeCutoffMap[key].first = thisRcut;
        // gTypeCutoffMap[key].third = pow(thisRcut + skinThickness_, 2);
        // sanity check
        
        if (userChoseCutoff_) {
          if (abs(GrCut[i][j] - userCutoff_) > 0.0001) {
            sprintf(painCave.errMsg,
                    "ForceMatrixDecomposition::createGtypeCutoffMap " 
                    "user-specified rCut (%lf) does not match computed group Cutoff\n", userCutoff_);
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();            
          }
        }
      }
    }
  }

  void ForceMatrixDecomposition::getGroupCutoffs(int &cg1, int &cg2, RealType &rcut, RealType &rcutsq, RealType &rlistsq) {
    int i, j;   
#ifdef IS_MPI
    i = groupRowToGtype[cg1];
    j = groupColToGtype[cg2];
#else
    i = groupToGtype[cg1];
    j = groupToGtype[cg2];
#endif    
    rcut = GrCut[i][j];
    rcutsq = GrCutSq[i][j];
    rlistsq = GrlistSq[i][j];
    return;
    //return gTypeCutoffMap[make_pair(i,j)];
  }

  int ForceMatrixDecomposition::getTopologicalDistance(int atom1, int atom2) {
    for (unsigned int j = 0; j < toposForAtom[atom1].size(); j++) {
      if (toposForAtom[atom1][j] == atom2) 
        return topoDist[atom1][j];
    }                                           
    return 0;
  }

  void ForceMatrixDecomposition::zeroWorkArrays() {
    pairwisePot = 0.0;
    embeddingPot = 0.0;
    excludedPot = 0.0;
    excludedSelfPot = 0.0;

#ifdef IS_MPI
    if (storageLayout_ & DataStorage::dslForce) {
      fill(atomRowData.force.begin(), atomRowData.force.end(), V3Zero);
      fill(atomColData.force.begin(), atomColData.force.end(), V3Zero);
    }

    if (storageLayout_ & DataStorage::dslTorque) {
      fill(atomRowData.torque.begin(), atomRowData.torque.end(), V3Zero);
      fill(atomColData.torque.begin(), atomColData.torque.end(), V3Zero);
    }
    
    fill(pot_row.begin(), pot_row.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));

    fill(pot_col.begin(), pot_col.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));   

    fill(expot_row.begin(), expot_row.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));

    fill(expot_col.begin(), expot_col.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));   

    if (storageLayout_ & DataStorage::dslParticlePot) {    
      fill(atomRowData.particlePot.begin(), atomRowData.particlePot.end(),
           0.0);
      fill(atomColData.particlePot.begin(), atomColData.particlePot.end(),
           0.0);
    }

    if (storageLayout_ & DataStorage::dslDensity) {      
      fill(atomRowData.density.begin(), atomRowData.density.end(), 0.0);
      fill(atomColData.density.begin(), atomColData.density.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslFunctional) {   
      fill(atomRowData.functional.begin(), atomRowData.functional.end(),
           0.0);
      fill(atomColData.functional.begin(), atomColData.functional.end(),
           0.0);
    }

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {      
      fill(atomRowData.functionalDerivative.begin(), 
           atomRowData.functionalDerivative.end(), 0.0);
      fill(atomColData.functionalDerivative.begin(), 
           atomColData.functionalDerivative.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslSkippedCharge) {      
      fill(atomRowData.skippedCharge.begin(), 
           atomRowData.skippedCharge.end(), 0.0);
      fill(atomColData.skippedCharge.begin(), 
           atomColData.skippedCharge.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslFlucQForce) {      
      fill(atomRowData.flucQFrc.begin(), 
           atomRowData.flucQFrc.end(), 0.0);
      fill(atomColData.flucQFrc.begin(), 
           atomColData.flucQFrc.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslElectricField) {    
      fill(atomRowData.electricField.begin(), 
           atomRowData.electricField.end(), V3Zero);
      fill(atomColData.electricField.begin(), 
           atomColData.electricField.end(), V3Zero);
    }

#endif
    // even in parallel, we need to zero out the local arrays:

    if (storageLayout_ & DataStorage::dslParticlePot) {      
      fill(snap_->atomData.particlePot.begin(), 
           snap_->atomData.particlePot.end(), 0.0);
    }
    
    if (storageLayout_ & DataStorage::dslDensity) {      
      fill(snap_->atomData.density.begin(), 
           snap_->atomData.density.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslFunctional) {
      fill(snap_->atomData.functional.begin(), 
           snap_->atomData.functional.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {      
      fill(snap_->atomData.functionalDerivative.begin(), 
           snap_->atomData.functionalDerivative.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslSkippedCharge) {      
      fill(snap_->atomData.skippedCharge.begin(), 
           snap_->atomData.skippedCharge.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslElectricField) {      
      fill(snap_->atomData.electricField.begin(), 
           snap_->atomData.electricField.end(), V3Zero);
    }
  }


  void ForceMatrixDecomposition::distributeData()  {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI
    
    // gather up the atomic positions
    AtomPlanVectorRow->gather(snap_->atomData.position, 
                              atomRowData.position);
    AtomPlanVectorColumn->gather(snap_->atomData.position, 
                                 atomColData.position);
    
    // gather up the cutoff group positions

    cgPlanVectorRow->gather(snap_->cgData.position, 
                            cgRowData.position);

    cgPlanVectorColumn->gather(snap_->cgData.position, 
                               cgColData.position);



    if (needVelocities_) {
      // gather up the atomic velocities
      AtomPlanVectorColumn->gather(snap_->atomData.velocity, 
                                   atomColData.velocity);
      
      cgPlanVectorColumn->gather(snap_->cgData.velocity, 
                                 cgColData.velocity);
    }

    
    // if needed, gather the atomic rotation matrices
    if (storageLayout_ & DataStorage::dslAmat) {
      AtomPlanMatrixRow->gather(snap_->atomData.aMat, 
                                atomRowData.aMat);
      AtomPlanMatrixColumn->gather(snap_->atomData.aMat, 
                                   atomColData.aMat);
    }

    // if needed, gather the atomic eletrostatic information
    if (storageLayout_ & DataStorage::dslDipole) {
      AtomPlanVectorRow->gather(snap_->atomData.dipole, 
                                atomRowData.dipole);
      AtomPlanVectorColumn->gather(snap_->atomData.dipole, 
                                   atomColData.dipole);
    }

    if (storageLayout_ & DataStorage::dslQuadrupole) {
      AtomPlanMatrixRow->gather(snap_->atomData.quadrupole, 
                                atomRowData.quadrupole);
      AtomPlanMatrixColumn->gather(snap_->atomData.quadrupole, 
                                   atomColData.quadrupole);
    }
        
    // if needed, gather the atomic fluctuating charge values
    if (storageLayout_ & DataStorage::dslFlucQPosition) {
      AtomPlanRealRow->gather(snap_->atomData.flucQPos, 
                              atomRowData.flucQPos);
      AtomPlanRealColumn->gather(snap_->atomData.flucQPos, 
                                 atomColData.flucQPos);
    }

#endif      
  }
  
  /* collects information obtained during the pre-pair loop onto local
   * data structures.
   */
  void ForceMatrixDecomposition::collectIntermediateData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI
    
    if (storageLayout_ & DataStorage::dslDensity) {
      
      AtomPlanRealRow->scatter(atomRowData.density, 
                               snap_->atomData.density);
      
      int n = snap_->atomData.density.size();
      vector<RealType> rho_tmp(n, 0.0);
      AtomPlanRealColumn->scatter(atomColData.density, rho_tmp);
      for (int i = 0; i < n; i++)
        snap_->atomData.density[i] += rho_tmp[i];
    }

    // this isn't necessary if we don't have polarizable atoms, but
    // we'll leave it here for now.
    if (storageLayout_ & DataStorage::dslElectricField) {
      
      AtomPlanVectorRow->scatter(atomRowData.electricField, 
                                 snap_->atomData.electricField);
      
      int n = snap_->atomData.electricField.size();
      vector<Vector3d> field_tmp(n, V3Zero);
      AtomPlanVectorColumn->scatter(atomColData.electricField, 
                                    field_tmp);
      for (int i = 0; i < n; i++)
        snap_->atomData.electricField[i] += field_tmp[i];
    }
#endif
  }

  /*
   * redistributes information obtained during the pre-pair loop out to 
   * row and column-indexed data structures
   */
  void ForceMatrixDecomposition::distributeIntermediateData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI
    if (storageLayout_ & DataStorage::dslFunctional) {
      AtomPlanRealRow->gather(snap_->atomData.functional, 
                              atomRowData.functional);
      AtomPlanRealColumn->gather(snap_->atomData.functional, 
                                 atomColData.functional);
    }
    
    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      AtomPlanRealRow->gather(snap_->atomData.functionalDerivative, 
                              atomRowData.functionalDerivative);
      AtomPlanRealColumn->gather(snap_->atomData.functionalDerivative, 
                                 atomColData.functionalDerivative);
    }
#endif
  }
  
  
  void ForceMatrixDecomposition::collectData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI    
    int n = snap_->atomData.force.size();
    vector<Vector3d> frc_tmp(n, V3Zero);
    
    AtomPlanVectorRow->scatter(atomRowData.force, frc_tmp);
    for (int i = 0; i < n; i++) {
      snap_->atomData.force[i] += frc_tmp[i];
      frc_tmp[i] = 0.0;
    }
    
    AtomPlanVectorColumn->scatter(atomColData.force, frc_tmp);
    for (int i = 0; i < n; i++) {
      snap_->atomData.force[i] += frc_tmp[i];
    }
        
    if (storageLayout_ & DataStorage::dslTorque) {

      int nt = snap_->atomData.torque.size();
      vector<Vector3d> trq_tmp(nt, V3Zero);

      AtomPlanVectorRow->scatter(atomRowData.torque, trq_tmp);
      for (int i = 0; i < nt; i++) {
        snap_->atomData.torque[i] += trq_tmp[i];
        trq_tmp[i] = 0.0;
      }
      
      AtomPlanVectorColumn->scatter(atomColData.torque, trq_tmp);
      for (int i = 0; i < nt; i++)
        snap_->atomData.torque[i] += trq_tmp[i];
    }

    if (storageLayout_ & DataStorage::dslSkippedCharge) {

      int ns = snap_->atomData.skippedCharge.size();
      vector<RealType> skch_tmp(ns, 0.0);

      AtomPlanRealRow->scatter(atomRowData.skippedCharge, skch_tmp);
      for (int i = 0; i < ns; i++) {
        snap_->atomData.skippedCharge[i] += skch_tmp[i];
        skch_tmp[i] = 0.0;
      }
      
      AtomPlanRealColumn->scatter(atomColData.skippedCharge, skch_tmp);
      for (int i = 0; i < ns; i++) 
        snap_->atomData.skippedCharge[i] += skch_tmp[i];
            
    }
    
    if (storageLayout_ & DataStorage::dslFlucQForce) {

      int nq = snap_->atomData.flucQFrc.size();
      vector<RealType> fqfrc_tmp(nq, 0.0);

      AtomPlanRealRow->scatter(atomRowData.flucQFrc, fqfrc_tmp);
      for (int i = 0; i < nq; i++) {
        snap_->atomData.flucQFrc[i] += fqfrc_tmp[i];
        fqfrc_tmp[i] = 0.0;
      }
      
      AtomPlanRealColumn->scatter(atomColData.flucQFrc, fqfrc_tmp);
      for (int i = 0; i < nq; i++) 
        snap_->atomData.flucQFrc[i] += fqfrc_tmp[i];
            
    }

    if (storageLayout_ & DataStorage::dslElectricField) {

      int nef = snap_->atomData.electricField.size();
      vector<Vector3d> efield_tmp(nef, V3Zero);

      AtomPlanVectorRow->scatter(atomRowData.electricField, efield_tmp);
      for (int i = 0; i < nef; i++) {
        snap_->atomData.electricField[i] += efield_tmp[i];
        efield_tmp[i] = 0.0;
      }
      
      AtomPlanVectorColumn->scatter(atomColData.electricField, efield_tmp);
      for (int i = 0; i < nef; i++)
        snap_->atomData.electricField[i] += efield_tmp[i];
    }


    nLocal_ = snap_->getNumberOfAtoms();

    vector<potVec> pot_temp(nLocal_, 
                            Vector<RealType, N_INTERACTION_FAMILIES> (0.0));
    vector<potVec> expot_temp(nLocal_, 
                              Vector<RealType, N_INTERACTION_FAMILIES> (0.0));

    // scatter/gather pot_row into the members of my column
          
    AtomPlanPotRow->scatter(pot_row, pot_temp);
    AtomPlanPotRow->scatter(expot_row, expot_temp);

    for (int ii = 0;  ii < pot_temp.size(); ii++ ) 
      pairwisePot += pot_temp[ii];

    for (int ii = 0;  ii < expot_temp.size(); ii++ ) 
      excludedPot += expot_temp[ii];
        
    if (storageLayout_ & DataStorage::dslParticlePot) {
      // This is the pairwise contribution to the particle pot.  The
      // embedding contribution is added in each of the low level
      // non-bonded routines.  In single processor, this is done in
      // unpackInteractionData, not in collectData.
      for (int ii = 0; ii < N_INTERACTION_FAMILIES; ii++) {
        for (int i = 0; i < nLocal_; i++) {
          // factor of two is because the total potential terms are divided
          // by 2 in parallel due to row/ column scatter       
          snap_->atomData.particlePot[i] += 2.0 * pot_temp[i](ii);
        }
      }
    }

    fill(pot_temp.begin(), pot_temp.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));
    fill(expot_temp.begin(), expot_temp.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));
      
    AtomPlanPotColumn->scatter(pot_col, pot_temp);    
    AtomPlanPotColumn->scatter(expot_col, expot_temp);    
    
    for (int ii = 0;  ii < pot_temp.size(); ii++ )
      pairwisePot += pot_temp[ii];    

    for (int ii = 0;  ii < expot_temp.size(); ii++ )
      excludedPot += expot_temp[ii];    

    if (storageLayout_ & DataStorage::dslParticlePot) {
      // This is the pairwise contribution to the particle pot.  The
      // embedding contribution is added in each of the low level
      // non-bonded routines.  In single processor, this is done in
      // unpackInteractionData, not in collectData.
      for (int ii = 0; ii < N_INTERACTION_FAMILIES; ii++) {
        for (int i = 0; i < nLocal_; i++) {
          // factor of two is because the total potential terms are divided
          // by 2 in parallel due to row/ column scatter       
          snap_->atomData.particlePot[i] += 2.0 * pot_temp[i](ii);
        }
      }
    }
    
    if (storageLayout_ & DataStorage::dslParticlePot) {
      int npp = snap_->atomData.particlePot.size();
      vector<RealType> ppot_temp(npp, 0.0);

      // This is the direct or embedding contribution to the particle
      // pot.
      
      AtomPlanRealRow->scatter(atomRowData.particlePot, ppot_temp);
      for (int i = 0; i < npp; i++) {
        snap_->atomData.particlePot[i] += ppot_temp[i];
      }

      fill(ppot_temp.begin(), ppot_temp.end(), 0.0);
      
      AtomPlanRealColumn->scatter(atomColData.particlePot, ppot_temp);
      for (int i = 0; i < npp; i++) {
        snap_->atomData.particlePot[i] += ppot_temp[i];
      }
    }

    for (int ii = 0; ii < N_INTERACTION_FAMILIES; ii++) {
      RealType ploc1 = pairwisePot[ii];
      RealType ploc2 = 0.0;
      MPI_Allreduce(&ploc1, &ploc2, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      pairwisePot[ii] = ploc2;
    }

    for (int ii = 0; ii < N_INTERACTION_FAMILIES; ii++) {
      RealType ploc1 = excludedPot[ii];
      RealType ploc2 = 0.0;
      MPI_Allreduce(&ploc1, &ploc2, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      excludedPot[ii] = ploc2;
    }

    // Here be dragons.
    MPI_Comm col = colComm.getComm();

    MPI_Allreduce(MPI_IN_PLACE, 
                  &snap_->frameData.conductiveHeatFlux[0], 3, 
                  MPI_REALTYPE, MPI_SUM, col);


#endif

  }

  /** 
   * Collects information obtained during the post-pair (and embedding
   * functional) loops onto local data structures.
   */
  void ForceMatrixDecomposition::collectSelfData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();

#ifdef IS_MPI
    for (int ii = 0; ii < N_INTERACTION_FAMILIES; ii++) {
      RealType ploc1 = embeddingPot[ii];
      RealType ploc2 = 0.0;
      MPI_Allreduce(&ploc1, &ploc2, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      embeddingPot[ii] = ploc2;
    }    
    for (int ii = 0; ii < N_INTERACTION_FAMILIES; ii++) {
      RealType ploc1 = excludedSelfPot[ii];
      RealType ploc2 = 0.0;
      MPI_Allreduce(&ploc1, &ploc2, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      excludedSelfPot[ii] = ploc2;
    }    
#endif
    
  }



  int& ForceMatrixDecomposition::getNAtomsInRow() {   
#ifdef IS_MPI
    return nAtomsInRow_;
#else
    return nLocal_;
#endif
  }

  /**
   * returns the list of atoms belonging to this group.  
   */
  vector<int>& ForceMatrixDecomposition::getAtomsInGroupRow(int cg1){
#ifdef IS_MPI
    return groupListRow_[cg1];
#else 
    return groupList_[cg1];
#endif
  }

  vector<int>& ForceMatrixDecomposition::getAtomsInGroupColumn(int cg2){
#ifdef IS_MPI
    return groupListCol_[cg2];
#else 
    return groupList_[cg2];
#endif
  }
  
  Vector3d ForceMatrixDecomposition::getIntergroupVector(int cg1, int cg2){
    Vector3d d;
    
#ifdef IS_MPI
    d = cgColData.position[cg2] - cgRowData.position[cg1];
#else
    d = snap_->cgData.position[cg2] - snap_->cgData.position[cg1];
#endif
    
    if (usePeriodicBoundaryConditions_) {
      snap_->wrapVector(d);
    }
    return d;    
  }

  Vector3d& ForceMatrixDecomposition::getGroupVelocityColumn(int cg2){
#ifdef IS_MPI
    return cgColData.velocity[cg2];
#else
    return snap_->cgData.velocity[cg2];
#endif
  }

  Vector3d& ForceMatrixDecomposition::getAtomVelocityColumn(int atom2){
#ifdef IS_MPI
    return atomColData.velocity[atom2];
#else
    return snap_->atomData.velocity[atom2];
#endif
  }


  Vector3d ForceMatrixDecomposition::getAtomToGroupVectorRow(int atom1, int cg1){

    Vector3d d;
    
#ifdef IS_MPI
    d = cgRowData.position[cg1] - atomRowData.position[atom1];
#else
    d = snap_->cgData.position[cg1] - snap_->atomData.position[atom1];
#endif
    if (usePeriodicBoundaryConditions_) {
      snap_->wrapVector(d);
    }
    return d;    
  }
  
  Vector3d ForceMatrixDecomposition::getAtomToGroupVectorColumn(int atom2, int cg2){
    Vector3d d;
    
#ifdef IS_MPI
    d = cgColData.position[cg2] - atomColData.position[atom2];
#else
    d = snap_->cgData.position[cg2] - snap_->atomData.position[atom2];
#endif
    if (usePeriodicBoundaryConditions_) {
      snap_->wrapVector(d);
    }
    return d;    
  }

  RealType& ForceMatrixDecomposition::getMassFactorRow(int atom1) {
#ifdef IS_MPI
    return massFactorsRow[atom1];
#else
    return massFactors[atom1];
#endif
  }

  RealType& ForceMatrixDecomposition::getMassFactorColumn(int atom2) {
#ifdef IS_MPI
    return massFactorsCol[atom2];
#else
    return massFactors[atom2];
#endif

  }
    
  Vector3d ForceMatrixDecomposition::getInteratomicVector(int atom1, int atom2){
    Vector3d d;
    
#ifdef IS_MPI
    d = atomColData.position[atom2] - atomRowData.position[atom1];
#else
    d = snap_->atomData.position[atom2] - snap_->atomData.position[atom1];
#endif
    if (usePeriodicBoundaryConditions_) {
      snap_->wrapVector(d);
    }
    return d;    
  }

  vector<int>& ForceMatrixDecomposition::getExcludesForAtom(int atom1) {
    return excludesForAtom[atom1];
  }

  /**
   * We need to exclude some overcounted interactions that result from
   * the parallel decomposition.
   */
  bool ForceMatrixDecomposition::skipAtomPair(int atom1, int atom2, int cg1, int cg2) {
    int unique_id_1, unique_id_2;
        
#ifdef IS_MPI
    // in MPI, we have to look up the unique IDs for each atom
    unique_id_1 = AtomRowToGlobal[atom1];
    unique_id_2 = AtomColToGlobal[atom2];
    // group1 = cgRowToGlobal[cg1];
    // group2 = cgColToGlobal[cg2];
#else
    unique_id_1 = AtomLocalToGlobal[atom1];
    unique_id_2 = AtomLocalToGlobal[atom2];
    int group1 = cgLocalToGlobal[cg1];
    int group2 = cgLocalToGlobal[cg2];
#endif   

    if (unique_id_1 == unique_id_2) return true;

#ifdef IS_MPI
    // this prevents us from doing the pair on multiple processors
    if (unique_id_1 < unique_id_2) {
      if ((unique_id_1 + unique_id_2) % 2 == 0) return true;
    } else {
      if ((unique_id_1 + unique_id_2) % 2 == 1) return true;
    }
#endif    

#ifndef IS_MPI
    if (group1 == group2) {
      if (unique_id_1 < unique_id_2) return true;
    }
#endif
    
    return false;
  }

  /**
   * We need to handle the interactions for atoms who are involved in
   * the same rigid body as well as some short range interactions
   * (bonds, bends, torsions) differently from other interactions.
   * We'll still visit the pairwise routines, but with a flag that
   * tells those routines to exclude the pair from direct long range
   * interactions.  Some indirect interactions (notably reaction
   * field) must still be handled for these pairs.
   */
  bool ForceMatrixDecomposition::excludeAtomPair(int atom1, int atom2) {

    // excludesForAtom was constructed to use row/column indices in the MPI
    // version, and to use local IDs in the non-MPI version:
    
    for (vector<int>::iterator i = excludesForAtom[atom1].begin();
         i != excludesForAtom[atom1].end(); ++i) {
      if ( (*i) == atom2 ) return true;
    }

    return false;
  }


  void ForceMatrixDecomposition::addForceToAtomRow(int atom1, Vector3d fg){
#ifdef IS_MPI
    atomRowData.force[atom1] += fg;
#else
    snap_->atomData.force[atom1] += fg;
#endif
  }

  void ForceMatrixDecomposition::addForceToAtomColumn(int atom2, Vector3d fg){
#ifdef IS_MPI
    atomColData.force[atom2] += fg;
#else
    snap_->atomData.force[atom2] += fg;
#endif
  }

    // filling interaction blocks with pointers
  void ForceMatrixDecomposition::fillInteractionData(InteractionData &idat, 
                                                     int atom1, int atom2) {

    idat.excluded = excludeAtomPair(atom1, atom2);
   
#ifdef IS_MPI
    //idat.atypes = make_pair( atypesRow[atom1], atypesCol[atom2]);
    idat.atid1 = identsRow[atom1];
    idat.atid2 = identsCol[atom2];

    if (regionsRow[atom1] >= 0 && regionsCol[atom2] >= 0) {
      idat.sameRegion = (regionsRow[atom1] == regionsCol[atom2]);
    } else {
      idat.sameRegion = false;
    }

    if (storageLayout_ & DataStorage::dslAmat) {
      idat.A1 = &(atomRowData.aMat[atom1]);
      idat.A2 = &(atomColData.aMat[atom2]);
    }
    
    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(atomRowData.torque[atom1]);
      idat.t2 = &(atomColData.torque[atom2]);
    }

    if (storageLayout_ & DataStorage::dslDipole) {
      idat.dipole1 = &(atomRowData.dipole[atom1]);
      idat.dipole2 = &(atomColData.dipole[atom2]);
    }

    if (storageLayout_ & DataStorage::dslQuadrupole) {
      idat.quadrupole1 = &(atomRowData.quadrupole[atom1]);
      idat.quadrupole2 = &(atomColData.quadrupole[atom2]);
    }

    if (storageLayout_ & DataStorage::dslDensity) {
      idat.rho1 = &(atomRowData.density[atom1]);
      idat.rho2 = &(atomColData.density[atom2]);
    }

    if (storageLayout_ & DataStorage::dslFunctional) {
      idat.frho1 = &(atomRowData.functional[atom1]);
      idat.frho2 = &(atomColData.functional[atom2]);
    }

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      idat.dfrho1 = &(atomRowData.functionalDerivative[atom1]);
      idat.dfrho2 = &(atomColData.functionalDerivative[atom2]);
    }

    if (storageLayout_ & DataStorage::dslParticlePot) {
      idat.particlePot1 = &(atomRowData.particlePot[atom1]);
      idat.particlePot2 = &(atomColData.particlePot[atom2]);
    }

    if (storageLayout_ & DataStorage::dslSkippedCharge) {              
      idat.skippedCharge1 = &(atomRowData.skippedCharge[atom1]);
      idat.skippedCharge2 = &(atomColData.skippedCharge[atom2]);
    }

    if (storageLayout_ & DataStorage::dslFlucQPosition) {              
      idat.flucQ1 = &(atomRowData.flucQPos[atom1]);
      idat.flucQ2 = &(atomColData.flucQPos[atom2]);
    }

#else
    
    //idat.atypes = make_pair( atypesLocal[atom1], atypesLocal[atom2]);
    idat.atid1 = idents[atom1];
    idat.atid2 = idents[atom2];

    if (regions[atom1] >= 0 && regions[atom2] >= 0) {
      idat.sameRegion = (regions[atom1] == regions[atom2]);
    } else {
      idat.sameRegion = false;
    }

    if (storageLayout_ & DataStorage::dslAmat) {
      idat.A1 = &(snap_->atomData.aMat[atom1]);
      idat.A2 = &(snap_->atomData.aMat[atom2]);
    }

    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(snap_->atomData.torque[atom1]);
      idat.t2 = &(snap_->atomData.torque[atom2]);
    }

    if (storageLayout_ & DataStorage::dslDipole) {
      idat.dipole1 = &(snap_->atomData.dipole[atom1]);
      idat.dipole2 = &(snap_->atomData.dipole[atom2]);
    }

    if (storageLayout_ & DataStorage::dslQuadrupole) {
      idat.quadrupole1 = &(snap_->atomData.quadrupole[atom1]);
      idat.quadrupole2 = &(snap_->atomData.quadrupole[atom2]);
    }

    if (storageLayout_ & DataStorage::dslDensity) {     
      idat.rho1 = &(snap_->atomData.density[atom1]);
      idat.rho2 = &(snap_->atomData.density[atom2]);
    }

    if (storageLayout_ & DataStorage::dslFunctional) {
      idat.frho1 = &(snap_->atomData.functional[atom1]);
      idat.frho2 = &(snap_->atomData.functional[atom2]);
    }

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      idat.dfrho1 = &(snap_->atomData.functionalDerivative[atom1]);
      idat.dfrho2 = &(snap_->atomData.functionalDerivative[atom2]);
    }

    if (storageLayout_ & DataStorage::dslParticlePot) {
      idat.particlePot1 = &(snap_->atomData.particlePot[atom1]);
      idat.particlePot2 = &(snap_->atomData.particlePot[atom2]);
    }

    if (storageLayout_ & DataStorage::dslSkippedCharge) {
      idat.skippedCharge1 = &(snap_->atomData.skippedCharge[atom1]);
      idat.skippedCharge2 = &(snap_->atomData.skippedCharge[atom2]);
    }

    if (storageLayout_ & DataStorage::dslFlucQPosition) {              
      idat.flucQ1 = &(snap_->atomData.flucQPos[atom1]);
      idat.flucQ2 = &(snap_->atomData.flucQPos[atom2]);
    }

#endif
  }

  
  void ForceMatrixDecomposition::unpackInteractionData(InteractionData &idat, int atom1, int atom2) {    
#ifdef IS_MPI
    pot_row[atom1] += RealType(0.5) *  *(idat.pot);
    pot_col[atom2] += RealType(0.5) *  *(idat.pot);
    expot_row[atom1] += RealType(0.5) *  *(idat.excludedPot);
    expot_col[atom2] += RealType(0.5) *  *(idat.excludedPot);

    atomRowData.force[atom1] += *(idat.f1);
    atomColData.force[atom2] -= *(idat.f1);

    if (storageLayout_ & DataStorage::dslFlucQForce) {              
      atomRowData.flucQFrc[atom1] -= *(idat.dVdFQ1);
      atomColData.flucQFrc[atom2] -= *(idat.dVdFQ2);
    }

    if (storageLayout_ & DataStorage::dslElectricField) {              
      atomRowData.electricField[atom1] += *(idat.eField1);
      atomColData.electricField[atom2] += *(idat.eField2);
    }

#else
    pairwisePot += *(idat.pot);
    excludedPot += *(idat.excludedPot);

    snap_->atomData.force[atom1] += *(idat.f1);
    snap_->atomData.force[atom2] -= *(idat.f1);

    if (idat.doParticlePot) {
      // This is the pairwise contribution to the particle pot.  The
      // embedding contribution is added in each of the low level
      // non-bonded routines.  In parallel, this calculation is done
      // in collectData, not in unpackInteractionData.
      snap_->atomData.particlePot[atom1] += *(idat.vpair) * *(idat.sw);
      snap_->atomData.particlePot[atom2] += *(idat.vpair) * *(idat.sw);
    }
    
    if (storageLayout_ & DataStorage::dslFlucQForce) {              
      snap_->atomData.flucQFrc[atom1] -= *(idat.dVdFQ1);
      snap_->atomData.flucQFrc[atom2] -= *(idat.dVdFQ2);
    }

    if (storageLayout_ & DataStorage::dslElectricField) {              
      snap_->atomData.electricField[atom1] += *(idat.eField1);
      snap_->atomData.electricField[atom2] += *(idat.eField2);
    }

#endif
    
  }

  /*
   * buildNeighborList
   *
   * first element of pair is row-indexed CutoffGroup
   * second element of pair is column-indexed CutoffGroup
   */
  void ForceMatrixDecomposition::buildNeighborList(vector<pair<int,int> >& neighborList) {
    
    neighborList.clear();
    groupCutoffs cuts;
    bool doAllPairs = false;

    RealType rList_ = (largestRcut_ + skinThickness_);
    RealType rcut, rcutsq, rlistsq;
    Snapshot* snap_ = sman_->getCurrentSnapshot();
    Mat3x3d box;
    Mat3x3d invBox;

    Vector3d rs, scaled, dr;
    Vector3i whichCell;
    int cellIndex;

#ifdef IS_MPI
    cellListRow_.clear();
    cellListCol_.clear();
#else
    cellList_.clear();
#endif
    
    if (!usePeriodicBoundaryConditions_) {
      box = snap_->getBoundingBox();
      invBox = snap_->getInvBoundingBox();
    } else {
      box = snap_->getHmat();
      invBox = snap_->getInvHmat();
    }
    
    Vector3d boxX = box.getColumn(0);
    Vector3d boxY = box.getColumn(1);
    Vector3d boxZ = box.getColumn(2);
    
    nCells_.x() = int( boxX.length() / rList_ );
    nCells_.y() = int( boxY.length() / rList_ );
    nCells_.z() = int( boxZ.length() / rList_ );
    
    // handle small boxes where the cell offsets can end up repeating cells
    
    if (nCells_.x() < 3) doAllPairs = true;
    if (nCells_.y() < 3) doAllPairs = true;
    if (nCells_.z() < 3) doAllPairs = true;
    
    int nCtot = nCells_.x() * nCells_.y() * nCells_.z();
    
#ifdef IS_MPI
    cellListRow_.resize(nCtot);
    cellListCol_.resize(nCtot);
#else
    cellList_.resize(nCtot);
#endif
    
    if (!doAllPairs) {
#ifdef IS_MPI
      
      for (int i = 0; i < nGroupsInRow_; i++) {
        rs = cgRowData.position[i];
        
        // scaled positions relative to the box vectors
        scaled = invBox * rs;
        
        // wrap the vector back into the unit box by subtracting integer box 
        // numbers
        for (int j = 0; j < 3; j++) {
          scaled[j] -= roundMe(scaled[j]);
          scaled[j] += 0.5;
          // Handle the special case when an object is exactly on the
          // boundary (a scaled coordinate of 1.0 is the same as
          // scaled coordinate of 0.0)
          if (scaled[j] >= 1.0) scaled[j] -= 1.0;
        }
        
        // find xyz-indices of cell that cutoffGroup is in.
        whichCell.x() = nCells_.x() * scaled.x();
        whichCell.y() = nCells_.y() * scaled.y();
        whichCell.z() = nCells_.z() * scaled.z();
        
        // find single index of this cell:
        cellIndex = Vlinear(whichCell, nCells_);
        
        // add this cutoff group to the list of groups in this cell;
        cellListRow_[cellIndex].push_back(i);
      }
      for (int i = 0; i < nGroupsInCol_; i++) {
        rs = cgColData.position[i];
        
        // scaled positions relative to the box vectors
        scaled = invBox * rs;
        
        // wrap the vector back into the unit box by subtracting integer box 
        // numbers
        for (int j = 0; j < 3; j++) {
          scaled[j] -= roundMe(scaled[j]);
          scaled[j] += 0.5;
          // Handle the special case when an object is exactly on the
          // boundary (a scaled coordinate of 1.0 is the same as
          // scaled coordinate of 0.0)
          if (scaled[j] >= 1.0) scaled[j] -= 1.0;
        }
        
        // find xyz-indices of cell that cutoffGroup is in.
        whichCell.x() = nCells_.x() * scaled.x();
        whichCell.y() = nCells_.y() * scaled.y();
        whichCell.z() = nCells_.z() * scaled.z();
        
        // find single index of this cell:
        cellIndex = Vlinear(whichCell, nCells_);
        
        // add this cutoff group to the list of groups in this cell;
        cellListCol_[cellIndex].push_back(i);
      }
      
#else
      for (int i = 0; i < nGroups_; i++) {
        rs = snap_->cgData.position[i];
        
        // scaled positions relative to the box vectors
        scaled = invBox * rs;
        
        // wrap the vector back into the unit box by subtracting integer box 
        // numbers
        for (int j = 0; j < 3; j++) {
          scaled[j] -= roundMe(scaled[j]);
          scaled[j] += 0.5;
          // Handle the special case when an object is exactly on the
          // boundary (a scaled coordinate of 1.0 is the same as
          // scaled coordinate of 0.0)
          if (scaled[j] >= 1.0) scaled[j] -= 1.0;
        }
        
        // find xyz-indices of cell that cutoffGroup is in.
        whichCell.x() = int(nCells_.x() * scaled.x());
        whichCell.y() = int(nCells_.y() * scaled.y());
        whichCell.z() = int(nCells_.z() * scaled.z());
        
        // find single index of this cell:
        cellIndex = Vlinear(whichCell, nCells_);
        
        // add this cutoff group to the list of groups in this cell;
        cellList_[cellIndex].push_back(i);
      }

#endif

      for (int m1z = 0; m1z < nCells_.z(); m1z++) {
        for (int m1y = 0; m1y < nCells_.y(); m1y++) {
          for (int m1x = 0; m1x < nCells_.x(); m1x++) {
            Vector3i m1v(m1x, m1y, m1z);
            int m1 = Vlinear(m1v, nCells_);
            
            for (vector<Vector3i>::iterator os = cellOffsets_.begin();
                 os != cellOffsets_.end(); ++os) {
              
              Vector3i m2v = m1v + (*os);
             

              if (m2v.x() >= nCells_.x()) {
                m2v.x() = 0;           
              } else if (m2v.x() < 0) {
                m2v.x() = nCells_.x() - 1; 
              }
              
              if (m2v.y() >= nCells_.y()) {
                m2v.y() = 0;           
              } else if (m2v.y() < 0) {
                m2v.y() = nCells_.y() - 1; 
              }
              
              if (m2v.z() >= nCells_.z()) {
                m2v.z() = 0;           
              } else if (m2v.z() < 0) {
                m2v.z() = nCells_.z() - 1; 
              }

              int m2 = Vlinear (m2v, nCells_);
              
#ifdef IS_MPI
              for (vector<int>::iterator j1 = cellListRow_[m1].begin(); 
                   j1 != cellListRow_[m1].end(); ++j1) {
                for (vector<int>::iterator j2 = cellListCol_[m2].begin(); 
                     j2 != cellListCol_[m2].end(); ++j2) {
                  
                  // In parallel, we need to visit *all* pairs of row
                  // & column indicies and will divide labor in the
                  // force evaluation later.
                  dr = cgColData.position[(*j2)] - cgRowData.position[(*j1)];
                  if (usePeriodicBoundaryConditions_) {
                    snap_->wrapVector(dr);
                  }
                  getGroupCutoffs( (*j1), (*j2), rcut, rcutsq, rlistsq );
                  if (dr.lengthSquare() < rlistsq) {
                    neighborList.push_back(make_pair((*j1), (*j2)));
                  }                  
                }
              }
#else
              for (vector<int>::iterator j1 = cellList_[m1].begin(); 
                   j1 != cellList_[m1].end(); ++j1) {
                for (vector<int>::iterator j2 = cellList_[m2].begin(); 
                     j2 != cellList_[m2].end(); ++j2) {
     
                  // Always do this if we're in different cells or if
                  // we're in the same cell and the global index of
                  // the j2 cutoff group is greater than or equal to
                  // the j1 cutoff group.  Note that Rappaport's code
                  // has a "less than" conditional here, but that
                  // deals with atom-by-atom computation.  OpenMD
                  // allows atoms within a single cutoff group to
                  // interact with each other.

                  if (m2 != m1 || (*j2) >= (*j1) ) {

                    dr = snap_->cgData.position[(*j2)] - snap_->cgData.position[(*j1)];
                    if (usePeriodicBoundaryConditions_) {
                      snap_->wrapVector(dr);
                    }
                    getGroupCutoffs( (*j1), (*j2), rcut, rcutsq, rlistsq );
                    if (dr.lengthSquare() < rlistsq) {
                      neighborList.push_back(make_pair((*j1), (*j2)));
                    }
                  }
                }
              }
#endif
            }
          }
        }
      }
    } else {
      // branch to do all cutoff group pairs
#ifdef IS_MPI
      for (int j1 = 0; j1 < nGroupsInRow_; j1++) {
        for (int j2 = 0; j2 < nGroupsInCol_; j2++) {    
          dr = cgColData.position[j2] - cgRowData.position[j1];
          if (usePeriodicBoundaryConditions_) {
            snap_->wrapVector(dr);
          }
          getGroupCutoffs( j1, j2, rcut, rcutsq, rlistsq);
          if (dr.lengthSquare() < rlistsq) {
            neighborList.push_back(make_pair(j1, j2));
          }
        }
      }      
#else
      // include all groups here.
      for (int j1 = 0; j1 < nGroups_; j1++) {
        // include self group interactions j2 == j1
        for (int j2 = j1; j2 < nGroups_; j2++) {
          dr = snap_->cgData.position[j2] - snap_->cgData.position[j1];
          if (usePeriodicBoundaryConditions_) {
            snap_->wrapVector(dr);
          }
          getGroupCutoffs( j1, j2, rcut, rcutsq, rlistsq );
          if (dr.lengthSquare() < rlistsq) {
            neighborList.push_back(make_pair(j1, j2));
          }
        }    
      }
#endif
    }
      
    // save the local cutoff group positions for the check that is
    // done on each loop:
    saved_CG_positions_.clear();
    for (int i = 0; i < nGroups_; i++)
      saved_CG_positions_.push_back(snap_->cgData.position[i]);
  }
} //end namespace OpenMD
