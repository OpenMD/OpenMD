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
#include "parallel/ForceMatrixDecomposition.hpp"
#include "math/SquareMatrix3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "brains/SnapshotManager.hpp"
#include "brains/PairList.hpp"

using namespace std;
namespace OpenMD {

  /**
   * distributeInitialData is essentially a copy of the older fortran 
   * SimulationSetup
   */
  
  void ForceMatrixDecomposition::distributeInitialData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
    ff_ = info_->getForceField();
    nLocal_ = snap_->getNumberOfAtoms();
    nGroups_ = snap_->getNumberOfCutoffGroups();

    // gather the information for atomtype IDs (atids):
    identsLocal = info_->getIdentArray();
    AtomLocalToGlobal = info_->getGlobalAtomIndices();
    cgLocalToGlobal = info_->getGlobalGroupIndices();
    vector<int> globalGroupMembership = info_->getGlobalGroupMembership();
    vector<RealType> massFactorsLocal = info_->getMassFactors();
    PairList excludes = info_->getExcludedInteractions();
    PairList oneTwo = info_->getOneTwoInteractions();
    PairList oneThree = info_->getOneThreeInteractions();
    PairList oneFour = info_->getOneFourInteractions();

#ifdef IS_MPI
 
    AtomCommIntRow = new Communicator<Row,int>(nLocal_);
    AtomCommRealRow = new Communicator<Row,RealType>(nLocal_);
    AtomCommVectorRow = new Communicator<Row,Vector3d>(nLocal_);
    AtomCommMatrixRow = new Communicator<Row,Mat3x3d>(nLocal_);
    AtomCommPotRow = new Communicator<Row,potVec>(nLocal_);

    AtomCommIntColumn = new Communicator<Column,int>(nLocal_);
    AtomCommRealColumn = new Communicator<Column,RealType>(nLocal_);
    AtomCommVectorColumn = new Communicator<Column,Vector3d>(nLocal_);
    AtomCommMatrixColumn = new Communicator<Column,Mat3x3d>(nLocal_);
    AtomCommPotColumn = new Communicator<Column,potVec>(nLocal_);

    cgCommIntRow = new Communicator<Row,int>(nGroups_);
    cgCommVectorRow = new Communicator<Row,Vector3d>(nGroups_);
    cgCommIntColumn = new Communicator<Column,int>(nGroups_);
    cgCommVectorColumn = new Communicator<Column,Vector3d>(nGroups_);

    nAtomsInRow_ = AtomCommIntRow->getSize();
    nAtomsInCol_ = AtomCommIntColumn->getSize();
    nGroupsInRow_ = cgCommIntRow->getSize();
    nGroupsInCol_ = cgCommIntColumn->getSize();

    // Modify the data storage objects with the correct layouts and sizes:
    atomRowData.resize(nAtomsInRow_);
    atomRowData.setStorageLayout(storageLayout_);
    atomColData.resize(nAtomsInCol_);
    atomColData.setStorageLayout(storageLayout_);
    cgRowData.resize(nGroupsInRow_);
    cgRowData.setStorageLayout(DataStorage::dslPosition);
    cgColData.resize(nGroupsInCol_);
    cgColData.setStorageLayout(DataStorage::dslPosition);
        
    identsRow.reserve(nAtomsInRow_);
    identsCol.reserve(nAtomsInCol_);
    
    AtomCommIntRow->gather(identsLocal, identsRow);
    AtomCommIntColumn->gather(identsLocal, identsCol);
    
    AtomCommIntRow->gather(AtomLocalToGlobal, AtomRowToGlobal);
    AtomCommIntColumn->gather(AtomLocalToGlobal, AtomColToGlobal);
    
    cgCommIntRow->gather(cgLocalToGlobal, cgRowToGlobal);
    cgCommIntColumn->gather(cgLocalToGlobal, cgColToGlobal);

    AtomCommRealRow->gather(massFactorsLocal, massFactorsRow);
    AtomCommRealColumn->gather(massFactorsLocal, massFactorsCol);

    groupListRow_.clear();
    groupListRow_.reserve(nGroupsInRow_);
    for (int i = 0; i < nGroupsInRow_; i++) {
      int gid = cgRowToGlobal[i];
      for (int j = 0; j < nAtomsInRow_; j++) {
        int aid = AtomRowToGlobal[j];
        if (globalGroupMembership[aid] == gid)
          groupListRow_[i].push_back(j);
      }      
    }

    groupListCol_.clear();
    groupListCol_.reserve(nGroupsInCol_);
    for (int i = 0; i < nGroupsInCol_; i++) {
      int gid = cgColToGlobal[i];
      for (int j = 0; j < nAtomsInCol_; j++) {
        int aid = AtomColToGlobal[j];
        if (globalGroupMembership[aid] == gid)
          groupListCol_[i].push_back(j);
      }      
    }

    skipsForRowAtom.clear();
    skipsForRowAtom.reserve(nAtomsInRow_);
    for (int i = 0; i < nAtomsInRow_; i++) {
      int iglob = AtomRowToGlobal[i];
      for (int j = 0; j < nAtomsInCol_; j++) {
        int jglob = AtomColToGlobal[j];        
        if (excludes.hasPair(iglob, jglob)) 
          skipsForRowAtom[i].push_back(j);       
      }      
    }

    toposForRowAtom.clear();
    toposForRowAtom.reserve(nAtomsInRow_);
    for (int i = 0; i < nAtomsInRow_; i++) {
      int iglob = AtomRowToGlobal[i];
      int nTopos = 0;
      for (int j = 0; j < nAtomsInCol_; j++) {
        int jglob = AtomColToGlobal[j];        
        if (oneTwo.hasPair(iglob, jglob)) {
          toposForRowAtom[i].push_back(j);
          topoDistRow[i][nTopos] = 1;
          nTopos++;
        }
        if (oneThree.hasPair(iglob, jglob)) {
          toposForRowAtom[i].push_back(j);
          topoDistRow[i][nTopos] = 2;
          nTopos++;
        }
        if (oneFour.hasPair(iglob, jglob)) {
          toposForRowAtom[i].push_back(j);
          topoDistRow[i][nTopos] = 3;
          nTopos++;
        }
      }      
    }

#endif

    groupList_.clear();
    groupList_.reserve(nGroups_);
    for (int i = 0; i < nGroups_; i++) {
      int gid = cgLocalToGlobal[i];
      for (int j = 0; j < nLocal_; j++) {
        int aid = AtomLocalToGlobal[j];
        if (globalGroupMembership[aid] == gid)
          groupList_[i].push_back(j);
      }      
    }

    skipsForLocalAtom.clear();
    skipsForLocalAtom.reserve(nLocal_);

    for (int i = 0; i < nLocal_; i++) {
      int iglob = AtomLocalToGlobal[i];
      for (int j = 0; j < nLocal_; j++) {
        int jglob = AtomLocalToGlobal[j];        
        if (excludes.hasPair(iglob, jglob)) 
          skipsForLocalAtom[i].push_back(j);       
      }      
    }

    toposForLocalAtom.clear();
    toposForLocalAtom.reserve(nLocal_);
    for (int i = 0; i < nLocal_; i++) {
      int iglob = AtomLocalToGlobal[i];
      int nTopos = 0;
      for (int j = 0; j < nLocal_; j++) {
        int jglob = AtomLocalToGlobal[j];        
        if (oneTwo.hasPair(iglob, jglob)) {
          toposForLocalAtom[i].push_back(j);
          topoDistLocal[i][nTopos] = 1;
          nTopos++;
        }
        if (oneThree.hasPair(iglob, jglob)) {
          toposForLocalAtom[i].push_back(j);
          topoDistLocal[i][nTopos] = 2;
          nTopos++;
        }
        if (oneFour.hasPair(iglob, jglob)) {
          toposForLocalAtom[i].push_back(j);
          topoDistLocal[i][nTopos] = 3;
          nTopos++;
        }
      }      
    }
  }
   
  void ForceMatrixDecomposition::zeroWorkArrays() {

    for (int j = 0; j < N_INTERACTION_FAMILIES; j++) {
      longRangePot_[j] = 0.0;
    }

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
    
    pot_local = Vector<RealType, N_INTERACTION_FAMILIES>(0.0);

    if (storageLayout_ & DataStorage::dslParticlePot) {    
      fill(atomRowData.particlePot.begin(), atomRowData.particlePot.end(), 0.0);
      fill(atomColData.particlePot.begin(), atomColData.particlePot.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslDensity) {      
      fill(atomRowData.density.begin(), atomRowData.density.end(), 0.0);
      fill(atomColData.density.begin(), atomColData.density.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslFunctional) {   
      fill(atomRowData.functional.begin(), atomRowData.functional.end(), 0.0);
      fill(atomColData.functional.begin(), atomColData.functional.end(), 0.0);
    }

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {      
      fill(atomRowData.functionalDerivative.begin(), 
           atomRowData.functionalDerivative.end(), 0.0);
      fill(atomColData.functionalDerivative.begin(), 
           atomColData.functionalDerivative.end(), 0.0);
    }

#else
    
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
#endif
    
  }


  void ForceMatrixDecomposition::distributeData()  {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI
    
    // gather up the atomic positions
    AtomCommVectorRow->gather(snap_->atomData.position, 
                              atomRowData.position);
    AtomCommVectorColumn->gather(snap_->atomData.position, 
                                 atomColData.position);
    
    // gather up the cutoff group positions
    cgCommVectorRow->gather(snap_->cgData.position, 
                            cgRowData.position);
    cgCommVectorColumn->gather(snap_->cgData.position, 
                               cgColData.position);
    
    // if needed, gather the atomic rotation matrices
    if (storageLayout_ & DataStorage::dslAmat) {
      AtomCommMatrixRow->gather(snap_->atomData.aMat, 
                                atomRowData.aMat);
      AtomCommMatrixColumn->gather(snap_->atomData.aMat, 
                                   atomColData.aMat);
    }
    
    // if needed, gather the atomic eletrostatic frames
    if (storageLayout_ & DataStorage::dslElectroFrame) {
      AtomCommMatrixRow->gather(snap_->atomData.electroFrame, 
                                atomRowData.electroFrame);
      AtomCommMatrixColumn->gather(snap_->atomData.electroFrame, 
                                   atomColData.electroFrame);
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
      
      AtomCommRealRow->scatter(atomRowData.density, 
                               snap_->atomData.density);
      
      int n = snap_->atomData.density.size();
      vector<RealType> rho_tmp(n, 0.0);
      AtomCommRealColumn->scatter(atomColData.density, rho_tmp);
      for (int i = 0; i < n; i++)
        snap_->atomData.density[i] += rho_tmp[i];
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
      AtomCommRealRow->gather(snap_->atomData.functional, 
                              atomRowData.functional);
      AtomCommRealColumn->gather(snap_->atomData.functional, 
                                 atomColData.functional);
    }
    
    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      AtomCommRealRow->gather(snap_->atomData.functionalDerivative, 
                              atomRowData.functionalDerivative);
      AtomCommRealColumn->gather(snap_->atomData.functionalDerivative, 
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
    
    AtomCommVectorRow->scatter(atomRowData.force, frc_tmp);
    for (int i = 0; i < n; i++) {
      snap_->atomData.force[i] += frc_tmp[i];
      frc_tmp[i] = 0.0;
    }
    
    AtomCommVectorColumn->scatter(atomColData.force, frc_tmp);
    for (int i = 0; i < n; i++)
      snap_->atomData.force[i] += frc_tmp[i];
    
    
    if (storageLayout_ & DataStorage::dslTorque) {

      int nt = snap_->atomData.force.size();
      vector<Vector3d> trq_tmp(nt, V3Zero);

      AtomCommVectorRow->scatter(atomRowData.torque, trq_tmp);
      for (int i = 0; i < n; i++) {
        snap_->atomData.torque[i] += trq_tmp[i];
        trq_tmp[i] = 0.0;
      }
      
      AtomCommVectorColumn->scatter(atomColData.torque, trq_tmp);
      for (int i = 0; i < n; i++)
        snap_->atomData.torque[i] += trq_tmp[i];
    }
    
    nLocal_ = snap_->getNumberOfAtoms();

    vector<potVec> pot_temp(nLocal_, 
                            Vector<RealType, N_INTERACTION_FAMILIES> (0.0));

    // scatter/gather pot_row into the members of my column
          
    AtomCommPotRow->scatter(pot_row, pot_temp);

    for (int ii = 0;  ii < pot_temp.size(); ii++ )
      pot_local += pot_temp[ii];
    
    fill(pot_temp.begin(), pot_temp.end(), 
         Vector<RealType, N_INTERACTION_FAMILIES> (0.0));
      
    AtomCommPotColumn->scatter(pot_col, pot_temp);    
    
    for (int ii = 0;  ii < pot_temp.size(); ii++ )
      pot_local += pot_temp[ii];
    
#endif
  }

  int ForceMatrixDecomposition::getNAtomsInRow() {   
#ifdef IS_MPI
    return nAtomsInRow_;
#else
    return nLocal_;
#endif
  }

  /**
   * returns the list of atoms belonging to this group.  
   */
  vector<int> ForceMatrixDecomposition::getAtomsInGroupRow(int cg1){
#ifdef IS_MPI
    return groupListRow_[cg1];
#else 
    return groupList_[cg1];
#endif
  }

  vector<int> ForceMatrixDecomposition::getAtomsInGroupColumn(int cg2){
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
    
    snap_->wrapVector(d);
    return d;    
  }


  Vector3d ForceMatrixDecomposition::getAtomToGroupVectorRow(int atom1, int cg1){

    Vector3d d;
    
#ifdef IS_MPI
    d = cgRowData.position[cg1] - atomRowData.position[atom1];
#else
    d = snap_->cgData.position[cg1] - snap_->atomData.position[atom1];
#endif

    snap_->wrapVector(d);
    return d;    
  }
  
  Vector3d ForceMatrixDecomposition::getAtomToGroupVectorColumn(int atom2, int cg2){
    Vector3d d;
    
#ifdef IS_MPI
    d = cgColData.position[cg2] - atomColData.position[atom2];
#else
    d = snap_->cgData.position[cg2] - snap_->atomData.position[atom2];
#endif
    
    snap_->wrapVector(d);
    return d;    
  }

  RealType ForceMatrixDecomposition::getMassFactorRow(int atom1) {
#ifdef IS_MPI
    return massFactorsRow[atom1];
#else
    return massFactorsLocal[atom1];
#endif
  }

  RealType ForceMatrixDecomposition::getMassFactorColumn(int atom2) {
#ifdef IS_MPI
    return massFactorsCol[atom2];
#else
    return massFactorsLocal[atom2];
#endif

  }
    
  Vector3d ForceMatrixDecomposition::getInteratomicVector(int atom1, int atom2){
    Vector3d d;
    
#ifdef IS_MPI
    d = atomColData.position[atom2] - atomRowData.position[atom1];
#else
    d = snap_->atomData.position[atom2] - snap_->atomData.position[atom1];
#endif

    snap_->wrapVector(d);
    return d;    
  }

  vector<int> ForceMatrixDecomposition::getSkipsForRowAtom(int atom1) {
#ifdef IS_MPI
    return skipsForRowAtom[atom1];
#else
    return skipsForLocalAtom[atom1];
#endif
  }

  /**
   * There are a number of reasons to skip a pair or a
   * particle. Mostly we do this to exclude atoms who are involved in
   * short range interactions (bonds, bends, torsions), but we also
   * need to exclude some overcounted interactions that result from
   * the parallel decomposition.
   */
  bool ForceMatrixDecomposition::skipAtomPair(int atom1, int atom2) {
    int unique_id_1, unique_id_2;

#ifdef IS_MPI
    // in MPI, we have to look up the unique IDs for each atom
    unique_id_1 = AtomRowToGlobal[atom1];
    unique_id_2 = AtomColToGlobal[atom2];

    // this situation should only arise in MPI simulations
    if (unique_id_1 == unique_id_2) return true;
    
    // this prevents us from doing the pair on multiple processors
    if (unique_id_1 < unique_id_2) {
      if ((unique_id_1 + unique_id_2) % 2 == 0) return true;
    } else {
      if ((unique_id_1 + unique_id_2) % 2 == 1) return true; 
    }
#else
    // in the normal loop, the atom numbers are unique
    unique_id_1 = atom1;
    unique_id_2 = atom2;
#endif
    
#ifdef IS_MPI
    for (vector<int>::iterator i = skipsForRowAtom[atom1].begin();
         i != skipsForRowAtom[atom1].end(); ++i) {
      if ( (*i) == unique_id_2 ) return true;
    }    
#else
    for (vector<int>::iterator i = skipsForLocalAtom[atom1].begin();
         i != skipsForLocalAtom[atom1].end(); ++i) {
      if ( (*i) == unique_id_2 ) return true;
    }    
#endif
  }

  int ForceMatrixDecomposition::getTopoDistance(int atom1, int atom2) {
    
#ifdef IS_MPI
    for (int i = 0; i < toposForRowAtom[atom1].size(); i++) {
      if ( toposForRowAtom[atom1][i] == atom2 ) return topoDistRow[atom1][i];
    }
#else
    for (int i = 0; i < toposForLocalAtom[atom1].size(); i++) {
      if ( toposForLocalAtom[atom1][i] == atom2 ) return topoDistLocal[atom1][i];
    }
#endif

    // zero is default for unconnected (i.e. normal) pair interactions
    return 0;
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
  InteractionData ForceMatrixDecomposition::fillInteractionData(int atom1, int atom2) {    
    InteractionData idat;

#ifdef IS_MPI
    
    idat.atypes = make_pair( ff_->getAtomType(identsRow[atom1]), 
                             ff_->getAtomType(identsCol[atom2]) );

    
    if (storageLayout_ & DataStorage::dslAmat) {
      idat.A1 = &(atomRowData.aMat[atom1]);
      idat.A2 = &(atomColData.aMat[atom2]);
    }
    
    if (storageLayout_ & DataStorage::dslElectroFrame) {
      idat.eFrame1 = &(atomRowData.electroFrame[atom1]);
      idat.eFrame2 = &(atomColData.electroFrame[atom2]);
    }

    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(atomRowData.torque[atom1]);
      idat.t2 = &(atomColData.torque[atom2]);
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

#else

    idat.atypes = make_pair( ff_->getAtomType(identsLocal[atom1]), 
                             ff_->getAtomType(identsLocal[atom2]) );

    if (storageLayout_ & DataStorage::dslAmat) {
      idat.A1 = &(snap_->atomData.aMat[atom1]);
      idat.A2 = &(snap_->atomData.aMat[atom2]);
    }

    if (storageLayout_ & DataStorage::dslElectroFrame) {
      idat.eFrame1 = &(snap_->atomData.electroFrame[atom1]);
      idat.eFrame2 = &(snap_->atomData.electroFrame[atom2]);
    }

    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(snap_->atomData.torque[atom1]);
      idat.t2 = &(snap_->atomData.torque[atom2]);
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

#endif
    return idat;
  }

  
  void ForceMatrixDecomposition::unpackInteractionData(InteractionData idat, int atom1, int atom2) {    
#ifdef IS_MPI
    pot_row[atom1] += 0.5 *  *(idat.pot);
    pot_col[atom2] += 0.5 *  *(idat.pot);

    atomRowData.force[atom1] += *(idat.f1);
    atomColData.force[atom2] -= *(idat.f1);
#else
    longRangePot_ += *(idat.pot);
    
    snap_->atomData.force[atom1] += *(idat.f1);
    snap_->atomData.force[atom2] -= *(idat.f1);
#endif

  }


  InteractionData ForceMatrixDecomposition::fillSkipData(int atom1, int atom2){

    InteractionData idat;
#ifdef IS_MPI
    idat.atypes = make_pair( ff_->getAtomType(identsRow[atom1]), 
                             ff_->getAtomType(identsCol[atom2]) );

    if (storageLayout_ & DataStorage::dslElectroFrame) {
      idat.eFrame1 = &(atomRowData.electroFrame[atom1]);
      idat.eFrame2 = &(atomColData.electroFrame[atom2]);
    }
    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(atomRowData.torque[atom1]);
      idat.t2 = &(atomColData.torque[atom2]);
    }
#else
    idat.atypes = make_pair( ff_->getAtomType(identsLocal[atom1]), 
                             ff_->getAtomType(identsLocal[atom2]) );

    if (storageLayout_ & DataStorage::dslElectroFrame) {
      idat.eFrame1 = &(snap_->atomData.electroFrame[atom1]);
      idat.eFrame2 = &(snap_->atomData.electroFrame[atom2]);
    }
    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(snap_->atomData.torque[atom1]);
      idat.t2 = &(snap_->atomData.torque[atom2]);
    }
#endif    
  }

  /*
   * buildNeighborList
   *
   * first element of pair is row-indexed CutoffGroup
   * second element of pair is column-indexed CutoffGroup
   */
  vector<pair<int, int> > ForceMatrixDecomposition::buildNeighborList() {
      
    vector<pair<int, int> > neighborList;
#ifdef IS_MPI
    cellListRow_.clear();
    cellListCol_.clear();
#else
    cellList_.clear();
#endif

    // dangerous to not do error checking.
    RealType rCut_;
 
    RealType rList_ = (rCut_ + skinThickness_);
    RealType rl2 = rList_ * rList_;
    Snapshot* snap_ = sman_->getCurrentSnapshot();
    Mat3x3d Hmat = snap_->getHmat();
    Vector3d Hx = Hmat.getColumn(0);
    Vector3d Hy = Hmat.getColumn(1);
    Vector3d Hz = Hmat.getColumn(2);

    nCells_.x() = (int) ( Hx.length() )/ rList_;
    nCells_.y() = (int) ( Hy.length() )/ rList_;
    nCells_.z() = (int) ( Hz.length() )/ rList_;

    Mat3x3d invHmat = snap_->getInvHmat();
    Vector3d rs, scaled, dr;
    Vector3i whichCell;
    int cellIndex;

#ifdef IS_MPI
    for (int i = 0; i < nGroupsInRow_; i++) {
      rs = cgRowData.position[i];
      // scaled positions relative to the box vectors
      scaled = invHmat * rs;
      // wrap the vector back into the unit box by subtracting integer box 
      // numbers
      for (int j = 0; j < 3; j++) 
        scaled[j] -= roundMe(scaled[j]);
     
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
      scaled = invHmat * rs;
      // wrap the vector back into the unit box by subtracting integer box 
      // numbers
      for (int j = 0; j < 3; j++) 
        scaled[j] -= roundMe(scaled[j]);

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
      scaled = invHmat * rs;
      // wrap the vector back into the unit box by subtracting integer box 
      // numbers
      for (int j = 0; j < 3; j++) 
        scaled[j] -= roundMe(scaled[j]);

      // find xyz-indices of cell that cutoffGroup is in.
      whichCell.x() = nCells_.x() * scaled.x();
      whichCell.y() = nCells_.y() * scaled.y();
      whichCell.z() = nCells_.z() * scaled.z();

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
                               
                // Always do this if we're in different cells or if
                // we're in the same cell and the global index of the
                // j2 cutoff group is less than the j1 cutoff group

                if (m2 != m1 || cgColToGlobal[(*j2)] < cgRowToGlobal[(*j1)]) {
                  dr = cgColData.position[(*j2)] - cgRowData.position[(*j1)];
                  snap_->wrapVector(dr);
                  if (dr.lengthSquare() < rl2) {
                    neighborList.push_back(make_pair((*j1), (*j2)));
                  }
                }
              }
            }
#else
            for (vector<int>::iterator j1 = cellList_[m1].begin(); 
                 j1 != cellList_[m1].end(); ++j1) {
              for (vector<int>::iterator j2 = cellList_[m2].begin(); 
                   j2 != cellList_[m2].end(); ++j2) {
                               
                // Always do this if we're in different cells or if
                // we're in the same cell and the global index of the
                // j2 cutoff group is less than the j1 cutoff group

                if (m2 != m1 || (*j2) < (*j1)) {
                  dr = snap_->cgData.position[(*j2)] - snap_->cgData.position[(*j1)];
                  snap_->wrapVector(dr);
                  if (dr.lengthSquare() < rl2) {
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

    // save the local cutoff group positions for the check that is
    // done on each loop:
    saved_CG_positions_.clear();
    for (int i = 0; i < nGroups_; i++)
      saved_CG_positions_.push_back(snap_->cgData.position[i]);

    return neighborList;
  }
} //end namespace OpenMD
