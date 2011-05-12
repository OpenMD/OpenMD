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

using namespace std;
namespace OpenMD {

  /**
   * distributeInitialData is essentially a copy of the older fortran 
   * SimulationSetup
   */
  
  void ForceMatrixDecomposition::distributeInitialData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI    
    int nLocal = snap_->getNumberOfAtoms();
    int nGroups = snap_->getNumberOfCutoffGroups();
    
    AtomCommIntRow = new Communicator<Row,int>(nLocal);
    AtomCommRealRow = new Communicator<Row,RealType>(nLocal);
    AtomCommVectorRow = new Communicator<Row,Vector3d>(nLocal);
    AtomCommMatrixRow = new Communicator<Row,Mat3x3d>(nLocal);

    AtomCommIntColumn = new Communicator<Column,int>(nLocal);
    AtomCommRealColumn = new Communicator<Column,RealType>(nLocal);
    AtomCommVectorColumn = new Communicator<Column,Vector3d>(nLocal);
    AtomCommMatrixColumn = new Communicator<Column,Mat3x3d>(nLocal);

    cgCommIntRow = new Communicator<Row,int>(nGroups);
    cgCommVectorRow = new Communicator<Row,Vector3d>(nGroups);
    cgCommIntColumn = new Communicator<Column,int>(nGroups);
    cgCommVectorColumn = new Communicator<Column,Vector3d>(nGroups);

    int nAtomsInRow = AtomCommIntRow->getSize();
    int nAtomsInCol = AtomCommIntColumn->getSize();
    int nGroupsInRow = cgCommIntRow->getSize();
    int nGroupsInCol = cgCommIntColumn->getSize();

    // Modify the data storage objects with the correct layouts and sizes:
    atomRowData.resize(nAtomsInRow);
    atomRowData.setStorageLayout(storageLayout_);
    atomColData.resize(nAtomsInCol);
    atomColData.setStorageLayout(storageLayout_);
    cgRowData.resize(nGroupsInRow);
    cgRowData.setStorageLayout(DataStorage::dslPosition);
    cgColData.resize(nGroupsInCol);
    cgColData.setStorageLayout(DataStorage::dslPosition);
    
    vector<vector<RealType> > pot_row(N_INTERACTION_FAMILIES, 
                                      vector<RealType> (nAtomsInRow, 0.0));
    vector<vector<RealType> > pot_col(N_INTERACTION_FAMILIES,
                                      vector<RealType> (nAtomsInCol, 0.0));


    vector<RealType> pot_local(N_INTERACTION_FAMILIES, 0.0);
    
    // gather the information for atomtype IDs (atids):
    vector<int> identsLocal = info_->getIdentArray();
    identsRow.reserve(nAtomsInRow);
    identsCol.reserve(nAtomsInCol);
    
    AtomCommIntRow->gather(identsLocal, identsRow);
    AtomCommIntColumn->gather(identsLocal, identsCol);
    
    AtomLocalToGlobal = info_->getGlobalAtomIndices();
    AtomCommIntRow->gather(AtomLocalToGlobal, AtomRowToGlobal);
    AtomCommIntColumn->gather(AtomLocalToGlobal, AtomColToGlobal);
    
    cgLocalToGlobal = info_->getGlobalGroupIndices();
    cgCommIntRow->gather(cgLocalToGlobal, cgRowToGlobal);
    cgCommIntColumn->gather(cgLocalToGlobal, cgColToGlobal);

    // still need:
    // topoDist
    // exclude
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
  
  void ForceMatrixDecomposition::collectIntermediateData() {
    snap_ = sman_->getCurrentSnapshot();
    storageLayout_ = sman_->getStorageLayout();
#ifdef IS_MPI
    
    if (storageLayout_ & DataStorage::dslDensity) {
      
      AtomCommRealRow->scatter(atomRowData.density, 
                               snap_->atomData.density);
      
      int n = snap_->atomData.density.size();
      std::vector<RealType> rho_tmp(n, 0.0);
      AtomCommRealColumn->scatter(atomColData.density, rho_tmp);
      for (int i = 0; i < n; i++)
        snap_->atomData.density[i] += rho_tmp[i];
    }
#endif
  }
  
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
    
    int nLocal = snap_->getNumberOfAtoms();

    vector<vector<RealType> > pot_temp(N_INTERACTION_FAMILIES, 
                                       vector<RealType> (nLocal, 0.0));
    
    for (int i = 0; i < N_INTERACTION_FAMILIES; i++) {
      AtomCommRealRow->scatter(pot_row[i], pot_temp[i]);
      for (int ii = 0;  ii < pot_temp[i].size(); ii++ ) {
        pot_local[i] += pot_temp[i][ii];
      }
    }
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

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      idat.dfrho1 = &(atomRowData.functionalDerivative[atom1]);
      idat.dfrho2 = &(atomColData.functionalDerivative[atom2]);
    }
#else
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

    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      idat.dfrho1 = &(snap_->atomData.functionalDerivative[atom1]);
      idat.dfrho2 = &(snap_->atomData.functionalDerivative[atom2]);
    }
#endif
    
  }
  InteractionData ForceMatrixDecomposition::fillSkipData(int atom1, int atom2){
    InteractionData idat;
    skippedCharge1
      skippedCharge2
      rij
      d
    electroMult
    sw
    f
#ifdef IS_MPI

    if (storageLayout_ & DataStorage::dslElectroFrame) {
      idat.eFrame1 = &(atomRowData.electroFrame[atom1]);
      idat.eFrame2 = &(atomColData.electroFrame[atom2]);
    }
    if (storageLayout_ & DataStorage::dslTorque) {
      idat.t1 = &(atomRowData.torque[atom1]);
      idat.t2 = &(atomColData.torque[atom2]);
    }

    
  }
  SelfData ForceMatrixDecomposition::fillSelfData(int atom1) {
  }


  /*
   * buildNeighborList
   *
   * first element of pair is row-indexed CutoffGroup
   * second element of pair is column-indexed CutoffGroup
   */
  vector<pair<int, int> >  buildNeighborList() {
    Vector3d dr, invWid, rs, shift;
    Vector3i cc, m1v, m2s;
    RealType rrNebr;
    int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;


    vector<pair<int, int> > neighborList;   
    Vector3i nCells;
    Vector3d invWid, r;

    rList_ = (rCut_ + skinThickness_);
    rl2 = rList_ * rList_;

    snap_ = sman_->getCurrentSnapshot();
    Mat3x3d Hmat = snap_->getHmat();
    Vector3d Hx = Hmat.getColumn(0);
    Vector3d Hy = Hmat.getColumn(1);
    Vector3d Hz = Hmat.getColumn(2);

    nCells.x() = (int) ( Hx.length() )/ rList_;
    nCells.y() = (int) ( Hy.length() )/ rList_;
    nCells.z() = (int) ( Hz.length() )/ rList_;

    for (i = 0; i < nGroupsInRow; i++) {
      rs = cgRowData.position[i];
      snap_->scaleVector(rs);     
    }
    

    VDiv (invWid, cells, region);
    for (n = nMol; n < nMol + cells.componentProduct(); n ++) cellList[n] = -1;
    for (n = 0; n < nMol; n ++) {
      VSAdd (rs, mol[n].r, 0.5, region);
      VMul (cc, rs, invWid);
      c = VLinear (cc, cells) + nMol;
      cellList[n] = cellList[c];
      cellList[c] = n;
    }
    nebrTabLen = 0;
    for (m1z = 0; m1z < cells.z(); m1z++) {
      for (m1y = 0; m1y < cells.y(); m1y++) {
        for (m1x = 0; m1x < cells.x(); m1x++) {
          Vector3i m1v(m1x, m1y, m1z);
          m1 = VLinear(m1v, cells) + nMol;
          for (offset = 0; offset < nOffset_; offset++) {
            m2v = m1v + cellOffsets_[offset];
            shift = V3Zero();

            if (m2v.x() >= cells.x) {
              m2v.x() = 0;           
              shift.x() = region.x();  
            } else if (m2v.x() < 0) {
              m2v.x() = cells.x() - 1; 
              shift.x() = - region.x();
            }

            if (m2v.y() >= cells.y()) {
              m2v.y() = 0;           
              shift.y() = region.y();  
            } else if (m2v.y() < 0) {
              m2v.y() = cells.y() - 1; 
              shift.y() = - region.y();
            }

            m2 = VLinear (m2v, cells) + nMol;
            for (j1 = cellList[m1]; j1 >= 0; j1 = cellList[j1]) {
              for (j2 = cellList[m2]; j2 >= 0; j2 = cellList[j2]) {
                if (m1 != m2 || j2 < j1) {
                  dr = mol[j1].r - mol[j2].r;
                  VSub (dr, mol[j1].r, mol[j2].r);
                  VVSub (dr, shift);
                  if (VLenSq (dr) < rrNebr) {
                    neighborList.push_back(make_pair(j1, j2));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  
} //end namespace OpenMD
