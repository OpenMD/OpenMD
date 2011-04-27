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
#ifdef IS_MPI    
    Snapshot* snap = sman_->getCurrentSnapshot();
    int nLocal = snap->getNumberOfAtoms();
    int nGroups = snap->getNumberOfCutoffGroups();

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
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    
    // gather up the atomic positions
    AtomCommVectorRow->gather(snap->atomData.position, 
                            snap->atomIData.position);
    AtomCommVectorColumn->gather(snap->atomData.position, 
                            snap->atomJData.position);
    
    // gather up the cutoff group positions
    cgCommVectorRow->gather(snap->cgData.position, 
                          snap->cgIData.position);
    cgCommVectorColumn->gather(snap->cgData.position, 
                          snap->cgJData.position);
    
    // if needed, gather the atomic rotation matrices
    if (snap->atomData.getStorageLayout() & DataStorage::dslAmat) {
      AtomCommMatrixRow->gather(snap->atomData.aMat, 
                              snap->atomIData.aMat);
      AtomCommMatrixColumn->gather(snap->atomData.aMat, 
                              snap->atomJData.aMat);
    }
    
    // if needed, gather the atomic eletrostatic frames
    if (snap->atomData.getStorageLayout() & DataStorage::dslElectroFrame) {
      AtomCommMatrixRow->gather(snap->atomData.electroFrame, 
                              snap->atomIData.electroFrame);
      AtomCommMatrixColumn->gather(snap->atomData.electroFrame, 
                              snap->atomJData.electroFrame);
    }
#endif      
  }
  
  void ForceMatrixDecomposition::collectIntermediateData() {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslDensity) {

      AtomCommRealRow->scatter(snap->atomIData.density, 
                             snap->atomData.density);

      int n = snap->atomData.density.size();
      std::vector<RealType> rho_tmp(n, 0.0);
      AtomCommRealColumn->scatter(snap->atomJData.density, rho_tmp);
      for (int i = 0; i < n; i++)
        snap->atomData.density[i] += rho_tmp[i];
    }
#endif
  }
  
  void ForceMatrixDecomposition::distributeIntermediateData() {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    if (snap->atomData.getStorageLayout() & DataStorage::dslFunctional) {
      AtomCommRealRow->gather(snap->atomData.functional, 
                            snap->atomIData.functional);
      AtomCommRealColumn->gather(snap->atomData.functional, 
                            snap->atomJData.functional);
    }
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslFunctionalDerivative) {
      AtomCommRealRow->gather(snap->atomData.functionalDerivative, 
                            snap->atomIData.functionalDerivative);
      AtomCommRealColumn->gather(snap->atomData.functionalDerivative, 
                            snap->atomJData.functionalDerivative);
    }
#endif
  }
  
  
  void ForceMatrixDecomposition::collectData() {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    
    int n = snap->atomData.force.size();
    vector<Vector3d> frc_tmp(n, V3Zero);
    
    AtomCommVectorRow->scatter(snap->atomIData.force, frc_tmp);
    for (int i = 0; i < n; i++) {
      snap->atomData.force[i] += frc_tmp[i];
      frc_tmp[i] = 0.0;
    }
    
    AtomCommVectorColumn->scatter(snap->atomJData.force, frc_tmp);
    for (int i = 0; i < n; i++)
      snap->atomData.force[i] += frc_tmp[i];
    
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslTorque) {

      int nt = snap->atomData.force.size();
      vector<Vector3d> trq_tmp(nt, V3Zero);

      AtomCommVectorRow->scatter(snap->atomIData.torque, trq_tmp);
      for (int i = 0; i < n; i++) {
        snap->atomData.torque[i] += trq_tmp[i];
        trq_tmp[i] = 0.0;
      }
      
      AtomCommVectorColumn->scatter(snap->atomJData.torque, trq_tmp);
      for (int i = 0; i < n; i++)
        snap->atomData.torque[i] += trq_tmp[i];
    }
    
    int nLocal = snap->getNumberOfAtoms();

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
  
} //end namespace OpenMD
