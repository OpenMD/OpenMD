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
#include "parallel/ForceDecomposition.hpp"
#include "math/SquareMatrix3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "brains/SnapshotManager.hpp"

using namespace std;
namespace OpenMD {

  /**
   * distributeInitialData is essentially a copy of the older fortran 
   * SimulationSetup
   */
  
  void ForceDecomposition::distributeInitialData() {
#ifdef IS_MPI    
    Snapshot* snap = sman_->getCurrentSnapshot();
    int nLocal = snap->getNumberOfAtoms();
    int nGroups = snap->getNumberOfCutoffGroups();

    AtomCommIntI = new Communicator<Row,int>(nLocal);
    AtomCommRealI = new Communicator<Row,RealType>(nLocal);
    AtomCommVectorI = new Communicator<Row,Vector3d>(nLocal);
    AtomCommMatrixI = new Communicator<Row,Mat3x3d>(nLocal);

    AtomCommIntJ = new Communicator<Column,int>(nLocal);
    AtomCommRealJ = new Communicator<Column,RealType>(nLocal);
    AtomCommVectorJ = new Communicator<Column,Vector3d>(nLocal);
    AtomCommMatrixJ = new Communicator<Column,Mat3x3d>(nLocal);

    cgCommIntI = new Communicator<Row,int>(nGroups);
    cgCommVectorI = new Communicator<Row,Vector3d>(nGroups);
    cgCommIntJ = new Communicator<Column,int>(nGroups);
    cgCommVectorJ = new Communicator<Column,Vector3d>(nGroups);

    int nAtomsInRow = AtomCommIntI->getSize();
    int nAtomsInCol = AtomCommIntJ->getSize();
    int nGroupsInRow = cgCommIntI->getSize();
    int nGroupsInCol = cgCommIntJ->getSize();

    vector<vector<RealType> > pot_row(N_INTERACTION_FAMILIES, 
                                      vector<RealType> (nAtomsInRow, 0.0));
    vector<vector<RealType> > pot_col(N_INTERACTION_FAMILIES,
                                      vector<RealType> (nAtomsInCol, 0.0));
    
    vector<RealType> pot_local(N_INTERACTION_FAMILIES, 0.0);

    // gather the information for atomtype IDs (atids):
    AtomCommIntI->gather(info_->getIdentArray(), identsRow);
    AtomCommIntJ->gather(info_->getIdentArray(), identsCol);

    AtomLocalToGlobal = info_->getLocalToGlobalAtomIndex();
    AtomCommIntI->gather(AtomLocalToGlobal, AtomRowToGlobal);
    AtomCommIntJ->gather(AtomLocalToGlobal, AtomColToGlobal);

    cgLocalToGlobal = info_->getLocalToGlobalCutoffGroupIndex();
    cgCommIntI->gather(cgLocalToGlobal, cgRowToGlobal);
    cgCommIntJ->gather(cgLocalToGlobal, cgColToGlobal);

      
      



    // still need:
    // topoDist
    // exclude
#endif
  }
    


  void ForceDecomposition::distributeData()  {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    
    // gather up the atomic positions
    AtomCommVectorI->gather(snap->atomData.position, 
                            snap->atomIData.position);
    AtomCommVectorJ->gather(snap->atomData.position, 
                            snap->atomJData.position);
    
    // gather up the cutoff group positions
    cgCommVectorI->gather(snap->cgData.position, 
                          snap->cgIData.position);
    cgCommVectorJ->gather(snap->cgData.position, 
                          snap->cgJData.position);
    
    // if needed, gather the atomic rotation matrices
    if (snap->atomData.getStorageLayout() & DataStorage::dslAmat) {
      AtomCommMatrixI->gather(snap->atomData.aMat, 
                              snap->atomIData.aMat);
      AtomCommMatrixJ->gather(snap->atomData.aMat, 
                              snap->atomJData.aMat);
    }
    
    // if needed, gather the atomic eletrostatic frames
    if (snap->atomData.getStorageLayout() & DataStorage::dslElectroFrame) {
      AtomCommMatrixI->gather(snap->atomData.electroFrame, 
                              snap->atomIData.electroFrame);
      AtomCommMatrixJ->gather(snap->atomData.electroFrame, 
                              snap->atomJData.electroFrame);
    }
#endif      
  }
  
  void ForceDecomposition::collectIntermediateData() {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslDensity) {

      AtomCommRealI->scatter(snap->atomIData.density, 
                             snap->atomData.density);

      int n = snap->atomData.density.size();
      std::vector<RealType> rho_tmp(n, 0.0);
      AtomCommRealJ->scatter(snap->atomJData.density, rho_tmp);
      for (int i = 0; i < n; i++)
        snap->atomData.density[i] += rho_tmp[i];
    }
#endif
  }
  
  void ForceDecomposition::distributeIntermediateData() {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    if (snap->atomData.getStorageLayout() & DataStorage::dslFunctional) {
      AtomCommRealI->gather(snap->atomData.functional, 
                            snap->atomIData.functional);
      AtomCommRealJ->gather(snap->atomData.functional, 
                            snap->atomJData.functional);
    }
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslFunctionalDerivative) {
      AtomCommRealI->gather(snap->atomData.functionalDerivative, 
                            snap->atomIData.functionalDerivative);
      AtomCommRealJ->gather(snap->atomData.functionalDerivative, 
                            snap->atomJData.functionalDerivative);
    }
#endif
  }
  
  
  void ForceDecomposition::collectData() {
#ifdef IS_MPI
    Snapshot* snap = sman_->getCurrentSnapshot();
    
    int n = snap->atomData.force.size();
    vector<Vector3d> frc_tmp(n, V3Zero);
    
    AtomCommVectorI->scatter(snap->atomIData.force, frc_tmp);
    for (int i = 0; i < n; i++) {
      snap->atomData.force[i] += frc_tmp[i];
      frc_tmp[i] = 0.0;
    }
    
    AtomCommVectorJ->scatter(snap->atomJData.force, frc_tmp);
    for (int i = 0; i < n; i++)
      snap->atomData.force[i] += frc_tmp[i];
    
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslTorque) {

      int nt = snap->atomData.force.size();
      vector<Vector3d> trq_tmp(nt, V3Zero);

      AtomCommVectorI->scatter(snap->atomIData.torque, trq_tmp);
      for (int i = 0; i < n; i++) {
        snap->atomData.torque[i] += trq_tmp[i];
        trq_tmp[i] = 0.0;
      }
      
      AtomCommVectorJ->scatter(snap->atomJData.torque, trq_tmp);
      for (int i = 0; i < n; i++)
        snap->atomData.torque[i] += trq_tmp[i];
    }
    
    int nLocal = snap->getNumberOfAtoms();

    vector<vector<RealType> > pot_temp(N_INTERACTION_FAMILIES, 
                                       vector<RealType> (nLocal, 0.0));
    
    for (int i = 0; i < N_INTERACTION_FAMILIES; i++) {
      AtomCommRealI->scatter(pot_row[i], pot_temp[i]);
      for (int ii = 0;  ii < pot_temp[i].size(); ii++ ) {
        pot_local[i] += pot_temp[i][ii];
      }
    }
#endif
  }
  
} //end namespace OpenMD
