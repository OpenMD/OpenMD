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
#include "parallel/Communicator.hpp"
#include "math/SquareMatrix3.hpp"

namespace OpenMD {

  void ForceDecomposition::distributeInitialData() {
#ifdef IS_MPI

    int nAtoms;
    int nGroups;

    AtomCommRealI = new Comm<I,RealType>(nAtoms);
    AtomCommVectorI = new Comm<I,Vector3d>(nAtoms);
    AtomCommMatrixI = new Comm<I,Mat3x3d>(nAtoms);

    AtomCommRealJ = new Comm<J,RealType>(nAtoms);
    AtomCommVectorJ = new Comm<J,Vector3d>(nAtoms);
    AtomCommMatrixJ = new Comm<J,Mat3x3d>(nAtoms);

    cgCommVectorI = new Comm<I,Vector3d>(nGroups);
    cgCommVectorJ = new Comm<J,Vector3d>(nGroups);
    // more to come
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
    // gather up the atomic positions
    
    if (snap->atomData.getStorageLayout() & DataStorage::dslDensity) {
      AtomCommRealI->scatter(snap->atomIData.density, 
                            snap->atomData.density);
      std::vector<RealType> rho_tmp;
      int n = snap->getNumberOfAtoms();
      rho_tmp.reserve( n );
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
#endif
  }
  
} //end namespace OpenMD
