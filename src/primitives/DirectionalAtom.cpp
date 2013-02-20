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
 
#include "primitives/DirectionalAtom.hpp"
#include "types/DirectionalAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"
namespace OpenMD {
  
  DirectionalAtom::DirectionalAtom(AtomType* dAtomType) 
    : Atom(dAtomType) {
    objType_= otDAtom;

    DirectionalAdapter da = DirectionalAdapter(dAtomType);
    I_ = da.getI();

    MultipoleAdapter ma = MultipoleAdapter(dAtomType);
    if (ma.isDipole()) {
      dipole_ = ma.getDipole();
    }
    if (ma.isQuadrupole()) {
      quadrupole_ = ma.getQuadrupole();
    }

    // Check if one of the diagonal inertia tensor of this directional
    // atom is zero:
    int nLinearAxis = 0;
    Mat3x3d inertiaTensor = getI();
    for (int i = 0; i < 3; i++) {    
      if (fabs(inertiaTensor(i, i)) < OpenMD::epsilon) {
        linear_ = true;
        linearAxis_ = i;
        ++ nLinearAxis;
      }
    }

    if (nLinearAxis > 1) {
      sprintf( painCave.errMsg,
               "Directional Atom warning.\n"
               "\tOpenMD found more than one axis in this directional atom with a vanishing \n"
               "\tmoment of inertia.");
      painCave.isFatal = 0;
      simError();
    }    
  }
  
  Mat3x3d DirectionalAtom::getI() {
    return I_;     
  }    
  
  void DirectionalAtom::setPrevA(const RotMat3x3d& a) {
    ((snapshotMan_->getPrevSnapshot())->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      RotMat3x3d atrans = a.transpose();
      
      if (atomType_->isDipole()) {
        ((snapshotMan_->getPrevSnapshot())->*storage_).dipole[localIndex_] = atrans * dipole_;
      }

      if (atomType_->isQuadrupole()) {
        ((snapshotMan_->getPrevSnapshot())->*storage_).quadrupole[localIndex_] = atrans * quadrupole_ * a;
      }
    }
  }
  
  
  void DirectionalAtom::setA(const RotMat3x3d& a) {
    ((snapshotMan_->getCurrentSnapshot())->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      RotMat3x3d atrans = a.transpose();

      if (atomType_->isDipole()) {
        ((snapshotMan_->getCurrentSnapshot())->*storage_).dipole[localIndex_] = atrans * dipole_;
      }

      if (atomType_->isQuadrupole()) {
        ((snapshotMan_->getCurrentSnapshot())->*storage_).quadrupole[localIndex_] = atrans * quadrupole_ * a;
      }
    }
   
  }    
  
  void DirectionalAtom::setA(const RotMat3x3d& a, int snapshotNo) {
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      RotMat3x3d atrans = a.transpose();
      
      if (atomType_->isDipole()) {
        ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).dipole[localIndex_] = atrans * dipole_;
      }

      if (atomType_->isQuadrupole()) {
        ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).quadrupole[localIndex_] = atrans * quadrupole_ * a;
      }
    }

  }    
  
  void DirectionalAtom::rotateBy(const RotMat3x3d& m) {
    setA(m *getA());
  }
  
  std::vector<RealType> DirectionalAtom::getGrad() {
    std::vector<RealType> grad(6, 0.0);
    Vector3d force;
    Vector3d torque;
    Vector3d myEuler;
    RealType phi, theta;
    // RealType psi;
    RealType cphi, sphi, ctheta, stheta;
    Vector3d ephi;
    Vector3d etheta;
    Vector3d epsi;
    
    force = getFrc();
    torque =getTrq();
    myEuler = getA().toEulerAngles();
    
    phi = myEuler[0];
    theta = myEuler[1];
    // psi = myEuler[2];
    
    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);
    
    // get unit vectors along the phi, theta and psi rotation axes
    
    ephi[0] = 0.0;
    ephi[1] = 0.0;
    ephi[2] = 1.0;
    
    //etheta[0] = -sphi;
    //etheta[1] =  cphi;
    //etheta[2] =  0.0;
    
    etheta[0] = cphi;
    etheta[1] = sphi;
    etheta[2] = 0.0;
    
    epsi[0] = stheta * cphi;
    epsi[1] = stheta * sphi;
    epsi[2] = ctheta;
    
    //gradient is equal to -force
    for (int j = 0 ; j<3; j++)
      grad[j] = -force[j];
    
    for (int j = 0; j < 3; j++ ) {      
      grad[3] -= torque[j]*ephi[j];
      grad[4] -= torque[j]*etheta[j];
      grad[5] -= torque[j]*epsi[j];      
    }
    
    return grad;
  }    
  
  void DirectionalAtom::accept(BaseVisitor* v) {
    v->visit(this);
  }
}

