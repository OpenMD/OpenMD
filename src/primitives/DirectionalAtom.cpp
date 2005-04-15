/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#include "primitives/DirectionalAtom.hpp"
#include "utils/simError.h"
namespace oopse {

  DirectionalAtom::DirectionalAtom(DirectionalAtomType* dAtomType) 
    : Atom(dAtomType){
      objType_= otDAtom;
      if (dAtomType->isMultipole()) {
        electroBodyFrame_ = dAtomType->getElectroBodyFrame();
      }

      //check if one of the diagonal inertia tensor of this directional atom  is zero
      int nLinearAxis = 0;
      Mat3x3d inertiaTensor = getI();
      for (int i = 0; i < 3; i++) {    
        if (fabs(inertiaTensor(i, i)) < oopse::epsilon) {
	  linear_ = true;
	  linearAxis_ = i;
	  ++ nLinearAxis;
        }
      }

      if (nLinearAxis > 1) {
        sprintf( painCave.errMsg,
		 "Directional Atom error.\n"
		 "\tOOPSE found more than one axis in this directional atom with a vanishing \n"
		 "\tmoment of inertia.");
        painCave.isFatal = 1;
        simError();
      }
      
    }

  Mat3x3d DirectionalAtom::getI() {
    return static_cast<DirectionalAtomType*>(getAtomType())->getI();
  }    

  void DirectionalAtom::setPrevA(const RotMat3x3d& a) {
    ((snapshotMan_->getPrevSnapshot())->*storage_).aMat[localIndex_] = a;
    if (atomType_->isMultipole()) {
      ((snapshotMan_->getPrevSnapshot())->*storage_).electroFrame[localIndex_] = a.transpose() * electroBodyFrame_;
    }
  }

      
  void DirectionalAtom::setA(const RotMat3x3d& a) {
    ((snapshotMan_->getCurrentSnapshot())->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      ((snapshotMan_->getCurrentSnapshot())->*storage_).electroFrame[localIndex_] = a.transpose() * electroBodyFrame_;
    }
  }    
    
  void DirectionalAtom::setA(const RotMat3x3d& a, int snapshotNo) {
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).aMat[localIndex_] = a;

    if (atomType_->isMultipole()) {
      ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).electroFrame[localIndex_] = a.transpose() * electroBodyFrame_;    
    }
  }    

  void DirectionalAtom::rotateBy(const RotMat3x3d& m) {
    setA(m *getA());
  }

  std::vector<double> DirectionalAtom::getGrad() {
    std::vector<double> grad(6, 0.0);
    Vector3d force;
    Vector3d torque;
    Vector3d myEuler;
    double phi, theta, psi;
    double cphi, sphi, ctheta, stheta;
    Vector3d ephi;
    Vector3d etheta;
    Vector3d epsi;

    force = getFrc();
    torque =getTrq();
    myEuler = getA().toEulerAngles();

    phi = myEuler[0];
    theta = myEuler[1];
    psi = myEuler[2];

    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);

    // get unit vectors along the phi, theta and psi rotation axes

    ephi[0] = 0.0;
    ephi[1] = 0.0;
    ephi[2] = 1.0;

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

      grad[3] += torque[j]*ephi[j];
      grad[4] += torque[j]*etheta[j];
      grad[5] += torque[j]*epsi[j];

    }
    
    return grad;
  }    

  void DirectionalAtom::accept(BaseVisitor* v) {
    v->visit(this);
  }    

}

