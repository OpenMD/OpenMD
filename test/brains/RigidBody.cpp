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

#include "primitives/RigidBody.hpp"

namespace OpenMD {

RigidBody::RigidBody() : StuntDouble(otRigidBody, &Snapshot::rigidbodyData){

}

void RigidBody::setPrevA(const RotMat3x3d& a) {
    ((snapshotMan_->getPrevSnapshot())->*storage_).aMat[localIndex_] = a;
    ((snapshotMan_->getPrevSnapshot())->*storage_).unitVector[localIndex_] = a.inverse() * sU_.getColumn(2);

    std::vector<Atom*>::iterator i;
    for (i = atoms_.begin(); i != atoms_.end(); ++i) {
        if ((*i)->isDirectional()) {
            (*i)->setPrevA(a * (*i)->getPrevA());
        }
    }

}

      
void RigidBody::setA(const RotMat3x3d& a) {
    ((snapshotMan_->getCurrentSnapshot())->*storage_).aMat[localIndex_] = a;
    ((snapshotMan_->getCurrentSnapshot())->*storage_).unitVector[localIndex_] = a.inverse() * sU_.getColumn(2);

    std::vector<Atom*>::iterator i;
    for (i = atoms_.begin(); i != atoms_.end(); ++i) {
        if ((*i)->isDirectional()) {
            (*i)->setA(a * (*i)->getA());
        }
    }
}    
    
void RigidBody::setA(const RotMat3x3d& a, int snapshotNo) {
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).aMat[localIndex_] = a;
    ((snapshotMan_->getSnapshot(snapshotNo))->*storage_).unitVector[localIndex_] = a.inverse() * sU_.getColumn(2);    

    std::vector<Atom*>::iterator i;
    for (i = atoms_.begin(); i != atoms_.end(); ++i) {
        if ((*i)->isDirectional()) {
            (*i)->setA(a * (*i)->getA(snapshotNo), snapshotNo);
        }
    }

}   

void  DirectionalAtom::setUnitFrameFromEuler(double phi, double theta, double psi) {
    sU_.setupRotMat(phi,theta,psi);
}

Mat3x3d RigidBody::getI() {
    return inertiaTensor_;
}    

std::vector<double> RigidBody::getGrad() {
    vector<double> grad(6, 0.0);
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

void RigidBody::accept(BaseVisitor* v) {
    v->visit(this);
}    

void  RigidBody::calcRefCoords() {
  /*
  double mtmp;
  vec3 apos;
  double refCOM[3];
  vec3 ptmp;
  double Itmp[3][3];
  double evals[3];
  double evects[3][3];
  double r, r2, len;

  // First, find the center of mass:
  
  mass = 0.0;
  for (j=0; j<3; j++)
    refCOM[j] = 0.0;
  
  for (i = 0; i < atoms_.size(); i++) {
    mtmp = atoms_[i]->getMass();
    mass += mtmp;

    apos = refCoords[i];
    
    for(j = 0; j < 3; j++) {
      refCOM[j] += apos[j]*mtmp;     
    }    
  }
  
  for(j = 0; j < 3; j++) 
    refCOM[j] /= mass;

// Next, move the origin of the reference coordinate system to the COM:

  for (i = 0; i < atoms_.size(); i++) {
    apos = refCoords[i];
    for (j=0; j < 3; j++) {
      apos[j] = apos[j] - refCOM[j];
    }
    refCoords[i] = apos;
  }

// Moment of Inertia calculation

  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++)
      Itmp[i][j] = 0.0;  
  
  for (it = 0; it < atoms_.size(); it++) {

    mtmp = atoms_[it]->getMass();
    ptmp = refCoords[it];
    r= norm3(ptmp.vec);
    r2 = r*r;
    
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        
        if (i==j) Itmp[i][j] += mtmp * r2;

        Itmp[i][j] -= mtmp * ptmp.vec[i]*ptmp.vec[j];
      }
    }
  }
  
  diagonalize3x3(Itmp, evals, sU);
  
  // zero out I and then fill the diagonals with the moments of inertia:

  n_linear_coords = 0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      I[i][j] = 0.0;  
    }
    I[i][i] = evals[i];

    if (fabs(evals[i]) < momIntTol) {
      is_linear = true;
      n_linear_coords++;
      linear_axis = i;
    }
  }

  if (n_linear_coords > 1) {
          sprintf( painCave.errMsg,
               "RigidBody error.\n"
               "\tOpenMD found more than one axis in this rigid body with a vanishing \n"
               "\tmoment of inertia.  This can happen in one of three ways:\n"
               "\t 1) Only one atom was specified, or \n"
               "\t 2) All atoms were specified at the same location, or\n"
               "\t 3) The programmers did something stupid.\n"
               "\tIt is silly to use a rigid body to describe this situation.  Be smarter.\n"
               );
      painCave.isFatal = 1;
      simError();
  }
  
  // renormalize column vectors:
  
  for (i=0; i < 3; i++) {
    len = 0.0;
    for (j = 0; j < 3; j++) {
      len += sU[i][j]*sU[i][j];
    }
    len = sqrt(len);
    for (j = 0; j < 3; j++) {
      sU[i][j] /= len;
    }
  }
  */
}

void  RigidBody::calcForcesAndTorques() {
    unsigned int i;
    unsigned int j;
    //Vector3d apos;
    Vector3d afrc;
    Vector3d atrq;
    Vector3d rpos;
    Vector3d frc;
    Vector3d trq;
    //Vector3d pos;

    zeroForces();
    
    //pos = getPos();
    frc = getFrc();
    trq = getTrq();

    for (i = 0; i < atoms_.size(); i++) {

        afrc = atoms_[i]->getFrc();

        //apos = atoms_[i]->getPos(apos);
        //rpos = apos - pos;
        rpos = refCoords_[i];
        
        frc += afrc;

        trq[0] += rpos[1]*afrc[2] - rpos[2]*afrc[1];
        trq[1] += rpos[2]*afrc[0] - rpos[0]*afrc[2];
        trq[2] += rpos[0]*afrc[1] - rpos[1]*afrc[0];

        // If the atom has a torque associated with it, then we also need to 
        // migrate the torques onto the center of mass:

        if (atoms_[i]->isDirectional()) {
            atrq = atoms_[i]->getTrq();
            trq += atrq;
        }
        
    }
    
    setFrc(frc);
    setTrq(trq);
    
}

void  RigidBody::updateAtoms() {
    unsigned int i;
    unsigned int j;
    Vector3d ref;
    Vector3d apos;
    DirectionalAtom* dAtom;
    Vector3d pos = getPos();
    RotMat3x3d A = getA();
    
    for (i = 0; i < atoms_.size(); i++) {
     
        ref = body2Lab(refCoords_[i]);

        apos = pos + ref;

        atoms_[i]->setPos(apos);

        if (atoms_[i]->isDirectional()) {
          
          dAtom = (DirectionalAtom *) atoms_[i];
          dAtom->rotateBy( A );      
        }

    }
  
}


bool RigidBody::getAtomPos(Vector3d& pos, unsigned int index) {
    if (index < atoms_.size()) {

        Vector3d ref = body2Lab(refCoords_[index]);
        pos = getPos() + ref;
        return true;
    } else {
        std::cerr << index << " is an invalid index, current rigid body contains " 
                      << atoms_.size() << "atoms" << std::endl;
        return false;
    }    
}

bool RigidBody::getAtomPos(Vector3d& pos, Atom* atom) {
    std::vector<Atom*>::iterator i;
    i = find(atoms_.begin(), atoms_.end(), atom);
    if (i != atoms_.end()) {
        //RigidBody class makes sure refCoords_ and atoms_ match each other 
        Vector3d ref = body2Lab(refCoords_[i - atoms_.begin()]);
        pos = getPos() + ref;
        return true;
    } else {
        std::cerr << "Atom " << atom->getGlobalIndex() 
                      <<" does not belong to Rigid body "<< getGlobalIndex() << std::endl; 
        return false;
    }
}
bool RigidBody::getAtomVel(Vector3d& vel, unsigned int index) {

    //velRot = $(A\cdot skew(I^{-1}j))^{T}refCoor$

    if (index < atoms_.size()) {

        Vector3d velRot;
        Mat3x3d skewMat;;
        Vector3d ref = refCoords_[index];
        Vector3d ji = getJ();
        Mat3x3d I =  getI();

        skewMat(0, 0) =0;
        skewMat(0, 1) = ji[2] /I(2, 2);
        skewMat(0, 2) = -ji[1] /I(1, 1);

        skewMat(1, 0) = -ji[2] /I(2, 2);
        skewMat(1, 1) = 0;
        skewMat(1, 2) = ji[0]/I(0, 0);

        skewMat(2, 0) =ji[1] /I(1, 1);
        skewMat(2, 1) = -ji[0]/I(0, 0);
        skewMat(2, 2) = 0;

        velRot = (getA() * skewMat).transpose() * ref;

        vel =getVel() + velRot;
        return true;
        
    } else {
        std::cerr << index << " is an invalid index, current rigid body contains " 
                      << atoms_.size() << "atoms" << std::endl;
        return false;
    }
}

bool RigidBody::getAtomVel(Vector3d& vel, Atom* atom) {

    std::vector<Atom*>::iterator i;
    i = find(atoms_.begin(), atoms_.end(), atom);
    if (i != atoms_.end()) {
        return getAtomVel(vel, i - atoms_.begin());
    } else {
        std::cerr << "Atom " << atom->getGlobalIndex() 
                      <<" does not belong to Rigid body "<< getGlobalIndex() << std::endl;    
        return false;
    }    
}

bool RigidBody::getAtomRefCoor(Vector3d& coor, unsigned int index) {
    if (index < atoms_.size()) {

        coor = refCoords_[index];
        return true;
    } else {
        std::cerr << index << " is an invalid index, current rigid body contains " 
                      << atoms_.size() << "atoms" << std::endl;
        return false;
    }

}

bool RigidBody::getAtomRefCoor(Vector3d& coor, Atom* atom) {
    std::vector<Atom*>::iterator i;
    i = find(atoms_.begin(), atoms_.end(), atom);
    if (i != atoms_.end()) {
        //RigidBody class makes sure refCoords_ and atoms_ match each other 
        coor = refCoords_[i - atoms_.begin()];
        return true;
    } else {
        std::cerr << "Atom " << atom->getGlobalIndex() 
                      <<" does not belong to Rigid body "<< getGlobalIndex() << std::endl;    
        return false;
    }

}

}

