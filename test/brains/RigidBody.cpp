/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

#include "primitives/RigidBody.hpp"

namespace oopse {

RigidBody::RigidBody() : objType_(otRigidBody), storage_(&Snapshot::rigidbodyData){

}

Mat3x3d RigidBody::getI() {
    return inertiaTensor_;
}    

void RigidBody::setI(Mat3x3d& I) {
    inertiaTensor_ = I;
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

}

