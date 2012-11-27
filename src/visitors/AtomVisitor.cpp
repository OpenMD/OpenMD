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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <cstring>
#include "visitors/AtomVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "types/MultipoleAdapter.hpp"
#include "types/GayBerneAdapter.hpp"

namespace OpenMD {
  void BaseAtomVisitor::visit(RigidBody *rb) {
    //vector<Atom*> myAtoms;
    //vector<Atom*>::iterator atomIter;

    //myAtoms = rb->getAtoms();

    //for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
    //  (*atomIter)->accept(this);
  }

  void BaseAtomVisitor::setVisited(Atom *atom) {
    GenericData *data;
    data = atom->getPropertyByName("VISITED");

    //if visited property is not existed, add it as new property
    if (data == NULL) {
      data = new GenericData();
      data->setID("VISITED");
      atom->addProperty(data);
    }
  }

  bool BaseAtomVisitor::isVisited(Atom *atom) {
    GenericData *data;
    data = atom->getPropertyByName("VISITED");
    return data == NULL ? false : true;
  }

  //------------------------------------------------------------------------//
	
  void DefaultAtomVisitor::visit(Atom *atom) {
    AtomData *atomData;
    AtomInfo *atomInfo;
    Vector3d  pos;
    Vector3d  vel;
    Vector3d  frc;
    Vector3d  u;

    if (isVisited(atom))
      return;
    
    atomInfo = new AtomInfo;
    
    atomData = new AtomData;
    atomData->setID("ATOMDATA");
    
    pos = atom->getPos();
    vel = atom->getVel();
    frc = atom->getFrc();
    atomInfo->atomTypeName = atom->getType();
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->vel[0] = vel[0];
    atomInfo->vel[1] = vel[1];
    atomInfo->vel[2] = vel[2];
    atomInfo->hasVelocity = true;
    atomInfo->frc[0] = frc[0];
    atomInfo->frc[1] = frc[1];
    atomInfo->frc[2] = frc[2];
    atomInfo->hasForce = true;
    atomInfo->vec[0] = 0.0;
    atomInfo->vec[1] = 0.0;
    atomInfo->vec[2] = 0.0;
    
    atomData->addAtomInfo(atomInfo);
    
    atom->addProperty(atomData);
    
    setVisited(atom);
  }
  
  void DefaultAtomVisitor::visit(DirectionalAtom *datom) {
    AtomData *atomData;
    AtomInfo *atomInfo;
    Vector3d  pos;
    Vector3d  vel;
    Vector3d  frc;
    Vector3d  u;

    if (isVisited(datom))
      return;
    
    pos = datom->getPos();
    vel = datom->getVel();
    frc = datom->getFrc();

    GayBerneAdapter gba = GayBerneAdapter(datom->getAtomType());
    MultipoleAdapter ma = MultipoleAdapter(datom->getAtomType());
    
    if (gba.isGayBerne()) {
      u = datom->getA().transpose()*V3Z;         
    } else if (ma.isDipole()) {
      u = datom->getDipole();
    } else if (ma.isQuadrupole()) {
      u = datom->getQuadrupole().getColumn(2);
    }
    atomData = new AtomData;
    atomData->setID("ATOMDATA");
    atomInfo = new AtomInfo;

    atomInfo->atomTypeName = datom->getType();
    atomInfo->pos[0] = pos[0];
    atomInfo->pos[1] = pos[1];
    atomInfo->pos[2] = pos[2];
    atomInfo->vel[0] = vel[0];
    atomInfo->vel[1] = vel[1];
    atomInfo->vel[2] = vel[2];
    atomInfo->hasVelocity = true;
    atomInfo->frc[0] = frc[0];
    atomInfo->frc[1] = frc[1];
    atomInfo->frc[2] = frc[2];
    atomInfo->hasForce = true;
    atomInfo->vec[0] = u[0];
    atomInfo->vec[1] = u[1];
    atomInfo->vec[2] = u[2];
    atomInfo->hasVector = true;

    atomData->addAtomInfo(atomInfo);

    datom->addProperty(atomData);

    setVisited(datom);
  }

  const std::string DefaultAtomVisitor::toString() {
    char   buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: copy atom infomation into atom data\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
  }
} //namespace OpenMD
