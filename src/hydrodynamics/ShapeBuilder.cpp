/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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
#include "hydrodynamics/ShapeBuilder.hpp"
#include "hydrodynamics/CompositeShape.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "types/GayBerneAdapter.hpp"
#include "types/LennardJonesAdapter.hpp"

namespace OpenMD {

  Shape* ShapeBuilder::createShape(StuntDouble* sd) {
    Shape* currShape = NULL;
    if (sd->isDirectionalAtom()) {
      currShape = internalCreateShape(static_cast<DirectionalAtom*>(sd));
    } else if (sd->isAtom()) {
      currShape = internalCreateShape(static_cast<Atom*>(sd));
    } else if (sd->isRigidBody()) {
      currShape = internalCreateShape(static_cast<RigidBody*>(sd));
    }
    return currShape;
  }

  Shape* ShapeBuilder::createShape(Molecule* mol) {
    Molecule::IntegrableObjectIterator j;
    StuntDouble* sd;

    if (mol->getNAtoms() == 1) {
       Shape* molShape = createShape(mol->getAtomAt(0));
       return molShape;
    } else {
      CompositeShape* molShape = new CompositeShape();
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j)) {
	Shape* currShape  = createShape(sd);
	if (currShape != NULL) {	      
	  molShape->addShape(currShape);
	}
      }
      molShape->setName( mol->getType() );
      return molShape;	   
    }
    return NULL;
  }

  Shape* ShapeBuilder::internalCreateShape(Atom* atom) {
    AtomType* atomType      = atom->getAtomType();
    Shape* currShape        = NULL;
    LennardJonesAdapter lja = LennardJonesAdapter(atomType);
    std::cerr << atom->getType() << "\n";
    if (lja.isLennardJones()) {
      currShape = new Sphere(atom->getPos(), lja.getSigma() / 2.0);
      currShape->setName( atom->getType() );
    } else {
      int obanum(0);
      std::vector<AtomType*> atChain = atomType->allYourBase();
      std::vector<AtomType*>::iterator i;
      for (i = atChain.begin(); i != atChain.end(); ++i) {
        obanum = etab.GetAtomicNum((*i)->getName().c_str());
        if (obanum != 0) {
          currShape = new Sphere(atom->getPos(), etab.GetVdwRad(obanum));
	  currShape->setName( atom->getType() );
          break;
        }
      }
      if (obanum == 0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                "Could not find atom type in default element.txt\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    }
    return currShape;
  }

  Shape* ShapeBuilder::internalCreateShape(DirectionalAtom* datom) {
    AtomType* atomType      = datom->getAtomType();
    Shape* currShape        = NULL;
    LennardJonesAdapter lja = LennardJonesAdapter(atomType);
    GayBerneAdapter gba     = GayBerneAdapter(atomType);
    if (gba.isGayBerne()) {
      currShape = new Ellipsoid(datom->getPos(), gba.getL() / 2.0,
                                gba.getD() / 2.0, datom->getA());
    } else if (lja.isLennardJones()) {
      currShape = new Sphere(datom->getPos(), lja.getSigma() / 2.0);
    } else {
      int obanum = etab.GetAtomicNum((datom->getType()).c_str());
      if (obanum != 0) {
        currShape = new Sphere(datom->getPos(), etab.GetVdwRad(obanum));
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                "Could not find atom type in default element.txt\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    }
    currShape->setName( datom->getType() );
    return currShape;
  }

  Shape* ShapeBuilder::internalCreateShape(RigidBody* rb) {
    std::vector<Atom*>::iterator ai;
    CompositeShape* compositeShape = new CompositeShape;
    Atom* atom;
    for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
      Shape* currShape = NULL;
      if (atom->isDirectionalAtom()) {
        currShape = internalCreateShape(static_cast<DirectionalAtom*>(atom));
      } else if (atom->isAtom()) {
        currShape = internalCreateShape(static_cast<Atom*>(atom));
      }
      if (currShape != NULL) {
	compositeShape->addShape(currShape);
	compositeShape->setName( rb->getType() );
      }
    }
    return compositeShape;
  }
}  // namespace OpenMD
