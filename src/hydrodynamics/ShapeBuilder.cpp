/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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
        Shape* currShape = createShape(sd);
        if (currShape != NULL) { molShape->addShape(currShape); }
      }
      molShape->setName(mol->getType());
      return molShape;
    }
    return NULL;
  }

  Shape* ShapeBuilder::internalCreateShape(Atom* atom) {
    AtomType* atomType      = atom->getAtomType();
    Shape* currShape        = NULL;
    LennardJonesAdapter lja = LennardJonesAdapter(atomType);
    if (lja.isLennardJones()) {
      currShape = new Sphere(atom->getPos(), lja.getSigma() / 2.0);
      currShape->setName(atom->getType());
    } else {
      int obanum(0);
      std::vector<AtomType*> atChain = atomType->allYourBase();
      std::vector<AtomType*>::iterator i;
      for (i = atChain.begin(); i != atChain.end(); ++i) {
        obanum = etab.GetAtomicNum((*i)->getName().c_str());
        if (obanum != 0) {
          currShape = new Sphere(atom->getPos(), etab.GetVdwRad(obanum));
          currShape->setName(atom->getType());
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
    currShape->setName(datom->getType());
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
        compositeShape->setName(rb->getType());
      }
    }
    return compositeShape;
  }
}  // namespace OpenMD
