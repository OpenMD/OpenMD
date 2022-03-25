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

#include "selection/HullFinder.hpp"

#include "math/AlphaHull.hpp"
#include "math/ConvexHull.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  HullFinder::HullFinder(SimInfo* info) : info_(info) {
    nObjects_.push_back(info_->getNGlobalAtoms() +
                        info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());
    nObjects_.push_back(info_->getNGlobalMolecules());

    stuntdoubles_.resize(nObjects_[STUNTDOUBLE]);
    bonds_.resize(nObjects_[BOND]);
    bends_.resize(nObjects_[BEND]);
    torsions_.resize(nObjects_[TORSION]);
    inversions_.resize(nObjects_[INVERSION]);
    molecules_.resize(nObjects_[MOLECULE]);

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    StuntDouble* sd;
    Atom* atom;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      molecules_[mol->getGlobalIndex()] = mol;

      // Hull is constructed from all known integrable objects.
      for (sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        localSites_.push_back(sd);
      }

      // selection can include atoms (which may be a subset of the IOs)
      for (atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        stuntdoubles_[atom->getGlobalIndex()] = atom;
      }

      // and rigid bodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL;
           rb = mol->nextRigidBody(rbIter)) {
        stuntdoubles_[rb->getGlobalIndex()] = rb;
      }

      // These others are going to be inferred from the objects on the hull:
      for (bond = mol->beginBond(bondIter); bond != NULL;
           bond = mol->nextBond(bondIter)) {
        bonds_[bond->getGlobalIndex()] = bond;
      }
      for (bend = mol->beginBend(bendIter); bend != NULL;
           bend = mol->nextBend(bendIter)) {
        bends_[bend->getGlobalIndex()] = bend;
      }
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL;
           torsion = mol->nextTorsion(torsionIter)) {
        torsions_[torsion->getGlobalIndex()] = torsion;
      }
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL;
           inversion = mol->nextInversion(inversionIter)) {
        inversions_[inversion->getGlobalIndex()] = inversion;
      }
    }
#ifdef HAVE_QHULL
    surfaceMesh_ = new ConvexHull();
#endif
  }

  AlphaHullFinder::AlphaHullFinder(SimInfo* info) : HullFinder(info) {
#ifdef HAVE_QHULL
    delete surfaceMesh_;
    surfaceMesh_ = new AlphaHull(0.0);
#endif
  }

  void AlphaHullFinder::setAlpha(RealType alpha) {
#ifdef HAVE_QHULL
    delete surfaceMesh_;
    surfaceMesh_ = new AlphaHull(alpha);
#endif
  }

  HullFinder::~HullFinder() {
#ifdef HAVE_QHULL
    delete surfaceMesh_;
#endif
  }

  SelectionSet HullFinder::findHull() {
    SelectionSet ssResult(nObjects_);
#ifdef HAVE_QHULL
    surfaceMesh_->computeHull(localSites_);

    std::vector<Triangle> sMesh = surfaceMesh_->getMesh();
    // Loop over the mesh faces
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;

    // This will work in parallel because the triangles returned by the mesh
    // have a NULL stuntDouble if this processor doesn't own

    for (face = sMesh.begin(); face != sMesh.end(); ++face) {
      Triangle thisTriangle               = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle.getVertices();
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex) {
        if ((*vertex) != NULL) {
          ssResult.bitsets_[STUNTDOUBLE].setBitOn((*vertex)->getGlobalIndex());
        }
      }
    }
    surfaceArea_ = surfaceMesh_->getArea();
#else
    sprintf(painCave.errMsg,
            "HullFinder : Hull calculation is not possible without libqhull.\n"
            "\tPlease rebuild OpenMD with qhull enabled.");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
#endif

    return ssResult;
  }

  SelectionSet HullFinder::findHull(int) {
    SelectionSet ssResult(nObjects_);
#ifdef HAVE_QHULL
    surfaceMesh_->computeHull(localSites_);
    std::vector<Triangle> sMesh = surfaceMesh_->getMesh();
    // Loop over the mesh faces
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;

    // This will work in parallel because the triangles returned by the mesh
    // have a NULL stuntDouble if this processor doesn't own the

    for (face = sMesh.begin(); face != sMesh.end(); ++face) {
      Triangle thisTriangle               = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle.getVertices();
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex) {
        if ((*vertex) != NULL) {
          ssResult.bitsets_[STUNTDOUBLE].setBitOn((*vertex)->getGlobalIndex());
        }
      }
    }
    surfaceArea_ = surfaceMesh_->getArea();
    volume_      = surfaceMesh_->getVolume();

#else
    sprintf(painCave.errMsg,
            "HullFinder : Hull calculation is not possible without libqhull.\n"
            "\tPlease rebuild OpenMD with qhull enabled.");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal  = 1;
    simError();
#endif

    return ssResult;
  }
}  // namespace OpenMD
