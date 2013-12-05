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

#include "selection/HullFinder.hpp"
#include "primitives/Molecule.hpp"
#include "math/ConvexHull.hpp"

namespace OpenMD {

  HullFinder::HullFinder(SimInfo* info) : info_(info) {

    nObjects_.push_back(info_->getNGlobalAtoms()+info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());

    stuntdoubles_.resize(nObjects_[STUNTDOUBLE]);
    bonds_.resize(nObjects_[BOND]);
    bends_.resize(nObjects_[BEND]);
    torsions_.resize(nObjects_[TORSION]);
    inversions_.resize(nObjects_[INVERSION]);
    
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
        
      // Hull is constructed from all known integrable objects.
      for (sd = mol->beginIntegrableObject(ioi);
           sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        localSites_.push_back(sd);
      }
      
      // selection can include atoms (which may be a subset of the IOs)
      for(atom = mol->beginAtom(ai); atom != NULL; 
          atom = mol->nextAtom(ai)) {
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

  HullFinder::~HullFinder() {
    delete surfaceMesh_;
  }

  SelectionSet HullFinder::findHull() {
    SelectionSet ssResult(nObjects_);
#ifdef HAVE_QHULL
    surfaceMesh_->computeHull(localSites_);
#else
    sprintf( painCave.errMsg,
             "HullFinder : Hull calculation is not possible without libqhull.\n"
             "\tPlease rebuild OpenMD with qhull enabled.");
    painCave.severity = OPENMD_ERROR;
    painCave.isFatal = 1;
    simError();
#endif
    
    std::vector<Triangle> sMesh = surfaceMesh_->getMesh();
    // Loop over the mesh faces
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;

    // This will work in parallel because the triangles returned by the mesh
    // have a NULL stuntDouble if this processor doesn't own the 
    
    for (face = sMesh.begin(); face != sMesh.end(); ++face) {
      Triangle thisTriangle = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle.getVertices();
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex) {
        if ((*vertex) != NULL) {
          ssResult.bitsets_[STUNTDOUBLE].setBitOn((*vertex)->getGlobalIndex());
        }
      }
    }
    return ssResult;
  }

  SelectionSet HullFinder::findHull(int frame) {
    SelectionSet ssResult(nObjects_);
#ifdef HAVE_QHULL
    surfaceMesh_->computeHull(localSites_);
#else
    sprintf( painCave.errMsg,
             "HullFinder : Hull calculation is not possible without libqhull.\n"
             "\tPlease rebuild OpenMD with qhull enabled.");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
#endif
    
    std::vector<Triangle> sMesh = surfaceMesh_->getMesh();
    // Loop over the mesh faces
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;

    // This will work in parallel because the triangles returned by the mesh
    // have a NULL stuntDouble if this processor doesn't own the 
    
    for (face = sMesh.begin(); face != sMesh.end(); ++face) {
      Triangle thisTriangle = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle.getVertices();
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex) {
        if ((*vertex) != NULL) {
          ssResult.bitsets_[STUNTDOUBLE].setBitOn((*vertex)->getGlobalIndex());
        }
      }
    }
    surfaceArea_ = surfaceMesh_->getArea();
    return ssResult;
  }

}
