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

#include "selection/DistanceFinder.hpp"
#include "primitives/Molecule.hpp"
namespace oopse {

DistanceFinder::DistanceFinder(SimInfo* info) : info_(info) {

    nStuntDoubles_ = info_->getNGlobalAtoms() + info_->getNGlobalRigidBodies();
    stuntdoubles_.resize(nStuntDoubles_);
    
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;

    
    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
        
        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            stuntdoubles_[atom->getGlobalIndex()] = atom;
        }

        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
            stuntdoubles_[rb->getGlobalIndex()] = rb;
        }
        
    }    

}

BitSet DistanceFinder::find(const BitSet& bs, double distance) {
    StuntDouble * center;
    Vector3d centerPos;
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    BitSet bsResult(nStuntDoubles_);
    assert(bsResult.size() == bs.size());
    
    for (int i = bs.nextOnBit(-1); i != -1; i = bs.nextOnBit(i)) {
        center = stuntdoubles_[i];
        centerPos = center->getPos();
        for (int j = 0; j < stuntdoubles_.size(); ++j) {
            Vector3d r =centerPos - stuntdoubles_[j]->getPos();
            currSnapshot->wrapVector(r);
            if (r.length() <= distance) {
                bsResult.setBitOn(j);
            }
        }
    }

    return bsResult;
}

}
