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

#include "applications/hydrodynamics/BeadModel.hpp"

namespace oopse {
bool BeadModel::createBeads(std::vector<BeadParam>& beads) {

    if (sd_->isAtom()) {
        createSingleBead(static_cast<Atom*>(sd_), beads);
    }
    else if (sd_->isRigidBody()) {
        RigidBody* rb = static_cast<RigidBody*>(sd_);
        std::vector<Atom*>::iterator ai; 
        Atom* atom;
        for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
            if (!createSingleBead(atom, beads))
                return false;
        }
    }

    return true;
}

bool BeadModel::createSingleBead(Atom* atom, std::vector<BeadParam>& beads) {
    AtomType* atomType = atom->getAtomType();

    if (atomType->isGayBerne()) {
        return false;
    } else if (atomType->isLennardJones()){
        GenericData* data = atomType->getPropertyByName("LennardJones");
        if (data != NULL) {
            LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);

            if (ljData != NULL) {
                LJParam ljParam = ljData->getData();
                BeadParam currBead;
                currBead.atomName = atom->getType();
                currBead.pos = atom->getPos();
                currBead.radius = ljParam.sigma/2.0;
                beads.push_back(currBead);
        } else {
            sprintf( painCave.errMsg,
            "Can not cast GenericData to LJParam\n");
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();          
            }       
        }
    }

    return true;
}

}
