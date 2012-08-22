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

#include "applications/hydrodynamics/BeadModel.hpp"
#include "types/LennardJonesAdapter.hpp"

namespace OpenMD {
  bool BeadModel::createBeads(std::vector<BeadParam>& beads) {
    
    if (sd_->isAtom()) {
      if (!createSingleBead(static_cast<Atom*>(sd_), beads)) {
        sprintf( painCave.errMsg,
                 "BeadModel::createBeads Error: GayBerne and other non-spheric atoms should use RoughShell model\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();    
        return false;
      }
    }
    else if (sd_->isRigidBody()) {
      RigidBody* rb = static_cast<RigidBody*>(sd_);
      std::vector<Atom*>::iterator ai; 
      Atom* atom;
      for (atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
        if (!createSingleBead(atom, beads)) {
          sprintf( painCave.errMsg,
                   "BeadModel::createBeads Error: GayBerne and other non-spheric atoms should use RoughShell model\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();    
          return false;
        }
      }
    }    
    return true;
  }
  
  bool BeadModel::createSingleBead(Atom* atom, std::vector<BeadParam>& beads) {
    AtomType* atomType = atom->getAtomType();
    LennardJonesAdapter lja = LennardJonesAdapter(atomType);
    if (atomType->isGayBerne()) {
      return false;
    } else if (lja.isLennardJones()){
      BeadParam currBead;
      currBead.atomName = atom->getType();
      currBead.pos = atom->getPos();
      currBead.radius = lja.getSigma()/2.0;
      beads.push_back(currBead);
    } else {
      int obanum = etab.GetAtomicNum((atom->getType()).c_str());
      if (obanum != 0) {
        BeadParam currBead;      
        currBead.atomName = atom->getType();
        currBead.pos = atom->getPos();        
        currBead.radius = etab.GetVdwRad(obanum);
        std::cout << "using rvdw = " << currBead.radius << " for atomic number " << obanum << "\n";
        beads.push_back(currBead);
      } else {
        sprintf( painCave.errMsg,
                 "Could not find atom type in default element.txt\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();          
      }
    }
    return true;    
  }
}
