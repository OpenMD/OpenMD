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
 
#ifndef PRIMITIVES_CUTOFFGROUP_HPP

#define PRIMITIVES_CUTOFFGROUP_HPP

#include "primitives/Atom.hpp"
#include "math/Vector3.hpp"

namespace oopse {
  class CutoffGroup {
  public:

    CutoffGroup() {
      haveTotalMass = false;
      totalMass = 0.0;
    }

    void addAtom(Atom *atom) {
      cutoffAtomList.push_back(atom);
    }

    Atom *beginAtom(std::vector<Atom *>::iterator & i) {
      i = cutoffAtomList.begin();
      return i != cutoffAtomList.end() ? *i : NULL;
    }

    Atom *nextAtom(std::vector<Atom *>::iterator & i) {
      i++;
      return i != cutoffAtomList.end() ? *i : NULL;
    }

    std::vector<Atom*> getAtoms() { return cutoffAtomList; }
    double getMass() {
      std::vector<Atom *>::iterator i;
      Atom * atom;
      double mass;

      if (!haveTotalMass) {
	totalMass = 0;

	for(atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
	  mass = atom->getMass();
	  totalMass += mass;
	}

	haveTotalMass = true;
      }

      return totalMass;
    }

    void getCOM(Vector3d & com) {
      std::vector<Atom *>::iterator i;
      Atom * atom;
      Vector3d pos;
      double mass;

      com[0] = 0;
      com[1] = 0;
      com[2] = 0;
      totalMass = getMass();

      if (cutoffAtomList.size() == 1) {
	com = beginAtom(i)->getPos();
      } else {
	for(atom = beginAtom(i); atom != NULL; atom = nextAtom(i)) {
	  mass = atom->getMass();
	  pos = atom->getPos();
	  com += pos * mass;
	}

	com /= totalMass;
      }
    }

    int getNumAtom() {
      return cutoffAtomList.size();
    }

    int getGlobalIndex() {
      return globalIndex;
    }

    void setGlobalIndex(int id) {
      this->globalIndex = id;
    }

  private:

    std::vector<Atom *>cutoffAtomList;
    bool haveTotalMass;
    double totalMass;
    int globalIndex;
  };

}      //end namespace oopse

#endif //PRIMITIVES_CUTOFFGROUP_HPP  
