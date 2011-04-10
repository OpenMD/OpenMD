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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#ifndef PARALLEL_DECOMPOSITION_HPP
#define PARALLEL_DECOMPOSITION_HPP

#include "brains/SimInfo.hpp"
#include "nonbonded/NonBondedInteraction.hpp"

using namespace std;
namespace OpenMD {
  
  /**
   * @class Decomposition 
   * Decomposition is an interface for passing out and collecting information
   * from many processors at various stages of the main non-bonded ForceLoop.
   *
   * The pairwise force calculation has an outer-running loop (the "I"
   * loop) and an inner-running loop (the "J" loop).  In parallel
   * decompositions, these loop over different groups of atoms on
   * different processors.  Between each set of computations on the
   * local processor, data must be exchanged among the processors.
   * This can happen at different times in the calculation:
   *
   *  distributeInitialData      (parallel communication - one time only)
   *  distributeData             (parallel communication - every ForceLoop)
   *
   *  loop iLoop over nLoops     (nLoops may be 1, 2, or until self consistent)
   *  |  loop over i
   *  |  | loop over j
   *  |  | | localComputation
   *  |  |  end
   *  |  end
   *  |  if (nLoops > 1):
   *  |  |   collectIntermediateData    (parallel communication)
   *  |  |   distributeIntermediateData (parallel communication)
   *  |  endif
   *  end
   * collectData                        (parallel communication)
   *
   * Decomposition provides the interface for ForceLoop to do the
   * communication steps and to iterate using the correct set of atoms
   * and cutoff groups.
   */
  class Decomposition {
  public:

    Decomposition(SimInfo* info) : info_(info) {}
    virtual ~Decomposition() {}
    
    virtual void distributeInitialData() = 0;
    virtual void distributeData() = 0;
    virtual void collectIntermediateData() = 0;
    virtual void distributeIntermediateData() = 0;
    virtual void collectData() = 0;

    // neighbor list routines
    virtual bool checkNeighborList() = 0;
    virtual vector<pair<int, int> >  buildNeighborList() = 0;

    // group bookkeeping
    virtual pair<int, int> getGroupTypes(int cg1, int cg2) = 0;

    // Group->atom bookkeeping
    virtual vector<int> getAtomsInGroupI(int cg1) = 0;
    virtual vector<int> getAtomsInGroupJ(int cg2) = 0;
    virtual Vector3d getAtomToGroupVectorI(int atom1, int cg1) = 0;
    virtual Vector3d getAtomToGroupVectorJ(int atom2, int cg2) = 0;
    virtual RealType getMfactI(int atom1) = 0;
    virtual RealType getMfactJ(int atom2) = 0;

    // spatial data
    virtual Vector3d getIntergroupVector(int cg1, int cg2) = 0;
    virtual Vector3d getInteratomicVector(int atom1, int atom2) = 0;
       
    // atom bookkeeping
    virtual vector<int> getAtomList() = 0;
    virtual vector<int> getSkipsForAtom(int atom1) = 0;
    virtual bool skipAtomPair(int atom1, int atom2) = 0;
    virtual void addForceToAtomI(int atom1, Vector3d fg) = 0;
    virtual void addForceToAtomJ(int atom2, Vector3d fg) = 0;

    // filling interaction blocks with pointers
    virtual InteractionData fillInteractionData(int atom1, int atom2) = 0;
    virtual InteractionData fillSkipData(int atom1, int atom2) = 0;
    virtual SelfData fillSelfData(int atom1) = 0;
    
  protected:
    SimInfo* info_;   
    map<pair<int, int>, int> topoDist; //< topoDist gives the
                                       //topological distance between
                                       //two atomic sites.  This
                                       //declaration is agnostic
                                       //regarding the parallel
                                       //decomposition.  The two
                                       //indices could be local or row
                                       //& column.  It will be up to
                                       //the specific decomposition
                                       //method to fill this.
    map<pair<int, int>, bool> exclude; //< exclude is the set of pairs
                                       //to leave out of non-bonded
                                       //force evaluations.  This
                                       //declaration is agnostic
                                       //regarding the parallel
                                       //decomposition.  The two
                                       //indices could be local or row
                                       //& column.  It will be up to
                                       //the specific decomposition
                                       //method to fill this.
  };    
}
#endif
