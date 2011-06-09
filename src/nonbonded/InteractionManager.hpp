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
 
#ifndef NONBONDED_INTERACTIONMANAGER_HPP
#define NONBONDED_INTERACTIONMANAGER_HPP

#include "brains/SimInfo.hpp"
#include "types/AtomType.hpp"
#include "nonbonded/LJ.hpp"
#include "nonbonded/GB.hpp"
#include "nonbonded/Sticky.hpp"
#include "nonbonded/EAM.hpp"
#include "nonbonded/SC.hpp"
#include "nonbonded/Morse.hpp"
#include "nonbonded/Electrostatic.hpp"
#include "nonbonded/MAW.hpp"
#include "nonbonded/SwitchingFunction.hpp"

using namespace std;

namespace OpenMD {

  /**
   * @class InteractionManager InteractionManager is responsible for
   * keeping track of the non-bonded interactions (C++)
   */
  class InteractionManager {

  public:
    InteractionManager();
    ~InteractionManager();
    void setSimInfo(SimInfo* info) {info_ = info;} 
    void initialize();

    // Fortran support routines

    void doPrePair(InteractionData idat);
    void doPreForce(SelfData sdat);
    void doPair(InteractionData idat);    
    void doSkipCorrection(InteractionData idat);
    void doSelfCorrection(SelfData sdat);
    RealType getSuggestedCutoffRadius(int *atid1);   
    RealType getSuggestedCutoffRadius(AtomType *atype);
    
  private:
    bool initialized_;

    void setupElectrostatics();

    SimInfo* info_;
    LJ* lj_;
    GB* gb_;
    Sticky* sticky_;
    EAM* eam_;
    SC* sc_;
    Morse* morse_;
    Electrostatic* electrostatic_;
    MAW* maw_;
    
    map<int, AtomType*> typeMap_;
    /**
     * Each pair of atom types can have multiple interactions, so the 
     * natural data structures are a map between the pair, and a set
     * of non-bonded interactions.
     */
    map<pair<AtomType*, AtomType*>, set<NonBondedInteraction*> > interactions_;    
  };
}
#endif
