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
  enum CutoffMethod {
    HARD,
    SWITCHED,
    SHIFTED_POTENTIAL,
    SHIFTED_FORCE
  };

  /**
   * @class InteractionManager InteractionManager is responsible for
   * keeping track of the non-bonded interactions (C++) and providing
   * an interface to the low-level loop (Fortran). 
   */
  class InteractionManager {

  public:
    static InteractionManager* Instance();
    static void setSimInfo(SimInfo* info) {info_ = info;} 
    static void initialize();

    // Fortran support routines

    static void doPrePair(InteractionData idat);
    static void doPreForce(SelfData sdat);
    static void doPair(InteractionData idat);    
    static void doSkipCorrection(InteractionData idat);
    static void doSelfCorrection(SelfData sdat);
    static RealType getSuggestedCutoffRadius(int *atid1);   
    static RealType getSuggestedCutoffRadius(AtomType *atype); 
    SwitchingFunction* getSwitchingFunction() {return switcher_;}
    
  private:
    virtual ~InteractionManager() { }
    // singleton pattern, prevent reconstruction
    InteractionManager() { }
    InteractionManager(InteractionManager const&) {};
    InteractionManager& operator=(InteractionManager const&) {};
    static InteractionManager* _instance; 

    static bool initialized_; 

    static void setupCutoffs();
    static void setupSwitching();
    static void setupElectrostatics();

    static SimInfo* info_;
    static LJ* lj_;
    static GB* gb_;
    static Sticky* sticky_;
    static EAM* eam_;
    static SC* sc_;
    static Morse* morse_;
    static Electrostatic* electrostatic_;
    static MAW* maw_;
    static SwitchingFunction* switcher_;

    static RealType rCut_;            /**< cutoff radius for non-bonded interactions */
    static RealType rSwitch_;         /**< inner radius of switching function */
    static CutoffMethod cutoffMethod_;/**< Cutoff Method for most non-bonded interactions */
    static SwitchingFunctionType sft_;/**< Type of switching function in use */

    static RealType vdwScale_[4];
    static RealType electrostaticScale_[4];
   
    static map<int, AtomType*> typeMap_;
    /**
     * Each pair of atom types can have multiple interactions, so the 
     * natural data structures are a map between the pair, and a set
     * of non-bonded interactions.
     */
    static map<pair<AtomType*, AtomType*>, set<NonBondedInteraction*> > interactions_;    
  };
}
#endif
