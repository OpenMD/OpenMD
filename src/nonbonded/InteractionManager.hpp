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
#include "UseTheForce/ForceField.hpp"
#include "nonbonded/LJ.hpp"
#include "nonbonded/GB.hpp"
#include "nonbonded/Sticky.hpp"
#include "nonbonded/EAM.hpp"
#include "nonbonded/SC.hpp"
#include "nonbonded/Morse.hpp"
#include "nonbonded/Electrostatic.hpp"

using namespace std;

namespace OpenMD {
  /**
   * @class InteractionManager InteractionManager is responsible for
   * keeping track of the non-bonded interactions (C++) and providing
   * an interface to the low-level loop (Fortran). 
   */
  class InteractionManager {

  public:
    static InteractionManager* Instance();
    static void setForceField(ForceField *ff) {forceField_ = ff;};    



    // Fortran support routines

    static void doPrePair(int *atid1, int *atid2, RealType *rij, RealType *rho_i_at_j, RealType *rho_j_at_i);
    static void doPreForce(int *atid, RealType *rho, RealType *frho, RealType *dfrhodrho);
    static void doPair(int *atid1, int *atid2, RealType *d, RealType *r, RealType *r2, RealType *rcut, RealType *sw, RealType *vdwMult,RealType *electroMult, RealType *pot, RealType *vpair, RealType *f1, RealType *eFrame1, RealType *eFrame2, RealType *A1, RealType *A2, RealType *t1, RealType *t2, RealType *rho1, RealType *rho2, RealType *dfrho1, RealType *dfrho2, RealType *fshift1, RealType *fshift2);    
    static void doSkipCorrection(int *atid1, int *atid2, RealType *d, RealType *r, RealType *skippedCharge1, RealType *skippedCharge2, RealType *sw, RealType *electroMult, RealType *pot, RealType *vpair, RealType *f1, RealType *eFrame1, RealType *eFrame2, RealType *t1, RealType *t2);
    static void doSelfCorrection(int *atid, RealType *eFrame, RealType *skippedCharge, RealType *pot, RealType *t);
    static RealType getSuggestedCutoffRadius(int *atid1);   
    
  private:
    virtual ~InteractionManager() { }
    // singleton pattern, prevent reconstruction
    InteractionManager() { }
    InteractionManager(InteractionManager const&) {};
    InteractionManager& operator=(InteractionManager const&) {};
    static InteractionManager* _instance; 

    static void initialize();
    static bool initialized_; 

    static ForceField* forceField_;
    static LJ* lj_;
    static GB* gb_;
    static Sticky* sticky_;
    static EAM* eam_;
    static SC* sc_;
    static Morse* morse_;
    static Electrostatic* electrostatic_;

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
