 /*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */
 
#ifndef NONBONDED_SC_HPP
#define NONBONDED_SC_HPP

#include "nonbonded/NonBondedInteraction.hpp"
#include "brains/ForceField.hpp"
#include "math/CubicSpline.hpp"
#include "types/SuttonChenAdapter.hpp"

namespace OpenMD {

  struct SCAtomData {
    RealType c;
    RealType m;
    RealType n;
    RealType alpha;
    RealType epsilon;
    RealType rCut;
  };
  
  struct SCInteractionData {
    RealType alpha;
    RealType epsilon;
    RealType m;
    RealType n;
    RealType rCut;
    RealType vCut;
    CubicSpline* V;
    CubicSpline* phi;
    bool explicitlySet;
  };
    
  class SC : public MetallicInteraction {
    
  public:    
    SC();
    ~SC();
    void setForceField(ForceField *ff) {forceField_ = ff;};
    void setSimulatedAtomTypes(set<AtomType*> &simtypes) {simTypes_ = simtypes; initialize();};
    void addType(AtomType* atomType);
    void addExplicitInteraction(AtomType* atype1, AtomType* atype2, RealType epsilon, RealType m, RealType n, RealType alpha);
    void calcDensity(InteractionData &idat);
    void calcFunctional(SelfData &sdat);
    void calcForce(InteractionData &idat);
    virtual string getName() {return name_;}
    virtual int getHash() { return SC_INTERACTION; }
    virtual RealType getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes);
   
  private:
    void initialize();
    RealType getAlpha(AtomType* atomType1, AtomType* atomType2);
    RealType getEpsilon(AtomType* atomType1, AtomType* atomType2);
    RealType getM(AtomType* atomType1, AtomType* atomType2);
    RealType getN(AtomType* atomType1, AtomType* atomType2);
    
    string name_;
    bool initialized_;
    set<int> SCtypes;         /**< The set of AtomType idents that are SC types */
    vector<int> SCtids;       /**< The mapping from AtomType ident -> SC type ident */
    vector<SCAtomData> SCdata; /**< The EAM atomic data indexed by SC type ident */
    vector<vector<SCInteractionData> > MixingMap;  /**< The mixing parameters between two SC types */
    int nSC_;

    ForceField* forceField_;
    set<AtomType*> simTypes_;

    int np_;
    
  };
}

                               
#endif
