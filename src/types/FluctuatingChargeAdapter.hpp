/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef TYPES_FLUCTUATINGCHARGEADAPTER_HPP
#define TYPES_FLUCTUATINGCHARGEADAPTER_HPP

#include "utils/GenericData.hpp"
#include "types/AtomType.hpp"

using namespace std;
namespace OpenMD {

  string const FQtypeID = "FlucQ";

  struct FluctuatingAtypeParameters {
    RealType chargeMass;
    bool hasMultipleMinima;
    RealType electronegativity; /** (relative) electronegativity */
    RealType hardness;     /** diagonal "Coulomb" Jii */
    int slaterN;           /** principal quantum number for Slater orbitals */
    RealType slaterZeta;   /** off-diagonal Slater exponent */
    RealType curvature;    /** multiple minima curvature for charge states */
    RealType coupling;     /** multiple minima EVB coupling for charge states */
    vector<pair<int, RealType> > diabaticStates; /** multiple minima diabats */
  };
  typedef SimpleTypeData<FluctuatingAtypeParameters*> FluctuatingAtypeData;   

  class FluctuatingChargeAdapter {
  public:
    FluctuatingChargeAdapter(AtomType* AT) { at_ = AT; };

    void makeFluctuatingCharge(RealType chargeMass, RealType electronegativity, RealType hardness, int slaterN, RealType slaterZeta);
    void makeFluctuatingCharge(RealType chargeMass, RealType electronegativity, RealType hardness, int slaterN);

    void makeFluctuatingCharge(RealType chargeMass, RealType curvature, RealType coupling, vector<pair<int, RealType> > diabaticStates);
    
    bool isFluctuatingCharge();
    bool hasMultipleMinima();
    RealType getChargeMass();
    RealType getElectronegativity();
    RealType getHardness();
    int getSlaterN();
    RealType getSlaterZeta();
    RealType getCoupling();
    RealType getCurvature();
    vector<pair<int, RealType> > getDiabaticStates();

  private:
    AtomType* at_;
    FluctuatingAtypeParameters* getFluctuatingChargeParam();    
  };
}
#endif
