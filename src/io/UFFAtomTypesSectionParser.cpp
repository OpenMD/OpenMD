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


#include "io/UFFAtomTypesSectionParser.hpp"
#include "types/UFFAdapter.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"

namespace OpenMD {
  
  UFFAtomTypesSectionParser::UFFAtomTypesSectionParser(ForceFieldOptions& options) :
    options_(options){
    setSectionName("UFFAtomTypes");
  }
  
  void UFFAtomTypesSectionParser::parseLine(ForceField& ff,const std::string& line,
                                            int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();    
    
    // in UFFAtomTypesSectionParser, a line at least contains 12 tokens:
    // Atom r1 theta0 x1 D1 zeta Z1 Vi Uj Xi Hard Radius
    
    if (nTokens < 12)  {
      sprintf(painCave.errMsg, "UFFAtomTypesSectionParser Error: "
              "Not enough tokens at line %d\n", lineNo);
      painCave.isFatal = 1;
      simError();                    
    } else {
      
      std::string atomTypeName = tokenizer.nextToken();    
      AtomType* atomType = ff.getAtomType(atomTypeName);
      
      if (atomType != NULL) {
        UFFAdapter uffa = UFFAdapter(atomType);
        
        RealType r1 = tokenizer.nextTokenAsDouble();    
        RealType theta0 = tokenizer.nextTokenAsDouble();
        RealType x1 = tokenizer.nextTokenAsDouble();    
        RealType D1 = tokenizer.nextTokenAsDouble();    
        RealType zeta = tokenizer.nextTokenAsDouble();  
        RealType Z1 = tokenizer.nextTokenAsDouble();    
        RealType Vi = tokenizer.nextTokenAsDouble();    
        RealType Uj = tokenizer.nextTokenAsDouble();    
        RealType Xi = tokenizer.nextTokenAsDouble();    
        RealType Hard = tokenizer.nextTokenAsDouble();  
        RealType Radius = tokenizer.nextTokenAsDouble();

        r1 *= options_.getDistanceUnitScaling();
        theta0 *= options_.getAngleUnitScaling();
        x1 *= options_.getDistanceUnitScaling();
        D1 *= options_.getEnergyUnitScaling();
        Z1 *= options_.getChargeUnitScaling();
        Vi *= options_.getEnergyUnitScaling();
        Uj *= options_.getEnergyUnitScaling();

        uffa.makeUFF(r1,  theta0,  x1,  D1, zeta,  Z1,  Vi,  Uj, Xi,  Hard,  Radius);

      } else {
        sprintf(painCave.errMsg, "UFFAtomTypesSectionParser Error: Atom Type [%s]"
                " has not been created yet\n", atomTypeName.c_str());
        painCave.isFatal = 1;
        simError();    
      }      
    }        
  }

  void UFFAtomTypesSectionParser::validateSection(ForceField& ff) {
    // ForceField::AtomTypeContainer* atomTypes = ff.getAtomTypes();
    // ForceField::AtomTypeContainer::MapTypeIterator i;
    // AtomType* at;

    // std::vector<AtomType*> uffTypes;
    
    // for (at = atomTypes->beginType(i); at != NULL; at = atomTypes->nextType(i)) {      

    //   UFFAdapter uffa = UFFAdapter(at);      
    //   if (uffa.isUFF())
    //     uffTypes.push_back(at);
    // }

    // for (std::size_t i = 0; i < uffTypes.size(); ++i) {
    //   RealType ri = uffTypes[i]->getR1();
    //   RealType chiI = uffTypes[i]->getXi();
    //   RealType Ra = uffTypes[i]->getX1();
    //   RealType ka = uffTypes[i]->getD1();
      
    //   for (std::size_t j = 0; j < uffTypes.size(); ++j) {
    //     RealType rj = uffTypes[j]->getR1();
    //     RealType chiJ = uffTypes[j]->getXi();
    //     RealType Rb = uffTypes[j]->getX1();
    //     RealType kb = uffTypes[j]->getD1();

    //     // Precompute the equilibrium bond distance
    //     // From equation 3
    //     rbo = -0.1332*(ri+rj)*log(bondorder);
    //     // From equation 4
    //     ren = ri*rj*(pow((sqrt(chiI) - sqrt(chiJ)),2.0)) / (chiI*ri + chiJ*rj);
    //     // From equation 2
    //     // NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
    //     // There is a typo in the published paper
    //     rij = ri + rj + rbo - ren;

    //     kab = sqrt(ka * kb);

    //     // ka now represents the xij in equation 20 -- the expected vdw distance
    //     kaSquared = (Ra * Rb);
    //     ka = sqrt(kaSquared);
      
  }
} //end namespace OpenMD

