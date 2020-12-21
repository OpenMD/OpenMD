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

#include "io/StickyAtomTypesSectionParser.hpp"
#include "types/StickyAdapter.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"
using namespace std;
namespace OpenMD {
  
  StickyAtomTypesSectionParser::StickyAtomTypesSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("StickyAtomTypes");
  }
  
  void StickyAtomTypesSectionParser::parseLine(ForceField& ff,
                                               const string& line, 
                                               int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();    
    RealType dus = options_.getDistanceUnitScaling();
    RealType eus = options_.getEnergyUnitScaling();

    //in AtomTypeSection, a line at least contains 8 tokens
    //atomTypeName and 7 different sticky parameters
    if (nTokens < 8)  {
      sprintf(painCave.errMsg, "StickyAtomTypesSectionParser Error: Not enough "
              "tokens at line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();                      
    } else {
      
      string atomTypeName = tokenizer.nextToken();    
      AtomType* atomType = ff.getAtomType(atomTypeName);
      
      if (atomType != NULL) {
        
        StickyAdapter sa = StickyAdapter(atomType);

        RealType w0 = tokenizer.nextTokenAsDouble();
        RealType v0 = eus * tokenizer.nextTokenAsDouble();
        RealType v0p = eus * tokenizer.nextTokenAsDouble();
        RealType rl = dus * tokenizer.nextTokenAsDouble();
        RealType ru = dus * tokenizer.nextTokenAsDouble();
        RealType rlp = dus * tokenizer.nextTokenAsDouble();
        RealType rup = dus * tokenizer.nextTokenAsDouble();   
        bool isPower = false;

        sa.makeSticky(w0, v0, v0p, rl, ru, rlp, rup, isPower);
        
      } else {
        sprintf(painCave.errMsg, "StickyAtomTypesSectionParser Error: "
                "Can not find matching AtomType %s\n",
                atomTypeName.c_str());
        painCave.isFatal = 1;
        simError();     
      }
    }         
  }    
} //end namespace OpenMD

