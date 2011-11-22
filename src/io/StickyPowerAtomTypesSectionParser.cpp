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

#include "io/StickyPowerAtomTypesSectionParser.hpp"
#include "types/AtomType.hpp"
#include "types/DirectionalAtomType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "utils/simError.h"
namespace OpenMD {
  
  StickyPowerAtomTypesSectionParser::StickyPowerAtomTypesSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("StickyPowerAtomTypes");
  }
  
  void StickyPowerAtomTypesSectionParser::parseLine(ForceField& ff,
                                               const std::string& line, 
                                               int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();    
    
    //in AtomTypeSection, a line at least contains 8 tokens
    //atomTypeName and 7 different sticky parameters
    if (nTokens < 8)  {
      sprintf(painCave.errMsg, 
              "StickyPowerAtomTypesSectionParser Error: Not enough tokens at "
              "line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();                      
    } else {
      
      std::string atomTypeName = tokenizer.nextToken();    
      AtomType* atomType = ff.getAtomType(atomTypeName);
      
      if (atomType != NULL) {
        DirectionalAtomType* dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);
        
        if (dAtomType != NULL) {
          StickyParam sticky; 
          sticky.w0 = tokenizer.nextTokenAsDouble();
          sticky.v0 = tokenizer.nextTokenAsDouble();
          sticky.v0p = tokenizer.nextTokenAsDouble();
          sticky.rl = tokenizer.nextTokenAsDouble();
          sticky.ru = tokenizer.nextTokenAsDouble();
          sticky.rlp = tokenizer.nextTokenAsDouble();
          sticky.rup = tokenizer.nextTokenAsDouble();   
          
          dAtomType->addProperty(new StickyParamGenericData("Sticky", sticky));
          dAtomType->setStickyPower();
        } else {
          sprintf(painCave.errMsg, 
                  "StickyPowerAtomTypesSectionParser Error: Not enough tokens "
                  "at line %d\n",
                  lineNo);
          painCave.isFatal = 1;
          simError();            
          std::cerr << "StickyPowerAtomTypesSectionParser Warning:" 
            << std::endl;
        }
      } else {
        sprintf(painCave.errMsg, 
                "StickyPowerAtomTypesSectionParser Error: Can not find matched "
                "AtomType %s\n",
                atomTypeName.c_str());
        painCave.isFatal = 1;
        simError();    
        
      }
      
    }    
    
  }
    
} //end namespace OpenMD



