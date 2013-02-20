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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "io/BaseAtomTypesSectionParser.hpp"
#include "types/AtomType.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"
namespace OpenMD {

  BaseAtomTypesSectionParser::BaseAtomTypesSectionParser() {
    setSectionName("BaseAtomTypes");
  }

  void BaseAtomTypesSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();    

    //in BaseAtomTypeSection, a line at least contains 2 tokens
    //atomTypeName and mass
    if (nTokens < 2)  {
      sprintf(painCave.errMsg, "BaseAtomTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
            
    } else {

      std::string baseAtomTypeName = tokenizer.nextToken();    
      AtomType* baseAtomType = ff.getAtomType(baseAtomTypeName);

      if (baseAtomType == NULL) {
      	baseAtomType = new AtomType();
      	int ident = ff.getNAtomType() + 1;
      	baseAtomType->setIdent(ident); 
      	baseAtomType->setName(baseAtomTypeName);
      	ff.addAtomType(baseAtomTypeName, baseAtomType);
      } 

   //    else {
//         sprintf(painCave.errMsg, "BaseAtomTypesSectionParser Error: Duplicate BaseAtomType at line %d\n",
//                 lineNo);
//         painCave.isFatal = 1;
//         simError();
//       }
        
      RealType mass = tokenizer.nextTokenAsDouble();              
      baseAtomType->setMass(mass);
      if (tokenizer.hasMoreTokens()) {
          RealType nelectron = tokenizer.nextTokenAsDouble();
          baseAtomType->addProperty(new DoubleGenericData("nelectron", nelectron));
      }               
    }    


  }

} //end namespace OpenMD
