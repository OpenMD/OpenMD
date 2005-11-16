/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 *
 *  OptionSectionParser.cpp
 *  OOPSE-2.0
 *
 *  Created by Charles F. Vardeman II on 11/15/05.
 *  @author  Charles F. Vardeman II 
 *  @version $Id: OptionSectionParser.cpp,v 1.2 2005-11-16 23:10:02 tim Exp $
 *
 */

#include "io/OptionSectionParser.hpp"
#include "types/AtomType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "utils/simError.h"
#include "utils/StringUtils.hpp"
namespace oopse {

  bool ForceFieldOptions::setData(const std::string& keyword, const std::string& value) {
      bool result;
      ParamMap::iterator i =parameters_.find(keyword);
      if (i != parameters_.end()) {
        if(isType<int>(value)){
          int ival = lexi_cast<int>(value);
          result = i->second->setData(ival);
        }      
        else if (isType<double>(value)){
          double dval = lexi_cast<double>(value);
          result = i->second->setData(dval);
        } else{
           result = i->second->setData(value);
        }
      } else {
        sprintf(painCave.errMsg,  "%s is an unrecognized keyword\n", keyword.c_str() );
	  painCave.isFatal = 0;
	  simError();        
      }

      return result;
  }
 
  OptionSectionParser::OptionSectionParser(ForceFieldOptions& options) : options_(options) {
    setSectionName("Options");        
  }
  
  void OptionSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    
    StringTokenizer tokenizer(line);
    
    if (tokenizer.countTokens() >= 2) {
      std::string optionName = tokenizer.nextToken();
      std::string optionValue = tokenizer.nextToken();
      
      options_.setData(optionName, optionValue);
      
    } else {
      sprintf(painCave.errMsg, "OptionSectionParser Error: Not enough tokens at line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();    
    }
    
  }

  void OptionSectionParser::validateSection() {
    options_.validateOptions();
  }

} //end namespace oopse  
