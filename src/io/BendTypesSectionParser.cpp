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

#include "utils/NumericConstant.hpp" 
#include "io/BendTypesSectionParser.hpp"
#include "types/HarmonicBendType.hpp"
#include "types/UreyBradleyBendType.hpp"
#include "types/CubicBendType.hpp"
#include "types/QuarticBendType.hpp"
#include "types/PolynomialBendType.hpp"
#include "types/CosineBendType.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"
namespace OpenMD {

  BendTypesSectionParser::BendTypesSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("BendTypes");

    stringToEnumMap_["Harmonic"] =  btHarmonic;       
    stringToEnumMap_["GhostBend"] =  btGhostBend;                
    stringToEnumMap_["UreyBradley"] =  btUreyBradley;                    
    stringToEnumMap_["Cubic"] = btCubic;
    stringToEnumMap_["Quartic"] = btQuartic;
    stringToEnumMap_["Polynomial"] = btPolynomial;    
    stringToEnumMap_["Cosine"] = btCosine;    
  }

  void BendTypesSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    BendType* bendType = NULL;

    int nTokens = tokenizer.countTokens();

    if (nTokens < 5) {
      sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
      return;
    }
    
    std::string at1 = tokenizer.nextToken();
    std::string at2 = tokenizer.nextToken();
    std::string at3 = tokenizer.nextToken();
    BendTypeEnum bt = getBendTypeEnum(tokenizer.nextToken());
    RealType theta0 = tokenizer.nextTokenAsDouble() / 180.0 * NumericConstant::PI; //convert to rad
    nTokens -= 5;

    //switch is a maintain nightmare
    switch(bt) {
            
    case btHarmonic :
            
      if (nTokens < 1) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {

	RealType ktheta = tokenizer.nextTokenAsDouble();
	bendType = new HarmonicBendType(theta0, ktheta);
      }
      break;
    case btGhostBend :
      if (nTokens < 1) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {
	RealType ktheta = tokenizer.nextTokenAsDouble();
	bendType = new HarmonicBendType(theta0, ktheta);                
      }
      break;            

    case btUreyBradley :
      if (nTokens < 3) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {
	RealType ktheta = tokenizer.nextTokenAsDouble();
	RealType s0 =  tokenizer.nextTokenAsDouble();
	RealType kub = tokenizer.nextTokenAsDouble();
	bendType = new UreyBradleyBendType(theta0, ktheta, s0, kub);                
      }
      break; 
            
    case btCubic :
      if (nTokens < 4) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {

	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
	bendType = new CubicBendType(theta0, k3, k2, k1, k0);
      }
      break;
            
    case btQuartic :
      if (nTokens < 5) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {

	RealType k4 = tokenizer.nextTokenAsDouble();
	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
	bendType = new QuarticBendType(theta0, k4, k3, k2, k1, k0);
      }
      break;

    case btPolynomial :
      if (nTokens < 2 || nTokens % 2 != 0) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {
	int nPairs = nTokens / 2;
	int power;
	RealType coefficient;
	PolynomialBendType* pbt = new PolynomialBendType(theta0);
                
	for (int i = 0; i < nPairs; ++i) {
	  power = tokenizer.nextTokenAsInt();
	  coefficient = tokenizer.nextTokenAsDouble();
	  pbt->setCoefficient(power, coefficient);
	}
      }
            
      break;

    case btCosine :
            
      if (nTokens < 1) {
	sprintf(painCave.errMsg, "BendTypesSectionParser Error: Not enough tokens at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      } else {

	RealType ktheta = tokenizer.nextTokenAsDouble();
	bendType = new CosineBendType(theta0, ktheta);
      }
      break;

    case btUnknown :
    default:
      sprintf(painCave.errMsg, "BendTypesSectionParser Error: Unknown Bond Type at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
      break;
            
    }

    if (bendType != NULL) {
      ff.addBendType(at1, at2, at3, bendType);
    }
  }

  BendTypesSectionParser::BendTypeEnum BendTypesSectionParser::getBendTypeEnum(const std::string& str) {
    std::map<std::string, BendTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? btUnknown : i->second;
  }

} //end namespace OpenMD


