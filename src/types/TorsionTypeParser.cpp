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

#include "utils/Constants.hpp" 
#include "types/TorsionTypeParser.hpp"
#include "types/CubicTorsionType.hpp"
#include "types/QuarticTorsionType.hpp"
#include "types/PolynomialTorsionType.hpp"
#include "types/CharmmTorsionType.hpp"
#include "types/OplsTorsionType.hpp"
#include "types/TrappeTorsionType.hpp"
#include "types/HarmonicTorsionType.hpp"
#include "utils/OpenMDException.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  TorsionTypeParser::TorsionTypeParser() : trans180_(true) {
    stringToEnumMap_["GhostTorsion"] = ttGhostTorsion;
    stringToEnumMap_["Cubic"] = ttCubic;
    stringToEnumMap_["Quartic"] = ttQuartic;
    stringToEnumMap_["Polynomial"] = ttPolynomial;
    stringToEnumMap_["Charmm"] =  ttCharmm;
    stringToEnumMap_["Opls"] =  ttOpls;
    stringToEnumMap_["Trappe"] =  ttTrappe;
    stringToEnumMap_["Harmonic"] =  ttHarmonic;
  }

  TorsionType* TorsionTypeParser::parseTypeAndPars(const std::string& type,
                                                   std::vector<RealType> pars) {
    
    std::string line(type);

    std::vector<RealType>::iterator it;
    for (it = pars.begin(); it != pars.end(); ++it) {
      line.append("\t");
      line.append( OpenMD::to_string(*it) );      
    }
    
    return parseLine( line );    
  }

  TorsionType* TorsionTypeParser::parseLine(const std::string& line) {
    
    StringTokenizer tokenizer(line);
    TorsionType* torsionType = NULL;
    int nTokens = tokenizer.countTokens();

    if (nTokens < 1) {
      throw OpenMDException("TorsionTypeParser: Not enough tokens");
    }
    
    TorsionTypeEnum tt = getTorsionTypeEnum(tokenizer.nextToken());

    nTokens -= 1;

    switch(tt) {

    case ttGhostTorsion:
      if (nTokens < 4) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {
        
	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
        
        if (trans180_)
          torsionType = new CubicTorsionType( k3, k2,  k1, k0);
        else 
          torsionType = new CubicTorsionType(-k3, k2, -k1, k0);
      }
      break;
            
    case ttCubic :
      if (nTokens < 4) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {

	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
        if (trans180_)
          torsionType = new CubicTorsionType( k3, k2,  k1, k0);
        else 
          torsionType = new CubicTorsionType(-k3, k2, -k1, k0);
      }
      break;
            
    case ttQuartic:
      if (nTokens < 5) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {

	RealType k4 = tokenizer.nextTokenAsDouble();
	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
        if (trans180_)
          torsionType = new QuarticTorsionType( k4,  k3, k2,  k1, k0);
        else 
          torsionType = new QuarticTorsionType( k4, -k3, k2, -k1, k0);
      }
      break;
        
    case ttPolynomial:
      if (nTokens < 2 || nTokens % 2 != 0) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {
	int nPairs = nTokens / 2;
	int power;
	RealType coefficient;
	torsionType = new PolynomialTorsionType();

        PolynomialTorsionType* ptt = dynamic_cast<PolynomialTorsionType*>(torsionType); 
                
	for (int i = 0; i < nPairs; ++i) {
	  power = tokenizer.nextTokenAsInt();
          bool isOdd = power % 2 == 0 ? false : true;

	  coefficient = tokenizer.nextTokenAsDouble();

          if (!trans180_ && isOdd)
            coefficient = -coefficient;

          ptt->setCoefficient(power, coefficient);
          
	}
      }
            
      break;
             
    case ttCharmm:
            
      if (nTokens < 3 || nTokens % 3 != 0) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {
	int nSets = nTokens / 3;
  
        std::vector<CharmmTorsionParameter> parameters;             
	for (int i = 0; i < nSets; ++i) {
          CharmmTorsionParameter currParam;
	  currParam.kchi = tokenizer.nextTokenAsDouble();
	  currParam.n = tokenizer.nextTokenAsInt();
	  currParam.delta = tokenizer.nextTokenAsDouble() / 180.0 * Constants::PI; //convert to rad

          bool isOdd = currParam.n % 2 == 0 ? false : true;
          if (!trans180_) {           
            currParam.delta = Constants::PI - currParam.delta;
            if (isOdd) currParam.kchi = -currParam.kchi;
          }

          parameters.push_back(currParam);
	}

	torsionType = new CharmmTorsionType(parameters);

      }

      break;

    case ttOpls:
            
      if (nTokens < 3) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {
	RealType v1 = tokenizer.nextTokenAsDouble();
	RealType v2 = tokenizer.nextTokenAsDouble();
	RealType v3 = tokenizer.nextTokenAsDouble();
        
	torsionType = new OplsTorsionType(v1, v2, v3, trans180_);
      }

      break;
            

    case ttTrappe:
            
      if (nTokens < 4) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {

	RealType c0 = tokenizer.nextTokenAsDouble();
	RealType c1 = tokenizer.nextTokenAsDouble();
	RealType c2 = tokenizer.nextTokenAsDouble();
	RealType c3 = tokenizer.nextTokenAsDouble();
	torsionType = new TrappeTorsionType(c0, c1, c2, c3, trans180_);
      }

      break;

    case ttHarmonic:
            
      if (nTokens < 2) {
        throw OpenMDException("TorsionTypeParser: Not enough tokens");
      } else {
        // Most torsions don't have specific angle information since
        // they are cosine polynomials.  This one is different,
        // however.  To match our other force field files
        // (particularly for bends): 
        //
        // phi0 should be read in degrees
        // d0 should be read in kcal / mol / degrees^2

        RealType degreesPerRadian = 180.0 / Constants::PI;

        // convert to radians
	RealType phi0 = tokenizer.nextTokenAsDouble() / degreesPerRadian;
        
        // convert to kcal / mol / radians^2
	RealType d0 = tokenizer.nextTokenAsDouble() * pow(degreesPerRadian,2);

        if (!trans180_)
          phi0 = Constants::PI - phi0;

	torsionType = new HarmonicTorsionType(d0, phi0);
      }

      break;
            
    case ttUnknown :
    default:
      throw OpenMDException("TorsionTypeParser: Unknown Torsion Type");
    }
    return torsionType;
  }

  TorsionTypeParser::TorsionTypeEnum TorsionTypeParser::getTorsionTypeEnum(const std::string& str) {
    std::map<std::string, TorsionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);
    
    return i == stringToEnumMap_.end() ? ttUnknown : i->second;
  }
  
}



