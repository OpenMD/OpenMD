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
 */
 
#include "io/TorsionTypesSectionParser.hpp"
#include "types/TorsionType.hpp"
#include "types/CubicTorsionType.hpp"
#include "types/QuarticTorsionType.hpp"
#include "types/PolynomialTorsionType.hpp"
#include "types/CharmmTorsionType.hpp"
#include "UseTheForce/ForceField.hpp"

namespace oopse {

TorsionTypesSectionParser::TorsionTypesSectionParser() {
    setSectionName("TorsionTypes");


    stringToEnumMap_["GhostTorsion"] == ttGhostTorsion;
    stringToEnumMap_["Cubic"] = ttCubic;
    stringToEnumMap_["Quartic"] = ttQuartic;
    stringToEnumMap_["Polynomial"] = ttPolynomial;
    stringToEnumMap_["Charmm"] =  ttCharmm;

}

void TorsionTypesSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    TorsionType* torsionType = NULL;

    int nTokens = tokenizer.countTokens();

    if (nTokens < 5) {

        return;
    }
    
    std::string at1 = tokenizer.nextToken();
    std::string at2 = tokenizer.nextToken();
    std::string at3 = tokenizer.nextToken();
    std::string at4 = tokenizer.nextToken();
    TorsionTypeEnum tt = getTorsionTypeEnum(tokenizer.nextToken());

    nTokens -= 5;

    switch(tt) {

        case TorsionTypesSectionParser::ttGhostTorsion:
            if (nTokens < 4) {

            } else {

                double k3 = tokenizer.nextTokenAsDouble();
                double k2 = tokenizer.nextTokenAsDouble();
                double k1 = tokenizer.nextTokenAsDouble();
                double k0 = tokenizer.nextTokenAsDouble();
                
                torsionType = new CubicTorsionType(k3, k2, k1, k0);
            }
            break;
            
        case TorsionTypesSectionParser::ttCubic :
            if (nTokens < 4) {

            } else {

                double k3 = tokenizer.nextTokenAsDouble();
                double k2 = tokenizer.nextTokenAsDouble();
                double k1 = tokenizer.nextTokenAsDouble();
                double k0 = tokenizer.nextTokenAsDouble();
                
                torsionType = new CubicTorsionType(k3, k2, k1, k0);
            }
            break;
            
        case TorsionTypesSectionParser::ttQuartic:
            if (nTokens < 5) {

            } else {

                double k4 = tokenizer.nextTokenAsDouble();
                double k3 = tokenizer.nextTokenAsDouble();
                double k2 = tokenizer.nextTokenAsDouble();
                double k1 = tokenizer.nextTokenAsDouble();
                double k0 = tokenizer.nextTokenAsDouble();
                
                torsionType = new QuarticTorsionType( k4, k3, k2, k1, k0);
            }
            break;

        
        case TorsionTypesSectionParser::ttPolynomial:
            if (nTokens < 2 || nTokens % 2 != 0) {

            } else {
                int nPairs = nTokens / 2;
                int power;
                double coefficient;
                PolynomialTorsionType* ptt = new PolynomialTorsionType();
                
                for (int i = 0; i < nPairs; ++i) {
                    power = tokenizer.nextTokenAsInt();
                    coefficient = tokenizer.nextTokenAsDouble();
                    ptt->setCoefficient(power, coefficient);
                }
            }
            
            break;
             
        case TorsionTypesSectionParser::ttCharmm:
            
            if (nTokens < 3 || nTokens % 3 != 0) {

            } else {
                int nSets = nTokens / 3;
  
                CharmmTorsionType* ctt = new CharmmTorsionType();
                
                for (int i = 0; i < nSets; ++i) {
                    double kchi = tokenizer.nextTokenAsDouble();
                    int n = tokenizer.nextTokenAsInt();
                    double delta = tokenizer.nextTokenAsDouble();
    
                    ctt->setCharmmTorsionParameter(kchi, n, delta);
                }
            }

            break;
            
        case TorsionTypesSectionParser::ttUnknown :
        default:

            break;
            
    }

    if (torsionType != NULL) {
        ff.addTorsionType(at1, at2, at3, at4, torsionType);
    }

}

TorsionTypesSectionParser::TorsionTypeEnum TorsionTypesSectionParser::getTorsionTypeEnum(const std::string& str) {
    std::map<std::string, TorsionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? TorsionTypesSectionParser::ttUnknown : i->second;
}

} //end namespace oopse



