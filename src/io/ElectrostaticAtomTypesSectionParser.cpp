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
 
#include "io/ElectrostaticAtomTypesSectionParser.hpp"
#include "UseTheForce/ForceField.hpp"
#include "utils/NumericConstant.hpp"
namespace oopse {

ElectrostaticAtomTypesSectionParser::ElectrostaticAtomTypesSectionParser() {
    setSectionName("ElectrostaticAtomTypes");
}

void ElectrostaticAtomTypesSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();    

    //in AtomTypeSection, a line at least contains 2 tokens
    //atomTypeName and biggest rank
    //for the time being, we only support up to quadrupole
    // "name" must match the name in the AtomTypes section
    // charge is given in units of electrons (1.61 x 10^-19 C)
    // Directionality for dipoles and quadrupoles must be given because the body-fixed
    // reference frame for directional atoms is determined by the *mass* distribution and
    // not by the charge distribution.  
    // Dipoles are given in units of Debye  
    // Quadrupoles are given in units of 
    // name 0 charge
    // name 1 charge |u| [theta phi psi]
    // name 2 charge |u| Qxx Qyy Qzz [theta phi psi]
    
    if (nTokens < 2)  {
        std::cerr << "ElectrostaticAtomTypesSectionParser Error: Not enought Tokens at line " << lineNo << std::endl;                  
    } else {

        std::string atomTypeName = tokenizer.nextToken();    
        int biggestRank = tokenizer.nextTokenAsInt();
        nTokens -=  2;
        
        AtomType* atomType = ff.getAtomType(atomTypeName);
        DirectionalAtomType* dAtomType;
        if (atomType != NULL) {
                       
            switch (biggestRank) {
                case 0 :
                    parseCharge(tokenizer, atomType);
                    break;

               case 1 :

                    dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);            
                    if (dAtomType == NULL) {
                        std::cerr << "ElectrostaticAtomTypesSectionParser Warning:" << std::endl;
                    }

                    parseCharge(tokenizer, dAtomType);
                    parseDipole(tokenizer, dAtomType);
                    parseElectroBodyFrame(tokenizer, dAtomType);
                    break;

               case 2:

                    dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);            
                    if (dAtomType == NULL) {
                        std::cerr << "ElectrostaticAtomTypesSectionParser Warning:" << std::endl;
                    }

                    parseCharge(tokenizer, dAtomType);
                    parseDipole(tokenizer, dAtomType);
                    parseQuadruple(tokenizer, dAtomType);
                    parseElectroBodyFrame(tokenizer, dAtomType);
                    break;
                    
               default :
                    break;

            }
            
        } else {
            std::cerr << "ElectrostaticAtomTypesSectionParser Error: Can not find matched AtomType " << atomTypeName
                          << "at line " << lineNo << std::endl;
        }
                       
    }    


}


void ElectrostaticAtomTypesSectionParser::parseCharge(StringTokenizer& tokenizer,
    AtomType* atomType) {

    double charge = tokenizer.nextTokenAsDouble();

    if (fabs(charge) > NumericConstant::epsilon) {
        atomType->addProperty(new DoubleGenericData("Charge", charge));
        atomType->setCharge();
    }
}
void ElectrostaticAtomTypesSectionParser::parseDipole(StringTokenizer& tokenizer, 
    DirectionalAtomType* dAtomType) {

    double dipole = tokenizer.nextTokenAsDouble();

    if (fabs(dipole) > NumericConstant::epsilon) {

        dAtomType->addProperty(new DoubleGenericData("Dipole", dipole));
        dAtomType->setDipole();
    }
}

void ElectrostaticAtomTypesSectionParser::parseQuadruple(StringTokenizer& tokenizer,
    DirectionalAtomType* dAtomType) {

    Vector3d Q;
    Q[0] = tokenizer.nextTokenAsDouble();
    Q[1] = tokenizer.nextTokenAsDouble();
    Q[2] = tokenizer.nextTokenAsDouble();

    if (fabs(Q[0]) > NumericConstant::epsilon && fabs(Q[1]) > NumericConstant::epsilon 
        && fabs(Q[2]) > NumericConstant::epsilon) {    
        
        dAtomType->addProperty(new Vector3dGenericData("Quadrupole", Q));
        dAtomType->setQuadrupole();
    }
}
void ElectrostaticAtomTypesSectionParser::parseElectroBodyFrame(StringTokenizer& tokenizer,
    DirectionalAtomType* dAtomType) {

    double phi;
    double theta;
    double psi;

    if (tokenizer.countTokens() >=3 ) {
        phi = tokenizer.nextTokenAsDouble()/180.0;
        theta = tokenizer.nextTokenAsDouble()/180.0;
        psi = tokenizer.nextTokenAsDouble()/180.0;
    } else {
        phi = 0.0;
        theta = 0.0;
        psi = 0.0;    
    }
    
    RotMat3x3d electroBodyFrame(phi, theta, psi);
    dAtomType->setElectroBodyFrame(electroBodyFrame);
        
}

} //end namespace oopse



