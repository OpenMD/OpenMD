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
 
#include "io/EAMAtomTypesSectionParser.hpp"
#include "types/AtomType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "utils/simError.h"
namespace oopse {

EAMAtomTypesSectionParser::EAMAtomTypesSectionParser() {
    setSectionName("EAMAtomTypes");
}

void EAMAtomTypesSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){

    StringTokenizer tokenizer(line);

    if (tokenizer.countTokens() >= 2) {
        std::string atomTypeName = tokenizer.nextToken();
        std::string potentialParamFile = tokenizer.nextToken();

        AtomType* atomType = ff.getAtomType(atomTypeName);
        if (atomType != NULL) {
            atomType->setEAM();                            
            parseEAMParamFile(ff, atomType, potentialParamFile, atomType->getIdent());                                                    
        } else {

        }
        
    } else {
      
    }
            

}

void EAMAtomTypesSectionParser::parseEAMParamFile(ForceField& ff, AtomType* atomType, 
    const std::string& potentialParamFile, int ident) {

    ifstrstream* ppfStream = ff.openForceFieldFile(potentialParamFile);
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string line;

    //skip first line
    ppfStream->getline(buffer, bufferSize);


    //The Second line contains atomic number, atomic mass, a lattice constant and lattic type
    int junk;
    double mass;
    double latticeConstant; 
    std::string lattice;
    if (ppfStream->getline(buffer, bufferSize)) {       
        StringTokenizer tokenizer1(buffer);
        
        if (tokenizer1.countTokens() >= 4) {
            junk = tokenizer1.nextTokenAsInt();
            mass = tokenizer1.nextTokenAsDouble();
            latticeConstant = tokenizer1.nextTokenAsDouble();
            lattice = tokenizer1.nextToken();
        }else {
            std::cerr << "Not enought tokens" << std::endl;
        }
    } else {

    }
    
    // The third line is nrho, drho, nr, dr and rcut
    EAMParam eamParam;
    eamParam.latticeConstant = latticeConstant;
    
    if (ppfStream->getline(buffer, bufferSize)) {
        StringTokenizer tokenizer2(buffer);

        if (tokenizer2.countTokens() >= 5){
            eamParam.nrho = tokenizer2.nextTokenAsInt();
            eamParam.drho = tokenizer2.nextTokenAsDouble();
            eamParam.nr = tokenizer2.nextTokenAsInt();
            eamParam.dr = tokenizer2.nextTokenAsDouble();
            eamParam.rcut = tokenizer2.nextTokenAsDouble();
        }else {
            std::cerr << "Not enought tokens" << std::endl;
        }
    } else {

    }

    parseEAMArray(*ppfStream, eamParam.Frhovals, eamParam.nrho);    
    parseEAMArray(*ppfStream, eamParam.rvals, eamParam.nr);
    parseEAMArray(*ppfStream, eamParam.rhovals, eamParam.nr);
    
    atomType->addProperty(new EAMParamGenericData("EAM", eamParam));
}

void EAMAtomTypesSectionParser::parseEAMArray(std::istream& input, 
    std::vector<double>& array, int num) {
    
    const int dataPerLine = 5;
    if (num % dataPerLine != 0) {

    }

    int nlinesToRead = num / dataPerLine;
    
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string line;
    int readLines = num/dataPerLine;
    int lineCount = 0;

    while(lineCount < nlinesToRead && input.getline(buffer, bufferSize) ){

        StringTokenizer tokenizer(buffer);
        if (tokenizer.countTokens() >= dataPerLine) {
            for (int i = 0; i < dataPerLine; ++i) {
                array.push_back(tokenizer.nextTokenAsDouble());
            }
        } else {

        }
        ++lineCount;
    }

    if (lineCount < nlinesToRead) {
        
    }
    
}


} //end namespace oopse

