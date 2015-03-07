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
 
#include "io/EAMAtomTypesSectionParser.hpp"
#include "types/AtomType.hpp"
#include "types/EAMAdapter.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  EAMAtomTypesSectionParser::EAMAtomTypesSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("EAMAtomTypes");
  }

  void EAMAtomTypesSectionParser::parseLine(ForceField& ff,
                                            const string& line, int lineNo){

    StringTokenizer tokenizer(line);

    if (tokenizer.countTokens() >= 2) {
      string atomTypeName = tokenizer.nextToken();
      string potentialParamFile = tokenizer.nextToken();

      AtomType* atomType = ff.getAtomType(atomTypeName);      
      if (atomType != NULL) {        

        EAMAdapter ea = EAMAdapter(atomType);
        parseEAMParamFile(ff, ea, potentialParamFile, atomType->getIdent());                                                    
      } else {
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: Can not find AtomType [%s]\n",
                atomTypeName.c_str());
	painCave.isFatal = 1;
	simError();  
      }
        
    } else {
      sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();    
    }
            
  }

  void EAMAtomTypesSectionParser::parseEAMParamFile(ForceField& ff, 
                                                    EAMAdapter ea, 
						    const string& potentialParamFile, 
                                                    int ident) {

    ifstrstream* ppfStream = ff.openForceFieldFile(potentialParamFile);
    const int bufferSize = 65535;
    char buffer[bufferSize];

    // skip first line
    ppfStream->getline(buffer, bufferSize);

    // The second line contains atomic number, atomic mass, a lattice
    // constant and lattice type
    RealType latticeConstant(0.0); 
    string lattice;

    // The third line is nrho, drho, nr, dr and rcut 
    int nrho(0);
    RealType drho(0.0);
    int nr(0);
    RealType dr(0.0);
    RealType rcut(0.0);
    vector<RealType> F;
    vector<RealType> Z;
    vector<RealType> rho;
      
    if (ppfStream->getline(buffer, bufferSize)) {       
      StringTokenizer tokenizer1(buffer);
        
      if (tokenizer1.countTokens() >= 4) {
	int junk = tokenizer1.nextTokenAsInt();
	RealType mass = tokenizer1.nextTokenAsDouble();
	latticeConstant = tokenizer1.nextTokenAsDouble();
	lattice = tokenizer1.nextToken();
      }else {
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: Not enough tokens\n");
	painCave.isFatal = 1;
	simError();  
      }
    }
  
    
    if (ppfStream->getline(buffer, bufferSize)) {
      StringTokenizer tokenizer2(buffer);
      
      if (tokenizer2.countTokens() >= 5){
	nrho = tokenizer2.nextTokenAsInt();
	drho = tokenizer2.nextTokenAsDouble();
        nr = tokenizer2.nextTokenAsInt();
	dr = tokenizer2.nextTokenAsDouble();
	rcut = tokenizer2.nextTokenAsDouble();
      }else {
        
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: Not enough tokens\n");
	painCave.isFatal = 1;
	simError();            
        
      }
    } 
    
    parseEAMArray(*ppfStream, F,   nrho);    
    parseEAMArray(*ppfStream, Z,   nr);
    parseEAMArray(*ppfStream, rho, nr);
    
    ea.makeEAM(latticeConstant, nrho, drho, nr, dr, rcut, Z, rho, F);

    delete ppfStream;
  }

  void EAMAtomTypesSectionParser::parseEAMArray(istream& input, 
						vector<RealType>& array, 
                                                int num) {
    
    const int dataPerLine = 5;
    if (num % dataPerLine != 0) {

    }

    int nlinesToRead = num / dataPerLine;
    
    const int bufferSize = 65535;
    char buffer[bufferSize];
    int lineCount = 0;

    while(lineCount < nlinesToRead && input.getline(buffer, bufferSize) ){
      
      StringTokenizer tokenizer(buffer);
      if (tokenizer.countTokens() >= dataPerLine) {
	for (int i = 0; i < dataPerLine; ++i) {
	  array.push_back(tokenizer.nextTokenAsDouble());
	}
      } else {
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: Not enough tokens\n");
	painCave.isFatal = 1;
	simError();  
      }
      ++lineCount;
    }

    if (lineCount < nlinesToRead) {
      sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: Not enough lines to read\n");
      painCave.isFatal = 1;
      simError();          
    }
    
  }


} //end namespace OpenMD

