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

  EAMAtomTypesSectionParser::EAMAtomTypesSectionParser(ForceFieldOptions& options) :
    options_(options){
    setSectionName("EAMAtomTypes");
  }

  void EAMAtomTypesSectionParser::parseLine(ForceField& ff,
                                            const string& line, int lineNo){

    eus_ = options_.getMetallicEnergyUnitScaling();
    dus_ = options_.getDistanceUnitScaling();            

    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();

    if (tokenizer.countTokens() >= 2) {
      string atomTypeName = tokenizer.nextToken();
      string eamParameterType = tokenizer.nextToken();
      nTokens -= 2;
      AtomType* atomType = ff.getAtomType(atomTypeName);      
      if (atomType != NULL) {        

        EAMAdapter ea = EAMAdapter(atomType);
        toUpper(eamParameterType);
        
        if (eamParameterType == "FUNCFL") {
          string potentialParamFile = tokenizer.nextToken();
          parseEAMfuncflFile(ff, ea, potentialParamFile, atomType->getIdent());
        } else if(eamParameterType == "ZHOU") {
          if (nTokens < 20) {
            sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                    "Not enough tokens at line %d\n", lineNo);
            painCave.isFatal = 1;
            simError();
          } else {
            std::string lattice = tokenizer.nextToken();
            toUpper(lattice);

            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
            RealType latticeConstant;
            // default to FCC if we don't specify HCP or BCC:
            if (lattice == "HCP") 
              latticeConstant = re;
            else if (lattice == "BCC")
              latticeConstant = 2.0 * re / sqrt(3.0);
            else
              latticeConstant = 2.0 * re / sqrt(2.0);
            
            RealType fe         = tokenizer.nextTokenAsDouble();
            RealType rhoe       = tokenizer.nextTokenAsDouble();
            RealType alpha      = tokenizer.nextTokenAsDouble();
            RealType beta       = tokenizer.nextTokenAsDouble();
            RealType A          = eus_ * tokenizer.nextTokenAsDouble();
            RealType B          = eus_ * tokenizer.nextTokenAsDouble();
            RealType kappa      = tokenizer.nextTokenAsDouble();
            RealType lambda     = tokenizer.nextTokenAsDouble();
            std::vector<RealType> Fn;
            Fn.push_back(eus_ * tokenizer.nextTokenAsDouble());
            Fn.push_back(eus_ * tokenizer.nextTokenAsDouble());
            Fn.push_back(eus_ * tokenizer.nextTokenAsDouble());
            Fn.push_back(eus_ * tokenizer.nextTokenAsDouble());
            std::vector<RealType> F;
            F.push_back(eus_ * tokenizer.nextTokenAsDouble());
            F.push_back(eus_ * tokenizer.nextTokenAsDouble());
            F.push_back(eus_ * tokenizer.nextTokenAsDouble());
            F.push_back(eus_ * tokenizer.nextTokenAsDouble());
            RealType eta        = tokenizer.nextTokenAsDouble();
            RealType Fe         = eus_ * tokenizer.nextTokenAsDouble();
            ea.makeEAM(latticeConstant, lattice, re, fe, rhoe, alpha, beta, A,
                       B, kappa, lambda, Fn, F, eta, Fe, false);
          }
        
        } else {
          sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: %s "
                  "is not a recognized EAM type\n", eamParameterType.c_str()); 
          painCave.isFatal = 1; 
          simError();     
        }
        
      } else {
        sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                "Can not find AtomType [%s]\n", atomTypeName.c_str());
        painCave.isFatal = 1;
        simError();  
      }
      
    } else {
      sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
              "Not enough tokens at line %d\n", lineNo);
      painCave.isFatal = 1;
      simError();    
    }    
  }

  void EAMAtomTypesSectionParser::parseEAMfuncflFile(ForceField& ff, 
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
    int atomicNumber;
    RealType atomicMass;
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
        atomicNumber = tokenizer1.nextTokenAsInt();
        atomicMass = tokenizer1.nextTokenAsDouble();
	latticeConstant = tokenizer1.nextTokenAsDouble() * dus_;
	lattice = tokenizer1.nextToken();
      }else {
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                "Not enough tokens\n");
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
	dr = tokenizer2.nextTokenAsDouble() * dus_;
	rcut = tokenizer2.nextTokenAsDouble() * dus_;
      } else {
        
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                "Not enough tokens\n");
	painCave.isFatal = 1;
	simError();
        
      }
    } 
    
    parseEAMArray(*ppfStream, F,   nrho);    
    parseEAMArray(*ppfStream, Z,   nr);
    parseEAMArray(*ppfStream, rho, nr);

    // Convert to kcal/mol using energy unit scaling in force field:
    std::transform(F.begin(), F.end(), F.begin(),
                   std::bind1st(std::multiplies<RealType>(), eus_));
    
    ea.makeEAM(latticeConstant, lattice, nrho, drho, nr, dr, rcut, Z, rho,
               F, true);
    
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
	sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                "Not enough tokens\n");
	painCave.isFatal = 1;
	simError();  
      }
      ++lineCount;
    }
    
    if (lineCount < nlinesToRead) {
      sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
              "Not enough lines to read\n");
      painCave.isFatal = 1;
      simError();          
    }
    
  }
  
  
} //end namespace OpenMD

