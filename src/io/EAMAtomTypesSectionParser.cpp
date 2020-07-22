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

#include "io/EAMAtomTypesSectionParser.hpp"
#include "types/AtomType.hpp"
#include "types/EAMAdapter.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

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
      std::string atomTypeName = tokenizer.nextToken();
      std::string eamParameterType = tokenizer.nextToken();
      nTokens -= 2;
      AtomType* atomType = ff.getAtomType(atomTypeName);
      if (atomType != NULL) {

        EAMAdapter ea = EAMAdapter(atomType);
        toUpper(eamParameterType);

        if (eamParameterType == "FUNCFL") {
          std::string funcflFile = tokenizer.nextToken();
          parseFuncflFile(ff, ea, funcflFile, atomType->getIdent());
        } else if(eamParameterType == "ZHOU" ||
                  eamParameterType == "ZHOU2001") {
          if (nTokens < 20) {
            sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                    "Not enough tokens at line %d\n", lineNo);
            painCave.isFatal = 1;
            simError();
          } else {

            std::string lattice = tokenizer.nextToken();
            toUpper(lattice);
            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
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

            ea.makeZhou2001(lattice, re, fe, rhoe, alpha, beta, A,
                            B, kappa, lambda, Fn, F, eta, Fe);

          }
        } else if(eamParameterType == "ZHOU2004") {
          if (nTokens < 23) {
            sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                    "Not enough tokens at line %d\n", lineNo);
            painCave.isFatal = 1;
            simError();
          } else {

            std::string lattice = tokenizer.nextToken();
            toUpper(lattice);
            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
            RealType fe         = tokenizer.nextTokenAsDouble();
            RealType rhoe       = tokenizer.nextTokenAsDouble();
            RealType rhos       = tokenizer.nextTokenAsDouble();
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
            RealType rhol       = tokenizer.nextTokenAsDouble();
            RealType rhoh       = tokenizer.nextTokenAsDouble();

            ea.makeZhou2004(lattice, re, fe, rhoe, rhos, alpha, beta,
                            A, B, kappa, lambda, Fn, F, eta, Fe, rhol, rhoh);
          }

        } else if(eamParameterType == "ZHOU2005") {
          if (nTokens < 22) {
            sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                    "Not enough tokens at line %d\n", lineNo);
            painCave.isFatal = 1;
            simError();
          } else {

            std::string lattice = tokenizer.nextToken();
            toUpper(lattice);
            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
            RealType fe         = tokenizer.nextTokenAsDouble();
            RealType rhoe       = tokenizer.nextTokenAsDouble();
            RealType rhos       = tokenizer.nextTokenAsDouble();
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
            RealType F3minus    = eus_ * tokenizer.nextTokenAsDouble();
            RealType F3plus     = eus_ * tokenizer.nextTokenAsDouble();
            RealType eta        = tokenizer.nextTokenAsDouble();
            RealType Fe         = eus_ * tokenizer.nextTokenAsDouble();

            ea.makeZhou2005(lattice, re, fe, rhoe, rhos, alpha, beta,
                            A, B, kappa, lambda, Fn, F, F3minus, F3plus,
                            eta, Fe);
          }
        } else if(eamParameterType == "ZHOU2005OXYGEN") {
          if (nTokens < 36) {
            sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                    "Not enough tokens at line %d\n", lineNo);
            painCave.isFatal = 1;
            simError();
          } else {
            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
            RealType fe         = tokenizer.nextTokenAsDouble();
            RealType alpha      = tokenizer.nextTokenAsDouble();
            RealType beta       = tokenizer.nextTokenAsDouble();
            RealType A          = eus_ * tokenizer.nextTokenAsDouble();
            RealType B          = eus_ * tokenizer.nextTokenAsDouble();
            RealType kappa      = tokenizer.nextTokenAsDouble();
            RealType lambda     = tokenizer.nextTokenAsDouble();
            RealType gamma      = tokenizer.nextTokenAsDouble();
            RealType nu         = tokenizer.nextTokenAsDouble();

            std::vector<RealType> OrhoLimits;
            OrhoLimits.push_back(tokenizer.nextTokenAsDouble());
            OrhoLimits.push_back(tokenizer.nextTokenAsDouble());
            OrhoLimits.push_back(tokenizer.nextTokenAsDouble());
            OrhoLimits.push_back(tokenizer.nextTokenAsDouble());
            OrhoLimits.push_back(tokenizer.nextTokenAsDouble());

            std::vector<RealType> OrhoE;
            OrhoE.push_back(tokenizer.nextTokenAsDouble());
            OrhoE.push_back(tokenizer.nextTokenAsDouble());
            OrhoE.push_back(tokenizer.nextTokenAsDouble());
            OrhoE.push_back(tokenizer.nextTokenAsDouble());
            OrhoE.push_back(tokenizer.nextTokenAsDouble());

            std::vector<std::vector<RealType> > OF;
            OF.resize(5);
            OF[0].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[0].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[0].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[0].push_back(eus_ * tokenizer.nextTokenAsDouble());

            OF[1].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[1].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[1].push_back(eus_ * tokenizer.nextTokenAsDouble());

            OF[2].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[2].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[2].push_back(eus_ * tokenizer.nextTokenAsDouble());

            OF[3].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[3].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[3].push_back(eus_ * tokenizer.nextTokenAsDouble());

            OF[4].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[4].push_back(eus_ * tokenizer.nextTokenAsDouble());
            OF[4].push_back(eus_ * tokenizer.nextTokenAsDouble());

            ea.makeZhou2005Oxygen(re, fe, alpha, beta, A, B, kappa, lambda,
                                  gamma, nu, OrhoLimits, OrhoE, OF);
          }
        } else if(eamParameterType == "ZHOUROSE") {
          if (nTokens < 10) {
              sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                      "Not enough tokens at line %d\n", lineNo);
              painCave.isFatal = 1;
            simError();
          } else {

            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
            RealType fe         = tokenizer.nextTokenAsDouble();
            RealType rhoe       = tokenizer.nextTokenAsDouble();
            RealType alpha      = tokenizer.nextTokenAsDouble();
            RealType beta       = tokenizer.nextTokenAsDouble();
            RealType A          = eus_ * tokenizer.nextTokenAsDouble();
            RealType B          = eus_ * tokenizer.nextTokenAsDouble();
            RealType kappa      = tokenizer.nextTokenAsDouble();
            RealType lambda     = tokenizer.nextTokenAsDouble();
            RealType F0         = eus_ * tokenizer.nextTokenAsDouble();

            ea.makeZhouRose(re, fe, rhoe, alpha, beta, A, B, kappa, lambda, F0);
            
          }          
        }else if(eamParameterType == "OXYFUNCFL") {
          if (nTokens < 9) {
              sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                      "Not enough tokens at line %d\n", lineNo);
              painCave.isFatal = 1;
            simError();
          } else {
            std::string funcflFile = tokenizer.nextToken();

            RealType re         = dus_ * tokenizer.nextTokenAsDouble();
            RealType fe         = tokenizer.nextTokenAsDouble();
            RealType alpha      = tokenizer.nextTokenAsDouble();
            RealType beta       = tokenizer.nextTokenAsDouble();
            RealType A          = eus_ * tokenizer.nextTokenAsDouble();
            RealType B          = eus_ * tokenizer.nextTokenAsDouble();
            RealType kappa      = tokenizer.nextTokenAsDouble();
            RealType lambda     = tokenizer.nextTokenAsDouble();
            
            ifstrstream* ppfStream = ff.openForceFieldFile(funcflFile);
            const int bufferSize = 65535;
            char buffer[bufferSize];

          // skip first line
            ppfStream->getline(buffer, bufferSize);

            std::vector<RealType> F;
    
            int nrho = 0;
            RealType drho = 0.0;
            if (ppfStream->getline(buffer, bufferSize)) {
              StringTokenizer tokenizer1(buffer);

              if (tokenizer1.countTokens() == 2) {
                nrho = tokenizer1.nextTokenAsInt();
                drho = tokenizer1.nextTokenAsDouble();
              }else {
                  sprintf(painCave.errMsg, "EAMAtomTypesSectionParser Error: "
                  "Not enough tokens\n");
                  painCave.isFatal = 1;
                  simError();
                  }
            }
    
            parseEAMArray(*ppfStream, F, nrho);
            delete ppfStream;

            // Convert to eV using energy unit scaling in force field:
            std::transform(F.begin(), F.end(), F.begin(),
                   std::bind(std::multiplies<RealType>(), eus_, placeholders::_1));

            ea.makeOxygenFuncfl(re, fe, alpha, beta, A, B, kappa, lambda, drho, nrho, F);
                              
    
            
    }         
    }else {
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

  void EAMAtomTypesSectionParser::parseFuncflFile(ForceField& ff,
                                                  EAMAdapter ea,
                                                  const string& funcflFile,
                                                  int ident) {

    ifstrstream* ppfStream = ff.openForceFieldFile(funcflFile);
    const int bufferSize = 65535;
    char buffer[bufferSize];

    // skip first line
    ppfStream->getline(buffer, bufferSize);

    // The second line contains atomic number, atomic mass, a lattice
    // constant and lattice type
    int atomicNumber;
    RealType atomicMass;
    RealType latticeConstant(0.0);
    std::string lattice;

    // The third line is nrho, drho, nr, dr and rcut
    int nrho(0);
    RealType drho(0.0);
    int nr(0);
    RealType dr(0.0);
    RealType rcut(0.0);
    std::vector<RealType> F;
    std::vector<RealType> Z;
    std::vector<RealType> rho;

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

    FIX_UNUSED(atomicNumber);
    FIX_UNUSED(atomicMass);

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
                   std::bind(std::multiplies<RealType>(), eus_, placeholders::_1));

    ea.makeFuncfl(latticeConstant, lattice, nrho, drho, nr, dr, rcut, Z, rho,
                  F);

    delete ppfStream;
  }

  void EAMAtomTypesSectionParser::parseEAMArray(istream& input,
						std::vector<RealType>& array,
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
