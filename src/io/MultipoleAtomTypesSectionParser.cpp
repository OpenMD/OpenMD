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
 
#include "io/MultipoleAtomTypesSectionParser.hpp"
#include "brains/ForceField.hpp"
#include "utils/NumericConstant.hpp"
#include "utils/simError.h"

namespace OpenMD {

  MultipoleAtomTypesSectionParser::MultipoleAtomTypesSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("MultipoleAtomTypes");
  }

  void MultipoleAtomTypesSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();    

    // name multipole_type theta phi psi 
    // "name" must match the name in the AtomTypes section
    // avaliable multipole type is d (dipole), s (split dipole) ,  q (quadrupoles), dq(dipole plus quadrupole)
    // and sq(split dipole plus quadrupole)
    // Directionality for dipoles and quadrupoles is given by three euler angles (phi, theta, psi),
    //because the body-fixed reference frame for directional atoms is determined by the *mass* 
    //distribution and not by the charge distribution.  
    // Dipoles are given in units of Debye  
    // Quadrupoles are given in units of 
    // examples:
    // name d phi theta psi dipole_moment
    // name s phi theta psi dipole_moment splitdipole_distance
    // name q phi theta psi Qxx Qyy Qzz
    // name dq phi theta psi dipole_moment Qxx Qyy Qzz
    // name sq phi theta psi dipole_moment splitdipole_distance Qxx Qyy Qzz
        
    if (nTokens < 5)  {
      sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
    } else {

      std::string atomTypeName = tokenizer.nextToken();    
      std::string multipoleType = tokenizer.nextToken();
      RealType phi = tokenizer.nextTokenAsDouble() * NumericConstant::PI /180.0;
      RealType theta = tokenizer.nextTokenAsDouble() * NumericConstant::PI /180.0;
      RealType psi = tokenizer.nextTokenAsDouble() * NumericConstant::PI /180.0;        
      nTokens -=  5;

      AtomType* atomType = ff.getAtomType(atomTypeName);
      if (atomType == NULL) {
	sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: Can not find matched AtomType[%s] at line %d\n",
		atomTypeName.c_str(), lineNo);
	painCave.isFatal = 1;
	simError();
      }
        
      MultipoleAdapter ma = MultipoleAdapter(atomType);

      RotMat3x3d electroBodyFrame(0.0);

      electroBodyFrame.setupRotMat(phi, theta, psi);

      RealType dipoleMoment(0);
      RealType splitDipoleDistance(0);
      Vector3d quadrupoleMoments(V3Zero);
      bool isDipole(false);
      bool isSplitDipole(false);
      bool isQuadrupole(false);
      
      if (multipoleType== "d") {
	parseDipole(tokenizer, dipoleMoment, lineNo);
        isDipole = true;
      } else if (multipoleType== "s") {
	parseSplitDipole(tokenizer, dipoleMoment, splitDipoleDistance, lineNo);
        isDipole = true;
        isSplitDipole = true;
      } else if (multipoleType== "q") {
	parseQuadrupole( tokenizer, quadrupoleMoments, lineNo);
        isQuadrupole = true;
      } else if (multipoleType== "dq") {
	parseDipole(tokenizer, dipoleMoment, lineNo);
        isDipole = true;
	parseQuadrupole( tokenizer, quadrupoleMoments, lineNo);
        isQuadrupole = true;
      } else if (multipoleType== "sq") {
	parseSplitDipole(tokenizer, dipoleMoment, splitDipoleDistance, lineNo);
        isDipole = true;
        isSplitDipole = true;
	parseQuadrupole( tokenizer, quadrupoleMoments, lineNo);
        isQuadrupole = true;
      } else {
	sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: unrecognized multiple type at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      }

      ma.makeMultipole(electroBodyFrame, dipoleMoment, splitDipoleDistance,
                       quadrupoleMoments, isDipole, isSplitDipole, 
                       isQuadrupole);
    }
  }

  void MultipoleAtomTypesSectionParser::parseDipole(StringTokenizer& tokenizer, 
						    RealType& dipoleMoment,
                                                    int lineNo) {

    if (tokenizer.hasMoreTokens()) {
      dipoleMoment = tokenizer.nextTokenAsDouble();    
    } else {
      sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
    }
  }

  void MultipoleAtomTypesSectionParser::parseSplitDipole(StringTokenizer& tokenizer, 
                                                         RealType& dipoleMoment,
                                                         RealType& splitDipoleDistance,
                                                         int lineNo) {

    if (tokenizer.hasMoreTokens()) {
      parseDipole(tokenizer, dipoleMoment, lineNo);    
      splitDipoleDistance = tokenizer.nextTokenAsDouble();
    } else {
      sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
    }
  }

  void MultipoleAtomTypesSectionParser::parseQuadrupole(StringTokenizer& tokenizer,
                                                        Vector3d& quadrupoleMoments,
                                                        int lineNo) {
    int nTokens = tokenizer.countTokens();   
    if (nTokens >= 3) {
      
      quadrupoleMoments[0] = tokenizer.nextTokenAsDouble();
      quadrupoleMoments[1] = tokenizer.nextTokenAsDouble();
      quadrupoleMoments[2] = tokenizer.nextTokenAsDouble();

      RealType trace =  quadrupoleMoments.sum();
      
      if (fabs(trace) > OpenMD::epsilon) {
	sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: the trace of quadrupole moments is not zero at line %d\n",
		lineNo);
	painCave.isFatal = 1;
	simError();
      }

    } else {
      sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
    }
  }


} //end namespace OpenMD



