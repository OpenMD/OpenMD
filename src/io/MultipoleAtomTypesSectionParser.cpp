/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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

#include "io/MultipoleAtomTypesSectionParser.hpp"

#include "brains/ForceField.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  MultipoleAtomTypesSectionParser::MultipoleAtomTypesSectionParser(
      ForceFieldOptions&) {
    setSectionName("MultipoleAtomTypes");
  }

  void MultipoleAtomTypesSectionParser::parseLine(ForceField& ff,
                                                  const std::string& line,
                                                  int lineNo) {
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();

    // name multipole_type theta phi psi
    // "name" must match the name in the AtomTypes section
    //
    // avaliable multipole types are:
    // d  (dipole)
    // q  (quadrupole)
    // dq (dipole plus quadrupole)
    //
    // Directionality for dipoles and quadrupoles is given by three
    // euler angles (phi, theta, psi), because the body-fixed
    // reference frame for directional atoms is determined by the
    // *mass* distribution and not by the charge distribution.
    //
    // Dipoles are given in units of Debye
    // Quadrupoles are given in units of esu centibarn
    //
    // Examples:
    //
    // name d phi theta psi dipole_moment
    // name q phi theta psi Qxx Qyy Qzz
    // name dq phi theta psi dipole_moment Qxx Qyy Qzz

    if (nTokens < 5) {
      sprintf(painCave.errMsg,
              "MultipoleAtomTypesSectionParser Error: Not enough tokens at "
              "line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();
    } else {
      std::string atomTypeName  = tokenizer.nextToken();
      std::string multipoleType = tokenizer.nextToken();
      RealType phi   = tokenizer.nextTokenAsDouble() * Constants::PI / 180.0;
      RealType theta = tokenizer.nextTokenAsDouble() * Constants::PI / 180.0;
      RealType psi   = tokenizer.nextTokenAsDouble() * Constants::PI / 180.0;

      AtomType* atomType = ff.getAtomType(atomTypeName);
      if (atomType == NULL) {
        sprintf(painCave.errMsg,
                "MultipoleAtomTypesSectionParser Error: Can not find matched "
                "AtomType[%s] "
                "at line %d\n",
                atomTypeName.c_str(), lineNo);
        painCave.isFatal = 1;
        simError();
      }

      MultipoleAdapter ma = MultipoleAdapter(atomType);

      RotMat3x3d eFrame(0.0);

      eFrame.setupRotMat(phi, theta, psi);

      RealType dipoleMoment(0);
      Vector3d dipole(V3Zero);
      Vector3d quadrupoleMoments(V3Zero);
      Mat3x3d quadrupole(0.0);

      bool isDipole(false);
      bool isQuadrupole(false);

      if (multipoleType == "d") {
        parseDipole(tokenizer, dipoleMoment, lineNo);
        isDipole = true;
      } else if (multipoleType == "s") {
        sprintf(painCave.errMsg,
                "MultipoleAtomTypesSectionParser Error: \n"
                "\tsplit dipoles (type s) have been deprecated (line: %d)\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else if (multipoleType == "q") {
        parseQuadrupole(tokenizer, quadrupoleMoments, lineNo);
        isQuadrupole = true;
      } else if (multipoleType == "dq") {
        parseDipole(tokenizer, dipoleMoment, lineNo);
        isDipole = true;
        parseQuadrupole(tokenizer, quadrupoleMoments, lineNo);
        isQuadrupole = true;
      } else if (multipoleType == "sq") {
        sprintf(painCave.errMsg,
                "MultipoleAtomTypesSectionParser Error: \n"
                "\tsplit dipole quadrupoles (type sq) have been deprecated "
                "(line: %d)\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        sprintf(painCave.errMsg,
                "MultipoleAtomTypesSectionParser Error: unrecognized multiple "
                "type at line "
                "%d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      }
      if (isDipole) dipole = dipoleMoment * eFrame.transpose() * V3Z;
      if (isQuadrupole) {
        quadrupole(0, 0) = quadrupoleMoments(0);
        quadrupole(1, 1) = quadrupoleMoments(1);
        quadrupole(2, 2) = quadrupoleMoments(2);
        quadrupole       = eFrame.transpose() * quadrupole * eFrame;
      }

      ma.makeMultipole(dipole, quadrupole, isDipole, isQuadrupole);
    }
  }

  void MultipoleAtomTypesSectionParser::parseDipole(StringTokenizer& tokenizer,
                                                    RealType& dipoleMoment,
                                                    int lineNo) {
    if (tokenizer.hasMoreTokens()) {
      dipoleMoment = tokenizer.nextTokenAsDouble();
    } else {
      sprintf(painCave.errMsg,
              "MultipoleAtomTypesSectionParser Error: Not enough tokens at "
              "line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();
    }
  }

  void MultipoleAtomTypesSectionParser::parseQuadrupole(
      StringTokenizer& tokenizer, Vector3d& quadrupoleMoments, int lineNo) {
    int nTokens = tokenizer.countTokens();
    if (nTokens >= 3) {
      quadrupoleMoments[0] = tokenizer.nextTokenAsDouble();
      quadrupoleMoments[1] = tokenizer.nextTokenAsDouble();
      quadrupoleMoments[2] = tokenizer.nextTokenAsDouble();

      // RealType trace =  quadrupoleMoments.sum();
      //
      // if (fabs(trace) > OpenMD::epsilon) {
      //   sprintf(painCave.errMsg, "MultipoleAtomTypesSectionParser Error: the
      //   trace of quadrupole moments is not zero at line %d\n", 	lineNo);
      //   painCave.isFatal = 1; simError();
      // }

    } else {
      sprintf(painCave.errMsg,
              "MultipoleAtomTypesSectionParser Error: Not enough tokens at "
              "line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();
    }
  }
}  // namespace OpenMD
