/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "io/ShapeAtomTypesSectionParser.hpp"

#include <memory>

#include "brains/ForceField.hpp"
#include "math/RealSphericalHarmonic.hpp"
#include "math/SquareMatrix3.hpp"
#include "types/AtomType.hpp"
#include "types/DirectionalAdapter.hpp"
#include "types/ShapeAtomType.hpp"
#include "utils/CaseConversion.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ShapeAtomTypesSectionParser::ShapeAtomTypesSectionParser(ForceFieldOptions&) {
    setSectionName("ShapeAtomTypes");
  }

  void ShapeAtomTypesSectionParser::parseLine(ForceField& ff,
                                              const std::string& line,
                                              int lineNo) {
    StringTokenizer tokenizer(line);

    if (tokenizer.countTokens() >= 2) {
      std::string shapeTypeName = tokenizer.nextToken();
      std::string shapeFile     = tokenizer.nextToken();

      AtomType* atomType = ff.getAtomType(shapeTypeName);
      if (atomType == NULL) {
        atomType  = new AtomType();
        int ident = ff.getNAtomType() + 1;
        atomType->setIdent(ident);
        atomType->setName(shapeTypeName);
        ff.addAtomType(shapeTypeName, atomType);
      }

      parseShapeFile(ff, shapeFile, atomType);

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ShapesAtomTypesSectionParser Error: "
               "Not enough tokens at line %d\n",
               lineNo);
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  }

  void ShapeAtomTypesSectionParser::parseShapeFile(ForceField&,
                                                   std::string& shapeFileName,
                                                   AtomType* at) {
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string token;

    Mat3x3d momInert;
    RealSphericalHarmonic* rsh;
    std::vector<RealSphericalHarmonic*> functionVector;
    ifstrstream shapeStream;
    std::string tempString;
    std::string ffPath;
    char* tempPath;

    tempPath = getenv("FORCE_PARAM_PATH");

    if (tempPath == NULL) {
      // convert a macro from compiler to a string in c++
      STR_DEFINE(ffPath, FRC_PATH);
    } else {
      ffPath = tempPath;
    }

    shapeStream.open(shapeFileName.c_str());

    if (!shapeStream.is_open()) {
      tempString = ffPath;
      tempString += "/";
      tempString += shapeFileName;
      shapeFileName = tempString;

      shapeStream.open(shapeFileName.c_str());

      if (!shapeStream.is_open()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Error opening the shape file:\n"
                 "\t%s\n"
                 "\tHave you tried setting the FORCE_PARAM_PATH environment "
                 "variable?\n",
                 shapeFileName.c_str());
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }
    }

    ShapeAtomType* st = new ShapeAtomType();

    // first parse the info. in the ShapeInfo section
    findBegin(shapeStream, "ShapeInfo");
    shapeStream.getline(buffer, bufferSize);

    // loop over the interior of the ShapeInfo section
    while (!shapeStream.eof()) {
      // toss comment lines
      if (buffer[0] != '!' && buffer[0] != '#') {
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo(buffer);
        // blank lines are ignored
        if (tokenInfo.countTokens() != 0) {
          if (tokenInfo.countTokens() < 5) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a ShapeInfo line in file: %s\n",
                     shapeFileName.c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          } else {
            tokenInfo.skipToken();
            at->setMass(tokenInfo.nextTokenAsDouble());
            DirectionalAdapter da = DirectionalAdapter(at);
            Mat3x3d I;
            I(0, 0) = tokenInfo.nextTokenAsDouble();
            I(1, 1) = tokenInfo.nextTokenAsDouble();
            I(2, 2) = tokenInfo.nextTokenAsDouble();
            da.makeDirectional(I);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }

    // now grab the contact functions
    findBegin(shapeStream, "ContactFunctions");
    functionVector.clear();

    shapeStream.getline(buffer, bufferSize);
    while (!shapeStream.eof()) {
      // toss comment lines
      if (buffer[0] != '!' && buffer[0] != '#') {
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo1(buffer);
        // blank lines are ignored
        if (tokenInfo1.countTokens() != 0) {
          if (tokenInfo1.countTokens() < 4) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a ContactFunctions line in file: %s\n",
                     shapeFileName.c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          } else {
            // read in a spherical harmonic function
            rsh = new RealSphericalHarmonic();
            rsh->setL(tokenInfo1.nextTokenAsInt());
            rsh->setM(tokenInfo1.nextTokenAsInt());
            token = tokenInfo1.nextToken();
            toLower(token);
            if (token == "sin")
              rsh->makeSinFunction();
            else
              rsh->makeCosFunction();
            rsh->setCoefficient(tokenInfo1.nextTokenAsDouble());

            functionVector.push_back(rsh);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }

    // pass contact functions to ShapeType
    st->setContactFuncs(functionVector);

    // now grab the range functions
    findBegin(shapeStream, "RangeFunctions");
    functionVector.clear();

    shapeStream.getline(buffer, bufferSize);
    while (!shapeStream.eof()) {
      // toss comment lines
      if (buffer[0] != '!' && buffer[0] != '#') {
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo2(buffer);
        // blank lines are ignored
        if (tokenInfo2.countTokens() != 0) {
          if (tokenInfo2.countTokens() < 4) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a RangeFunctions line in file: %s\n",
                     shapeFileName.c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          } else {
            // read in a spherical harmonic function
            rsh = new RealSphericalHarmonic();
            rsh->setL(tokenInfo2.nextTokenAsInt());
            rsh->setM(tokenInfo2.nextTokenAsInt());
            token = tokenInfo2.nextToken();
            toLower(token);
            if (token == "sin")
              rsh->makeSinFunction();
            else
              rsh->makeCosFunction();
            rsh->setCoefficient(tokenInfo2.nextTokenAsDouble());

            functionVector.push_back(rsh);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }

    // pass range functions to ShapeType
    st->setRangeFuncs(functionVector);

    // finally grab the strength functions
    findBegin(shapeStream, "StrengthFunctions");
    functionVector.clear();

    shapeStream.getline(buffer, bufferSize);
    while (!shapeStream.eof()) {
      // toss comment lines
      if (buffer[0] != '!' && buffer[0] != '#') {
        // end marks section completion
        if (isEndLine(buffer)) break;
        StringTokenizer tokenInfo3(buffer);
        // blank lines are ignored
        if (tokenInfo3.countTokens() != 0) {
          if (tokenInfo3.countTokens() < 4) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "ShapesAtomTypesSectionParser Error: Not enough "
                     "information on a StrengthFunctions line in file: %s\n",
                     shapeFileName.c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal  = 1;
            simError();
          } else {
            // read in a spherical harmonic function
            rsh = new RealSphericalHarmonic();
            rsh->setL(tokenInfo3.nextTokenAsInt());
            rsh->setM(tokenInfo3.nextTokenAsInt());
            token = tokenInfo3.nextToken();
            toLower(token);
            if (token == "sin")
              rsh->makeSinFunction();
            else
              rsh->makeCosFunction();
            rsh->setCoefficient(tokenInfo3.nextTokenAsDouble());

            functionVector.push_back(rsh);
          }
        }
      }
      shapeStream.getline(buffer, bufferSize);
    }

    // pass strength functions to ShapeType
    st->setStrengthFuncs(functionVector);
    at->addProperty(
        std::shared_ptr<GenericData>(new ShapeAtypeData("Shape", st)));
    //  delete shapeStream;
  }
}  // namespace OpenMD
