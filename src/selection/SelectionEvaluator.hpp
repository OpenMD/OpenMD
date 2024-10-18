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

#ifndef SELECTION_SELECTIONEVALUATOR_HPP
#define SELECTION_SELECTIONEVALUATOR_HPP

#include <any>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "brains/SimInfo.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/DistanceFinder.hpp"
#include "selection/HullFinder.hpp"
#include "selection/IndexFinder.hpp"
#include "selection/NameFinder.hpp"
#include "selection/SelectionCompiler.hpp"
#include "selection/SelectionSet.hpp"
#include "selection/SelectionToken.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  /**
   * @class SelectionEvaluator SelectionEvaluator.hpp
   * "selection/SelectionEvaluator"
   * @brief Evalute the tokens compiled by SelectionCompiler and return a
   * OpenMDBitSet
   */
  class SelectionEvaluator {
  public:
    SelectionEvaluator(SimInfo* info);

    bool loadScriptString(const std::string& script);
    bool loadScriptFile(const std::string& filename);

    SelectionSet evaluate();
    SelectionSet evaluate(int frame);

    /**
     * Tests if the result from evaluation of script is dynamic.
     */
    bool isDynamic() { return isDynamic_; }

    bool hadRuntimeError() const { return error; }

    std::string getErrorMessage() const { return errorMessage; }

    int getLinenumber() { return linenumbers[pc]; }

    std::string getLine() {
      std::size_t ichBegin = lineIndices[pc];
      std::size_t ichEnd;
      if ((ichEnd = script.find('\r', ichBegin)) == std::string::npos &&
          (ichEnd = script.find('\n', ichBegin)) == std::string::npos) {
        ichEnd = script.size();
      }
      return script.substr(ichBegin, ichEnd);
    }
    bool hasSurfaceArea() { return hasSurfaceArea_; }
    RealType getSurfaceArea() {
      if (hasSurfaceArea_) {
        return surfaceArea_;
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SelectionEvaluator Error: %s\n", "No Surface Area For You!");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
        return 0.0;
      }
    }
    bool hasVolume() { return hasVolume_; }
    RealType getVolume() {
      if (hasVolume_) {
        return volume_;
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "SelectionEvaluator Error: %s\n", "No Volume For You!");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
        return 0.0;
      }
    }

  private:
    void clearState();

    bool loadScript(const std::string& filename, const std::string& script);

    bool loadScriptFileInternal(const std::string& filename);

    SelectionSet createSelectionSets();
    void clearDefinitionsAndLoadPredefined();

    void define();
    void select(SelectionSet& bs);
    void select(SelectionSet& bs, int frame);
    void predefine(const std::string& script);

    void instructionDispatchLoop(SelectionSet& bs);
    void instructionDispatchLoop(SelectionSet& bs, int frame);

    void withinInstruction(const Token& instruction, SelectionSet& bs);
    void withinInstruction(const Token& instruction, SelectionSet& bs,
                           int frame);

    SelectionSet alphaHullInstruction(const Token& instruction);
    SelectionSet alphaHullInstruction(const Token& instruction, int frame);

    SelectionSet allInstruction();

    SelectionSet comparatorInstruction(const Token& instruction);
    SelectionSet comparatorInstruction(const Token& instruction, int frame);
    void compareProperty(StuntDouble* sd, SelectionSet& bs, int property,
                         int comparator, float comparisonValue);
    void compareProperty(StuntDouble* sd, SelectionSet& bs, int property,
                         int comparator, float comparisonValue, int frame);
    void compareProperty(Molecule* mol, SelectionSet& bs, int property,
                         int comparator, float comparisonValue);
    void compareProperty(Molecule* mol, SelectionSet& bs, int property,
                         int comparator, float comparisonValue, int frame);
    SelectionSet nameInstruction(const std::string& name);
    SelectionSet indexInstruction(const std::any& value);
    SelectionSet expression(const std::vector<Token>& tokens, int pc);
    SelectionSet expression(const std::vector<Token>& tokens, int pc,
                            int frame);

    SelectionSet lookupValue(const std::string& variable);

    SelectionSet hull();
    SelectionSet hull(int frame);

    void evalError(const std::string& message) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SelectionEvaluator Error: %s\n", message.c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    void unrecognizedCommand(const Token& token) {
      evalError("unrecognized command:" +
                std::any_cast<std::string>(token.value));
    }

    void unrecognizedExpression() { evalError("unrecognized expression"); }

    void unrecognizedAtomProperty(int) {
      evalError("unrecognized atom property");
    }

    void unrecognizedMoleculeProperty(int) {
      evalError("unrecognized molecule property");
    }

    void unrecognizedIdentifier(const std::string& identifier) {
      evalError("unrecognized identifier:" + identifier);
    }

    void invalidIndexRange(std::pair<int, int> range) {
      evalError("invalid index range: [" + toString(range.first) + ", " +
                toString(range.second) + ")");
    }

    void invalidIndex(int index) {
      evalError("invalid index : " + toString(index));
    }

    bool containDynamicToken(const std::vector<Token>& tokens);

    RealType getCharge(Atom* atom);
    RealType getCharge(Atom* atom, int frame);

    SelectionCompiler compiler;

    // const static int scriptLevelMax = 10;
    // int scriptLevel;

    // Context stack[scriptLevelMax];

    std::string filename;
    std::string script;
    std::vector<int> linenumbers;
    std::vector<int> lineIndices;
    std::vector<std::vector<Token>> aatoken;
    unsigned int pc;  // program counter

    bool error {false};
    std::string errorMessage;

    std::vector<Token> statement;
    int statementLength;

    SimInfo* info {nullptr};
    NameFinder nameFinder;
    DistanceFinder distanceFinder;
    HullFinder hullFinder;
    AlphaHullFinder alphaHullFinder;
    IndexFinder indexFinder;
    vector<int> nObjects;

    using VariablesType = std::map<std::string, std::any>;
    VariablesType variables;

    bool isDynamic_ {false};
    bool isLoaded_ {false};
    bool hasSurfaceArea_ {false};
    RealType surfaceArea_;
    bool hasVolume_ {false};
    RealType volume_;
  };
}  // namespace OpenMD

#endif
