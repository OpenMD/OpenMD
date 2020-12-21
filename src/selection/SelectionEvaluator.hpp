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

#ifndef SELECTION_SELECTIONEVALUATOR_HPP
#define SELECTION_SELECTIONEVALUATOR_HPP

#include <map>
#include <string>
#include <vector>
#include <fstream>

#include "brains/SimInfo.hpp"
#include "selection/SelectionToken.hpp"
#include "selection/SelectionCompiler.hpp"
#include "selection/NameFinder.hpp"
#include "selection/DistanceFinder.hpp"
#include "selection/HullFinder.hpp"
#include "selection/IndexFinder.hpp"
#include "selection/SelectionSet.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/StringUtils.hpp"
namespace OpenMD {


  /**
   * @class SelectionEvaluator SelectionEvaluator.hpp "selection/SelectionEvaluator"
   * @brief Evalute the tokens compiled by SelectionCompiler and return a OpenMDBitSet
   */
  class SelectionEvaluator{
  public:

    SelectionEvaluator(SimInfo* info);

    bool loadScriptString(const std::string& script);
    bool loadScriptFile(const std::string& filename);
        
    SelectionSet evaluate();
    SelectionSet evaluate(int frame);
        
    /**
     * Tests if the result from evaluation of script is dynamic.
     */         
    bool isDynamic() {
      return isDynamic_;
    }

    bool hadRuntimeError() const{
      return error;
    }

    std::string getErrorMessage() const {
      return errorMessage;
    }


    int getLinenumber() {
      return linenumbers[pc];
    }

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
        sprintf( painCave.errMsg,
                 "SelectionEvaluator Error: %s\n", "No Surface Area For You!");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
        return 0.0;
      }
    }
    bool hasVolume() { return hasVolume_; }
    RealType getVolume() { 
      if (hasVolume_) {
        return volume_;
      } else {
        sprintf( painCave.errMsg,
                 "SelectionEvaluator Error: %s\n", "No Volume For You!");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
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
    void withinInstruction(const Token& instruction, SelectionSet& bs, int frame);

    SelectionSet alphaHullInstruction(const Token& instruction);
    SelectionSet alphaHullInstruction(const Token& instruction, int frame);
    
    SelectionSet allInstruction();
        
    SelectionSet comparatorInstruction(const Token& instruction); 
    SelectionSet comparatorInstruction(const Token& instruction, int frame); 
    void compareProperty(StuntDouble* sd, SelectionSet& bs, int property, int comparator, float comparisonValue);
    void compareProperty(StuntDouble* sd, SelectionSet& bs, int property, int comparator, float comparisonValue, int frame);
    void compareProperty(Molecule* mol, SelectionSet& bs, int property, int comparator, float comparisonValue);
    void compareProperty(Molecule* mol, SelectionSet& bs, int property, int comparator, float comparisonValue, int frame);
    SelectionSet nameInstruction(const std::string& name);
    SelectionSet indexInstruction(const boost::any& value);
    SelectionSet expression(const std::vector<Token>& tokens, int pc);
    SelectionSet expression(const std::vector<Token>& tokens, int pc, int frame);

    SelectionSet lookupValue(const std::string& variable);

    SelectionSet hull();
    SelectionSet hull(int frame);
        
    void evalError(const std::string& message) {
      sprintf( painCave.errMsg,
               "SelectionEvaluator Error: %s\n", message.c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }

    void unrecognizedCommand(const Token& token) {
      evalError("unrecognized command:" + boost::any_cast<std::string>(token.value));
    }        

    void unrecognizedExpression() {
      evalError("unrecognized expression");
    }

    void unrecognizedAtomProperty(int property){
      evalError("unrecognized atom property");
    }
    
    void unrecognizedMoleculeProperty(int property){
      evalError("unrecognized molecule property");
    }

    void unrecognizedIdentifier(const std::string& identifier) {
      evalError("unrecognized identifier:" + identifier);
    }    

    void invalidIndexRange(std::pair<int, int> range) {
      evalError("invalid index range: [" + toString(range.first) + ", " + toString(range.second) + ")");
    }

    void invalidIndex(int index) {
      evalError("invalid index : " + toString(index) );
    }

        
    bool containDynamicToken(const std::vector<Token>& tokens);

    RealType getCharge(Atom* atom);
    RealType getCharge(Atom* atom, int frame);
        
    SelectionCompiler compiler;

    //const static int scriptLevelMax = 10;
    //int scriptLevel;

    //Context stack[scriptLevelMax];

    std::string filename;
    std::string script;
    std::vector<int> linenumbers;
    std::vector<int> lineIndices;
    std::vector<std::vector<Token> > aatoken;
    unsigned int pc; // program counter

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

    typedef std::map<std::string, boost::any > VariablesType;
    VariablesType variables;

    bool isDynamic_ {false};
    bool isLoaded_ {false};
    bool hasSurfaceArea_ {false};
    RealType surfaceArea_;
    bool hasVolume_ {false};
    RealType volume_;        
  };
}

#endif
