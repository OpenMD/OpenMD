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

#ifndef SELECTION_SELECTIONEVALUATOR_HPP
#define SELECTION_SELECTIONEVALUATOR_HPP

#include <map>
#include <string>
#include <vector>

#include "brains/SimInfo.hpp"
#include "selection/Token.hpp"
#include "selection/SelectionCompiler.hpp"
#include "utils/BitSet.hpp"

namespace oopse {

class Context {
    public:
        
        void clear() {
            linenumbers.clear();
            lineIndices.clear();
            aatoken.clear();
        }
        
        std::string filename;
        std::string script;
        std::vector<int> linenumbers;
        std::vector<int> lineIndices;
        std::vector<std::vector<Token> > aatoken;
        int pc;
};

/**
 * @class SelectionEvaluator SelectionEvaluator.hpp "selection/SelectionEvaluator"
 * @brief Evalute the tokens compiled by SelectionCompiler and return a BitSet
 */
class SelectionEvaluator{
    public:

        SelectionEvaluator(SimInfo* info, const std::string& script);


        BitSet evaluate();
        
        BitSet evaluate(int frameNo);

        /**
         * Tests if the result from evaluation of script is static.
         * @return true if the result keeps the same even the frame change, otherwise return false
         */         
        bool isStatic();

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
            int ichBegin = lineIndices[pc];
            int ichEnd;
            if ((ichEnd = script.find('\r', ichBegin)) == std::string::npos &&
                (ichEnd = script.find('\n', ichBegin)) == std::string::npos) {
                ichEnd = script.size();
            }            
            return script.substr(ichBegin, ichEnd);
        }
  
    private:

        void clearState();
        
        bool loadScript(const std::string& filename, const std::string& script); 
        bool loadScriptString(const std::string& script);
        bool loadScriptFileInternal(const std::string& filename);
        bool loadScriptFile(const std::string& filename);

        void clearDefinitionsAndLoadPredefined();
         
        void define();
        void select();
        void predefine();

        void instructionDispatchLoop();

        void withinInstruction(const Token& instruction, BitSet& , BitSet&);

        BitSet expression(std::vector<Token>& tokens, int pc);


        void evalError(const std::string& message);

        void unrecognizedCommand(const Token& token) {
            evalError("unrecognized command:" + token.value);
        }        
        
        SelectionCompiler compiler;

        const static int scriptLevelMax = 10;
        int scriptLevel;

        Context stack[scriptLevelMax];

        std::string filename;
        std::string script;
        std::vector<int> linenumbers;
        std::vector<int> lineIndices;
        std::vector<std::vector<Token> > aatoken;
        int pc; // program counter

        bool error;
        std::string errorMessage;

        std::vector<Token> statement;
        int statementLength;

        SimInfo* info;

        std::map<std::string, std::vector<Token> > variables;
};

}
#endif
