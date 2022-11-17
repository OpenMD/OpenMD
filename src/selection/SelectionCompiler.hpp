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

#ifndef SELECTION_SELECTIONCOMPILER_HPP
#define SELECTION_SELECTIONCOMPILER_HPP

#include <any>
#include <iostream>
#include <string>
#include <vector>

#include "brains/SimInfo.hpp"
#include "selection/SelectionToken.hpp"
#include "selection/TokenMap.hpp"

namespace OpenMD {

  /**
   * @class SelectionCompiler SelectionCompiler.hpp
   "selection/SelectionCompiler.hpp"
   * @brief compile a selection script to tokens
   * @todo document
   * <pre>

   expression       :: = clauseOr

   clauseOr         ::= clauseAnd {OR clauseAnd}*

   clauseAnd        ::= clauseNot {AND clauseNot}*

   clauseNot        ::= NOT clauseNot | clausePrimitive

   clausePrimitive  ::= clauseComparator |
   clauseWithin |
   clauseAlphaHull |
   clauseName |
   none | all |
   ( clauseOr )

   clauseComparator ::= atomproperty comparatorop integer

   clauseWithin     ::= WITHIN ( clauseDistance , expression )

   clauseDistance   ::= integer | decimal

   clauseName::= *|string{.string{.string}}


   * </pre>
   */
  class SelectionCompiler {
  public:
    bool compile(const std::string& filename, const std::string& script);

    std::vector<int> getLineNumbers() { return lineNumbers; }

    std::vector<int> getLineIndices() { return lineIndices; }

    std::vector<std::vector<Token>> getAatokenCompiled() {
      return aatokenCompiled;
    }

    std::string getErrorMessage() {
      std::string strError = errorMessage;
      strError += " : " + errorLine + "\n";

      if (!filename.empty()) { strError += filename; }

      return strError;
    }

  private:
    bool internalCompile();

    bool lookingAtLeadingWhitespace();
    // bool lookingAtComment();
    bool lookingAtEndOfLine();
    bool lookingAtEndOfStatement();
    bool lookingAtString();
    bool lookingAtDecimal(bool allowNegative);
    bool lookingAtInteger(bool allowNegative);
    bool lookingAtLookupToken();
    bool lookingAtSpecialString();

    std::string getUnescapedStringLiteral();
    int getHexitValue(char ch);

    bool compileCommand(const std::vector<Token>& ltoken);
    bool compileExpression();
    bool compileExpression(int itoken);

    bool clauseOr();
    bool clauseAnd();
    bool clauseNot();
    bool clausePrimitive();
    bool clauseWithin();
    bool clauseAlphaHull();
    bool clauseComparator();
    bool clauseChemObjName();
    bool clauseIndex();
    Token tokenNext();
    std::any valuePeek();
    int tokPeek();

    bool addTokenToPostfix(const Token& token);
    bool isNameValid(const std::string& name);

    bool compileError(const std::string& errorMsg) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SelectionCompiler Error: %s\n", errorMsg.c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();

      error              = true;
      this->errorMessage = errorMsg;
      return false;
    }

    bool commandExpected() { return compileError("command expected"); }

    bool invalidExpressionToken(const std::string& ident) {
      return compileError("invalid expression token:" + ident);
    }

    bool unrecognizedToken() { return compileError("unrecognized token"); }

    bool badArgumentCount() { return compileError("bad argument count"); }

    bool endOfExpressionExpected() {
      return compileError("end of expression expected");
    }

    bool leftParenthesisExpected() {
      return compileError("left parenthesis expected");
    }

    bool rightParenthesisExpected() {
      return compileError("right parenthesis expected");
    }

    bool commaExpected() { return compileError("comma expected"); }

    bool unrecognizedExpressionToken() {
      std::any tmp = valuePeek();
      std::string tokenStr;

      try {
        tokenStr = std::any_cast<std::string>(tmp);
      } catch (const std::bad_any_cast&) {
        return compileError("any_cast error");
      }

      return compileError("unrecognized expression token:" + tokenStr);
    }

    bool comparisonOperatorExpected() {
      return compileError("comparison operator expected");
    }

    bool numberExpected() { return compileError("number expected"); }

    bool numberOrKeywordExpected() {
      return compileError("number or keyword expected");
    }

    std::string filename;
    std::string script;

    std::vector<int> lineNumbers;
    std::vector<int> lineIndices;
    std::vector<std::vector<Token>> aatokenCompiled;

    bool error;
    std::string errorMessage;
    std::string errorLine;

    int cchScript;
    short lineCurrent;

    int ichToken;
    int cchToken;
    std::vector<Token> atokenCommand;

    int ichCurrentCommand;

    std::vector<Token> ltokenPostfix;
    std::vector<Token> atokenInfix;
    std::size_t itokenInfix;

    // std::vector<Token> compiledTokens_;
  };
}  // namespace OpenMD

#endif
