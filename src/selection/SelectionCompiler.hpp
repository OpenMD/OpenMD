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

#ifndef SELECTION_SELECTIONCOMPILER_HPP
#define SELECTION_SELECTIONCOMPILER_HPP
#include <iostream>
#include <string>
#include <vector>

#include "selection/SelectionToken.hpp"
#include "selection/TokenMap.hpp"
#include "brains/SimInfo.hpp"

namespace OpenMD {


  /**
   * @class SelectionCompiler SelectionCompiler.hpp "selection/SelectionCompiler.hpp"
   * @brief compile a selection script to tokens
   * @todo document
   * <pre>

   expression       :: = clauseOr

   clauseOr         ::= clauseAnd {OR clauseAnd}*

   clauseAnd        ::= clauseNot {AND clauseNot}*

   clauseNot        ::= NOT clauseNot | clausePrimitive

   clausePrimitive  ::= clauseComparator |
   clauseWithin |
   clauseName |
   none | all |
   ( clauseOr )

   clauseComparator ::= atomproperty comparatorop integer

   clauseWithin     ::= WITHIN ( clauseDistance , expression )

   clauseDistance   ::= integer | decimal
        
   clauseName::= *|string{.string{.string}}


   * </pre>
   */
  class SelectionCompiler{
  public:
    bool compile(const std::string& filename, const std::string& script );
        

    std::vector<int> getLineNumbers() {
      return lineNumbers;
    }

    std::vector<int> getLineIndices() {
      return lineIndices;
    }

    std::vector<std::vector<Token> > getAatokenCompiled() {
      return aatokenCompiled;
    }

    std::string getErrorMessage() {
      std::string strError = errorMessage;
      strError += " : " + errorLine + "\n";

      if (!filename.empty()) {
	strError += filename;
      }

      return strError;
    }

        
  private:

    bool internalCompile();


    bool lookingAtLeadingWhitespace();
    //bool lookingAtComment();
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
    bool clauseComparator();
    bool clauseChemObjName();        
    bool clauseIndex();
    Token tokenNext();
    boost::any valuePeek();
    int tokPeek();

    bool addTokenToPostfix(const Token& token);
    bool isNameValid(const std::string& name);

    bool compileError(const std::string& errorMsg) {

      sprintf( painCave.errMsg,
               "SelectionCompiler Error: %s\n", errorMsg.c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();

      error = true;
      this->errorMessage = errorMsg;
      return false;
    }
        
    bool commandExpected() {
      return compileError("command expected");
    }

    bool invalidExpressionToken(const std::string& ident) {
      return compileError("invalid expression token:" + ident);
    }

    bool unrecognizedToken() {
      return compileError("unrecognized token");
    }

    bool badArgumentCount() {
      return compileError("bad argument count");
    }

    bool endOfExpressionExpected() {
      return compileError("end of expression expected");
    }

    bool leftParenthesisExpected() {
      return compileError("left parenthesis expected");
    }

    bool rightParenthesisExpected() {
      return compileError("right parenthesis expected");
    }

    bool commaExpected() {
      return compileError("comma expected");
    }

    bool unrecognizedExpressionToken() {
      boost::any tmp = valuePeek();
      std::string tokenStr;

      try {
	tokenStr = boost::any_cast<std::string>(tmp);                
      } catch(const boost::bad_any_cast &) {
	return compileError("any_cast error");
      }
            
      return compileError("unrecognized expression token:" + tokenStr);
    }

    bool comparisonOperatorExpected() {
      return compileError("comparison operator expected");
    }

    bool numberExpected() {
      return compileError("number expected");
    }        
        
    bool numberOrKeywordExpected() {
      return compileError("number or keyword expected");
    }        
        
    std::string filename;
    std::string script;

    std::vector<int> lineNumbers;
    std::vector<int> lineIndices;
    std::vector<std::vector<Token> >aatokenCompiled;

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
    int itokenInfix;

    //std::vector<Token> compiledTokens_;
  };

}
#endif

