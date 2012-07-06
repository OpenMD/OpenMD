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

#include "selection/SelectionCompiler.hpp"
#include "utils/StringUtils.hpp"
namespace OpenMD {

  bool SelectionCompiler::compile(const std::string& filename, 
                                  const std::string& script) {

    this->filename = filename;
    this->script = script;
    lineNumbers.clear();
    lineIndices.clear();
    aatokenCompiled.clear();
        
    if (internalCompile()) {
      return true;
    }
    
    int icharEnd;
    if ((icharEnd = script.find('\r', ichCurrentCommand)) == std::string::npos &&
        (icharEnd = script.find('\n', ichCurrentCommand)) == std::string::npos) {
      icharEnd = script.size();
    }
    errorLine = script.substr(ichCurrentCommand, icharEnd);
    return false;
  }

  bool SelectionCompiler::internalCompile(){

    cchScript = script.size();
    ichToken = 0;
    lineCurrent = 1;

    error = false;

    //std::vector<Token> lltoken;
    aatokenCompiled.clear();
    std::vector<Token> ltoken;

    Token tokenCommand;
    int tokCommand = Token::nada;

    for ( ; true; ichToken += cchToken) {
      if (lookingAtLeadingWhitespace())
	continue;
      //if (lookingAtComment())
      //    continue;
      bool endOfLine = lookingAtEndOfLine();
      if (endOfLine || lookingAtEndOfStatement()) {
	if (tokCommand != Token::nada) {
	  if (! compileCommand(ltoken)) {
	    return false;
	  }
	  aatokenCompiled.push_back(atokenCommand);
	  lineNumbers.push_back(lineCurrent);
	  lineIndices.push_back(ichCurrentCommand);
	  ltoken.clear();
	  tokCommand = Token::nada;
	}
            
	if (ichToken < cchScript) {
	  if (endOfLine)
	    ++lineCurrent;
	  continue;
	}
	break;
      }

      if (tokCommand != Token::nada) {
	if (lookingAtString()) {
	  std::string str = getUnescapedStringLiteral();
	  ltoken.push_back(Token(Token::string, str));
	  continue;
	}
	//if ((tokCommand & Token::specialstring) != 0 &&
	//    lookingAtSpecialString()) {
	//    std::string str = script.substr(ichToken, ichToken + cchToken);
	//    ltoken.push_back(Token(Token::string, str));
	//    continue;
	//}
	if (lookingAtDecimal((tokCommand & Token::negnums) != 0)) {
	  float value = lexi_cast<float>(script.substr(ichToken, cchToken));        
	  ltoken.push_back(Token(Token::decimal, boost::any(value)));
	  continue;
	}
	if (lookingAtInteger((tokCommand & Token::negnums) != 0)) {

	  int val = lexi_cast<int>(script.substr(ichToken, cchToken));
	  ltoken.push_back(Token(Token::integer,   boost::any(val)));
	  continue;
	}
      }
      
      if (lookingAtLookupToken()) {
	std::string ident = script.substr(ichToken, cchToken);
	Token token;            
	Token* pToken = TokenMap::getInstance()->getToken(ident);
	if (pToken != NULL) {
	  token = *pToken;
	} else {
	  token = Token(Token::identifier, ident);
	}
            
	int tok = token.tok;
            
	switch (tokCommand) {
	case Token::nada:
	  ichCurrentCommand = ichToken;
	  //tokenCommand = token;
	  tokCommand = tok;
	  if ((tokCommand & Token::command) == 0)
	    return commandExpected();
	  break;

	case Token::define:
	  if (ltoken.size() == 1) {
	    // we are looking at the variable name
	    if (tok != Token::identifier &&
		(tok & Token::predefinedset) != Token::predefinedset)
	      return invalidExpressionToken(ident);
	  } else {
	    // we are looking at the expression
	    if (tok != Token::identifier && 
		(tok & (Token::expression | Token::predefinedset)) == 0)
	      return invalidExpressionToken(ident);
	  }
                    
	  break;

	case Token::select:
	  if (tok != Token::identifier && (tok & Token::expression) == 0)
	    return invalidExpressionToken(ident);
	  break;
	}
	ltoken.push_back(token);
	continue;
      }

      if (ltoken.size() == 0) {
	return commandExpected();
      }
        
      return unrecognizedToken();
    }

    return true;
  }


  bool SelectionCompiler::lookingAtLeadingWhitespace() {

    int ichT = ichToken;
    while (ichT < cchScript && std::isspace(script[ichT])) {
      ++ichT;
    }
    cchToken = ichT - ichToken;
    return cchToken > 0;
  }

  bool SelectionCompiler::lookingAtEndOfLine() {
    if (ichToken == cchScript)
      return true;
    int ichT = ichToken;
    char ch = script[ichT];
    if (ch == '\r') {
      ++ichT;
      if (ichT < cchScript && script[ichT] == '\n')
	++ichT;
    } else if (ch == '\n') {
      ++ichT;
    } else {
      return false;
    }
    cchToken = ichT - ichToken;
    return true;
  }

  bool SelectionCompiler::lookingAtEndOfStatement() {
    if (ichToken == cchScript || script[ichToken] != ';')
      return false;
    cchToken = 1;
    return true;
  }

  bool SelectionCompiler::lookingAtString() {
    if (ichToken == cchScript)
      return false;
    if (script[ichToken] != '"')
      return false;
    // remove support for single quote
    // in order to use it in atom expressions
    //    char chFirst = script.charAt(ichToken);
    //    if (chFirst != '"' && chFirst != '\'')
    //      return false;
    int ichT = ichToken + 1;
    //    while (ichT < cchScript && script.charAt(ichT++) != chFirst)
    char ch;
    bool previousCharBackslash = false;
    while (ichT < cchScript) {
      ch = script[ichT++];
      if (ch == '"' && !previousCharBackslash)
        break;
      previousCharBackslash = ch == '\\' ? !previousCharBackslash : false;
    }
    cchToken = ichT - ichToken;

    return true;
  }

  
  std::string SelectionCompiler::getUnescapedStringLiteral() {
    /** @todo */
    std::string sb(cchToken - 2, ' ');
    
    int ichMax = ichToken + cchToken - 1;
    int ich = ichToken + 1;

    while (ich < ichMax) {
      char ch = script[ich++];
      if (ch == '\\' && ich < ichMax) {
	ch = script[ich++];
	switch (ch) {
	case 'b':
	  ch = '\b';
	  break;
	case 'n':
	  ch = '\n';
	  break;
	case 't':
	  ch = '\t';
	  break;
	case 'r':
	  ch = '\r';
	  // fall into
	case '"':
	case '\\':
	case '\'':
	  break;
	case 'x':
	case 'u':
	  int digitCount = ch == 'x' ? 2 : 4;
	  if (ich < ichMax) {
	    int unicode = 0;
	    for (int k = digitCount; --k >= 0 && ich < ichMax; ) {
	      char chT = script[ich];
	      int hexit = getHexitValue(chT);
	      if (hexit < 0)
		break;
	      unicode <<= 4;
	      unicode += hexit;
	      ++ich;
	    }
	    ch = (char)unicode;
	  }
	}
      }
      sb.append(1, ch);
    }

    return sb;
  }

  int SelectionCompiler::getHexitValue(char ch) {
    if (ch >= '0' && ch <= '9')
      return ch - '0';
    else if (ch >= 'a' && ch <= 'f')
      return 10 + ch - 'a';
    else if (ch >= 'A' && ch <= 'F')
      return 10 + ch - 'A';
    else
      return -1;
  }

  bool SelectionCompiler::lookingAtSpecialString() {
    int ichT = ichToken;
    char ch = script[ichT];
    while (ichT < cchScript && ch != ';' && ch != '\r' && ch != '\n') {
      ++ichT;
    }
    cchToken = ichT - ichToken;
    return cchToken > 0;
  }

  bool SelectionCompiler::lookingAtDecimal(bool allowNegative) {
    if (ichToken == cchScript) {
      return false;
    }
    
    int ichT = ichToken;
    if (script[ichT] == '-') {
      ++ichT;
    }
    bool digitSeen = false;
    char ch = 'X';
    while (ichT < cchScript && std::isdigit(ch = script[ichT])) {
      ++ichT;
      digitSeen = true;
    }

    if (ichT == cchScript || ch != '.') {
      return false;
    }

    // to support DMPC.1, let's check the character before the dot
    if (ch == '.' && (ichT > 0) && std::isalpha(script[ichT - 1])) {
      return false;
    }

    ++ichT;
    while (ichT < cchScript && std::isdigit(script[ichT])) {
      ++ichT;
      digitSeen = true;
    }
    cchToken = ichT - ichToken;
    return digitSeen;
  }

  bool SelectionCompiler::lookingAtInteger(bool allowNegative) {
    if (ichToken == cchScript) {
      return false;
    }
    int ichT = ichToken;
    if (allowNegative && script[ichToken] == '-') {
      ++ichT;
    }
    int ichBeginDigits = ichT;
    while (ichT < cchScript && std::isdigit(script[ichT])) {
      ++ichT;
    }
    if (ichBeginDigits == ichT) {
      return false;
    }
    cchToken = ichT - ichToken;
    return true;
  }

  bool SelectionCompiler::lookingAtLookupToken() {
    if (ichToken == cchScript) {
      return false;
    }

    int ichT = ichToken;
    char ch;
    switch (ch = script[ichT++]) {
    case '(':
    case ')':
    case ',':
    case '[':
    case ']':
      break;
    case '&':
    case '|':
      if (ichT < cchScript && script[ichT] == ch) {
	++ichT;
      }
      break;
    case '<':
    case '=':
    case '>':
      if (ichT < cchScript && ((ch = script[ichT]) == '<' || ch == '=' || ch == '>')) {
	++ichT;
      }
      break;
    case '/':
    case '!':
      if (ichT < cchScript && script[ichT] == '=') {
	++ichT;
      }
      break;
    default:
      if ((ch < 'a' || ch > 'z') && (ch < 'A' && ch > 'Z') && ch != '_') {
	return false;
      }
    case '*':
    case '?': // include question marks in identifier for atom expressions
      while (ichT < cchScript && !std::isspace(ch = script[ichT]) && 
	     (std::isalpha(ch) ||std::isdigit(ch) || ch == '_' || ch == '.' || ch == '*' || ch == '?' || ch == '+' || ch == '-' || ch == '[' || ch == ']') ){

	++ichT;
      }
      break;
    }

    cchToken = ichT - ichToken;

    return true;
  }

  bool SelectionCompiler::compileCommand(const std::vector<Token>& ltoken) {
    const Token& tokenCommand = ltoken[0];
    int tokCommand = tokenCommand.tok;

    atokenCommand = ltoken;
    if ((tokCommand & Token::expressionCommand) != 0 && !compileExpression()) {
      return false;
    }
    
    return true;
  }

  bool SelectionCompiler::compileExpression() {
    /** todo */
    unsigned int i = 1;
    int tokCommand = atokenCommand[0].tok;
    if (tokCommand == Token::define) {
      i = 2;
    } else if ((tokCommand & Token::embeddedExpression) != 0) {
      // look for the open parenthesis
      while (i < atokenCommand.size() &&
	     atokenCommand[i].tok != Token::leftparen)
        ++i;
    }

    if (i >= atokenCommand.size()) {
      return true;
    }
    return compileExpression(i);
  }

                  
  bool SelectionCompiler::addTokenToPostfix(const Token& token) {
    ltokenPostfix.push_back(token);
    return true;
  }

  bool SelectionCompiler::compileExpression(int itoken) {
    ltokenPostfix.clear();
    for (int i = 0; i < itoken; ++i) {
      addTokenToPostfix(atokenCommand[i]);
    }
    
    atokenInfix = atokenCommand;
    itokenInfix = itoken;

    addTokenToPostfix(Token::tokenExpressionBegin);
    if (!clauseOr()) {
      return false;
    }
    
    addTokenToPostfix(Token::tokenExpressionEnd);
    if (itokenInfix != atokenInfix.size()) {
      return endOfExpressionExpected();
    }

    atokenCommand = ltokenPostfix;
    return true;
  }

  Token SelectionCompiler::tokenNext() {
    if (itokenInfix == atokenInfix.size()) {
      return Token();
    }
    return atokenInfix[itokenInfix++];
  }

  boost::any SelectionCompiler::valuePeek() {
    if (itokenInfix == atokenInfix.size()) {
      return boost::any();
    } else {
      return atokenInfix[itokenInfix].value;
    }
  }

  int SelectionCompiler::tokPeek() {
    if (itokenInfix == atokenInfix.size()) {
      return 0;
    }else {
      return atokenInfix[itokenInfix].tok;
    }
  }

  bool SelectionCompiler::clauseOr() {
    if (!clauseAnd()) {
      return false;
    }
    
    while (tokPeek() == Token::opOr) {
      Token tokenOr = tokenNext();
      if (!clauseAnd()) {
	return false;
      }
      addTokenToPostfix(tokenOr);
    }
    return true;
  }

  bool SelectionCompiler::clauseAnd() {
    if (!clauseNot()) {
      return false;
    }

    while (tokPeek() == Token::opAnd) {
      Token tokenAnd = tokenNext();
      if (!clauseNot()) {
	return false;
      }
      addTokenToPostfix(tokenAnd);
    }
    return true;
  }

  bool SelectionCompiler::clauseNot() {
    if (tokPeek() == Token::opNot) {
      Token tokenNot = tokenNext();
      if (!clauseNot()) {
	return false;
      }
      return addTokenToPostfix(tokenNot);
    }
    return clausePrimitive();
  }

  bool SelectionCompiler::clausePrimitive() {
    int tok = tokPeek();
    switch (tok) {
    case Token::within:
      return clauseWithin();

    case Token::asterisk:
    case Token::identifier:
      return clauseChemObjName();

    case Token::integer :
      return clauseIndex();
    default:
      if ((tok & Token::atomproperty) == Token::atomproperty) {
	return clauseComparator();
      }
      if ((tok & Token::predefinedset) != Token::predefinedset) {
	break;
      }
      // fall into the code and below and just add the token
    case Token::all:
    case Token::none:
    case Token::hull:
      return addTokenToPostfix(tokenNext());
    case Token::leftparen:
      tokenNext();
      if (!clauseOr()) {
	return false;
      }
      if (tokenNext().tok != Token::rightparen) {
	return rightParenthesisExpected();
      }
      return true;
    }
    return unrecognizedExpressionToken();
  }

  bool SelectionCompiler::clauseComparator() {
    Token tokenAtomProperty = tokenNext();
    Token tokenComparator = tokenNext();
    if ((tokenComparator.tok & Token::comparator) == 0) {
      return comparisonOperatorExpected();
    }

    Token tokenValue = tokenNext();
    if (tokenValue.tok != Token::integer && tokenValue.tok != Token::decimal) {
      return numberExpected();
    }
    
    float val;
    if (tokenValue.value.type() == typeid(int)) {
      val = boost::any_cast<int>(tokenValue.value);
    } else if (tokenValue.value.type() == typeid(float)) {
      val = boost::any_cast<float>(tokenValue.value);
    } else {
      return false;
    }

    boost::any floatVal;
    floatVal = val;
    return addTokenToPostfix(Token(tokenComparator.tok,
				   tokenAtomProperty.tok, floatVal));
  }

  bool SelectionCompiler::clauseWithin() {
    tokenNext();                             // WITHIN
    if (tokenNext().tok != Token::leftparen) {  // (
      return leftParenthesisExpected();
    }
    
    boost::any distance;
    Token tokenDistance = tokenNext();       // distance
    switch(tokenDistance.tok) {
    case Token::integer:
    case Token::decimal:
      distance = tokenDistance.value;
      break;
    default:
      return numberOrKeywordExpected();
    }

    if (tokenNext().tok != Token::opOr) {       // ,
      return commaExpected();
    }
    
    if (! clauseOr()) {                        // *expression*
      return false;
    }
    
    if (tokenNext().tok != Token::rightparen) { // )T
      return rightParenthesisExpected();
    }
    
    return addTokenToPostfix(Token(Token::within, distance));
  }

  bool SelectionCompiler::clauseChemObjName() {
    Token token = tokenNext();
    if (token.tok == Token::identifier && token.value.type() == typeid(std::string)) {

      std::string name = boost::any_cast<std::string>(token.value);
      if (isNameValid(name)) {
	return addTokenToPostfix(Token(Token::name, name));
      } else {
	return compileError("invalid name: " + name);
      }
    } 

    return false;
        
  }

  bool SelectionCompiler::isNameValid(const std::string& name) {
    int nbracket = 0;
    int ndot = 0;
    for (unsigned int i = 0 ; i < name.size(); ++i) {
      switch(name[i]) {

      case '[' :
	++nbracket;
	break;
      case ']' :
	--nbracket;
	break;
      case '.' :
	++ndot;
	break;       
      }
    }

    //only allow 3 dots at most
    return (ndot <=3 && nbracket == 0) ? true : false;
  }

  bool SelectionCompiler::clauseIndex(){
    Token token = tokenNext();
    if (token.tok == Token::integer) {
      int index = boost::any_cast<int>(token.value);
      int tok = tokPeek();
      std::cout << "Token::to is " << Token::to << ", tok = " << tok << std::endl;
      if (tok == Token::to) {
	tokenNext();
	tok = tokPeek();
	if (tok != Token::integer) {
	  return numberExpected();
	}
            
	boost::any intVal = tokenNext().value;
	int first = index;
	if (intVal.type() != typeid(int)){
	  return false;
	}
	int second = boost::any_cast<int>(intVal);

	return addTokenToPostfix(Token(Token::index, boost::any(std::make_pair(first, second))));
            
      }else {
	return addTokenToPostfix(Token(Token::index, boost::any(index)));
      }
    } else {
      return numberExpected();
    }
  }

}
