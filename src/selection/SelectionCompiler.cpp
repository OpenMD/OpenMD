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

#include "selection/SelectionCompiler.hpp"
namespace oopse {

bool SelectionCompiler::compile(const std::string& filename, const std::string& script) {

    this->filename = filename;
    this->script = script;
    lineNumbers.clear();
    lineIndices.clear();
    aatokenCompiled.clear();
        
    if (internalcompile()) {
        return true;
    }
    
    int icharEnd;
    if ((icharEnd = script.find('\r', ichCurrentCommand)) == std::string:npos &&
        (icharEnd = script.find('\n', ichCurrentCommand)) == std::string:npos) {
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

    std::vector<Token> lltoken;
    std::vector<Token> ltoken;

    //Token tokenCommand = null;
    int tokCommand = Token.nada;

    for ( ; true; ichToken += cchToken) {
        if (lookingAtLeadingWhitespace())
            continue;
        if (lookingAtComment())
            continue;
        boolean endOfLine = lookingAtEndOfLine();
        if (endOfLine || lookingAtEndOfStatement()) {
            if (tokCommand != Token.nada) {
                if (! compileCommand(ltoken)) {
                    return false;
                }
            lltoken.push_back(atokenCommand);
            /** @todo*/
            int iCommand = lltoken.size();
            lineNumbers[iCommand] = lineCurrent;
            lineIndices[iCommand] = (short) ichCurrentCommand;
            ltoken.clear();
            tokCommand = Token.nada;
            }
            
            if (ichToken < cchScript) {
                if (endOfLine)
                    ++lineCurrent;
              continue;
            }
            break;
        }

        if (tokCommand != Token.nada) {
            if (lookingAtString()) {
                std::string str = getUnescapedStringLiteral();
                ltoken.push_back(Token(Token.string, str));
                continue;
            }
            if ((tokCommand & Token.specialstring) != 0 &&
                lookingAtSpecialString()) {
                std::string str = script.substr(ichToken, ichToken + cchToken);
                ltoken.push_back(Token(Token.string, str));
                continue;
            }
            if (lookingAtDecimal((tokCommand & Token.negnums) != 0)) {
                float value = lexi_cast<float>((script.substr(ichToken, ichToken + cchToken));          
                ltoken.push_back(Token(Token.decimal, new Float(value)));/**@todo*/
                continue;
            }
            if (lookingAtInteger((tokCommand & Token.negnums) != 0)) {
                std::string intString = script.substr(ichToken, ichToken + cchToken);
                int val = lexi_cast<int>(intString);
                ltoken.push_back(new Token(Token.integer, val, intString));/**@todo*/
                continue;
            }
        }
      
        if (lookingAtLookupToken()) {
            std::string ident = script.subst(ichToken, ichToken + cchToken);

            /**@todo */
            Token token = (Token) Token.map.get(ident);
            if (token == NULL) {

            }
            Token token(Token.identifier, ident);

            
            int tok = token.tok;
            
            switch (tokCommand) {
                case Token.nada:
                    ichCurrentCommand = ichToken;
                    //tokenCommand = token;
                    tokCommand = tok;
                    if ((tokCommand & Token.command) == 0)
                    return commandExpected();
                    break;

                case Token.define:
                    if (ltoken.size() == 1) {
                        // we are looking at the variable name
                        if (tok != Token.identifier &&
                        (tok & Token.predefinedset) != Token.predefinedset)
                        return invalidExpressionToken(ident);
                    } else {
                    // we are looking at the expression
                    if (tok != Token.identifier && tok != Token.set &&
                        (tok & (Token.expression | Token.predefinedset)) == 0)
                        return invalidExpressionToken(ident);
                    }
                    
                    break;

                case Token.select:
                    if (tok != Token.identifier && (tok & Token.expression) == 0)
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

    aatokenCompiled.push_back(lltoken);
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
    boolean previousCharBackslash = false;
    while (ichT < cchScript) {
      ch = script.[ichT++];
      if (ch == '"' && !previousCharBackslash)
        break;
      previousCharBackslash = ch == '\\' ? !previousCharBackslash : false;
    }
    cchToken = ichT - ichToken;
    return true;
  }

  
std::string SelectionCompiler::getUnescapedStringLiteral() {
    StringBuffer sb = new StringBuffer(cchToken - 2);
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
        sb.append(ch);
    }

    return "" + sb;
}

static int SelectionCompiler::getHexitValue(char ch) {
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

bool SelectionCompiler::lookingAtDecimal(boolean allowNegative) {
    if (ichToken == cchScript) {
        return false;
    }
    
    int ichT = ichToken;
    if (script[ichT] == '-') {
        ++ichT;
    }
    boolean digitSeen = false;
    char ch = 'X';
    while (ichT < cchScript && std::isdigit(ch = script[ichT])) {
        ++ichT;
        digitSeen = true;
    }

    if (ichT == cchScript || ch != '.') {
        return false;
    }

    // to support 1.ca, let's check the character after the dot
    // to determine if it is an alpha
    if (ch == '.' && (ichT + 1 < cchScript) && std::isalpha(script[ichT + 1])) {
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

bool SelectionCompiler::lookingAtInteger(boolean allowNegative) {
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
        case '*':
        case '-':
        case '[':
        case ']':
        case '+':
        case ':':
        case '@':
        case '.':
        case '%':
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
        case '?': // include question marks in identifier for atom expressions
            while (ichT < cchScript && (std::isalpha(ch = script[ichT]) ||std::isdigit(ch) ||
                ch == '_' || ch == '?') ||(ch == '^' && ichT > ichToken && std::isdigit(script[ichT - 1]))) {
                // hack for insertion codes embedded in an atom expression :-(
                // select c3^a
                ++ichT;
            }
        break;
    }
    cchToken = ichT - ichToken;
    return true;
}

bool SelectionCompiler::compileCommand(Vector ltoken) {
    /** todo */
    Token tokenCommand = (Token)ltoken.firstElement();
    int tokCommand = tokenCommand.tok;
        
    atokenCommand = new Token[ltoken.size()];
    ltoken.copyInto(atokenCommand);
    if ((tokCommand & Token.expressionCommand) != 0 && !compileExpression()) {
        return false;
    }
    return true;
}

bool SelectionCompiler::compileExpression() {
    /** todo */
    int i = 1;
    int tokCommand = atokenCommand[0].tok;
    if (tokCommand == Token.define)
      i = 2;
    else if ((tokCommand & Token.embeddedExpression) != 0) {
      // look for the open parenthesis
      while (i < atokenCommand.length &&
             atokenCommand[i].tok != Token.leftparen)
        ++i;
    }
    if (i >= atokenCommand.length)
      return true;
    return compileExpression(i);
  }

                  
bool SelectionCompiler::addTokenToPostfix(Token token) {
    ltokenPostfix.push_back(token);
    return true;
}

bool SelectionCompiler::compileExpression(int itoken) {
    ltokenPostfix = new Vector();
    for (int i = 0; i < itoken; ++i)
        addTokenToPostfix(atokenCommand[i]);

    atokenInfix = atokenCommand;
    itokenInfix = itoken;

    addTokenToPostfix(Token.tokenExpressionBegin);
    if (!clauseOr()) {
        return false;
    }
    
    addTokenToPostfix(Token.tokenExpressionEnd);
    if (itokenInfix != atokenInfix.length) {
        return endOfExpressionExpected();
    }

    atokenCommand = ltokenPostfix;
    return true;
}

Token SelectionCompiler::tokenNext() {
if (itokenInfix == atokenInfix.length)
return null;
return atokenInfix[itokenInfix++];
}

Object SelectionCompiler::valuePeek() {
    if (itokenInfix == atokenInfix.length) {
        return null;
    } else {
        return atokenInfix[itokenInfix].value;
    }
}

int SelectionCompiler::tokPeek() {
    if (itokenInfix == atokenInfix.length) {
        return 0;
    }else {
        return atokenInfix[itokenInfix].tok;
    }
}

bool SelectionCompiler::clauseOr() {
    if (!clauseAnd()) {
        return false;
    }
    
    while (tokPeek() == Token.opOr) {
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

    while (tokPeek() == Token.opAnd) {
        Token tokenAnd = tokenNext();
        if (!clauseNot()) {
            return false;
        }
        addTokenToPostfix(tokenAnd);
    }
    return true;
}

bool SelectionCompiler::clauseNot() {
    if (tokPeek() == Token.opNot) {
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
        case Token.within:
            return clauseWithin();
        case Token.hyphen: // selecting a negative residue spec
        case Token.integer:
        case Token.seqcode:
        case Token.asterisk:
        case Token.leftsquare:
        case Token.identifier:
        case Token.x:
        case Token.y:
        case Token.z:
        case Token.colon:
            return clauseResidueSpec();
        default:
            if ((tok & Token.atomproperty) == Token.atomproperty) {
                return clauseComparator();
            }
            if ((tok & Token.predefinedset) != Token.predefinedset) {
                break;
            }
            // fall into the code and below and just add the token
        case Token.all:
        case Token.none:
            return addTokenToPostfix(tokenNext());
        case Token.leftparen:
            tokenNext();
            if (!clauseOr()) {
                return false;
            }
            if (tokenNext().tok != Token.rightparen) {
                return rightParenthesisExpected();
            }
            return true;
    }
    return unrecognizedExpressionToken();
}

bool SelectionCompiler::clauseComparator() {
    Token tokenAtomProperty = tokenNext();
    Token tokenComparator = tokenNext();
    if ((tokenComparator.tok & Token.comparator) == 0) {
        return comparisonOperatorExpected();
    }

    Token tokenValue = tokenNext();
    if (tokenValue.tok != Token.integer) {
        return integerExpected();
    }
    int val = tokenValue.intValue;
    // note that a comparator instruction is a complicated instruction
    // int intValue is the tok of the property you are comparing
    // the value against which you are comparing is stored as an Integer
    // in the object value
    return addTokenToPostfix(new Token(tokenComparator.tok,
                       tokenAtomProperty.tok,
                       new Integer(val)));
}

bool SelectionCompiler::clauseWithin() {
    tokenNext();                             // WITHIN
    if (tokenNext().tok != Token.leftparen) {  // (
        return leftParenthesisExpected();
    }
    
    Object distance;
    Token tokenDistance = tokenNext();       // distance
    switch(tokenDistance.tok) {
        case Token.integer:
            distance = new Float((tokenDistance.intValue * 4) / 1000f);
            break;
        case Token.decimal:
            distance = tokenDistance.value;
            break;
        default:
            return numberOrKeywordExpected();
    }

    if (tokenNext().tok != Token.opOr) {       // ,
        return commaExpected();
    }
    
    if (! clauseOr()) {                        // *expression*
        return false;
    }
    
    if (tokenNext().tok != Token.rightparen) { // )T
        return rightParenthesisExpected();
    }
    
    return addTokenToPostfix(new Token(Token.within, distance));
}

bool SelectionCompiler:: clauseChemObject() {
}

bool SelectionCompiler:: clauseMolecule() {
}

bool SelectionCompiler:: clauseMolName() {
}

bool SelectionCompiler:: clauseMolIndex() {
}

bool SelectionCompiler:: clauseName() {
}

bool SelectionCompiler:: clauseIndex() {
}

bool SelectionCompiler:: clauseStuntDoubleName() {
}

bool SelectionCompiler:: clauseStuntDoubleIndex() {
}

}
