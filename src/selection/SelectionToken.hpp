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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#ifndef SELECTION_TOKEN_HPP
#define SELECTION_TOKEN_HPP

#include <map>
#include <string>
#include <utility>
#include "utils/any.hpp"

namespace OpenMD {


  /**
   * @class Token
   * @todo document
   * @note translate from jmol
   */
  class Token {

  public:


    int tok;
    boost::any value;
    int intValue;

    Token() { tok = unknown;}

    Token(const Token& token) {
      *this = token;
    }
        
    Token(int MyTok, int myIntValue, const boost::any& myValue) {
      this->tok = MyTok;
      this->intValue = myIntValue;
      this->value = myValue;
    }

    Token(int myTok, int myIntValue) {
      this->tok = myTok;
      this->intValue = myIntValue;
    }

    Token(int myTok) {
      this->tok = myTok;
    }

    Token(int myTok,  const boost::any& myValue) {
      this->tok = myTok;
      this->value = myValue;
    }
        
    const static int nada              =  0;
    const static int identifier        =  1;
    const static int integer           =  2;
    const static int decimal           =  3;
    const static int string            =  4;
    const static int unknown           =  5;
    const static int keyword           =  6;
    const static int whitespace        =  7;
    const static int comment           =  8;
    const static int endofline         = 9;
    const static int endofstatement    = 10;

    const static int command           = (1 <<  8);
    const static int expressionCommand = (1 <<  9); // expression command
    const static int expression        = (1 << 10); /// expression term

    // generally, the minus sign is used to denote atom ranges
    // this property is used for the few commands which allow negative integers
    const static int negnums      = (1 << 11);

    //expression involves coordinates which will change every frame, such as withins
    const static int dynamic        = (1 << 12);

    // every property is also valid in an expression context
    const static int atomproperty      = (1 << 12) | expression | negnums;
    // every predefined is also valid in an expression context
    const static int comparator        = (1 << 13) | expression;
    const static int predefinedset     = (1 << 14) | expression;
    const  static int embeddedExpression= (1 << 15); // embedded expression
    const static int index = (1 << 16) | expression;
    // rasmol commands
    const static int define       = command | expressionCommand | 1;
    const static int select       = command | expressionCommand | 2 ;

    //predefine
    //const static int selected    = predefinedset |0;

    // atom expression operators
    const static int leftparen    = expression |  0;
    const static int rightparen   = expression |  1;
    const static int to           = expression |  2;
    const static int opAnd        = expression |  3;
    const static int opOr         = expression |  4;
    const static int opNot        = expression |  5;
    const static int within       = expression | dynamic | 6;
    const static int asterisk     = expression |  7;
    const static int dot          = expression | 8;
    const static int all          = expression | 9 ; 
    const static int none         = expression | 10;
    const static int name         = expression | 11;
    const static int hull         = expression | dynamic | 12;

    // miguel 2005 01 01
    // these are used to demark the beginning and end of expressions
    // they do not exist in the source code, but are emitted by the
    // expression compiler
    const static int expressionBegin = expression | 100;
    const static int expressionEnd   = expression | 101;

    const static int mass         = atomproperty | 0;
    const static int charge       = atomproperty | dynamic | 1;
    const static int x            = atomproperty | dynamic | 2;
    const static int y            = atomproperty | dynamic | 3;
    const static int z            = atomproperty | dynamic | 4;
    const static int r            = atomproperty | dynamic | 5;
    const static int wrappedX     = atomproperty | dynamic | 6;
    const static int wrappedY     = atomproperty | dynamic | 7;
    const static int wrappedZ     = atomproperty | dynamic | 8;
        
    const static int opGT         = comparator |  dynamic | 0;
    const static int opGE         = comparator |  dynamic | 1;
    const static int opLE         = comparator |  dynamic | 2;
    const static int opLT         = comparator |  dynamic | 3;
    const static int opEQ         = comparator |  dynamic | 4;
    const static int opNE         = comparator |  dynamic | 5;
 
    static Token tokenExpressionBegin;
    static Token tokenExpressionEnd;

  };


}

#endif 
