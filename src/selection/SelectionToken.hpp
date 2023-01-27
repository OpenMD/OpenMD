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

#ifndef SELECTION_TOKEN_HPP
#define SELECTION_TOKEN_HPP

#include <any>

namespace OpenMD {

  /**
   * @class Token
   * @todo document
   * @note translate from jmol
   */
  class Token {
  public:
    int tok {unknown};
    int intValue {};
    std::any value {};

    Token() = default;
    explicit Token(int myTok) : tok {myTok} {}
    Token(int myTok, int myIntValue) : tok {myTok}, intValue {myIntValue} {}
    Token(int myTok, const std::any& myValue) : tok {myTok}, value {myValue} {}
    Token(int MyTok, int myIntValue, const std::any& myValue) :
        tok {MyTok}, intValue {myIntValue}, value {myValue} {}

    static constexpr int nada           = 0;
    static constexpr int identifier     = 1;
    static constexpr int integer        = 2;
    static constexpr int decimal        = 3;
    static constexpr int string         = 4;
    static constexpr int unknown        = 5;
    static constexpr int keyword        = 6;
    static constexpr int whitespace     = 7;
    static constexpr int comment        = 8;
    static constexpr int endofline      = 9;
    static constexpr int endofstatement = 10;

    static constexpr int command           = (1 << 8);
    static constexpr int expressionCommand = (1 << 9);   // expression command
    static constexpr int expression        = (1 << 10);  /// expression term

    // generally, the minus sign is used to denote atom ranges
    // this property is used for the few commands which allow negative integers
    static constexpr int negnums = (1 << 11);

    // expression involves coordinates which will change every frame, such as
    // withins
    static constexpr int dynamic = (1 << 12);

    // every property is also valid in an expression context
    static constexpr int atomproperty = (1 << 12) | expression | negnums;
    // every predefined is also valid in an expression context
    static constexpr int comparator         = (1 << 13) | expression;
    static constexpr int predefinedset      = (1 << 14) | expression;
    static constexpr int embeddedExpression = (1 << 15);  // embedded expression
    static constexpr int index              = (1 << 16) | expression;
    // rasmol commands
    static constexpr int define = command | expressionCommand | 1;
    static constexpr int select = command | expressionCommand | 2;

    // predefine
    // static constexpr int selected    = predefinedset |0;

    // atom expression operators
    static constexpr int leftparen  = expression | 0;
    static constexpr int rightparen = expression | 1;
    static constexpr int to         = expression | 2;
    static constexpr int opAnd      = expression | 3;
    static constexpr int opOr       = expression | 4;
    static constexpr int opNot      = expression | 5;
    static constexpr int within     = expression | dynamic | 6;
    static constexpr int asterisk   = expression | 7;
    static constexpr int dot        = expression | 8;
    static constexpr int all        = expression | 9;
    static constexpr int none       = expression | 10;
    static constexpr int name       = expression | 11;
    static constexpr int hull       = expression | dynamic | 12;
    static constexpr int alphahull  = expression | dynamic | 13;

    // miguel 2005 01 01
    // these are used to demark the beginning and end of expressions
    // they do not exist in the source code, but are emitted by the
    // expression compiler
    static constexpr int expressionBegin = expression | 100;
    static constexpr int expressionEnd   = expression | 101;

    static constexpr int mass     = atomproperty | 0;
    static constexpr int charge   = atomproperty | dynamic | 1;
    static constexpr int x        = atomproperty | dynamic | 2;
    static constexpr int y        = atomproperty | dynamic | 3;
    static constexpr int z        = atomproperty | dynamic | 4;
    static constexpr int r        = atomproperty | dynamic | 5;
    static constexpr int wrappedX = atomproperty | dynamic | 6;
    static constexpr int wrappedY = atomproperty | dynamic | 7;
    static constexpr int wrappedZ = atomproperty | dynamic | 8;
    static constexpr int atomno   = atomproperty | 9;

    static constexpr int opGT = comparator | dynamic | 0;
    static constexpr int opGE = comparator | dynamic | 1;
    static constexpr int opLE = comparator | dynamic | 2;
    static constexpr int opLT = comparator | dynamic | 3;
    static constexpr int opEQ = comparator | dynamic | 4;
    static constexpr int opNE = comparator | dynamic | 5;

    static Token tokenExpressionBegin;
    static Token tokenExpressionEnd;
  };
}  // namespace OpenMD

#endif
