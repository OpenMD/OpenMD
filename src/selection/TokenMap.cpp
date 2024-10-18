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

#include "selection/TokenMap.hpp"

namespace OpenMD {

  TokenMap::TokenMap() {
    tokenMap_.insert(TokenMapType::value_type(
        "define", Token(Token::define, std::string("define"))));
    tokenMap_.insert(TokenMapType::value_type(
        "select", Token(Token::select, std::string("select"))));
    // tokenMap_.insert(TokenMapType::value_type("selected",
    // Token(Token::selected, std::string("selected"))));

    // expressions
    tokenMap_.insert(TokenMapType::value_type(
        "(", Token(Token::leftparen, std::string("("))));
    tokenMap_.insert(TokenMapType::value_type(
        ")", Token(Token::rightparen, std::string(")"))));

    tokenMap_.insert(TokenMapType::value_type(
        "and", Token(Token::opAnd, std::string("and"))));
    tokenMap_.insert(
        TokenMapType::value_type("&", Token(Token::opAnd, std::string("and"))));
    tokenMap_.insert(TokenMapType::value_type(
        "&&", Token(Token::opAnd, std::string("and"))));

    tokenMap_.insert(
        TokenMapType::value_type("or", Token(Token::opOr, std::string("or"))));
    tokenMap_.insert(
        TokenMapType::value_type(",", Token(Token::opOr, std::string("or"))));
    tokenMap_.insert(
        TokenMapType::value_type("|", Token(Token::opOr, std::string("or"))));
    tokenMap_.insert(
        TokenMapType::value_type("||", Token(Token::opOr, std::string("or"))));

    tokenMap_.insert(TokenMapType::value_type(
        "not", Token(Token::opNot, std::string("not"))));
    tokenMap_.insert(
        TokenMapType::value_type("!", Token(Token::opNot, std::string("not"))));

    tokenMap_.insert(
        TokenMapType::value_type("<", Token(Token::opLT, std::string("<"))));
    tokenMap_.insert(
        TokenMapType::value_type("<=", Token(Token::opLE, std::string("<="))));
    tokenMap_.insert(
        TokenMapType::value_type(">=", Token(Token::opGE, std::string(">="))));
    tokenMap_.insert(
        TokenMapType::value_type(">", Token(Token::opGT, std::string(">"))));
    tokenMap_.insert(
        TokenMapType::value_type("==", Token(Token::opEQ, std::string("=="))));
    tokenMap_.insert(
        TokenMapType::value_type("!=", Token(Token::opNE, std::string("!="))));
    tokenMap_.insert(TokenMapType::value_type(
        "within", Token(Token::within, std::string("within"))));
    tokenMap_.insert(
        TokenMapType::value_type(".", Token(Token::dot, std::string("."))));
    tokenMap_.insert(TokenMapType::value_type(
        "mass", Token(Token::mass, std::string("mass"))));
    tokenMap_.insert(TokenMapType::value_type(
        "charge", Token(Token::charge, std::string("charge"))));
    tokenMap_.insert(TokenMapType::value_type(
        "hull", Token(Token::hull, std::string("hull"))));
    tokenMap_.insert(TokenMapType::value_type(
        "alphahull", Token(Token::alphahull, std::string("alphahull"))));
    tokenMap_.insert(
        TokenMapType::value_type("x", Token(Token::x, std::string("x"))));
    tokenMap_.insert(
        TokenMapType::value_type("y", Token(Token::y, std::string("y"))));
    tokenMap_.insert(
        TokenMapType::value_type("z", Token(Token::z, std::string("z"))));
    tokenMap_.insert(TokenMapType::value_type(
        "wrappedx", Token(Token::wrappedX, std::string("wrappedx"))));
    tokenMap_.insert(TokenMapType::value_type(
        "wrappedy", Token(Token::wrappedY, std::string("wrappedy"))));
    tokenMap_.insert(TokenMapType::value_type(
        "wrappedz", Token(Token::wrappedZ, std::string("wrappedz"))));
    tokenMap_.insert(TokenMapType::value_type(
        "atomno", Token(Token::atomno, std::string("atomno"))));
    tokenMap_.insert(
        TokenMapType::value_type("r", Token(Token::r, std::string("r"))));
    tokenMap_.insert(
        TokenMapType::value_type("to", Token(Token::to, std::string("to"))));

    tokenMap_.insert(
        TokenMapType::value_type("all", Token(Token::all, std::string("all"))));
    tokenMap_.insert(TokenMapType::value_type(
        "none", Token(Token::none, std::string("none"))));
  }

  Token* TokenMap::getToken(const std::string& ident) {
    std::map<std::string, Token>::iterator i = tokenMap_.find(ident);

    return i != tokenMap_.end() ? &(i->second) : NULL;
  }
}  // namespace OpenMD
