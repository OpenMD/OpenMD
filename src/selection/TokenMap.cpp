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

#include "selection/TokenMap.hpp"

namespace OpenMD {

  TokenMap* TokenMap::instance_ = NULL;

  TokenMap::TokenMap() {
    tokenMap_.insert(TokenMapType::value_type("define", Token(Token::define, std::string("define"))));
    tokenMap_.insert(TokenMapType::value_type("select", Token(Token::select, std::string("select")))); 
    //tokenMap_.insert(TokenMapType::value_type("selected", Token(Token::selected, std::string("selected")))); 
    
    //expressions
    tokenMap_.insert(TokenMapType::value_type("(", Token(Token::leftparen, std::string("("))));
    tokenMap_.insert(TokenMapType::value_type(")", Token(Token::rightparen, std::string(")"))));

    tokenMap_.insert(TokenMapType::value_type("and", Token(Token::opAnd, std::string("and"))));
    tokenMap_.insert(TokenMapType::value_type("&", Token(Token::opAnd, std::string("and"))));
    tokenMap_.insert(TokenMapType::value_type("&&", Token(Token::opAnd, std::string("and"))));

    tokenMap_.insert(TokenMapType::value_type("or", Token(Token::opOr, std::string("or"))));
    tokenMap_.insert(TokenMapType::value_type(",", Token(Token::opOr, std::string("or"))));
    tokenMap_.insert(TokenMapType::value_type("|", Token(Token::opOr, std::string("or"))));
    tokenMap_.insert(TokenMapType::value_type("||", Token(Token::opOr, std::string("or"))));

    tokenMap_.insert(TokenMapType::value_type("not", Token(Token::opNot, std::string("not"))));
    tokenMap_.insert(TokenMapType::value_type("!", Token(Token::opNot, std::string("not"))));
    
    tokenMap_.insert(TokenMapType::value_type("<", Token(Token::opLT, std::string("<"))));
    tokenMap_.insert(TokenMapType::value_type("<=", Token(Token::opLE, std::string("<="))));
    tokenMap_.insert(TokenMapType::value_type(">=", Token(Token::opGE, std::string(">="))));
    tokenMap_.insert(TokenMapType::value_type(">", Token(Token::opGT, std::string(">"))));
    tokenMap_.insert(TokenMapType::value_type("==", Token(Token::opEQ, std::string("=="))));
    tokenMap_.insert(TokenMapType::value_type("!=", Token(Token::opNE, std::string("!="))));
    tokenMap_.insert(TokenMapType::value_type("within", Token(Token::within, std::string("within"))));
    tokenMap_.insert(TokenMapType::value_type(".", Token(Token::dot, std::string("."))));
    tokenMap_.insert(TokenMapType::value_type("mass", Token(Token::mass, std::string("mass"))));
    tokenMap_.insert(TokenMapType::value_type("charge", Token(Token::charge, std::string("charge"))));
    tokenMap_.insert(TokenMapType::value_type("hull", Token(Token::hull, std::string("hull"))));
    tokenMap_.insert(TokenMapType::value_type("x", Token(Token::x, std::string("x"))));
    tokenMap_.insert(TokenMapType::value_type("y", Token(Token::y, std::string("y"))));
    tokenMap_.insert(TokenMapType::value_type("z", Token(Token::z, std::string("z"))));
    tokenMap_.insert(TokenMapType::value_type("wrappedx", Token(Token::wrappedX, std::string("wrappedx"))));
    tokenMap_.insert(TokenMapType::value_type("wrappedy", Token(Token::wrappedY, std::string("wrappedy"))));
    tokenMap_.insert(TokenMapType::value_type("wrappedz", Token(Token::wrappedZ, std::string("wrappedz"))));
    tokenMap_.insert(TokenMapType::value_type("r", Token(Token::r, std::string("r"))));
    tokenMap_.insert(TokenMapType::value_type("to", Token(Token::to, std::string("to"))));
    
    tokenMap_.insert(TokenMapType::value_type("all", Token(Token::all, std::string("all"))));
    tokenMap_.insert(TokenMapType::value_type("none", Token(Token::none, std::string("none"))));
  }

  Token* TokenMap::getToken(const std::string& ident) {
    std::map<std::string, Token>::iterator i = tokenMap_.find(ident);

    return i != tokenMap_.end() ? &(i->second) : NULL;
  }

  TokenMap::~TokenMap() {
    tokenMap_.clear();
  }
}
