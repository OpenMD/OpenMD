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

#ifndef SELECTION_TOKEN_HPP
#define SELECTION_TOKEN_HPP


namespace oopse {


/**
 * @class Token
 * @todo document
 */
class Token {

    public:

        int tok;
        int intValue;
        
        const static int nada              =  0;
        const static int identifier        =  1;
        const static int integer           =  2;
        const static int decimal           =  3;
        const static int string            =  4;
        const static int seqcode           =  5;
        const static int unknown           =  6;
        const static int keyword           =  7;
        const static int whitespace        =  8;
        const static int comment           =  9;
        const static int endofline         = 10;
        const static int endofstatement    = 11;

        const static int command           = (1 <<  8);
        const static int expressionCommand = (1 <<  9); // expression command
        const static int expression        = (1 << 15); /// expression term

        // every property is also valid in an expression context
        const static int atomproperty      = (1 << 16) | expression;
        // every predefined is also valid in an expression context
        const static int comparator        = (1 << 17) | expression;
        const static int predefinedset     = (1 << 18) | expression;
        // generally, the minus sign is used to denote atom ranges
        // this property is used for the few commands which allow negative integers
        const static int negnums      = (1 << 21);
        // for some commands the 'set' is optional

        // These are unrelated
        const static int varArgCount       = (1 << 4);
        const static int onDefault1        = (1 << 5) | 1;

        // rasmol commands
        const static int define       = command | expressionCommand |1;
        const static int select       = command |expressionCommand |2 ;
        const static int all          = expression | 11 ; 

        // atom expression operators
        const static int leftparen    = expression |  0;
        const static int rightparen   = expression |  1;
        const static int hyphen       = expression |  2;
        const static int opAnd        = expression |  3;
        const static int opOr         = expression |  4;
        const static int opNot        = expression |  5;
        const static int within       = expression |  6;
        const static int plus         = expression |  7;
        const static int pick         = expression |  8;
        const static int asterisk     = expression |  9;
        const static int dot          = expression | 11;
        const static int leftsquare   = expression | 12;
        const static int rightsquare  = expression | 13;
        const static int colon        = expression | 14;
        const static int slash        = expression | 15;

        // miguel 2005 01 01
        // these are used to demark the beginning and end of expressions
        // they do not exist in the source code, but are emitted by the
        // expression compiler
        const static int expressionBegin = expression | 100;
        const static int expressionEnd   = expression | 101;

        const static int atomno       = atomproperty | 0;
        const static int elemno       = atomproperty | 1;
        const static int resno        = atomproperty | 2;
        const static int radius       = atomproperty | 3 ; 
        const static int _bondedcount = atomproperty | 6;
        const static int _groupID     = atomproperty | 7;
        const static int _atomID      = atomproperty | 8;
        const static int _structure   = atomproperty | 9;
        const static int occupancy    = atomproperty | 10;
        const static int polymerLength= atomproperty | 11;

        const static int opGT         = comparator |  0;
        const static int opGE         = comparator |  1;
        const static int opLE         = comparator |  2;
        const static int opLT         = comparator |  3;
        const static int opEQ         = comparator |  4;
        const static int opNE         = comparator |  5;
 
        const static int x            =  expression |2;
        const static int y            =  expression | 3;
        const static int z            =  expression | 4;
        const static int none      =  expression | 5;
 
 
        const static Token tokenAll(all, "all");
        const static Token tokenAnd(opAnd, "and");
        const static Token tokenElemno(elemno, "elemno");
        const static Token tokenExpressionBegin(expressionBegin, "expressionBegin");
        const static Token tokenExpressionEnd(expressionEnd, "expressionEnd");


        const static String[] comparatorNames = {">", ">=", "<=", "<", "=", "!="};
        const static String[] atomPropertyNames = {
        "atomno", "elemno", "resno", "radius", "temperature", "model",
        "_bondedcount", "_groupID", "_atomID", "_structure"};

        /*
        Note that the RasMol scripting language is case-insensitive.
        So, the compiler turns all identifiers to lower-case before
        looking up in the hash table. 
        Therefore, the left column of this array *must* be lower-case
        */

        const static Object[] arrayPairs  = {
        // commands 
        "define",            new Token(define,   varArgCount, "define"), 
        "select",            new Token(select,   varArgCount, "select"), 
        // atom expressions
        "(",            new Token(leftparen, "("),
        ")",            new Token(rightparen, ")"),
        "-",            new Token(hyphen, "-"),
        "and",          tokenAnd,
        "&",            null,
        "&&",           null,
        "or",           new Token(opOr, "or"),
        ",",            null,
        "|",            null,
        "||",            null,
        "not",          new Token(opNot, "not"),
        "!",            null,
        "<",            new Token(opLT, "<"),
        "<=",           new Token(opLE, "<="),
        ">=",           new Token(opGE, ">="),
        ">",            new Token(opGT, ">="),
        "==",           new Token(opEQ, "=="),
        "=",            null,
        "!=",           new Token(opNE, "!="),
        "<>",           null,
        "/=",           null,
        "within",       new Token(within, "within"),
        "+",            new Token(plus, "+"),
        "pick",         new Token(pick, "pick"),
        ".",            new Token(dot, "."),
        "[",            new Token(leftsquare,  "["),
        "]",            new Token(rightsquare, "]"),
        ":",            new Token(colon, ":"),
        "/",            new Token(slash, "/"),

        "atomno",       new Token(atomno, "atomno"),
        "elemno",       tokenElemno,
        "_e",           tokenElemno,
        "resno",        new Token(resno, "resno"),
        "temperature",  new Token(temperature, "temperature"),
        "relativetemperature",  null,
        "_bondedcount", new Token(_bondedcount, "_bondedcount"),
        "_groupID",     new Token(_groupID, "_groupID"),
        "_g",           null,
        "_atomID",      new Token(_atomID, "_atomID"),
        "_a",           null,
        "_structure",   new Token(_structure, "_structure"),
        "occupancy",    new Token(occupancy, "occupancy"),
        "polymerlength",new Token(polymerLength, "polymerlength"),
 
        "x",            new Token(x, "x"),
        "y",            new Token(y, "y"),
        "z",            new Token(z, "z"),
        "*",            new Token(asterisk, "*"),
        "all",          tokenAll,
        "none",         new Token(none, "none"),
 
         };


  static Hashtable map = new Hashtable();
  static {
    Token tokenLast = null;
    String stringThis;
    Token tokenThis;
    for (int i = 0; i + 1 < arrayPairs.length; i += 2) {
      stringThis = (String) arrayPairs[i];
      tokenThis = (Token) arrayPairs[i + 1];
      if (tokenThis == null)
        tokenThis = tokenLast;
      if (map.get(stringThis) != null)
        System.out.println("duplicate token definition:" + stringThis);
      map.put(stringThis, tokenThis);
      tokenLast = tokenThis;
    }
  }
  

};

}

#endif 