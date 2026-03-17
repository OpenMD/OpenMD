// ANTLR v4 grammar for OpenMD Parser
// Single source grammar for C++ and Python targets.
//
// To generate C++ parser:
//   antlr4 -Dlanguage=Cpp -visitor OMD.g4
// 
// The @header and @parser/@lexer::members sections are C++-specific. To make
// a python grammer, strip these blocks and put the result in ../python/OMD.g4
//
// Then to generate Python parser:
//   antlr4 -Dlanguage=Python3 -visitor OMD.g4

grammar OMD;

// ============================================================================
// C++ TARGET SPECIFIC ACTIONS
// ============================================================================

@header {
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"
}

@parser::members {
    FilenameObserver* observer;

    void setObserver(FilenameObserver* osv) {
        observer = osv;
    }
}

@lexer::members {
    int deferredLineCount = 0;
    FilenameObserver* observer;

    void setObserver(FilenameObserver* osv) {
        observer = osv;
    }

    void initDeferredLineCount() {
        deferredLineCount = 0;
    }

    void deferredNewline() {
        deferredLineCount++;
    }
}

// ============================================================================
// PARSER RULES
// ============================================================================

omdfile
    : statement* EOF
    ;

statement
    : assignment
    | componentblock
    | moleculeblock
    | fragmentblock
    | zconstraintblock
    | restraintblock
    | flucqblock
    | rnemdblock
    | lightblock
    | minimizerblock
    ;

assignment
    : ID ASSIGNEQUAL constant SEMICOLON
    ;

constant
    : intConst
    | floatConst
    | vectorConst
    | TRUE
    | FALSE
    | ID
    | StringLiteral
    ;

componentblock
    : COMPONENT LCURLY assignment* RCURLY
    ;

zconstraintblock
    : ZCONSTRAINT LCURLY assignment* RCURLY
    ;

restraintblock
    : RESTRAINT LCURLY assignment* RCURLY
    ;

flucqblock
    : FLUCQ LCURLY assignment* RCURLY
    ;

rnemdblock
    : RNEMD LCURLY assignment* RCURLY
    ;

lightblock
    : LIGHT LCURLY assignment* RCURLY
    ;

minimizerblock
    : MINIMIZER LCURLY assignment* RCURLY
    ;

moleculeblock
    : MOLECULE LCURLY moleculestatement* RCURLY
    ;

moleculestatement
    : assignment
    | atomblock
    | bondblock
    | bendblock
    | torsionblock
    | inversionblock
    | rigidbodyblock
    | cutoffgroupblock
    | constraintblock
    | sequencestring
    ;

atomblock
    : ATOM LBRACKET intConst RBRACKET LCURLY atomstatement* RCURLY
    ;

atomstatement
    : assignment
    | POSITION LPAREN doubleNumberTuple RPAREN SEMICOLON
    | ORIENTATION LPAREN doubleNumberTuple RPAREN SEMICOLON
    | CHARGE LPAREN floatConst RPAREN SEMICOLON
    ;

bondblock
    : BOND (LBRACKET intConst RBRACKET)? LCURLY bondstatement* RCURLY
    ;

bondstatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    | FIXED LPAREN floatConst RPAREN SEMICOLON
    | HARMONIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | CUBIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | QUARTIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | POLYNOMIAL LPAREN doubleNumberTuple RPAREN SEMICOLON
    | MORSE LPAREN doubleNumberTuple RPAREN SEMICOLON
    ;

bendblock
    : BEND (LBRACKET intConst RBRACKET)? LCURLY bendstatement* RCURLY
    ;

bendstatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    | HARMONIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | GHOSTBEND LPAREN doubleNumberTuple RPAREN SEMICOLON
    | UREYBRADLEY LPAREN doubleNumberTuple RPAREN SEMICOLON
    | CUBIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | QUARTIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | POLYNOMIAL LPAREN doubleNumberTuple RPAREN SEMICOLON
    | COSINE LPAREN doubleNumberTuple RPAREN SEMICOLON
    ;

torsionblock
    : TORSION (LBRACKET intConst RBRACKET)? LCURLY torsionstatement* RCURLY
    ;

torsionstatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    | GHOSTTORSION LPAREN doubleNumberTuple RPAREN SEMICOLON
    | CUBIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | QUARTIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | POLYNOMIAL LPAREN doubleNumberTuple RPAREN SEMICOLON
    | CHARMM LPAREN doubleNumberTuple RPAREN SEMICOLON
    | OPLS LPAREN doubleNumberTuple RPAREN SEMICOLON
    | TRAPPE LPAREN doubleNumberTuple RPAREN SEMICOLON
    | HARMONIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    ;

inversionblock
    : INVERSION (LBRACKET intConst RBRACKET)? LCURLY inversionstatement* RCURLY
    ;

inversionstatement
    : assignment
    | CENTER LPAREN intConst RPAREN SEMICOLON
    | SATELLITES LPAREN inttuple RPAREN SEMICOLON
    | AMBERIMPROPER LPAREN doubleNumberTuple RPAREN SEMICOLON
    | IMPROPERCOSINE LPAREN doubleNumberTuple RPAREN SEMICOLON
    | HARMONIC LPAREN doubleNumberTuple RPAREN SEMICOLON
    | CENTRALATOMHEIGHT LPAREN doubleNumberTuple RPAREN SEMICOLON
    | DREIDING LPAREN doubleNumberTuple RPAREN SEMICOLON
    ;

rigidbodyblock
    : RIGIDBODY LBRACKET intConst RBRACKET LCURLY rigidbodystatement* RCURLY
    ;

rigidbodystatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    ;

cutoffgroupblock
    : CUTOFFGROUP (LBRACKET intConst RBRACKET)? LCURLY cutoffgroupstatement* RCURLY
    ;

cutoffgroupstatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    ;

nodesblock
    : NODES (LBRACKET intConst RBRACKET)? LCURLY nodesstatement* RCURLY
    ;

nodesstatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    ;

fragmentblock
    : FRAGMENT LCURLY fragmentstatement* RCURLY
    ;

fragmentstatement
    : assignment
    | atomblock
    | bondblock
    | bendblock
    | torsionblock
    | inversionblock
    | rigidbodyblock
    | cutoffgroupblock
    | constraintblock
    | nodesblock
    ;

constraintblock
    : CONSTRAINT (LBRACKET intConst RBRACKET)? LCURLY constraintstatement* RCURLY
    ;

constraintstatement
    : assignment
    | MEMBERS LPAREN inttuple RPAREN SEMICOLON
    ;

sequencestring
    : SEQUENCE ASSIGNEQUAL constant SEMICOLON
    ;

doubleNumberTuple
    : doubleNumber (COMMA doubleNumber)*
    ;

inttuple
    : intConst (COMMA intConst)*
    ;

intConst
    : NUM_INT
    | NUM_LONG
    ;

doubleNumber
    : intConst
    | floatConst
    ;

floatConst
    : NUM_FLOAT
    | NUM_DOUBLE
    ;

vectorConst
    : LPAREN doubleNumberTuple RPAREN
    ;

// ============================================================================
// LEXER RULES
// ============================================================================

// Boolean keywords — must come before ID
TRUE  : 'true' ;
FALSE : 'false' ;

// Block and structure keywords
COMPONENT       : 'component' ;
MOLECULE        : 'molecule' ;
ZCONSTRAINT     : 'zconstraint' ;
RESTRAINT       : 'restraint' ;
ATOM            : 'atom' ;
BOND            : 'bond' ;
BEND            : 'bend' ;
TORSION         : 'torsion' ;
INVERSION       : 'inversion' ;
RIGIDBODY       : 'rigidBody' ;
CUTOFFGROUP     : 'cutoffGroup' ;
CONSTRAINT      : 'constraint' ;
DISTANCE        : 'distance' ;
FRAGMENT        : 'fragment' ;
SEQUENCE        : 'sequence' ;
MEMBERS         : 'members' ;
CENTER          : 'center' ;
SATELLITES      : 'satellites' ;
POSITION        : 'position' ;
ORIENTATION     : 'orientation' ;
FLUCQ           : 'flucQ' ;
RNEMD           : 'RNEMD' ;
LIGHT           : 'light' ;
MINIMIZER       : 'minimizer' ;

// Potential type keywords
FIXED           : 'Fixed' ;
HARMONIC        : 'Harmonic' ;
CUBIC           : 'Cubic' ;
QUARTIC         : 'Quartic' ;
POLYNOMIAL      : 'Polynomial' ;
MORSE           : 'Morse' ;
GHOSTBEND       : 'GhostBend' ;
UREYBRADLEY     : 'UreyBradley' ;
COSINE          : 'Cosine' ;
GHOSTTORSION    : 'GhostTorsion' ;
CHARMM          : 'Charmm' ;
OPLS            : 'Opls' ;
TRAPPE          : 'Trappe' ;
AMBERIMPROPER   : 'AmberImproper' ;
IMPROPERCOSINE  : 'ImproperCosine' ;
CENTRALATOMHEIGHT : 'CentralAtomHeight' ;
DREIDING        : 'Dreiding' ;
CHARGE          : 'charge' ;
NODES           : 'nodes' ;

// Operators
ASSIGNEQUAL     : '=' ;
COLON           : ':' ;
COMMA           : ',' ;
QUESTIONMARK    : '?' ;
SEMICOLON       : ';' ;
DOT             : '.' ;

LPAREN          : '(' ;
RPAREN          : ')' ;
LBRACKET        : '[' ;
RBRACKET        : ']' ;
LCURLY          : '{' ;
RCURLY          : '}' ;

// Numeric literals
// NUM_LONG must precede NUM_INT so the mandatory [lL] suffix is captured
// before NUM_INT can match without it.
NUM_LONG
    : [+\-]? (
        '0' [xX] HEX_DIGIT+                     // Hexadecimal
      | '0' [0-7]+                              // Octal
      | [1-9] [0-9]*                            // Decimal
      | '0'                                     // Zero
      ) [lL]
    ;

NUM_INT
    : [+\-]? (
        '0' [xX] HEX_DIGIT+                     // Hexadecimal
      | '0' [0-7]+                              // Octal
      | [1-9] [0-9]*                            // Decimal
      | '0'                                     // Zero
      )
    ;

NUM_FLOAT
    : [+\-]? (
        [0-9]+ '.' [0-9]* EXPONENT? [fF]
      | '.' [0-9]+ EXPONENT? [fF]
      | [0-9]+ EXPONENT [fF]
      | [0-9]+ [fF]
      )
    ;

NUM_DOUBLE
    : [+\-]? (
        [0-9]+ '.' [0-9]* EXPONENT? [dD]?
      | '.' [0-9]+ EXPONENT? [dD]?
      | [0-9]+ EXPONENT [dD]?
      | [0-9]+ [dD]
      )
    ;

// Character and string literals
CharLiteral
    : '\'' (Escape | ~[']) '\''
    ;

StringLiteral
    : '"' (Escape | '\\' ('\r\n' | '\r' | '\n') | ~["\r\n\\])* '"'
    ;

// Identifiers — must come after all keywords so keywords are not lexed as ID
ID
    : [a-zA-Z_][a-zA-Z_0-9]*
    ;

// Whitespace
Whitespace
    : [ \t\f]+                   -> skip
    ;

Newline
    : ('\r\n' | '\r' | '\n')     -> skip
    ;

// Line continuation
LineContinuation
    : '\\' ('\r\n' | '\r' | '\n')  -> skip
    ;

// Comments
Comment
    : '/*' .*? '*/'              -> skip
    ;

CPPComment
    : '//' ~[\r\n]*              -> skip
    ;

// Preprocessor directives (includes handled externally before parsing)
PREPROC_DIRECTIVE
    : '#' ~[\r\n]*               -> skip
    ;

// Fragment rules
fragment
HEX_DIGIT
    : [0-9A-Fa-f]
    ;

fragment
EXPONENT
    : [eEdD] [+\-]? [0-9]+
    ;

fragment
Escape
    : '\\' (
        [abfnrtv"'\\?]
      | [0-3] [0-7]? [0-7]?
      | [4-7] [0-7]?
      | 'x' [0-9a-fA-F]+
      )
    ;

fragment
EndOfLine
    : '\r\n'
    | '\r'
    | '\n'
    ;
