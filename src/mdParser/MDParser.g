header
{

#include "antlr/CharScanner.hpp"
#include "utils/StringUtils.hpp"
#include "mdParser/FilenameObserver.hpp"
}

options
  {
  language = "Cpp";
  }
                     
class MDParser extends Parser;

options
{
        k = 3;
        exportVocab = MD;
        buildAST = true;
        codeGenMakeSwitchThreshold = 2;
        codeGenBitsetTestThreshold = 3;
        
}

tokens
{
  COMPONENT   = "component";
  MOLECULE    = "molecule";
  ZCONSTRAINT = "zconstraint";
  ATOM        = "atom";
  BOND        = "bond";
  BEND        = "bend";
  TORSION     = "torsion";
  RIGIDBODY   = "rigidBody";
  CUTOFFGROUP = "cutoffGroup";
  FRAGMENT    = "fragment";
  MEMBERS     = "members";
  POSITION    = "position";
  ORIENTATION = "orientation";
  ENDBLOCK;
}

{
    // Suppport C++-style single-line comments?
}

mdfile  : (statement)*
        ;

statement : assignment
          | componentblock
          | moleculeblock
          | zconstraintblock
          ;
            
assignment  : ID ASSIGNEQUAL^ constant SEMICOLON!
            ;
            
constant    : signedNumber
            | ID
            | StringLiteral
            ;

componentblock  : COMPONENT^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                ;
    
zconstraintblock  : ZCONSTRAINT^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                  ;
  
moleculeblock : MOLECULE^ LCURLY! (moleculestatement)*  RCURLY {#RCURLY->setType(ENDBLOCK);}
              ;

moleculestatement : assignment
                  | atomblock
                  | bondblock
                  | bendblock
                  | torsionblock
                  | rigidbodyblock
                  | cutoffgroupblock
                  | fragmentblock
                  ;

atomblock : ATOM^ LBRACKET! intConst RBRACKET! LCURLY! (atomstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
          ;

atomstatement : assignment
              | POSITION^ LPAREN! signedNumberTuple RPAREN! SEMICOLON!
              | ORIENTATION^  LPAREN! signedNumberTuple RPAREN! SEMICOLON!
              ;

                      
bondblock : BOND^ (LBRACKET! intConst! RBRACKET!)?  LCURLY!(bondstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
          ;

bondstatement : assignment
              | MEMBERS^ LPAREN! inttuple RPAREN! SEMICOLON!
              ;

bendblock : BEND^ (LBRACKET! intConst! RBRACKET!)? LCURLY!  (bendstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
          ;

bendstatement : assignment
              | MEMBERS^ LPAREN! inttuple RPAREN! SEMICOLON!
              ;

torsionblock  : TORSION^ (LBRACKET! intConst! RBRACKET!)?  LCURLY!(torsionstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
          ;

torsionstatement  : assignment
              | MEMBERS^ LPAREN! inttuple RPAREN! SEMICOLON!
              ;

rigidbodyblock  : RIGIDBODY^  LBRACKET! intConst RBRACKET! LCURLY!(rigidbodystatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                ;

rigidbodystatement  : assignment
              | MEMBERS^ LPAREN!  inttuple  RPAREN! SEMICOLON!
              ;

cutoffgroupblock  : CUTOFFGROUP^ (LBRACKET! intConst! RBRACKET!)? LCURLY! (cutoffgroupstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                  ;

cutoffgroupstatement  : assignment
              | MEMBERS^ LPAREN! inttuple RPAREN! SEMICOLON!
              ;

fragmentblock : FRAGMENT^ LBRACKET! intConst RBRACKET! LCURLY! (fragmentstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
              ;

fragmentstatement : assignment
              ;


              
signedNumberTuple   : signedNumber (COMMA! signedNumber)* 
              ;
                          
inttuple      : intConst (COMMA! intConst)*
              ;
              
protected
intConst
        :  OCTALINT | DECIMALINT | HEXADECIMALINT
        ;

protected
signedNumber  : (PLUS! | MINUS^)? 
                (intConst | floatConst)
              ;
              
protected
floatConst
        : 
          FLOATONE | FLOATTWO
        ;



class MDLexer extends Lexer;

options
  {
  k = 3;
  exportVocab = MD;
  testLiterals = false;
  }

tokens {
    DOT;
}

{


	int deferredLineCount;
	FilenameObserver* observer;
	
  public:
	void setObserver(FilenameObserver* osv) {observer = osv;}
  void initDeferredLineCount() { deferredLineCount = 0;}
  void deferredNewline() { 
        deferredLineCount++;
  }

	
  virtual void newline() { 
		for (;deferredLineCount>0;deferredLineCount--) {
			CharScanner::newline();
		}
    CharScanner::newline();
  }
    
}


// Operators:

ASSIGNEQUAL     : '=' ;
COLON           : ':' ;
COMMA           : ',' ;
QUESTIONMARK    : '?' ;
SEMICOLON       : ';' ;

LPAREN          : '(' ;
RPAREN          : ')' ;
LBRACKET        : '[' ;
RBRACKET        : ']' ;
LCURLY          : '{' ;
RCURLY          : '}' ;

PLUS            : '+' ;
MINUS           : '-' ;

/*
EQUAL           : "==" ;
NOTEQUAL        : "!=" ;
LESSTHANOREQUALTO     : "<=" ;
LESSTHAN              : "<" ;
GREATERTHANOREQUALTO  : ">=" ;
GREATERTHAN           : ">" ;

DIVIDE          : '/' ;
DIVIDEEQUAL     : "/=" ;
PLUS            : '+' ;
PLUSEQUAL       : "+=" ;
PLUSPLUS        : "++" ;
MINUS           : '-' ;
MINUSEQUAL      : "-=" ;
MINUSMINUS      : "--" ;
STAR            : '*' ;
TIMESEQUAL      : "*=" ;
MOD             : '%' ;
MODEQUAL        : "%=" ;
SHIFTRIGHT      : ">>" ;
SHIFTRIGHTEQUAL : ">>=" ;
SHIFTLEFT       : "<<" ;
SHIFTLEFTEQUAL  : "<<=" ;

AND            : "&&" ;
NOT            : '!' ;
OR             : "||" ;

AMPERSAND       : '&' ;
BITWISEANDEQUAL : "&=" ;
TILDE           : '~' ;
BITWISEOR       : '|' ;
BITWISEOREQUAL  : "|=" ;
BITWISEXOR      : '^' ;
BITWISEXOREQUAL : "^=" ;
*/


Whitespace  
  : 
    ( // whitespace ignored
      (' ' |'\t' | '\f')
    | // handle newlines
      ( '\r' '\n' // MS
      | '\r'    // Mac
      | '\n'    // Unix 
      ) { newline(); }
    | // handle continuation lines
      ( '\\' '\r' '\n'  // MS
      | '\\' '\r'   // Mac
      | '\\' '\n'   // Unix 
      ) {printf("CPP_parser.g continuation line detected\n");
        deferredNewline();}
    ) 
    {_ttype = ANTLR_USE_NAMESPACE(antlr)Token::SKIP;}
  ;

Comment  
  : 
    "/*"   
    ( {LA(2) != '/'}? '*'
    | EndOfLine {deferredNewline();}
    | ~('*'| '\r' | '\n')
    )*
    "*/" {_ttype = ANTLR_USE_NAMESPACE(antlr)Token::SKIP;}
  ;

CPPComment
  : 
    "//" (~('\n' | '\r'))* EndOfLine
    {_ttype = ANTLR_USE_NAMESPACE(antlr)Token::SKIP; newline();}                     
  ;

PREPROC_DIRECTIVE
  options{paraphrase = "a line directive";}
  : 
    '#' LineDirective
    {_ttype = ANTLR_USE_NAMESPACE(antlr)Token::SKIP; newline();} 
  ;

protected 
LineDirective
  :
  {
   deferredLineCount = 0;
  }
    ("line")?  // this would be for if the directive started "#line"
    (Space)+
    n:Decimal { setLine(oopse::lexi_cast<int>(n->getText()) - 1); } 
    (Space)+
    (sl:StringLiteral) {std::string filename = sl->getText().substr(1,sl->getText().length()-2); observer->notify(filename);} 
    ((Space)+ Decimal)* // To support cpp flags (GNU)
    EndOfLine
  ;

protected  
Space
  : 
    (' '|'\t'|'\f')
  ;


// Literals:

/*
 * Note that we do NOT handle tri-graphs nor multi-byte sequences.
 */

/*
 * Note that we can't have empty character constants (even though we
 * can have empty strings :-).
 */
CharLiteral
  : 
    '\'' (Escape | ~('\'')) '\''
  ;

/*
 * Can't have raw imbedded newlines in string constants.  Strict reading of
 * the standard gives odd dichotomy between newlines & carriage returns.
 * Go figure.
 */
StringLiteral
  : 
    '"'
    ( Escape
    | 
      ( "\\\r\n"   // MS 
      | "\\\r"     // MAC
      | "\\\n"     // Unix
      ) {deferredNewline();}
    | 
      ~('"'|'\r'|'\n'|'\\')
    )*
    '"'
  ;

protected
EndOfLine
  : 
    ( options{generateAmbigWarnings = false;}:
      "\r\n"  // MS
    | '\r'    // Mac
    | '\n'    // Unix
    )
  ;

/*
 * Handle the various escape sequences.
 *
 * Note carefully that these numeric escape *sequences* are *not* of the
 * same form as the C language numeric *constants*.
 *
 * There is no such thing as a binary numeric escape sequence.
 *
 * Octal escape sequences are either 1, 2, or 3 octal digits exactly.
 *
 * There is no such thing as a decimal escape sequence.
 *
 * Hexadecimal escape sequences are begun with a leading \x and continue
 * until a non-hexadecimal character is found.
 *
 * No real handling of tri-graph sequences, yet.
 */

protected
Escape  
  : 
    '\\'
    ( options{warnWhenFollowAmbig=false;}:
      'a'
    | 'b'
    | 'f'
    | 'n'
    | 'r'
    | 't'
    | 'v'
    | '"'
    | '\''
    | '\\'
    | '?'
    | ('0'..'3') (options{warnWhenFollowAmbig=false;}: Digit (options{warnWhenFollowAmbig=false;}: Digit)? )?
    | ('4'..'7') (options{warnWhenFollowAmbig=false;}: Digit)?
    | 'x' (options{warnWhenFollowAmbig=false;}: Digit | 'a'..'f' | 'A'..'F')+
    )
  ;

// Numeric Constants: 

protected
Digit
  : 
    '0'..'9'
  ;

protected
Decimal
  : 
    ('0'..'9')+
  ;

protected
LongSuffix
  : 'l'
  | 'L'
  ;

protected
UnsignedSuffix
  : 'u'
  | 'U'
  ;

protected
FloatSuffix
  : 'f'
  | 'F'
  ;

protected
Exponent
  : 
    ('e'|'E'|'d'|'D') ('+'|'-')? (Digit)+
  ;

protected
Vocabulary
  : 
    '\3'..'\377'
  ;

Number
  : 
    ( (Digit)+ ('.' | 'e' | 'E' | 'd' | 'D' ) )=> 
    (Digit)+
    ( '.' (Digit)* (Exponent)? {_ttype = FLOATONE;} //Zuo 3/12/01
    | Exponent                 {_ttype = FLOATTWO;} //Zuo 3/12/01
    )                          //{_ttype = DoubleDoubleConst;}
    (FloatSuffix               //{_ttype = FloatDoubleConst;}
    |LongSuffix                //{_ttype = LongDoubleConst;}
    )?
  | 
    '.'                        {_ttype = DOT;}
    ( (Digit)+ (Exponent)?   {_ttype = FLOATONE;} //Zuo 3/12/01
                                   //{_ttype = DoubleDoubleConst;}
      (FloatSuffix           //{_ttype = FloatDoubleConst;}
      |LongSuffix            //{_ttype = LongDoubleConst;}
      )?
    )?
  | 
    '0' ('0'..'7')*            //{_ttype = IntOctalConst;}
    (LongSuffix                //{_ttype = LongOctalConst;}
    |UnsignedSuffix            //{_ttype = UnsignedOctalConst;}
    )*                         {_ttype = OCTALINT;}
  | 
    '1'..'9' (Digit)*          //{_ttype = IntIntConst;}
    (LongSuffix                //{_ttype = LongIntConst;}
    |UnsignedSuffix            //{_ttype = UnsignedIntConst;}
    )*                         {_ttype = DECIMALINT;}  
  | 
    '0' ('x' | 'X') ('a'..'f' | 'A'..'F' | Digit)+
                                   //{_ttype = IntHexConst;}
    (LongSuffix                //{_ttype = LongHexConst;}
    |UnsignedSuffix            //{_ttype = UnsignedHexConst;}
    )*                         {_ttype = HEXADECIMALINT;}   
  ;

ID
  options {testLiterals = true;}
  : 
    ('a'..'z'|'A'..'Z'|'_')
    ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*
  ;
