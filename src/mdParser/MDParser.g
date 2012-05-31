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
  RESTRAINT   = "restraint";
  ATOM        = "atom";
  BOND        = "bond";
  BEND        = "bend";
  TORSION     = "torsion";
  INVERSION   = "inversion";
  RIGIDBODY   = "rigidBody";
  CUTOFFGROUP = "cutoffGroup";
  FRAGMENT    = "fragment";
  MEMBERS     = "members";
  CENTER      = "center";
  POSITION    = "position";
  ORIENTATION = "orientation";
  FLUCQ       = "flucQ";
  RNEMD       = "RNEMD";
  ENDBLOCK;
}


mdfile  : (statement)*
        ;

statement : assignment
    | componentblock
    | moleculeblock
    | zconstraintblock
    | restraintblock
    | flucqblock
    | rnemdblock
    ;

assignment  : ID ASSIGNEQUAL^ constant SEMICOLON!
            ;
            
constant    : intConst
						| floatConst
            | ID
            | StringLiteral
            ;

componentblock  : COMPONENT^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                ;
    
zconstraintblock  : ZCONSTRAINT^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                  ;

restraintblock  : RESTRAINT^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
                  ;

flucqblock  : FLUCQ^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
    ;

rnemdblock  : RNEMD^ LCURLY! (assignment)* RCURLY {#RCURLY->setType(ENDBLOCK);}
    ;
  
moleculeblock : MOLECULE^ LCURLY! (moleculestatement)*  RCURLY {#RCURLY->setType(ENDBLOCK);}
              ;

moleculestatement : assignment
                  | atomblock
                  | bondblock
                  | bendblock
                  | torsionblock
                  | inversionblock
                  | rigidbodyblock
                  | cutoffgroupblock
                  | fragmentblock
                  ;

atomblock : ATOM^ LBRACKET! intConst RBRACKET! LCURLY! (atomstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
          ;

atomstatement : assignment
              | POSITION^ LPAREN! doubleNumberTuple RPAREN! SEMICOLON!
              | ORIENTATION^  LPAREN! doubleNumberTuple RPAREN! SEMICOLON!
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

inversionblock  : INVERSION^ (LBRACKET! intConst! RBRACKET!)?  LCURLY!(inversionstatement)* RCURLY {#RCURLY->setType(ENDBLOCK);}
          ;

inversionstatement  : assignment
              | CENTER^ LPAREN! intConst RPAREN! SEMICOLON!
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


              
doubleNumberTuple   : doubleNumber (COMMA! doubleNumber)* 
              ;
                          
inttuple      : intConst (COMMA! intConst)*
              ;
              
protected
intConst
        :  NUM_INT | NUM_LONG
        ;

protected
doubleNumber  :  
                (intConst | floatConst)
              ;
              
protected
floatConst
        : 
          NUM_FLOAT | NUM_DOUBLE
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
    n:Decimal { setLine(OpenMD::lexi_cast<int>(n->getText()) - 1); } 
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


protected
Vocabulary
  : 
    '\3'..'\377'
  ;


ID
  options {testLiterals = true;}
  : 
    ('a'..'z'|'A'..'Z'|'_')
    ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*
  ;


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

// hexadecimal digit (again, note it's protected!)
protected
HEX_DIGIT
	:	('0'..'9'|'A'..'F'|'a'..'f')
	;


// a numeric literal
NUM_INT
	{
		bool isDecimal = false;
		ANTLR_USE_NAMESPACE(antlr)RefToken t = ANTLR_USE_NAMESPACE(antlr)nullToken;
	}
    : ('+' | '-')?
    (
      '.' {_ttype = DOT;}
            (	('0'..'9')+ (EXPONENT)? (f1:FLOAT_SUFFIX {t=f1;})?
            {
					if ( t &&
						  (t->getText().find('f') != ANTLR_USE_NAMESPACE(std)string::npos ||
							t->getText().find('F') != ANTLR_USE_NAMESPACE(std)string::npos ) ) {
						_ttype = NUM_FLOAT;
					}
					else {
						_ttype = NUM_DOUBLE; // assume double
					}
				}
            )?

	|	(	'0' {isDecimal = true;} // special case for just '0'
			(	('x'|'X')
				(											// hex
					// the 'e'|'E' and float suffix stuff look
					// like hex digits, hence the (...)+ doesn't
					// know when to stop: ambig.  ANTLR resolves
					// it correctly by matching immediately.  It
					// is therefor ok to hush warning.
					options {
						warnWhenFollowAmbig=false;
					}
				:	HEX_DIGIT
				)+
			|	//float or double with leading zero
				(('0'..'9')+ ('.'|EXPONENT|FLOAT_SUFFIX)) => ('0'..'9')+
			|	('0'..'7')+									// octal
			)?
		|	('1'..'9') ('0'..'9')*  {isDecimal=true;}		// non-zero decimal
		)
		(	('l'|'L') { _ttype = NUM_LONG; }

		// only check to see if it's a float if looks like decimal so far
		|	{isDecimal}?
            (   '.' ('0'..'9')* (EXPONENT)? (f2:FLOAT_SUFFIX {t=f2;})?
            |   EXPONENT (f3:FLOAT_SUFFIX {t=f3;})?
            |   f4:FLOAT_SUFFIX {t=f4;}
            )
            {
					if ( t &&
						  (t->getText().find('f') != ANTLR_USE_NAMESPACE(std)string::npos ||
							t->getText().find('F') != ANTLR_USE_NAMESPACE(std)string::npos ) ) {
						_ttype = NUM_FLOAT;
					}
					else {
						_ttype = NUM_DOUBLE; // assume double
					}
				}
        )?
  )
	;

// a couple protected methods to assist in matching floating point numbers
protected
EXPONENT
	:	('e'|'E'|'d'|'D') ('+'|'-')? ('0'..'9')+
	;

protected
FLOAT_SUFFIX
	:	'f'|'F'|'d'|'D'
	;
