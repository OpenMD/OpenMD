#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.7 (20160304): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

#ifndef CUSTOM_API
# define CUSTOM_API
#endif

#ifdef __cplusplus
struct CUSTOM_API MDTreeParserTokenTypes {
#endif
	enum {
		EOF_ = 1,
		COMPONENT = 4,
		MOLECULE = 5,
		ZCONSTRAINT = 6,
		RESTRAINT = 7,
		ATOM = 8,
		BOND = 9,
		BEND = 10,
		TORSION = 11,
		INVERSION = 12,
		RIGIDBODY = 13,
		CUTOFFGROUP = 14,
		CONSTRAINT = 15,
		DISTANCE = 16,
		FRAGMENT = 17,
		MEMBERS = 18,
		CENTER = 19,
		SATELLITES = 20,
		POSITION = 21,
		ORIENTATION = 22,
		FLUCQ = 23,
		RNEMD = 24,
		MINIMIZER = 25,
		FIXED = 26,
		HARMONIC = 27,
		CUBIC = 28,
		QUARTIC = 29,
		POLYNOMIAL = 30,
		MORSE = 31,
		GHOSTBEND = 32,
		UREYBRADLEY = 33,
		COSINE = 34,
		GHOSTTORSION = 35,
		CHARMM = 36,
		OPLS = 37,
		TRAPPE = 38,
		AMBERIMPROPER = 39,
		IMPROPERCOSINE = 40,
		CENTRALATOMHEIGHT = 41,
		DREIDING = 42,
		CHARGE = 43,
		ENDBLOCK = 44,
		ID = 45,
		ASSIGNEQUAL = 46,
		SEMICOLON = 47,
		StringLiteral = 48,
		LCURLY = 49,
		RCURLY = 50,
		LBRACKET = 51,
		RBRACKET = 52,
		LPAREN = 53,
		RPAREN = 54,
		COMMA = 55,
		NUM_INT = 56,
		NUM_LONG = 57,
		NUM_FLOAT = 58,
		NUM_DOUBLE = 59,
		DOT = 60,
		COLON = 61,
		QUESTIONMARK = 62,
		Whitespace = 63,
		Comment = 64,
		CPPComment = 65,
		PREPROC_DIRECTIVE = 66,
		LineDirective = 67,
		Space = 68,
		CharLiteral = 69,
		EndOfLine = 70,
		Escape = 71,
		Vocabulary = 72,
		Digit = 73,
		Decimal = 74,
		HEX_DIGIT = 75,
		EXPONENT = 76,
		FLOAT_SUFFIX = 77,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
