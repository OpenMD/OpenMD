#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.7 (2006-11-01): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

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
		SEQUENCE = 18,
		MEMBERS = 19,
		CENTER = 20,
		SATELLITES = 21,
		POSITION = 22,
		ORIENTATION = 23,
		FLUCQ = 24,
		RNEMD = 25,
		LIGHT = 26,
		MINIMIZER = 27,
		FIXED = 28,
		HARMONIC = 29,
		CUBIC = 30,
		QUARTIC = 31,
		POLYNOMIAL = 32,
		MORSE = 33,
		GHOSTBEND = 34,
		UREYBRADLEY = 35,
		COSINE = 36,
		GHOSTTORSION = 37,
		CHARMM = 38,
		OPLS = 39,
		TRAPPE = 40,
		AMBERIMPROPER = 41,
		IMPROPERCOSINE = 42,
		CENTRALATOMHEIGHT = 43,
		DREIDING = 44,
		CHARGE = 45,
		NODES = 46,
		ENDBLOCK = 47,
		ID = 48,
		ASSIGNEQUAL = 49,
		SEMICOLON = 50,
		StringLiteral = 51,
		LCURLY = 52,
		RCURLY = 53,
		LBRACKET = 54,
		RBRACKET = 55,
		LPAREN = 56,
		RPAREN = 57,
		COMMA = 58,
		NUM_INT = 59,
		NUM_LONG = 60,
		NUM_FLOAT = 61,
		NUM_DOUBLE = 62,
		DOT = 63,
		COLON = 64,
		QUESTIONMARK = 65,
		Whitespace = 66,
		Comment = 67,
		CPPComment = 68,
		PREPROC_DIRECTIVE = 69,
		LineDirective = 70,
		Space = 71,
		CharLiteral = 72,
		EndOfLine = 73,
		Escape = 74,
		Vocabulary = 75,
		Digit = 76,
		Decimal = 77,
		HEX_DIGIT = 78,
		EXPONENT = 79,
		FLOAT_SUFFIX = 80,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
