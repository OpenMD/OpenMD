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
		MINIMIZER = 26,
		FIXED = 27,
		HARMONIC = 28,
		CUBIC = 29,
		QUARTIC = 30,
		POLYNOMIAL = 31,
		MORSE = 32,
		GHOSTBEND = 33,
		UREYBRADLEY = 34,
		COSINE = 35,
		GHOSTTORSION = 36,
		CHARMM = 37,
		OPLS = 38,
		TRAPPE = 39,
		AMBERIMPROPER = 40,
		IMPROPERCOSINE = 41,
		CENTRALATOMHEIGHT = 42,
		DREIDING = 43,
		CHARGE = 44,
		NODES = 45,
		ENDBLOCK = 46,
		ID = 47,
		ASSIGNEQUAL = 48,
		SEMICOLON = 49,
		StringLiteral = 50,
		LCURLY = 51,
		RCURLY = 52,
		LBRACKET = 53,
		RBRACKET = 54,
		LPAREN = 55,
		RPAREN = 56,
		COMMA = 57,
		NUM_INT = 58,
		NUM_LONG = 59,
		NUM_FLOAT = 60,
		NUM_DOUBLE = 61,
		DOT = 62,
		COLON = 63,
		QUESTIONMARK = 64,
		Whitespace = 65,
		Comment = 66,
		CPPComment = 67,
		PREPROC_DIRECTIVE = 68,
		LineDirective = 69,
		Space = 70,
		CharLiteral = 71,
		EndOfLine = 72,
		Escape = 73,
		Vocabulary = 74,
		Digit = 75,
		Decimal = 76,
		HEX_DIGIT = 77,
		EXPONENT = 78,
		FLOAT_SUFFIX = 79,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
