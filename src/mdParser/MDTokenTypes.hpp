#ifndef INC_MDTokenTypes_hpp_
#define INC_MDTokenTypes_hpp_

/* $ANTLR 2.7.7 (20171128): "MDParser.g" -> "MDTokenTypes.hpp"$ */

#ifndef CUSTOM_API
# define CUSTOM_API
#endif

#ifdef __cplusplus
struct CUSTOM_API MDTokenTypes {
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
		ANALYZER = 26,
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
		ENDBLOCK = 45,
		ID = 46,
		ASSIGNEQUAL = 47,
		SEMICOLON = 48,
		StringLiteral = 49,
		LCURLY = 50,
		RCURLY = 51,
		LBRACKET = 52,
		RBRACKET = 53,
		LPAREN = 54,
		RPAREN = 55,
		COMMA = 56,
		NUM_INT = 57,
		NUM_LONG = 58,
		NUM_FLOAT = 59,
		NUM_DOUBLE = 60,
		DOT = 61,
		COLON = 62,
		QUESTIONMARK = 63,
		Whitespace = 64,
		Comment = 65,
		CPPComment = 66,
		PREPROC_DIRECTIVE = 67,
		LineDirective = 68,
		Space = 69,
		CharLiteral = 70,
		EndOfLine = 71,
		Escape = 72,
		Vocabulary = 73,
		Digit = 74,
		Decimal = 75,
		HEX_DIGIT = 76,
		EXPONENT = 77,
		FLOAT_SUFFIX = 78,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTokenTypes_hpp_*/
