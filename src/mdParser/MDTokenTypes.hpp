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
		ANALYSIS = 8,
		ATOM = 9,
		BOND = 10,
		BEND = 11,
		TORSION = 12,
		INVERSION = 13,
		RIGIDBODY = 14,
		CUTOFFGROUP = 15,
		CONSTRAINT = 16,
		DISTANCE = 17,
		FRAGMENT = 18,
		MEMBERS = 19,
		CENTER = 20,
		SATELLITES = 21,
		POSITION = 22,
		ORIENTATION = 23,
		FLUCQ = 24,
		RNEMD = 25,
		MINIMIZER = 26,
		ANALYZER = 27,
		NUDGEDELASTICBAND = 28,
		FIXED = 29,
		HARMONIC = 30,
		CUBIC = 31,
		QUARTIC = 32,
		POLYNOMIAL = 33,
		MORSE = 34,
		GHOSTBEND = 35,
		UREYBRADLEY = 36,
		COSINE = 37,
		GHOSTTORSION = 38,
		CHARMM = 39,
		OPLS = 40,
		TRAPPE = 41,
		AMBERIMPROPER = 42,
		IMPROPERCOSINE = 43,
		CENTRALATOMHEIGHT = 44,
		DREIDING = 45,
		CHARGE = 46,
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
#endif /*INC_MDTokenTypes_hpp_*/
