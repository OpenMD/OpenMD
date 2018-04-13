#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.7 (20171128): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

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
		ANALYZER = 26,
		NUDGEDELASTICBAND = 27,
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
