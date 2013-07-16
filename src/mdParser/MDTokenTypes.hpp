#ifndef INC_MDTokenTypes_hpp_
#define INC_MDTokenTypes_hpp_

/* $ANTLR 2.7.7 (20121118): "MDParser.g" -> "MDTokenTypes.hpp"$ */

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
		FRAGMENT = 15,
		MEMBERS = 16,
		CENTER = 17,
		SATELLITES = 18,
		POSITION = 19,
		ORIENTATION = 20,
		FLUCQ = 21,
		RNEMD = 22,
		MINIMIZER = 23,
		ENDBLOCK = 24,
		ID = 25,
		ASSIGNEQUAL = 26,
		SEMICOLON = 27,
		StringLiteral = 28,
		LCURLY = 29,
		RCURLY = 30,
		LBRACKET = 31,
		RBRACKET = 32,
		LPAREN = 33,
		RPAREN = 34,
		COMMA = 35,
		NUM_INT = 36,
		NUM_LONG = 37,
		NUM_FLOAT = 38,
		NUM_DOUBLE = 39,
		DOT = 40,
		COLON = 41,
		QUESTIONMARK = 42,
		Whitespace = 43,
		Comment = 44,
		CPPComment = 45,
		PREPROC_DIRECTIVE = 46,
		LineDirective = 47,
		Space = 48,
		CharLiteral = 49,
		EndOfLine = 50,
		Escape = 51,
		Vocabulary = 52,
		Digit = 53,
		Decimal = 54,
		HEX_DIGIT = 55,
		EXPONENT = 56,
		FLOAT_SUFFIX = 57,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTokenTypes_hpp_*/
