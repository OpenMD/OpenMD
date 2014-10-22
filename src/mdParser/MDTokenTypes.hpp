#ifndef INC_MDTokenTypes_hpp_
#define INC_MDTokenTypes_hpp_

/* $ANTLR 2.7.7 (20131114): "MDParser.g" -> "MDTokenTypes.hpp"$ */

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
		ENDBLOCK = 26,
		ID = 27,
		ASSIGNEQUAL = 28,
		SEMICOLON = 29,
		StringLiteral = 30,
		LCURLY = 31,
		RCURLY = 32,
		LBRACKET = 33,
		RBRACKET = 34,
		LPAREN = 35,
		RPAREN = 36,
		COMMA = 37,
		NUM_INT = 38,
		NUM_LONG = 39,
		NUM_FLOAT = 40,
		NUM_DOUBLE = 41,
		DOT = 42,
		COLON = 43,
		QUESTIONMARK = 44,
		Whitespace = 45,
		Comment = 46,
		CPPComment = 47,
		PREPROC_DIRECTIVE = 48,
		LineDirective = 49,
		Space = 50,
		CharLiteral = 51,
		EndOfLine = 52,
		Escape = 53,
		Vocabulary = 54,
		Digit = 55,
		Decimal = 56,
		HEX_DIGIT = 57,
		EXPONENT = 58,
		FLOAT_SUFFIX = 59,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTokenTypes_hpp_*/
