#ifndef INC_MDTokenTypes_hpp_
#define INC_MDTokenTypes_hpp_

/* $ANTLR 2.7.7 (20120725): "MDParser.g" -> "MDTokenTypes.hpp"$ */

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
		POSITION = 18,
		ORIENTATION = 19,
		FLUCQ = 20,
		RNEMD = 21,
		MINIMIZER = 22,
		ENDBLOCK = 23,
		ID = 24,
		ASSIGNEQUAL = 25,
		SEMICOLON = 26,
		StringLiteral = 27,
		LCURLY = 28,
		RCURLY = 29,
		LBRACKET = 30,
		RBRACKET = 31,
		LPAREN = 32,
		RPAREN = 33,
		COMMA = 34,
		NUM_INT = 35,
		NUM_LONG = 36,
		NUM_FLOAT = 37,
		NUM_DOUBLE = 38,
		DOT = 39,
		COLON = 40,
		QUESTIONMARK = 41,
		Whitespace = 42,
		Comment = 43,
		CPPComment = 44,
		PREPROC_DIRECTIVE = 45,
		LineDirective = 46,
		Space = 47,
		CharLiteral = 48,
		EndOfLine = 49,
		Escape = 50,
		Vocabulary = 51,
		Digit = 52,
		Decimal = 53,
		HEX_DIGIT = 54,
		EXPONENT = 55,
		FLOAT_SUFFIX = 56,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTokenTypes_hpp_*/
