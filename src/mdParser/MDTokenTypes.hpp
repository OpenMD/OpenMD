#ifndef INC_MDTokenTypes_hpp_
#define INC_MDTokenTypes_hpp_

/* $ANTLR 2.7.5 (20050406): "MDParser.g" -> "MDTokenTypes.hpp"$ */

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
		ATOM = 7,
		BOND = 8,
		BEND = 9,
		TORSION = 10,
		RIGIDBODY = 11,
		CUTOFFGROUP = 12,
		FRAGMENT = 13,
		MEMBERS = 14,
		POSITION = 15,
		ORIENTATION = 16,
		ENDBLOCK = 17,
		ID = 18,
		ASSIGNEQUAL = 19,
		SEMICOLON = 20,
		StringLiteral = 21,
		LCURLY = 22,
		RCURLY = 23,
		LBRACKET = 24,
		RBRACKET = 25,
		LPAREN = 26,
		RPAREN = 27,
		COMMA = 28,
		NUM_INT = 29,
		NUM_LONG = 30,
		NUM_FLOAT = 31,
		NUM_DOUBLE = 32,
		DOT = 33,
		COLON = 34,
		QUESTIONMARK = 35,
		Whitespace = 36,
		Comment = 37,
		CPPComment = 38,
		PREPROC_DIRECTIVE = 39,
		LineDirective = 40,
		Space = 41,
		CharLiteral = 42,
		EndOfLine = 43,
		Escape = 44,
		Vocabulary = 45,
		Digit = 46,
		Decimal = 47,
		HEX_DIGIT = 48,
		EXPONENT = 49,
		FLOAT_SUFFIX = 50,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTokenTypes_hpp_*/
