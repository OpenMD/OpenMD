#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.7 (20101020): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

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
		FRAGMENT = 15,
		MEMBERS = 16,
		CENTER = 17,
		POSITION = 18,
		ORIENTATION = 19,
		ENDBLOCK = 20,
		ID = 21,
		ASSIGNEQUAL = 22,
		SEMICOLON = 23,
		StringLiteral = 24,
		LCURLY = 25,
		RCURLY = 26,
		LBRACKET = 27,
		RBRACKET = 28,
		LPAREN = 29,
		RPAREN = 30,
		COMMA = 31,
		NUM_INT = 32,
		NUM_LONG = 33,
		NUM_FLOAT = 34,
		NUM_DOUBLE = 35,
		DOT = 36,
		COLON = 37,
		QUESTIONMARK = 38,
		Whitespace = 39,
		Comment = 40,
		CPPComment = 41,
		PREPROC_DIRECTIVE = 42,
		LineDirective = 43,
		Space = 44,
		CharLiteral = 45,
		EndOfLine = 46,
		Escape = 47,
		Vocabulary = 48,
		Digit = 49,
		Decimal = 50,
		HEX_DIGIT = 51,
		EXPONENT = 52,
		FLOAT_SUFFIX = 53,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
