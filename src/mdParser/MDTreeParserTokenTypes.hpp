#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.7 (20080702): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

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
		ATOM = 7,
		BOND = 8,
		BEND = 9,
		TORSION = 10,
		INVERSION = 11,
		RIGIDBODY = 12,
		CUTOFFGROUP = 13,
		FRAGMENT = 14,
		MEMBERS = 15,
		CENTER = 16,
		POSITION = 17,
		ORIENTATION = 18,
		ENDBLOCK = 19,
		ID = 20,
		ASSIGNEQUAL = 21,
		SEMICOLON = 22,
		StringLiteral = 23,
		LCURLY = 24,
		RCURLY = 25,
		LBRACKET = 26,
		RBRACKET = 27,
		LPAREN = 28,
		RPAREN = 29,
		COMMA = 30,
		NUM_INT = 31,
		NUM_LONG = 32,
		NUM_FLOAT = 33,
		NUM_DOUBLE = 34,
		DOT = 35,
		COLON = 36,
		QUESTIONMARK = 37,
		Whitespace = 38,
		Comment = 39,
		CPPComment = 40,
		PREPROC_DIRECTIVE = 41,
		LineDirective = 42,
		Space = 43,
		CharLiteral = 44,
		EndOfLine = 45,
		Escape = 46,
		Vocabulary = 47,
		Digit = 48,
		Decimal = 49,
		HEX_DIGIT = 50,
		EXPONENT = 51,
		FLOAT_SUFFIX = 52,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
