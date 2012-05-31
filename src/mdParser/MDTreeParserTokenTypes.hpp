#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.7 (20110725): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

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
		FLUCQ = 20,
		RNEMD = 21,
		ENDBLOCK = 22,
		ID = 23,
		ASSIGNEQUAL = 24,
		SEMICOLON = 25,
		StringLiteral = 26,
		LCURLY = 27,
		RCURLY = 28,
		LBRACKET = 29,
		RBRACKET = 30,
		LPAREN = 31,
		RPAREN = 32,
		COMMA = 33,
		NUM_INT = 34,
		NUM_LONG = 35,
		NUM_FLOAT = 36,
		NUM_DOUBLE = 37,
		DOT = 38,
		COLON = 39,
		QUESTIONMARK = 40,
		Whitespace = 41,
		Comment = 42,
		CPPComment = 43,
		PREPROC_DIRECTIVE = 44,
		LineDirective = 45,
		Space = 46,
		CharLiteral = 47,
		EndOfLine = 48,
		Escape = 49,
		Vocabulary = 50,
		Digit = 51,
		Decimal = 52,
		HEX_DIGIT = 53,
		EXPONENT = 54,
		FLOAT_SUFFIX = 55,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
