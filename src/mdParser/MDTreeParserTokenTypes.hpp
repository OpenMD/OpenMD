#ifndef INC_MDTreeParserTokenTypes_hpp_
#define INC_MDTreeParserTokenTypes_hpp_

/* $ANTLR 2.7.5 (20050406): "MDTreeParser.g" -> "MDTreeParserTokenTypes.hpp"$ */

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
		OCTALINT = 29,
		DECIMALINT = 30,
		HEXADECIMALINT = 31,
		FLOATONE = 32,
		FLOATTWO = 33,
		DOT = 34,
		COLON = 35,
		QUESTIONMARK = 36,
		Whitespace = 37,
		Comment = 38,
		CPPComment = 39,
		PREPROC_DIRECTIVE = 40,
		LineDirective = 41,
		Space = 42,
		CharLiteral = 43,
		EndOfLine = 44,
		Escape = 45,
		Digit = 46,
		Decimal = 47,
		LongSuffix = 48,
		UnsignedSuffix = 49,
		FloatSuffix = 50,
		Exponent = 51,
		Vocabulary = 52,
		Number = 53,
		MINUS = 54,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
