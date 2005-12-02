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
		PLUS = 32,
		MINUS = 33,
		FLOATONE = 34,
		FLOATTWO = 35,
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
		Digit = 48,
		Decimal = 49,
		LongSuffix = 50,
		UnsignedSuffix = 51,
		FloatSuffix = 52,
		Exponent = 53,
		Vocabulary = 54,
		Number = 55,
		NULL_TREE_LOOKAHEAD = 3
	};
#ifdef __cplusplus
};
#endif
#endif /*INC_MDTreeParserTokenTypes_hpp_*/
