#ifndef INC_MDParser_hpp_
#define INC_MDParser_hpp_

#include <antlr/config.hpp>
/* $ANTLR 2.7.7 (20160304): "MDParser.g" -> "MDParser.hpp"$ */
#include <antlr/TokenStream.hpp>
#include <antlr/TokenBuffer.hpp>
#include "MDTokenTypes.hpp"
#include <antlr/LLkParser.hpp>

#line 2 "MDParser.g"


#include "antlr/CharScanner.hpp"
#include "utils/StringUtils.hpp"
#include "mdParser/FilenameObserver.hpp"

#line 19 "MDParser.hpp"
class CUSTOM_API MDParser : public ANTLR_USE_NAMESPACE(antlr)LLkParser, public MDTokenTypes
{
#line 1 "MDParser.g"
#line 23 "MDParser.hpp"
public:
	void initializeASTFactory( ANTLR_USE_NAMESPACE(antlr)ASTFactory& factory );
protected:
	MDParser(ANTLR_USE_NAMESPACE(antlr)TokenBuffer& tokenBuf, int k);
public:
	MDParser(ANTLR_USE_NAMESPACE(antlr)TokenBuffer& tokenBuf);
protected:
	MDParser(ANTLR_USE_NAMESPACE(antlr)TokenStream& lexer, int k);
public:
	MDParser(ANTLR_USE_NAMESPACE(antlr)TokenStream& lexer);
	MDParser(const ANTLR_USE_NAMESPACE(antlr)ParserSharedInputState& state);
	int getNumTokens() const
	{
		return MDParser::NUM_TOKENS;
	}
	const char* getTokenName( int type ) const
	{
		if( type > getNumTokens() ) return 0;
		return MDParser::tokenNames[type];
	}
	const char* const* getTokenNames() const
	{
		return MDParser::tokenNames;
	}
	public: void mdfile();
	public: void statement();
	public: void assignment();
	public: void componentblock();
	public: void moleculeblock();
	public: void zconstraintblock();
	public: void restraintblock();
	public: void flucqblock();
	public: void rnemdblock();
	public: void minimizerblock();
	public: void constant();
	protected: void intConst();
	protected: void floatConst();
	protected: void vectorConst();
	public: void moleculestatement();
	public: void atomblock();
	public: void bondblock();
	public: void bendblock();
	public: void torsionblock();
	public: void inversionblock();
	public: void rigidbodyblock();
	public: void cutoffgroupblock();
	public: void fragmentblock();
	public: void constraintblock();
	public: void atomstatement();
	public: void doubleNumberTuple();
	public: void bondstatement();
	public: void inttuple();
	public: void bendstatement();
	public: void torsionstatement();
	public: void inversionstatement();
	public: void rigidbodystatement();
	public: void cutoffgroupstatement();
	public: void fragmentstatement();
	public: void constraintstatement();
	protected: void doubleNumber();
public:
	ANTLR_USE_NAMESPACE(antlr)RefAST getAST()
	{
		return returnAST;
	}
	
protected:
	ANTLR_USE_NAMESPACE(antlr)RefAST returnAST;
private:
	static const char* tokenNames[];
#ifndef NO_STATIC_CONSTS
	static const int NUM_TOKENS = 78;
#else
	enum {
		NUM_TOKENS = 78
	};
#endif
	
	static const unsigned long _tokenSet_0_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_0;
	static const unsigned long _tokenSet_1_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_1;
	static const unsigned long _tokenSet_2_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_2;
	static const unsigned long _tokenSet_3_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_3;
	static const unsigned long _tokenSet_4_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_4;
	static const unsigned long _tokenSet_5_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_5;
	static const unsigned long _tokenSet_6_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_6;
	static const unsigned long _tokenSet_7_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_7;
	static const unsigned long _tokenSet_8_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_8;
	static const unsigned long _tokenSet_9_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_9;
	static const unsigned long _tokenSet_10_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_10;
	static const unsigned long _tokenSet_11_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_11;
	static const unsigned long _tokenSet_12_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_12;
	static const unsigned long _tokenSet_13_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_13;
	static const unsigned long _tokenSet_14_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_14;
	static const unsigned long _tokenSet_15_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_15;
	static const unsigned long _tokenSet_16_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_16;
	static const unsigned long _tokenSet_17_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_17;
	static const unsigned long _tokenSet_18_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_18;
	static const unsigned long _tokenSet_19_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_19;
	static const unsigned long _tokenSet_20_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_20;
	static const unsigned long _tokenSet_21_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_21;
	static const unsigned long _tokenSet_22_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_22;
};

#endif /*INC_MDParser_hpp_*/
