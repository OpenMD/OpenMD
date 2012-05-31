#ifndef INC_MDTreeParser_hpp_
#define INC_MDTreeParser_hpp_

#include <antlr/config.hpp>
#include "MDTreeParserTokenTypes.hpp"
/* $ANTLR 2.7.7 (20110725): "MDTreeParser.g" -> "MDTreeParser.hpp"$ */
#include <antlr/TreeParser.hpp>

#line 2 "MDTreeParser.g"

#include <stack>
#include "io/Globals.hpp"
#include "utils/StringUtils.hpp"
using namespace std;
using namespace OpenMD;

#line 18 "MDTreeParser.hpp"
class CUSTOM_API MDTreeParser : public ANTLR_USE_NAMESPACE(antlr)TreeParser, public MDTreeParserTokenTypes
{
#line 21 "MDTreeParser.g"

  public:
    Globals* walkTree(ANTLR_USE_NAMESPACE(antlr)RefAST tree)
    {
      currConf = new Globals;
      blockStack.push(currConf);
      mdfile(tree);
      return currConf;
    }
  private:
    Globals* currConf;
    stack<DataHolder*> blockStack;    
#line 22 "MDTreeParser.hpp"
public:
	MDTreeParser();
	static void initializeASTFactory( ANTLR_USE_NAMESPACE(antlr)ASTFactory& factory );
	int getNumTokens() const
	{
		return MDTreeParser::NUM_TOKENS;
	}
	const char* getTokenName( int type ) const
	{
		if( type > getNumTokens() ) return 0;
		return MDTreeParser::tokenNames[type];
	}
	const char* const* getTokenNames() const
	{
		return MDTreeParser::tokenNames;
	}
	public: void mdfile(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void statement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void assignment(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void componentblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void moleculeblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void zconstraintblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void restraintblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void flucqblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void rnemdblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void constant(ANTLR_USE_NAMESPACE(antlr)RefAST _t,
		ANTLR_USE_NAMESPACE(antlr)RefAST id
	);
	protected: int  intConst(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	protected: RealType  floatConst(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void moleculestatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void atomblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void bondblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void bendblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void torsionblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void inversionblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void rigidbodyblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void cutoffgroupblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void fragmentblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void atomstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: vector<RealType>  doubleNumberTuple(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void bondstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: vector<int>  inttuple(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void bendstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void torsionstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void inversionstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void rigidbodystatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void cutoffgroupstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	public: void fragmentstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
	protected: RealType  doubleNumber(ANTLR_USE_NAMESPACE(antlr)RefAST _t);
public:
	ANTLR_USE_NAMESPACE(antlr)RefAST getAST()
	{
		return returnAST;
	}
	
protected:
	ANTLR_USE_NAMESPACE(antlr)RefAST returnAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST _retTree;
private:
	static const char* tokenNames[];
#ifndef NO_STATIC_CONSTS
	static const int NUM_TOKENS = 56;
#else
	enum {
		NUM_TOKENS = 56
	};
#endif
	
	static const unsigned long _tokenSet_0_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_0;
	static const unsigned long _tokenSet_1_data_[];
	static const ANTLR_USE_NAMESPACE(antlr)BitSet _tokenSet_1;
};

#endif /*INC_MDTreeParser_hpp_*/
