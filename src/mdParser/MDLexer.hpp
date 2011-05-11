#ifndef INC_MDLexer_hpp_
#define INC_MDLexer_hpp_

#include <antlr/config.hpp>
/* $ANTLR 2.7.7 (20101020): "MDParser.g" -> "MDLexer.hpp"$ */
#include <antlr/CommonToken.hpp>
#include <antlr/InputBuffer.hpp>
#include <antlr/BitSet.hpp>
#include "MDTokenTypes.hpp"
#include <antlr/CharScanner.hpp>
#line 2 "MDParser.g"


#include "antlr/CharScanner.hpp"
#include "utils/StringUtils.hpp"
#include "mdParser/FilenameObserver.hpp"

#line 19 "MDLexer.hpp"
class CUSTOM_API MDLexer : public ANTLR_USE_NAMESPACE(antlr)CharScanner, public MDTokenTypes
{
#line 186 "MDParser.g"



	int deferredLineCount;
	FilenameObserver* observer;
	
  public:
	void setObserver(FilenameObserver* osv) {observer = osv;}
  void initDeferredLineCount() { deferredLineCount = 0;}
  void deferredNewline() { 
        deferredLineCount++;
  }

	
  virtual void newline() { 
		for (;deferredLineCount>0;deferredLineCount--) {
			CharScanner::newline();
		}
    CharScanner::newline();
  }
    
#line 23 "MDLexer.hpp"
private:
	void initLiterals();
public:
	bool getCaseSensitiveLiterals() const
	{
		return true;
	}
public:
	MDLexer(ANTLR_USE_NAMESPACE(std)istream& in);
	MDLexer(ANTLR_USE_NAMESPACE(antlr)InputBuffer& ib);
	MDLexer(const ANTLR_USE_NAMESPACE(antlr)LexerSharedInputState& state);
	ANTLR_USE_NAMESPACE(antlr)RefToken nextToken();
	public: void mASSIGNEQUAL(bool _createToken);
	public: void mCOLON(bool _createToken);
	public: void mCOMMA(bool _createToken);
	public: void mQUESTIONMARK(bool _createToken);
	public: void mSEMICOLON(bool _createToken);
	public: void mLPAREN(bool _createToken);
	public: void mRPAREN(bool _createToken);
	public: void mLBRACKET(bool _createToken);
	public: void mRBRACKET(bool _createToken);
	public: void mLCURLY(bool _createToken);
	public: void mRCURLY(bool _createToken);
	public: void mWhitespace(bool _createToken);
	public: void mComment(bool _createToken);
	protected: void mEndOfLine(bool _createToken);
	public: void mCPPComment(bool _createToken);
	public: void mPREPROC_DIRECTIVE(bool _createToken);
	protected: void mLineDirective(bool _createToken);
	protected: void mSpace(bool _createToken);
	protected: void mDecimal(bool _createToken);
	public: void mStringLiteral(bool _createToken);
	public: void mCharLiteral(bool _createToken);
	protected: void mEscape(bool _createToken);
	protected: void mDigit(bool _createToken);
	protected: void mVocabulary(bool _createToken);
	public: void mID(bool _createToken);
	protected: void mHEX_DIGIT(bool _createToken);
	public: void mNUM_INT(bool _createToken);
	protected: void mEXPONENT(bool _createToken);
	protected: void mFLOAT_SUFFIX(bool _createToken);
private:
	
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
};

#endif /*INC_MDLexer_hpp_*/
