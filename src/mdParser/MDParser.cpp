/* $ANTLR 2.7.7 (20171128): "MDParser.g" -> "MDParser.cpp"$ */
#include "MDParser.hpp"
#include <antlr/NoViableAltException.hpp>
#include <antlr/SemanticException.hpp>
#include <antlr/ASTFactory.hpp>
#line 1 "MDParser.g"
#line 8 "MDParser.cpp"
MDParser::MDParser(ANTLR_USE_NAMESPACE(antlr)TokenBuffer& tokenBuf, int k)
: ANTLR_USE_NAMESPACE(antlr)LLkParser(tokenBuf,k)
{
}

MDParser::MDParser(ANTLR_USE_NAMESPACE(antlr)TokenBuffer& tokenBuf)
: ANTLR_USE_NAMESPACE(antlr)LLkParser(tokenBuf,3)
{
}

MDParser::MDParser(ANTLR_USE_NAMESPACE(antlr)TokenStream& lexer, int k)
: ANTLR_USE_NAMESPACE(antlr)LLkParser(lexer,k)
{
}

MDParser::MDParser(ANTLR_USE_NAMESPACE(antlr)TokenStream& lexer)
: ANTLR_USE_NAMESPACE(antlr)LLkParser(lexer,3)
{
}

MDParser::MDParser(const ANTLR_USE_NAMESPACE(antlr)ParserSharedInputState& state)
: ANTLR_USE_NAMESPACE(antlr)LLkParser(state,3)
{
}

void MDParser::mdfile() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST mdfile_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_0.member(LA(1)))) {
				statement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop3;
			}
			
		}
		_loop3:;
		} // ( ... )*
		mdfile_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_1);
	}
	returnAST = mdfile_AST;
}

void MDParser::statement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST statement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case COMPONENT:
		{
			componentblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case MOLECULE:
		{
			moleculeblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case ZCONSTRAINT:
		{
			zconstraintblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case RESTRAINT:
		{
			restraintblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case FLUCQ:
		{
			flucqblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case RNEMD:
		{
			rnemdblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case MINIMIZER:
		{
			minimizerblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		case ANALYZER:
		{
			analyzerblock();
			astFactory->addASTChild( currentAST, returnAST );
			statement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = statement_AST;
}

void MDParser::assignment() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST assignment_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp1_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp1_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp1_AST);
		match(ID);
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp2_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp2_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp2_AST);
		match(ASSIGNEQUAL);
		constant();
		astFactory->addASTChild( currentAST, returnAST );
		match(SEMICOLON);
		assignment_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_3);
	}
	returnAST = assignment_AST;
}

void MDParser::componentblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST componentblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp4_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp4_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp4_AST);
		match(COMPONENT);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop9;
			}
			
		}
		_loop9:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp6_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp6_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp6_AST);
		match(RCURLY);
#line 97 "MDParser.g"
		tmp6_AST->setType(ENDBLOCK);
#line 201 "MDParser.cpp"
		componentblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = componentblock_AST;
}

void MDParser::moleculeblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculeblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp7_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp7_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp7_AST);
		match(MOLECULE);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_4.member(LA(1)))) {
				moleculestatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop30;
			}
			
		}
		_loop30:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp9_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp9_AST);
		match(RCURLY);
#line 118 "MDParser.g"
		tmp9_AST->setType(ENDBLOCK);
#line 241 "MDParser.cpp"
		moleculeblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = moleculeblock_AST;
}

void MDParser::zconstraintblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST zconstraintblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp10_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp10_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp10_AST);
		match(ZCONSTRAINT);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop12;
			}
			
		}
		_loop12:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp12_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp12_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp12_AST);
		match(RCURLY);
#line 100 "MDParser.g"
		tmp12_AST->setType(ENDBLOCK);
#line 281 "MDParser.cpp"
		zconstraintblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = zconstraintblock_AST;
}

void MDParser::restraintblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST restraintblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp13_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp13_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp13_AST);
		match(RESTRAINT);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop15;
			}
			
		}
		_loop15:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp15_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp15_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp15_AST);
		match(RCURLY);
#line 103 "MDParser.g"
		tmp15_AST->setType(ENDBLOCK);
#line 321 "MDParser.cpp"
		restraintblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = restraintblock_AST;
}

void MDParser::flucqblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST flucqblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp16_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp16_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp16_AST);
		match(FLUCQ);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop18;
			}
			
		}
		_loop18:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp18_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp18_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp18_AST);
		match(RCURLY);
#line 106 "MDParser.g"
		tmp18_AST->setType(ENDBLOCK);
#line 361 "MDParser.cpp"
		flucqblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = flucqblock_AST;
}

void MDParser::rnemdblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST rnemdblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp19_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp19_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp19_AST);
		match(RNEMD);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop21;
			}
			
		}
		_loop21:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp21_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp21_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp21_AST);
		match(RCURLY);
#line 109 "MDParser.g"
		tmp21_AST->setType(ENDBLOCK);
#line 401 "MDParser.cpp"
		rnemdblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = rnemdblock_AST;
}

void MDParser::minimizerblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST minimizerblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp22_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp22_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp22_AST);
		match(MINIMIZER);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop24;
			}
			
		}
		_loop24:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp24_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp24_AST);
		match(RCURLY);
#line 112 "MDParser.g"
		tmp24_AST->setType(ENDBLOCK);
#line 441 "MDParser.cpp"
		minimizerblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = minimizerblock_AST;
}

void MDParser::analyzerblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST analyzerblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp25_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp25_AST);
		match(ANALYZER);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				assignment();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop27;
			}
			
		}
		_loop27:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp27_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp27_AST);
		match(RCURLY);
#line 115 "MDParser.g"
		tmp27_AST->setType(ENDBLOCK);
#line 481 "MDParser.cpp"
		analyzerblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = analyzerblock_AST;
}

void MDParser::constant() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST constant_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case NUM_INT:
		case NUM_LONG:
		{
			intConst();
			astFactory->addASTChild( currentAST, returnAST );
			constant_AST = currentAST.root;
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			floatConst();
			astFactory->addASTChild( currentAST, returnAST );
			constant_AST = currentAST.root;
			break;
		}
		case LPAREN:
		{
			vectorConst();
			astFactory->addASTChild( currentAST, returnAST );
			constant_AST = currentAST.root;
			break;
		}
		case ID:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp28_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp28_AST);
			match(ID);
			constant_AST = currentAST.root;
			break;
		}
		case StringLiteral:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp29_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp29_AST);
			match(StringLiteral);
			constant_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_5);
	}
	returnAST = constant_AST;
}

void MDParser::intConst() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST intConst_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case NUM_INT:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp30_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp30_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp30_AST);
			match(NUM_INT);
			intConst_AST = currentAST.root;
			break;
		}
		case NUM_LONG:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp31_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp31_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp31_AST);
			match(NUM_LONG);
			intConst_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_6);
	}
	returnAST = intConst_AST;
}

void MDParser::floatConst() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST floatConst_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case NUM_FLOAT:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp32_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp32_AST);
			match(NUM_FLOAT);
			floatConst_AST = currentAST.root;
			break;
		}
		case NUM_DOUBLE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp33_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp33_AST);
			match(NUM_DOUBLE);
			floatConst_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_7);
	}
	returnAST = floatConst_AST;
}

void MDParser::vectorConst() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST vectorConst_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp34_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp34_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp34_AST);
		match(LPAREN);
		doubleNumberTuple();
		astFactory->addASTChild( currentAST, returnAST );
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp35_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp35_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp35_AST);
		match(RPAREN);
		vectorConst_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_5);
	}
	returnAST = vectorConst_AST;
}

void MDParser::moleculestatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculestatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case ATOM:
		{
			atomblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case BOND:
		{
			bondblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case BEND:
		{
			bendblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case TORSION:
		{
			torsionblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case INVERSION:
		{
			inversionblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case RIGIDBODY:
		{
			rigidbodyblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case CUTOFFGROUP:
		{
			cutoffgroupblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case FRAGMENT:
		{
			fragmentblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		case CONSTRAINT:
		{
			constraintblock();
			astFactory->addASTChild( currentAST, returnAST );
			moleculestatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = moleculestatement_AST;
}

void MDParser::atomblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST atomblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp36_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp36_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp36_AST);
		match(ATOM);
		match(LBRACKET);
		intConst();
		astFactory->addASTChild( currentAST, returnAST );
		match(RBRACKET);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_9.member(LA(1)))) {
				atomstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop34;
			}
			
		}
		_loop34:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp40_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp40_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp40_AST);
		match(RCURLY);
#line 133 "MDParser.g"
		tmp40_AST->setType(ENDBLOCK);
#line 777 "MDParser.cpp"
		atomblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = atomblock_AST;
}

void MDParser::bondblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST bondblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp41_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp41_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp41_AST);
		match(BOND);
		{
		switch ( LA(1)) {
		case LBRACKET:
		{
			match(LBRACKET);
			intConst();
			match(RBRACKET);
			break;
		}
		case LCURLY:
		{
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_10.member(LA(1)))) {
				bondstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop39;
			}
			
		}
		_loop39:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp45_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp45_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp45_AST);
		match(RCURLY);
#line 143 "MDParser.g"
		tmp45_AST->setType(ENDBLOCK);
#line 836 "MDParser.cpp"
		bondblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = bondblock_AST;
}

void MDParser::bendblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST bendblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp46_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp46_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp46_AST);
		match(BEND);
		{
		switch ( LA(1)) {
		case LBRACKET:
		{
			match(LBRACKET);
			intConst();
			match(RBRACKET);
			break;
		}
		case LCURLY:
		{
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_11.member(LA(1)))) {
				bendstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop44;
			}
			
		}
		_loop44:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp50_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp50_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp50_AST);
		match(RCURLY);
#line 156 "MDParser.g"
		tmp50_AST->setType(ENDBLOCK);
#line 895 "MDParser.cpp"
		bendblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = bendblock_AST;
}

void MDParser::torsionblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST torsionblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp51_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp51_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp51_AST);
		match(TORSION);
		{
		switch ( LA(1)) {
		case LBRACKET:
		{
			match(LBRACKET);
			intConst();
			match(RBRACKET);
			break;
		}
		case LCURLY:
		{
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_12.member(LA(1)))) {
				torsionstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop49;
			}
			
		}
		_loop49:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp55_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp55_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp55_AST);
		match(RCURLY);
#line 170 "MDParser.g"
		tmp55_AST->setType(ENDBLOCK);
#line 954 "MDParser.cpp"
		torsionblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = torsionblock_AST;
}

void MDParser::inversionblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST inversionblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp56_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp56_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp56_AST);
		match(INVERSION);
		{
		switch ( LA(1)) {
		case LBRACKET:
		{
			match(LBRACKET);
			intConst();
			match(RBRACKET);
			break;
		}
		case LCURLY:
		{
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((_tokenSet_13.member(LA(1)))) {
				inversionstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop54;
			}
			
		}
		_loop54:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp60_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp60_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp60_AST);
		match(RCURLY);
#line 185 "MDParser.g"
		tmp60_AST->setType(ENDBLOCK);
#line 1013 "MDParser.cpp"
		inversionblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = inversionblock_AST;
}

void MDParser::rigidbodyblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST rigidbodyblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp61_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp61_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp61_AST);
		match(RIGIDBODY);
		match(LBRACKET);
		intConst();
		astFactory->addASTChild( currentAST, returnAST );
		match(RBRACKET);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == MEMBERS || LA(1) == ID)) {
				rigidbodystatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop58;
			}
			
		}
		_loop58:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp65_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp65_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp65_AST);
		match(RCURLY);
#line 198 "MDParser.g"
		tmp65_AST->setType(ENDBLOCK);
#line 1057 "MDParser.cpp"
		rigidbodyblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = rigidbodyblock_AST;
}

void MDParser::cutoffgroupblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST cutoffgroupblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp66_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp66_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp66_AST);
		match(CUTOFFGROUP);
		{
		switch ( LA(1)) {
		case LBRACKET:
		{
			match(LBRACKET);
			intConst();
			match(RBRACKET);
			break;
		}
		case LCURLY:
		{
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == MEMBERS || LA(1) == ID)) {
				cutoffgroupstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop63;
			}
			
		}
		_loop63:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp70_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp70_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp70_AST);
		match(RCURLY);
#line 205 "MDParser.g"
		tmp70_AST->setType(ENDBLOCK);
#line 1116 "MDParser.cpp"
		cutoffgroupblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = cutoffgroupblock_AST;
}

void MDParser::fragmentblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp71_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp71_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp71_AST);
		match(FRAGMENT);
		match(LBRACKET);
		intConst();
		astFactory->addASTChild( currentAST, returnAST );
		match(RBRACKET);
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == ID)) {
				fragmentstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop67;
			}
			
		}
		_loop67:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp75_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp75_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp75_AST);
		match(RCURLY);
#line 212 "MDParser.g"
		tmp75_AST->setType(ENDBLOCK);
#line 1160 "MDParser.cpp"
		fragmentblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = fragmentblock_AST;
}

void MDParser::constraintblock() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST constraintblock_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp76_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp76_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp76_AST);
		match(CONSTRAINT);
		{
		switch ( LA(1)) {
		case LBRACKET:
		{
			match(LBRACKET);
			intConst();
			match(RBRACKET);
			break;
		}
		case LCURLY:
		{
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		match(LCURLY);
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == MEMBERS || LA(1) == ID)) {
				constraintstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop72;
			}
			
		}
		_loop72:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp80_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp80_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp80_AST);
		match(RCURLY);
#line 218 "MDParser.g"
		tmp80_AST->setType(ENDBLOCK);
#line 1219 "MDParser.cpp"
		constraintblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = constraintblock_AST;
}

void MDParser::atomstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST atomstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			atomstatement_AST = currentAST.root;
			break;
		}
		case POSITION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp81_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp81_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp81_AST);
			match(POSITION);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			atomstatement_AST = currentAST.root;
			break;
		}
		case ORIENTATION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp85_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp85_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp85_AST);
			match(ORIENTATION);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			atomstatement_AST = currentAST.root;
			break;
		}
		case CHARGE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp89_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp89_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp89_AST);
			match(CHARGE);
			match(LPAREN);
			floatConst();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			atomstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_14);
	}
	returnAST = atomstatement_AST;
}

void MDParser::doubleNumberTuple() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST doubleNumberTuple_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		doubleNumber();
		astFactory->addASTChild( currentAST, returnAST );
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == COMMA)) {
				match(COMMA);
				doubleNumber();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop76;
			}
			
		}
		_loop76:;
		} // ( ... )*
		doubleNumberTuple_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_15);
	}
	returnAST = doubleNumberTuple_AST;
}

void MDParser::bondstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST bondstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			bondstatement_AST = currentAST.root;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp94_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp94_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp94_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		case FIXED:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp98_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp98_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp98_AST);
			match(FIXED);
			match(LPAREN);
			floatConst();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp102_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp102_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp102_AST);
			match(HARMONIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		case CUBIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp106_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp106_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp106_AST);
			match(CUBIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		case QUARTIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp110_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp110_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp110_AST);
			match(QUARTIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		case POLYNOMIAL:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp114_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp114_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp114_AST);
			match(POLYNOMIAL);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		case MORSE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp118_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp118_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp118_AST);
			match(MORSE);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bondstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_16);
	}
	returnAST = bondstatement_AST;
}

void MDParser::inttuple() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST inttuple_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		intConst();
		astFactory->addASTChild( currentAST, returnAST );
		{ // ( ... )*
		for (;;) {
			if ((LA(1) == COMMA)) {
				match(COMMA);
				intConst();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop79;
			}
			
		}
		_loop79:;
		} // ( ... )*
		inttuple_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_15);
	}
	returnAST = inttuple_AST;
}

void MDParser::bendstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST bendstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			bendstatement_AST = currentAST.root;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp123_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp123_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp123_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp127_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp127_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp127_AST);
			match(HARMONIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case GHOSTBEND:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp131_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp131_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp131_AST);
			match(GHOSTBEND);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case UREYBRADLEY:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp135_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp135_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp135_AST);
			match(UREYBRADLEY);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case CUBIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp139_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp139_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp139_AST);
			match(CUBIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case QUARTIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp143_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp143_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp143_AST);
			match(QUARTIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case POLYNOMIAL:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp147_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp147_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp147_AST);
			match(POLYNOMIAL);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		case COSINE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp151_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp151_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp151_AST);
			match(COSINE);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			bendstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_17);
	}
	returnAST = bendstatement_AST;
}

void MDParser::torsionstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST torsionstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			torsionstatement_AST = currentAST.root;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp155_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp155_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp155_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case GHOSTTORSION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp159_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp159_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp159_AST);
			match(GHOSTTORSION);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case CUBIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp163_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp163_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp163_AST);
			match(CUBIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case QUARTIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp167_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp167_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp167_AST);
			match(QUARTIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case POLYNOMIAL:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp171_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp171_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp171_AST);
			match(POLYNOMIAL);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case CHARMM:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp175_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp175_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp175_AST);
			match(CHARMM);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case OPLS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp179_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp179_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp179_AST);
			match(OPLS);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case TRAPPE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp183_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp183_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp183_AST);
			match(TRAPPE);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp187_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp187_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp187_AST);
			match(HARMONIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			torsionstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_18);
	}
	returnAST = torsionstatement_AST;
}

void MDParser::inversionstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST inversionstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			inversionstatement_AST = currentAST.root;
			break;
		}
		case CENTER:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp191_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp191_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp191_AST);
			match(CENTER);
			match(LPAREN);
			intConst();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		case SATELLITES:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp195_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp195_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp195_AST);
			match(SATELLITES);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		case AMBERIMPROPER:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp199_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp199_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp199_AST);
			match(AMBERIMPROPER);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		case IMPROPERCOSINE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp203_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp203_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp203_AST);
			match(IMPROPERCOSINE);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp207_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp207_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp207_AST);
			match(HARMONIC);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		case CENTRALATOMHEIGHT:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp211_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp211_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp211_AST);
			match(CENTRALATOMHEIGHT);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		case DREIDING:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp215_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp215_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp215_AST);
			match(DREIDING);
			match(LPAREN);
			doubleNumberTuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			inversionstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_19);
	}
	returnAST = inversionstatement_AST;
}

void MDParser::rigidbodystatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST rigidbodystatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			rigidbodystatement_AST = currentAST.root;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp219_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp219_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp219_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			rigidbodystatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_20);
	}
	returnAST = rigidbodystatement_AST;
}

void MDParser::cutoffgroupstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST cutoffgroupstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			cutoffgroupstatement_AST = currentAST.root;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp223_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp223_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp223_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			cutoffgroupstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_20);
	}
	returnAST = cutoffgroupstatement_AST;
}

void MDParser::fragmentstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		assignment();
		astFactory->addASTChild( currentAST, returnAST );
		fragmentstatement_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_21);
	}
	returnAST = fragmentstatement_AST;
}

void MDParser::constraintstatement() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST constraintstatement_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		switch ( LA(1)) {
		case ID:
		{
			assignment();
			astFactory->addASTChild( currentAST, returnAST );
			constraintstatement_AST = currentAST.root;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp227_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp227_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp227_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
			astFactory->addASTChild( currentAST, returnAST );
			match(RPAREN);
			match(SEMICOLON);
			constraintstatement_AST = currentAST.root;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_20);
	}
	returnAST = constraintstatement_AST;
}

void MDParser::doubleNumber() {
	returnAST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)ASTPair currentAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST doubleNumber_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		{
		switch ( LA(1)) {
		case NUM_INT:
		case NUM_LONG:
		{
			intConst();
			astFactory->addASTChild( currentAST, returnAST );
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			floatConst();
			astFactory->addASTChild( currentAST, returnAST );
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(LT(1), getFilename());
		}
		}
		}
		doubleNumber_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_22);
	}
	returnAST = doubleNumber_AST;
}

void MDParser::initializeASTFactory( ANTLR_USE_NAMESPACE(antlr)ASTFactory& factory )
{
	factory.setMaxNodeType(78);
}
const char* MDParser::tokenNames[] = {
	"<0>",
	"EOF",
	"<2>",
	"NULL_TREE_LOOKAHEAD",
	"\"component\"",
	"\"molecule\"",
	"\"zconstraint\"",
	"\"restraint\"",
	"\"atom\"",
	"\"bond\"",
	"\"bend\"",
	"\"torsion\"",
	"\"inversion\"",
	"\"rigidBody\"",
	"\"cutoffGroup\"",
	"\"constraint\"",
	"\"distance\"",
	"\"fragment\"",
	"\"members\"",
	"\"center\"",
	"\"satellites\"",
	"\"position\"",
	"\"orientation\"",
	"\"flucQ\"",
	"\"RNEMD\"",
	"\"minimizer\"",
	"\"analyzer\"",
	"\"Fixed\"",
	"\"Harmonic\"",
	"\"Cubic\"",
	"\"Quartic\"",
	"\"Polynomial\"",
	"\"Morse\"",
	"\"GhostBend\"",
	"\"UreyBradley\"",
	"\"Cosine\"",
	"\"GhostTorsion\"",
	"\"Charmm\"",
	"\"Opls\"",
	"\"Trappe\"",
	"\"AmberImproper\"",
	"\"ImproperCosine\"",
	"\"CentralAtomHeight\"",
	"\"Dreiding\"",
	"\"charge\"",
	"ENDBLOCK",
	"ID",
	"ASSIGNEQUAL",
	"SEMICOLON",
	"StringLiteral",
	"LCURLY",
	"RCURLY",
	"LBRACKET",
	"RBRACKET",
	"LPAREN",
	"RPAREN",
	"COMMA",
	"NUM_INT",
	"NUM_LONG",
	"NUM_FLOAT",
	"NUM_DOUBLE",
	"DOT",
	"COLON",
	"QUESTIONMARK",
	"Whitespace",
	"Comment",
	"CPPComment",
	"a line directive",
	"LineDirective",
	"Space",
	"CharLiteral",
	"EndOfLine",
	"Escape",
	"Vocabulary",
	"Digit",
	"Decimal",
	"HEX_DIGIT",
	"EXPONENT",
	"FLOAT_SUFFIX",
	0
};

const unsigned long MDParser::_tokenSet_0_data_[] = { 125829360UL, 16384UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" "minimizer" 
// "analyzer" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDParser::_tokenSet_1_data_[] = { 2UL, 0UL, 0UL, 0UL };
// EOF 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDParser::_tokenSet_2_data_[] = { 125829362UL, 16384UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" 
// "minimizer" "analyzer" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_2(_tokenSet_2_data_,4);
const unsigned long MDParser::_tokenSet_3_data_[] = { 4294901746UL, 548863UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "atom" "bond" "bend" 
// "torsion" "inversion" "rigidBody" "cutoffGroup" "constraint" "fragment" 
// "members" "center" "satellites" "position" "orientation" "flucQ" "RNEMD" 
// "minimizer" "analyzer" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" 
// "Morse" "GhostBend" "UreyBradley" "Cosine" "GhostTorsion" "Charmm" "Opls" 
// "Trappe" "AmberImproper" "ImproperCosine" "CentralAtomHeight" "Dreiding" 
// "charge" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_3(_tokenSet_3_data_,4);
const unsigned long MDParser::_tokenSet_4_data_[] = { 196352UL, 16384UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "constraint" "fragment" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_4(_tokenSet_4_data_,4);
const unsigned long MDParser::_tokenSet_5_data_[] = { 0UL, 65536UL, 0UL, 0UL };
// SEMICOLON 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_5(_tokenSet_5_data_,4);
const unsigned long MDParser::_tokenSet_6_data_[] = { 0UL, 27328512UL, 0UL, 0UL };
// SEMICOLON RBRACKET RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_6(_tokenSet_6_data_,4);
const unsigned long MDParser::_tokenSet_7_data_[] = { 0UL, 25231360UL, 0UL, 0UL };
// SEMICOLON RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_7(_tokenSet_7_data_,4);
const unsigned long MDParser::_tokenSet_8_data_[] = { 196352UL, 540672UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "constraint" "fragment" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_8(_tokenSet_8_data_,4);
const unsigned long MDParser::_tokenSet_9_data_[] = { 6291456UL, 20480UL, 0UL, 0UL };
// "position" "orientation" "charge" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_9(_tokenSet_9_data_,4);
const unsigned long MDParser::_tokenSet_10_data_[] = { 4161011712UL, 16385UL, 0UL, 0UL };
// "members" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" "Morse" 
// ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_10(_tokenSet_10_data_,4);
const unsigned long MDParser::_tokenSet_11_data_[] = { 4026793984UL, 16398UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostBend" "UreyBradley" 
// "Cosine" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_11(_tokenSet_11_data_,4);
const unsigned long MDParser::_tokenSet_12_data_[] = { 4026793984UL, 16624UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostTorsion" "Charmm" 
// "Opls" "Trappe" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_12(_tokenSet_12_data_,4);
const unsigned long MDParser::_tokenSet_13_data_[] = { 270008320UL, 20224UL, 0UL, 0UL };
// "center" "satellites" "Harmonic" "AmberImproper" "ImproperCosine" "CentralAtomHeight" 
// "Dreiding" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_13(_tokenSet_13_data_,4);
const unsigned long MDParser::_tokenSet_14_data_[] = { 6291456UL, 544768UL, 0UL, 0UL };
// "position" "orientation" "charge" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_14(_tokenSet_14_data_,4);
const unsigned long MDParser::_tokenSet_15_data_[] = { 0UL, 8388608UL, 0UL, 0UL };
// RPAREN 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_15(_tokenSet_15_data_,4);
const unsigned long MDParser::_tokenSet_16_data_[] = { 4161011712UL, 540673UL, 0UL, 0UL };
// "members" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" "Morse" 
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_16(_tokenSet_16_data_,4);
const unsigned long MDParser::_tokenSet_17_data_[] = { 4026793984UL, 540686UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostBend" "UreyBradley" 
// "Cosine" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_17(_tokenSet_17_data_,4);
const unsigned long MDParser::_tokenSet_18_data_[] = { 4026793984UL, 540912UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostTorsion" "Charmm" 
// "Opls" "Trappe" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_18(_tokenSet_18_data_,4);
const unsigned long MDParser::_tokenSet_19_data_[] = { 270008320UL, 544512UL, 0UL, 0UL };
// "center" "satellites" "Harmonic" "AmberImproper" "ImproperCosine" "CentralAtomHeight" 
// "Dreiding" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_19(_tokenSet_19_data_,4);
const unsigned long MDParser::_tokenSet_20_data_[] = { 262144UL, 540672UL, 0UL, 0UL };
// "members" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_20(_tokenSet_20_data_,4);
const unsigned long MDParser::_tokenSet_21_data_[] = { 0UL, 540672UL, 0UL, 0UL };
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_21(_tokenSet_21_data_,4);
const unsigned long MDParser::_tokenSet_22_data_[] = { 0UL, 25165824UL, 0UL, 0UL };
// RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_22(_tokenSet_22_data_,4);


