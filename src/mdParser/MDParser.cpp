/* $ANTLR 2.7.7 (20121118): "MDParser.g" -> "MDParser.cpp"$ */
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
#line 75 "MDParser.g"
		tmp6_AST->setType(ENDBLOCK);
#line 194 "MDParser.cpp"
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
				goto _loop27;
			}
			
		}
		_loop27:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp9_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp9_AST);
		match(RCURLY);
#line 93 "MDParser.g"
		tmp9_AST->setType(ENDBLOCK);
#line 234 "MDParser.cpp"
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
#line 78 "MDParser.g"
		tmp12_AST->setType(ENDBLOCK);
#line 274 "MDParser.cpp"
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
#line 81 "MDParser.g"
		tmp15_AST->setType(ENDBLOCK);
#line 314 "MDParser.cpp"
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
#line 84 "MDParser.g"
		tmp18_AST->setType(ENDBLOCK);
#line 354 "MDParser.cpp"
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
#line 87 "MDParser.g"
		tmp21_AST->setType(ENDBLOCK);
#line 394 "MDParser.cpp"
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
#line 90 "MDParser.g"
		tmp24_AST->setType(ENDBLOCK);
#line 434 "MDParser.cpp"
		minimizerblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = minimizerblock_AST;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp25_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp25_AST);
			match(ID);
			constant_AST = currentAST.root;
			break;
		}
		case StringLiteral:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp26_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp26_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp26_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp27_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp27_AST);
			match(NUM_INT);
			intConst_AST = currentAST.root;
			break;
		}
		case NUM_LONG:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp28_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp28_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp29_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp29_AST);
			match(NUM_FLOAT);
			floatConst_AST = currentAST.root;
			break;
		}
		case NUM_DOUBLE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp30_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp30_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp30_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp31_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp31_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp31_AST);
		match(LPAREN);
		doubleNumber();
		astFactory->addASTChild( currentAST, returnAST );
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp32_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp32_AST);
		match(COMMA);
		doubleNumber();
		astFactory->addASTChild( currentAST, returnAST );
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp33_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp33_AST);
		match(COMMA);
		doubleNumber();
		astFactory->addASTChild( currentAST, returnAST );
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp34_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp34_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp34_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp35_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp35_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp35_AST);
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
				goto _loop31;
			}
			
		}
		_loop31:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp39_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp39_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp39_AST);
		match(RCURLY);
#line 107 "MDParser.g"
		tmp39_AST->setType(ENDBLOCK);
#line 735 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp40_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp40_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp40_AST);
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
			if ((LA(1) == MEMBERS || LA(1) == ID)) {
				bondstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop36;
			}
			
		}
		_loop36:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp44_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp44_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp44_AST);
		match(RCURLY);
#line 116 "MDParser.g"
		tmp44_AST->setType(ENDBLOCK);
#line 794 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp45_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp45_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp45_AST);
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
			if ((LA(1) == MEMBERS || LA(1) == ID)) {
				bendstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop41;
			}
			
		}
		_loop41:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp49_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp49_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp49_AST);
		match(RCURLY);
#line 123 "MDParser.g"
		tmp49_AST->setType(ENDBLOCK);
#line 853 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp50_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp50_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp50_AST);
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
			if ((LA(1) == MEMBERS || LA(1) == ID)) {
				torsionstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop46;
			}
			
		}
		_loop46:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp54_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp54_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp54_AST);
		match(RCURLY);
#line 130 "MDParser.g"
		tmp54_AST->setType(ENDBLOCK);
#line 912 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp55_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp55_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp55_AST);
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
			if ((_tokenSet_10.member(LA(1)))) {
				inversionstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop51;
			}
			
		}
		_loop51:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp59_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp59_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp59_AST);
		match(RCURLY);
#line 137 "MDParser.g"
		tmp59_AST->setType(ENDBLOCK);
#line 971 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp60_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp60_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp60_AST);
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
				goto _loop55;
			}
			
		}
		_loop55:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp64_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp64_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp64_AST);
		match(RCURLY);
#line 145 "MDParser.g"
		tmp64_AST->setType(ENDBLOCK);
#line 1015 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp65_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp65_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp65_AST);
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
				goto _loop60;
			}
			
		}
		_loop60:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp69_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp69_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp69_AST);
		match(RCURLY);
#line 152 "MDParser.g"
		tmp69_AST->setType(ENDBLOCK);
#line 1074 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp70_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp70_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp70_AST);
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
				goto _loop64;
			}
			
		}
		_loop64:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp74_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp74_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp74_AST);
		match(RCURLY);
#line 159 "MDParser.g"
		tmp74_AST->setType(ENDBLOCK);
#line 1118 "MDParser.cpp"
		fragmentblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_8);
	}
	returnAST = fragmentblock_AST;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp75_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp75_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp75_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp79_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp79_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp79_AST);
			match(ORIENTATION);
			match(LPAREN);
			doubleNumberTuple();
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
		recover(ex,_tokenSet_11);
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
				goto _loop68;
			}
			
		}
		_loop68:;
		} // ( ... )*
		doubleNumberTuple_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_12);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp84_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp84_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp84_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
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
		recover(ex,_tokenSet_13);
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
				goto _loop71;
			}
			
		}
		_loop71:;
		} // ( ... )*
		inttuple_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_12);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp89_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp89_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp89_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
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
		recover(ex,_tokenSet_13);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp93_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp93_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp93_AST);
			match(MEMBERS);
			match(LPAREN);
			inttuple();
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
		recover(ex,_tokenSet_13);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp97_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp97_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp97_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp101_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp101_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp101_AST);
			match(SATELLITES);
			match(LPAREN);
			inttuple();
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
		recover(ex,_tokenSet_14);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp105_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp105_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp105_AST);
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
		recover(ex,_tokenSet_13);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp109_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp109_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp109_AST);
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
		recover(ex,_tokenSet_13);
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
		recover(ex,_tokenSet_15);
	}
	returnAST = fragmentstatement_AST;
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
		recover(ex,_tokenSet_16);
	}
	returnAST = doubleNumber_AST;
}

void MDParser::initializeASTFactory( ANTLR_USE_NAMESPACE(antlr)ASTFactory& factory )
{
	factory.setMaxNodeType(57);
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
	"\"fragment\"",
	"\"members\"",
	"\"center\"",
	"\"satellites\"",
	"\"position\"",
	"\"orientation\"",
	"\"flucQ\"",
	"\"RNEMD\"",
	"\"minimizer\"",
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

const unsigned long MDParser::_tokenSet_0_data_[] = { 48234736UL, 0UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" "minimizer" 
// ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDParser::_tokenSet_1_data_[] = { 2UL, 0UL, 0UL, 0UL };
// EOF 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDParser::_tokenSet_2_data_[] = { 48234738UL, 0UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" 
// "minimizer" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_2(_tokenSet_2_data_,4);
const unsigned long MDParser::_tokenSet_3_data_[] = { 1124073458UL, 0UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "atom" "bond" "bend" 
// "torsion" "inversion" "rigidBody" "cutoffGroup" "fragment" "members" 
// "center" "satellites" "position" "orientation" "flucQ" "RNEMD" "minimizer" 
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_3(_tokenSet_3_data_,4);
const unsigned long MDParser::_tokenSet_4_data_[] = { 33619712UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_4(_tokenSet_4_data_,4);
const unsigned long MDParser::_tokenSet_5_data_[] = { 134217728UL, 0UL, 0UL, 0UL };
// SEMICOLON 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_5(_tokenSet_5_data_,4);
const unsigned long MDParser::_tokenSet_6_data_[] = { 134217728UL, 13UL, 0UL, 0UL };
// SEMICOLON RBRACKET RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_6(_tokenSet_6_data_,4);
const unsigned long MDParser::_tokenSet_7_data_[] = { 134217728UL, 12UL, 0UL, 0UL };
// SEMICOLON RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_7(_tokenSet_7_data_,4);
const unsigned long MDParser::_tokenSet_8_data_[] = { 1107361536UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_8(_tokenSet_8_data_,4);
const unsigned long MDParser::_tokenSet_9_data_[] = { 35127296UL, 0UL, 0UL, 0UL };
// "position" "orientation" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_9(_tokenSet_9_data_,4);
const unsigned long MDParser::_tokenSet_10_data_[] = { 33947648UL, 0UL, 0UL, 0UL };
// "center" "satellites" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_10(_tokenSet_10_data_,4);
const unsigned long MDParser::_tokenSet_11_data_[] = { 1108869120UL, 0UL, 0UL, 0UL };
// "position" "orientation" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_11(_tokenSet_11_data_,4);
const unsigned long MDParser::_tokenSet_12_data_[] = { 0UL, 4UL, 0UL, 0UL };
// RPAREN 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_12(_tokenSet_12_data_,4);
const unsigned long MDParser::_tokenSet_13_data_[] = { 1107361792UL, 0UL, 0UL, 0UL };
// "members" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_13(_tokenSet_13_data_,4);
const unsigned long MDParser::_tokenSet_14_data_[] = { 1107689472UL, 0UL, 0UL, 0UL };
// "center" "satellites" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_14(_tokenSet_14_data_,4);
const unsigned long MDParser::_tokenSet_15_data_[] = { 1107296256UL, 0UL, 0UL, 0UL };
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_15(_tokenSet_15_data_,4);
const unsigned long MDParser::_tokenSet_16_data_[] = { 0UL, 12UL, 0UL, 0UL };
// RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_16(_tokenSet_16_data_,4);


