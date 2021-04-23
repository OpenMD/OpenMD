/* $ANTLR 2.7.7 (20160304): "MDParser.g" -> "MDParser.cpp"$ */
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
#line 95 "MDParser.g"
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
#line 113 "MDParser.g"
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
#line 98 "MDParser.g"
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
#line 101 "MDParser.g"
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
#line 104 "MDParser.g"
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
#line 107 "MDParser.g"
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
#line 110 "MDParser.g"
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
		doubleNumberTuple();
		astFactory->addASTChild( currentAST, returnAST );
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp32_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp32_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp33_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp33_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp37_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp37_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp37_AST);
		match(RCURLY);
#line 128 "MDParser.g"
		tmp37_AST->setType(ENDBLOCK);
#line 730 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp38_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp38_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp38_AST);
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
				goto _loop36;
			}
			
		}
		_loop36:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp42_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp42_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp42_AST);
		match(RCURLY);
#line 138 "MDParser.g"
		tmp42_AST->setType(ENDBLOCK);
#line 789 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp43_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp43_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp43_AST);
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
				goto _loop41;
			}
			
		}
		_loop41:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp47_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp47_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp47_AST);
		match(RCURLY);
#line 151 "MDParser.g"
		tmp47_AST->setType(ENDBLOCK);
#line 848 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp48_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp48_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp48_AST);
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
				goto _loop46;
			}
			
		}
		_loop46:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp52_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp52_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp52_AST);
		match(RCURLY);
#line 165 "MDParser.g"
		tmp52_AST->setType(ENDBLOCK);
#line 907 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp53_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp53_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp53_AST);
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
				goto _loop51;
			}
			
		}
		_loop51:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp57_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp57_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp57_AST);
		match(RCURLY);
#line 180 "MDParser.g"
		tmp57_AST->setType(ENDBLOCK);
#line 966 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp58_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp58_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp58_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp62_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp62_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp62_AST);
		match(RCURLY);
#line 193 "MDParser.g"
		tmp62_AST->setType(ENDBLOCK);
#line 1010 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp63_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp63_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp63_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp67_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp67_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp67_AST);
		match(RCURLY);
#line 200 "MDParser.g"
		tmp67_AST->setType(ENDBLOCK);
#line 1069 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp68_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp68_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp68_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp72_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp72_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp72_AST);
		match(RCURLY);
#line 207 "MDParser.g"
		tmp72_AST->setType(ENDBLOCK);
#line 1113 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp73_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp73_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp73_AST);
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
				goto _loop69;
			}
			
		}
		_loop69:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp77_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp77_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp77_AST);
		match(RCURLY);
#line 213 "MDParser.g"
		tmp77_AST->setType(ENDBLOCK);
#line 1172 "MDParser.cpp"
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp78_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp78_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp78_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp82_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp82_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp82_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp86_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp86_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp86_AST);
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
				goto _loop73;
			}
			
		}
		_loop73:;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp91_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp91_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp91_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp95_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp95_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp95_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp99_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp99_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp99_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp103_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp103_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp103_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp107_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp107_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp107_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp111_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp111_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp111_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp115_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp115_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp115_AST);
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
				goto _loop76;
			}
			
		}
		_loop76:;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp120_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp120_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp120_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp124_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp124_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp124_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp128_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp128_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp128_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp132_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp132_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp132_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp136_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp136_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp136_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp140_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp140_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp140_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp144_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp144_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp144_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp148_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp148_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp148_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp152_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp152_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp152_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp156_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp156_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp156_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp160_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp160_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp160_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp164_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp164_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp164_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp168_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp168_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp168_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp172_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp172_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp172_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp176_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp176_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp176_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp180_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp180_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp180_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp184_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp184_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp184_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp188_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp188_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp188_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp192_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp192_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp192_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp196_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp196_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp196_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp200_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp200_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp200_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp204_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp204_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp204_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp208_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp208_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp208_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp212_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp212_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp212_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp216_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp216_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp216_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp220_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp220_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp220_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp224_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp224_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp224_AST);
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
	factory.setMaxNodeType(77);
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

const unsigned long MDParser::_tokenSet_0_data_[] = { 58720496UL, 8192UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" "minimizer" 
// ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDParser::_tokenSet_1_data_[] = { 2UL, 0UL, 0UL, 0UL };
// EOF 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDParser::_tokenSet_2_data_[] = { 58720498UL, 8192UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" 
// "minimizer" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_2(_tokenSet_2_data_,4);
const unsigned long MDParser::_tokenSet_3_data_[] = { 4294901746UL, 274431UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "atom" "bond" "bend" 
// "torsion" "inversion" "rigidBody" "cutoffGroup" "constraint" "fragment" 
// "members" "center" "satellites" "position" "orientation" "flucQ" "RNEMD" 
// "minimizer" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" "Morse" 
// "GhostBend" "UreyBradley" "Cosine" "GhostTorsion" "Charmm" "Opls" "Trappe" 
// "AmberImproper" "ImproperCosine" "CentralAtomHeight" "Dreiding" "charge" 
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_3(_tokenSet_3_data_,4);
const unsigned long MDParser::_tokenSet_4_data_[] = { 196352UL, 8192UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "constraint" "fragment" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_4(_tokenSet_4_data_,4);
const unsigned long MDParser::_tokenSet_5_data_[] = { 0UL, 32768UL, 0UL, 0UL };
// SEMICOLON 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_5(_tokenSet_5_data_,4);
const unsigned long MDParser::_tokenSet_6_data_[] = { 0UL, 13664256UL, 0UL, 0UL };
// SEMICOLON RBRACKET RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_6(_tokenSet_6_data_,4);
const unsigned long MDParser::_tokenSet_7_data_[] = { 0UL, 12615680UL, 0UL, 0UL };
// SEMICOLON RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_7(_tokenSet_7_data_,4);
const unsigned long MDParser::_tokenSet_8_data_[] = { 196352UL, 270336UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "constraint" "fragment" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_8(_tokenSet_8_data_,4);
const unsigned long MDParser::_tokenSet_9_data_[] = { 6291456UL, 10240UL, 0UL, 0UL };
// "position" "orientation" "charge" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_9(_tokenSet_9_data_,4);
const unsigned long MDParser::_tokenSet_10_data_[] = { 4228120576UL, 8192UL, 0UL, 0UL };
// "members" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" "Morse" 
// ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_10(_tokenSet_10_data_,4);
const unsigned long MDParser::_tokenSet_11_data_[] = { 2013528064UL, 8199UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostBend" "UreyBradley" 
// "Cosine" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_11(_tokenSet_11_data_,4);
const unsigned long MDParser::_tokenSet_12_data_[] = { 2013528064UL, 8312UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostTorsion" "Charmm" 
// "Opls" "Trappe" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_12(_tokenSet_12_data_,4);
const unsigned long MDParser::_tokenSet_13_data_[] = { 135790592UL, 10112UL, 0UL, 0UL };
// "center" "satellites" "Harmonic" "AmberImproper" "ImproperCosine" "CentralAtomHeight" 
// "Dreiding" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_13(_tokenSet_13_data_,4);
const unsigned long MDParser::_tokenSet_14_data_[] = { 6291456UL, 272384UL, 0UL, 0UL };
// "position" "orientation" "charge" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_14(_tokenSet_14_data_,4);
const unsigned long MDParser::_tokenSet_15_data_[] = { 0UL, 4194304UL, 0UL, 0UL };
// RPAREN 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_15(_tokenSet_15_data_,4);
const unsigned long MDParser::_tokenSet_16_data_[] = { 4228120576UL, 270336UL, 0UL, 0UL };
// "members" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" "Morse" 
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_16(_tokenSet_16_data_,4);
const unsigned long MDParser::_tokenSet_17_data_[] = { 2013528064UL, 270343UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostBend" "UreyBradley" 
// "Cosine" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_17(_tokenSet_17_data_,4);
const unsigned long MDParser::_tokenSet_18_data_[] = { 2013528064UL, 270456UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostTorsion" "Charmm" 
// "Opls" "Trappe" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_18(_tokenSet_18_data_,4);
const unsigned long MDParser::_tokenSet_19_data_[] = { 135790592UL, 272256UL, 0UL, 0UL };
// "center" "satellites" "Harmonic" "AmberImproper" "ImproperCosine" "CentralAtomHeight" 
// "Dreiding" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_19(_tokenSet_19_data_,4);
const unsigned long MDParser::_tokenSet_20_data_[] = { 262144UL, 270336UL, 0UL, 0UL };
// "members" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_20(_tokenSet_20_data_,4);
const unsigned long MDParser::_tokenSet_21_data_[] = { 0UL, 270336UL, 0UL, 0UL };
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_21(_tokenSet_21_data_,4);
const unsigned long MDParser::_tokenSet_22_data_[] = { 0UL, 12582912UL, 0UL, 0UL };
// RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_22(_tokenSet_22_data_,4);


