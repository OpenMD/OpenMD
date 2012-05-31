/* $ANTLR 2.7.7 (20110725): "MDParser.g" -> "MDParser.cpp"$ */
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
#line 71 "MDParser.g"
		tmp6_AST->setType(ENDBLOCK);
#line 187 "MDParser.cpp"
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
				goto _loop24;
			}
			
		}
		_loop24:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp9_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp9_AST);
		match(RCURLY);
#line 86 "MDParser.g"
		tmp9_AST->setType(ENDBLOCK);
#line 227 "MDParser.cpp"
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
#line 74 "MDParser.g"
		tmp12_AST->setType(ENDBLOCK);
#line 267 "MDParser.cpp"
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
#line 77 "MDParser.g"
		tmp15_AST->setType(ENDBLOCK);
#line 307 "MDParser.cpp"
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
#line 80 "MDParser.g"
		tmp18_AST->setType(ENDBLOCK);
#line 347 "MDParser.cpp"
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
#line 83 "MDParser.g"
		tmp21_AST->setType(ENDBLOCK);
#line 387 "MDParser.cpp"
		rnemdblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = rnemdblock_AST;
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
		case ID:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp22_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp22_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp22_AST);
			match(ID);
			constant_AST = currentAST.root;
			break;
		}
		case StringLiteral:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp23_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp23_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp23_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp24_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp24_AST);
			match(NUM_INT);
			intConst_AST = currentAST.root;
			break;
		}
		case NUM_LONG:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp25_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp25_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp26_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp26_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp26_AST);
			match(NUM_FLOAT);
			floatConst_AST = currentAST.root;
			break;
		}
		case NUM_DOUBLE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp27_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp27_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp28_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp28_AST);
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
				goto _loop28;
			}
			
		}
		_loop28:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp32_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp32_AST);
		match(RCURLY);
#line 100 "MDParser.g"
		tmp32_AST->setType(ENDBLOCK);
#line 644 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp33_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp33_AST);
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
				goto _loop33;
			}
			
		}
		_loop33:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp37_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp37_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp37_AST);
		match(RCURLY);
#line 109 "MDParser.g"
		tmp37_AST->setType(ENDBLOCK);
#line 703 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp38_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp38_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp38_AST);
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
				goto _loop38;
			}
			
		}
		_loop38:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp42_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp42_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp42_AST);
		match(RCURLY);
#line 116 "MDParser.g"
		tmp42_AST->setType(ENDBLOCK);
#line 762 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp43_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp43_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp43_AST);
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
				goto _loop43;
			}
			
		}
		_loop43:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp47_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp47_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp47_AST);
		match(RCURLY);
#line 123 "MDParser.g"
		tmp47_AST->setType(ENDBLOCK);
#line 821 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp48_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp48_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp48_AST);
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
			if ((LA(1) == CENTER || LA(1) == ID)) {
				inversionstatement();
				astFactory->addASTChild( currentAST, returnAST );
			}
			else {
				goto _loop48;
			}
			
		}
		_loop48:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp52_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp52_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp52_AST);
		match(RCURLY);
#line 130 "MDParser.g"
		tmp52_AST->setType(ENDBLOCK);
#line 880 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp53_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp53_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp53_AST);
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
				goto _loop52;
			}
			
		}
		_loop52:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp57_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp57_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp57_AST);
		match(RCURLY);
#line 137 "MDParser.g"
		tmp57_AST->setType(ENDBLOCK);
#line 924 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp58_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp58_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp58_AST);
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
				goto _loop57;
			}
			
		}
		_loop57:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp62_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp62_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp62_AST);
		match(RCURLY);
#line 144 "MDParser.g"
		tmp62_AST->setType(ENDBLOCK);
#line 983 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp63_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp63_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp63_AST);
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
				goto _loop61;
			}
			
		}
		_loop61:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp67_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp67_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp67_AST);
		match(RCURLY);
#line 151 "MDParser.g"
		tmp67_AST->setType(ENDBLOCK);
#line 1027 "MDParser.cpp"
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp68_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp68_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp68_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp72_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp72_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp72_AST);
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
		recover(ex,_tokenSet_10);
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
				goto _loop65;
			}
			
		}
		_loop65:;
		} // ( ... )*
		doubleNumberTuple_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_11);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp77_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp77_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp77_AST);
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
		recover(ex,_tokenSet_12);
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
				goto _loop68;
			}
			
		}
		_loop68:;
		} // ( ... )*
		inttuple_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_11);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp82_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp82_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp82_AST);
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
		recover(ex,_tokenSet_12);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp86_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp86_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp86_AST);
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
		recover(ex,_tokenSet_12);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp90_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp90_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp90_AST);
			match(CENTER);
			match(LPAREN);
			intConst();
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
		recover(ex,_tokenSet_13);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp94_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp94_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp94_AST);
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
		recover(ex,_tokenSet_12);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp98_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp98_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp98_AST);
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
		recover(ex,_tokenSet_12);
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
		recover(ex,_tokenSet_14);
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
		recover(ex,_tokenSet_15);
	}
	returnAST = doubleNumber_AST;
}

void MDParser::initializeASTFactory( ANTLR_USE_NAMESPACE(antlr)ASTFactory& factory )
{
	factory.setMaxNodeType(55);
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
	"\"position\"",
	"\"orientation\"",
	"\"flucQ\"",
	"\"RNEMD\"",
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

const unsigned long MDParser::_tokenSet_0_data_[] = { 11534576UL, 0UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDParser::_tokenSet_1_data_[] = { 2UL, 0UL, 0UL, 0UL };
// EOF 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDParser::_tokenSet_2_data_[] = { 11534578UL, 0UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" 
// ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_2(_tokenSet_2_data_,4);
const unsigned long MDParser::_tokenSet_3_data_[] = { 281018354UL, 0UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "restraint" "atom" "bond" "bend" 
// "torsion" "inversion" "rigidBody" "cutoffGroup" "fragment" "members" 
// "center" "position" "orientation" "flucQ" "RNEMD" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_3(_tokenSet_3_data_,4);
const unsigned long MDParser::_tokenSet_4_data_[] = { 8453888UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_4(_tokenSet_4_data_,4);
const unsigned long MDParser::_tokenSet_5_data_[] = { 33554432UL, 0UL, 0UL, 0UL };
// SEMICOLON 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_5(_tokenSet_5_data_,4);
const unsigned long MDParser::_tokenSet_6_data_[] = { 1107296256UL, 3UL, 0UL, 0UL };
// SEMICOLON RBRACKET RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_6(_tokenSet_6_data_,4);
const unsigned long MDParser::_tokenSet_7_data_[] = { 33554432UL, 3UL, 0UL, 0UL };
// SEMICOLON RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_7(_tokenSet_7_data_,4);
const unsigned long MDParser::_tokenSet_8_data_[] = { 276889344UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_8(_tokenSet_8_data_,4);
const unsigned long MDParser::_tokenSet_9_data_[] = { 9175040UL, 0UL, 0UL, 0UL };
// "position" "orientation" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_9(_tokenSet_9_data_,4);
const unsigned long MDParser::_tokenSet_10_data_[] = { 277610496UL, 0UL, 0UL, 0UL };
// "position" "orientation" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_10(_tokenSet_10_data_,4);
const unsigned long MDParser::_tokenSet_11_data_[] = { 0UL, 1UL, 0UL, 0UL };
// RPAREN 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_11(_tokenSet_11_data_,4);
const unsigned long MDParser::_tokenSet_12_data_[] = { 276889600UL, 0UL, 0UL, 0UL };
// "members" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_12(_tokenSet_12_data_,4);
const unsigned long MDParser::_tokenSet_13_data_[] = { 276955136UL, 0UL, 0UL, 0UL };
// "center" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_13(_tokenSet_13_data_,4);
const unsigned long MDParser::_tokenSet_14_data_[] = { 276824064UL, 0UL, 0UL, 0UL };
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_14(_tokenSet_14_data_,4);
const unsigned long MDParser::_tokenSet_15_data_[] = { 0UL, 3UL, 0UL, 0UL };
// RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_15(_tokenSet_15_data_,4);


