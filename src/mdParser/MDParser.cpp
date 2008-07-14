/* $ANTLR 2.7.7 (20080702): "MDParser.g" -> "MDParser.cpp"$ */
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
#line 65 "MDParser.g"
		tmp6_AST->setType(ENDBLOCK);
#line 166 "MDParser.cpp"
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
				goto _loop15;
			}
			
		}
		_loop15:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp9_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp9_AST);
		match(RCURLY);
#line 71 "MDParser.g"
		tmp9_AST->setType(ENDBLOCK);
#line 206 "MDParser.cpp"
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
#line 68 "MDParser.g"
		tmp12_AST->setType(ENDBLOCK);
#line 246 "MDParser.cpp"
		zconstraintblock_AST = currentAST.root;
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		recover(ex,_tokenSet_2);
	}
	returnAST = zconstraintblock_AST;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp13_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp13_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp13_AST);
			match(ID);
			constant_AST = currentAST.root;
			break;
		}
		case StringLiteral:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp14_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp14_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp14_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp15_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp15_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp15_AST);
			match(NUM_INT);
			intConst_AST = currentAST.root;
			break;
		}
		case NUM_LONG:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp16_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp16_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp16_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp17_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp17_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp17_AST);
			match(NUM_FLOAT);
			floatConst_AST = currentAST.root;
			break;
		}
		case NUM_DOUBLE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp18_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp18_AST = astFactory->create(LT(1));
			astFactory->addASTChild(currentAST, tmp18_AST);
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp19_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp19_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp19_AST);
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
				goto _loop19;
			}
			
		}
		_loop19:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp23_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp23_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp23_AST);
		match(RCURLY);
#line 85 "MDParser.g"
		tmp23_AST->setType(ENDBLOCK);
#line 503 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp24_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp24_AST);
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
				goto _loop24;
			}
			
		}
		_loop24:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp28_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp28_AST);
		match(RCURLY);
#line 94 "MDParser.g"
		tmp28_AST->setType(ENDBLOCK);
#line 562 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp29_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp29_AST);
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
				goto _loop29;
			}
			
		}
		_loop29:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp33_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp33_AST);
		match(RCURLY);
#line 101 "MDParser.g"
		tmp33_AST->setType(ENDBLOCK);
#line 621 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp34_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp34_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp34_AST);
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
				goto _loop34;
			}
			
		}
		_loop34:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp38_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp38_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp38_AST);
		match(RCURLY);
#line 108 "MDParser.g"
		tmp38_AST->setType(ENDBLOCK);
#line 680 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp39_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp39_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp39_AST);
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
				goto _loop39;
			}
			
		}
		_loop39:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp43_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp43_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp43_AST);
		match(RCURLY);
#line 115 "MDParser.g"
		tmp43_AST->setType(ENDBLOCK);
#line 739 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp44_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp44_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp44_AST);
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
				goto _loop43;
			}
			
		}
		_loop43:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp48_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp48_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp48_AST);
		match(RCURLY);
#line 122 "MDParser.g"
		tmp48_AST->setType(ENDBLOCK);
#line 783 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp49_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp49_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp49_AST);
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
				goto _loop48;
			}
			
		}
		_loop48:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp53_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp53_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp53_AST);
		match(RCURLY);
#line 129 "MDParser.g"
		tmp53_AST->setType(ENDBLOCK);
#line 842 "MDParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp54_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp54_AST = astFactory->create(LT(1));
		astFactory->makeASTRoot(currentAST, tmp54_AST);
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
				goto _loop52;
			}
			
		}
		_loop52:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp58_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
		tmp58_AST = astFactory->create(LT(1));
		astFactory->addASTChild(currentAST, tmp58_AST);
		match(RCURLY);
#line 136 "MDParser.g"
		tmp58_AST->setType(ENDBLOCK);
#line 886 "MDParser.cpp"
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp59_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp59_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp59_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp63_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp63_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp63_AST);
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
				goto _loop56;
			}
			
		}
		_loop56:;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp68_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp68_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp68_AST);
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
				goto _loop59;
			}
			
		}
		_loop59:;
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp73_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp73_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp73_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp77_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp77_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp77_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp81_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp81_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp81_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp85_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp85_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp85_AST);
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
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp89_AST = ANTLR_USE_NAMESPACE(antlr)nullAST;
			tmp89_AST = astFactory->create(LT(1));
			astFactory->makeASTRoot(currentAST, tmp89_AST);
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
	factory.setMaxNodeType(52);
}
const char* MDParser::tokenNames[] = {
	"<0>",
	"EOF",
	"<2>",
	"NULL_TREE_LOOKAHEAD",
	"\"component\"",
	"\"molecule\"",
	"\"zconstraint\"",
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

const unsigned long MDParser::_tokenSet_0_data_[] = { 1048688UL, 0UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDParser::_tokenSet_1_data_[] = { 2UL, 0UL, 0UL, 0UL };
// EOF 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDParser::_tokenSet_2_data_[] = { 1048690UL, 0UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_2(_tokenSet_2_data_,4);
const unsigned long MDParser::_tokenSet_3_data_[] = { 35127282UL, 0UL, 0UL, 0UL };
// EOF "component" "molecule" "zconstraint" "atom" "bond" "bend" "torsion" 
// "inversion" "rigidBody" "cutoffGroup" "fragment" "members" "center" 
// "position" "orientation" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_3(_tokenSet_3_data_,4);
const unsigned long MDParser::_tokenSet_4_data_[] = { 1081216UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_4(_tokenSet_4_data_,4);
const unsigned long MDParser::_tokenSet_5_data_[] = { 4194304UL, 0UL, 0UL, 0UL };
// SEMICOLON 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_5(_tokenSet_5_data_,4);
const unsigned long MDParser::_tokenSet_6_data_[] = { 1749024768UL, 0UL, 0UL, 0UL };
// SEMICOLON RBRACKET RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_6(_tokenSet_6_data_,4);
const unsigned long MDParser::_tokenSet_7_data_[] = { 1614807040UL, 0UL, 0UL, 0UL };
// SEMICOLON RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_7(_tokenSet_7_data_,4);
const unsigned long MDParser::_tokenSet_8_data_[] = { 34635648UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_8(_tokenSet_8_data_,4);
const unsigned long MDParser::_tokenSet_9_data_[] = { 1441792UL, 0UL, 0UL, 0UL };
// "position" "orientation" ID 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_9(_tokenSet_9_data_,4);
const unsigned long MDParser::_tokenSet_10_data_[] = { 34996224UL, 0UL, 0UL, 0UL };
// "position" "orientation" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_10(_tokenSet_10_data_,4);
const unsigned long MDParser::_tokenSet_11_data_[] = { 536870912UL, 0UL, 0UL, 0UL };
// RPAREN 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_11(_tokenSet_11_data_,4);
const unsigned long MDParser::_tokenSet_12_data_[] = { 34635776UL, 0UL, 0UL, 0UL };
// "members" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_12(_tokenSet_12_data_,4);
const unsigned long MDParser::_tokenSet_13_data_[] = { 34668544UL, 0UL, 0UL, 0UL };
// "center" ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_13(_tokenSet_13_data_,4);
const unsigned long MDParser::_tokenSet_14_data_[] = { 34603008UL, 0UL, 0UL, 0UL };
// ID RCURLY 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_14(_tokenSet_14_data_,4);
const unsigned long MDParser::_tokenSet_15_data_[] = { 1610612736UL, 0UL, 0UL, 0UL };
// RPAREN COMMA 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDParser::_tokenSet_15(_tokenSet_15_data_,4);


