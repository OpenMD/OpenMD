/* $ANTLR 2.7.7 (2006-11-01): "MDTreeParser.g" -> "MDTreeParser.cpp"$ */
#include "MDTreeParser.hpp"
#include <antlr/Token.hpp>
#include <antlr/AST.hpp>
#include <antlr/NoViableAltException.hpp>
#include <antlr/MismatchedTokenException.hpp>
#include <antlr/SemanticException.hpp>
#include <antlr/BitSet.hpp>
#line 1 "MDTreeParser.g"
#line 11 "MDTreeParser.cpp"
MDTreeParser::MDTreeParser()
	: ANTLR_USE_NAMESPACE(antlr)TreeParser() {
}

void MDTreeParser::mdfile(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST mdfile_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_0.member(_t->getType()))) {
				statement(_t);
				_t = _retTree;
			}
			else {
				goto _loop3;
			}
			
		}
		_loop3:;
		} // ( ... )*
#line 34 "MDTreeParser.g"
		blockStack.top()->validate(); blockStack.pop();
#line 37 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::statement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST statement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case COMPONENT:
		{
			componentblock(_t);
			_t = _retTree;
			break;
		}
		case MOLECULE:
		{
			moleculeblock(_t);
			_t = _retTree;
			break;
		}
		case FRAGMENT:
		{
			fragmentblock(_t);
			_t = _retTree;
			break;
		}
		case ZCONSTRAINT:
		{
			zconstraintblock(_t);
			_t = _retTree;
			break;
		}
		case RESTRAINT:
		{
			restraintblock(_t);
			_t = _retTree;
			break;
		}
		case FLUCQ:
		{
			flucqblock(_t);
			_t = _retTree;
			break;
		}
		case RNEMD:
		{
			rnemdblock(_t);
			_t = _retTree;
			break;
		}
		case MINIMIZER:
		{
			minimizerblock(_t);
			_t = _retTree;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::assignment(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST assignment_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST id = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t6 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp1_AST_in = _t;
		match(_t,ASSIGNEQUAL);
		_t = _t->getFirstChild();
		id = _t;
		match(_t,ID);
		_t = _t->getNextSibling();
		constant(_t,id);
		_t = _retTree;
		_t = __t6;
		_t = _t->getNextSibling();
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::componentblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST componentblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t10 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp2_AST_in = _t;
		match(_t,COMPONENT);
		_t = _t->getFirstChild();
#line 71 "MDTreeParser.g"
		Component* currComponet = new Component(); blockStack.push(currComponet);
#line 157 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop12;
			}
			
		}
		_loop12:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp3_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t10;
		_t = _t->getNextSibling();
#line 73 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addComponent(currComponet);
#line 180 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::moleculeblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculeblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t34 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp4_AST_in = _t;
		match(_t,MOLECULE);
		_t = _t->getFirstChild();
#line 102 "MDTreeParser.g"
		MoleculeStamp* currMoleculeStamp = new MoleculeStamp(); blockStack.push(currMoleculeStamp);
#line 200 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_1.member(_t->getType()))) {
				moleculestatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop36;
			}
			
		}
		_loop36:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp5_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t34;
		_t = _t->getNextSibling();
#line 104 "MDTreeParser.g"
		blockStack.top()->validate(); blockStack.pop(); currConf->addMoleculeStamp(currMoleculeStamp);
#line 223 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::fragmentblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t41 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp6_AST_in = _t;
		match(_t,FRAGMENT);
		_t = _t->getFirstChild();
#line 122 "MDTreeParser.g"
		FragmentStamp* currFragmentStamp = new FragmentStamp(); blockStack.push(currFragmentStamp);
#line 243 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_2.member(_t->getType()))) {
				fragmentstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop43;
			}
			
		}
		_loop43:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp7_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t41;
		_t = _t->getNextSibling();
#line 124 "MDTreeParser.g"
		blockStack.top()->validate(); blockStack.pop(); currConf->addFragmentStamp(currFragmentStamp);
#line 266 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::zconstraintblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST zconstraintblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t14 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp8_AST_in = _t;
		match(_t,ZCONSTRAINT);
		_t = _t->getFirstChild();
#line 76 "MDTreeParser.g"
		ZConsStamp* currZConsStamp = new ZConsStamp(); blockStack.push(currZConsStamp);
#line 286 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop16;
			}
			
		}
		_loop16:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t14;
		_t = _t->getNextSibling();
#line 78 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addZConsStamp(currZConsStamp);
#line 309 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::restraintblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST restraintblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t18 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp10_AST_in = _t;
		match(_t,RESTRAINT);
		_t = _t->getFirstChild();
#line 81 "MDTreeParser.g"
		RestraintStamp* currRestraintStamp = new RestraintStamp(); blockStack.push(currRestraintStamp);
#line 329 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop20;
			}
			
		}
		_loop20:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp11_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t18;
		_t = _t->getNextSibling();
#line 83 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addRestraintStamp(currRestraintStamp);
#line 352 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::flucqblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST flucqblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t22 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp12_AST_in = _t;
		match(_t,FLUCQ);
		_t = _t->getFirstChild();
#line 86 "MDTreeParser.g"
		FluctuatingChargeParameters* flucQpars = new FluctuatingChargeParameters(); blockStack.push(flucQpars);
#line 372 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop24;
			}
			
		}
		_loop24:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp13_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t22;
		_t = _t->getNextSibling();
#line 88 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addFluctuatingChargeParameters(flucQpars);
#line 395 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::rnemdblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST rnemdblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t26 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp14_AST_in = _t;
		match(_t,RNEMD);
		_t = _t->getFirstChild();
#line 91 "MDTreeParser.g"
		RNEMD::RNEMDParameters* rnemdPars = new RNEMD::RNEMDParameters(); blockStack.push(rnemdPars);
#line 415 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop28;
			}
			
		}
		_loop28:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp15_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t26;
		_t = _t->getNextSibling();
#line 93 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addRNEMDParameters(rnemdPars);
#line 438 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::minimizerblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST minimizerblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t30 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp16_AST_in = _t;
		match(_t,MINIMIZER);
		_t = _t->getFirstChild();
#line 96 "MDTreeParser.g"
		MinimizerParameters* minimizerPars = new MinimizerParameters(); blockStack.push(minimizerPars);
#line 458 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop32;
			}
			
		}
		_loop32:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp17_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t30;
		_t = _t->getNextSibling();
#line 98 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addMinimizerParameters(minimizerPars);
#line 481 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::constant(ANTLR_USE_NAMESPACE(antlr)RefAST _t,
	ANTLR_USE_NAMESPACE(antlr)RefAST id
) {
	ANTLR_USE_NAMESPACE(antlr)RefAST constant_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST str1 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST str2 = ANTLR_USE_NAMESPACE(antlr)nullAST;
#line 51 "MDTreeParser.g"
	
	int ival;
	RealType dval;
	std::vector<RealType> dvec;
	
#line 503 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case NUM_INT:
		case NUM_LONG:
		{
			ival=intConst(_t);
			_t = _retTree;
#line 57 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), ival);
#line 516 "MDTreeParser.cpp"
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			dval=floatConst(_t);
			_t = _retTree;
#line 58 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), dval);
#line 526 "MDTreeParser.cpp"
			break;
		}
		case ID:
		{
			str1 = _t;
			match(_t,ID);
			_t = _t->getNextSibling();
#line 59 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), str1->getText());
#line 536 "MDTreeParser.cpp"
			break;
		}
		case StringLiteral:
		{
			str2 = _t;
			match(_t,StringLiteral);
			_t = _t->getNextSibling();
#line 60 "MDTreeParser.g"
			std::string s =  str2->getText();
			s = s.substr(1, s.length()-2);
			blockStack.top()->assign(id->getText(),s);
			
#line 549 "MDTreeParser.cpp"
			break;
		}
		case LPAREN:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t8 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp18_AST_in = _t;
			match(_t,LPAREN);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp19_AST_in = _t;
			match(_t,RPAREN);
			_t = _t->getNextSibling();
			_t = __t8;
			_t = _t->getNextSibling();
#line 65 "MDTreeParser.g"
			
			blockStack.top()->assign(id->getText(), dvec);
			
#line 569 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

int  MDTreeParser::intConst(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 376 "MDTreeParser.g"
	int ival;
#line 589 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST intConst_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST i1 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST i2 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case NUM_INT:
		{
			i1 = _t;
			match(_t,NUM_INT);
			_t = _t->getNextSibling();
#line 377 "MDTreeParser.g"
			ival = lexi_cast<int>(i1->getText());
#line 605 "MDTreeParser.cpp"
			break;
		}
		case NUM_LONG:
		{
			i2 = _t;
			match(_t,NUM_LONG);
			_t = _t->getNextSibling();
#line 378 "MDTreeParser.g"
			ival = lexi_cast<int>(i2->getText());
#line 615 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
	return ival;
}

RealType  MDTreeParser::floatConst(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 388 "MDTreeParser.g"
	RealType dval;
#line 636 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST floatConst_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST d1 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST d2 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case NUM_FLOAT:
		{
			d1 = _t;
			match(_t,NUM_FLOAT);
			_t = _t->getNextSibling();
#line 389 "MDTreeParser.g"
			dval = lexi_cast<RealType>(d1->getText());
#line 652 "MDTreeParser.cpp"
			break;
		}
		case NUM_DOUBLE:
		{
			d2 = _t;
			match(_t,NUM_DOUBLE);
			_t = _t->getNextSibling();
#line 390 "MDTreeParser.g"
			dval = lexi_cast<RealType>(d2->getText());
#line 662 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
	return dval;
}

vector<RealType>  MDTreeParser::doubleNumberTuple(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 360 "MDTreeParser.g"
	vector<RealType> dvec;
#line 683 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST doubleNumberTuple_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 360 "MDTreeParser.g"
	
	RealType dval;
	
#line 689 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt130=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if (((_t->getType() >= NUM_INT && _t->getType() <= NUM_DOUBLE))) {
				dval=doubleNumber(_t);
				_t = _retTree;
#line 364 "MDTreeParser.g"
				dvec.push_back(dval);
#line 702 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt130>=1 ) { goto _loop130; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt130++;
		}
		_loop130:;
		}  // ( ... )+
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
	return dvec;
}

void MDTreeParser::moleculestatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculestatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ATOM:
		{
			atomblock(_t);
			_t = _retTree;
			break;
		}
		case BOND:
		{
			bondblock(_t);
			_t = _retTree;
			break;
		}
		case BEND:
		{
			bendblock(_t);
			_t = _retTree;
			break;
		}
		case TORSION:
		{
			torsionblock(_t);
			_t = _retTree;
			break;
		}
		case INVERSION:
		{
			inversionblock(_t);
			_t = _retTree;
			break;
		}
		case RIGIDBODY:
		{
			rigidbodyblock(_t);
			_t = _retTree;
			break;
		}
		case CUTOFFGROUP:
		{
			cutoffgroupblock(_t);
			_t = _retTree;
			break;
		}
		case CONSTRAINT:
		{
			constraintblock(_t);
			_t = _retTree;
			break;
		}
		default:
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else if ((_t->getType() == ASSIGNEQUAL)) {
				sequencestring(_t);
				_t = _retTree;
			}
		else {
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::atomblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST atomblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 139 "MDTreeParser.g"
	
	int index;
	
#line 807 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t46 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp20_AST_in = _t;
		match(_t,ATOM);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 143 "MDTreeParser.g"
		AtomStamp* currAtomStamp = new AtomStamp(index); blockStack.push(currAtomStamp);
#line 818 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_3.member(_t->getType()))) {
				atomstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop48;
			}
			
		}
		_loop48:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp21_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t46;
		_t = _t->getNextSibling();
#line 145 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addAtomStamp(currAtomStamp); 
		
#line 846 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::bondblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bondblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t54 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp22_AST_in = _t;
		match(_t,BOND);
		_t = _t->getFirstChild();
#line 168 "MDTreeParser.g"
		BondStamp* currBondStamp = new BondStamp(); blockStack.push(currBondStamp);
#line 866 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_4.member(_t->getType()))) {
				bondstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop56;
			}
			
		}
		_loop56:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp23_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t54;
		_t = _t->getNextSibling();
#line 170 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBondStamp(currBondStamp); 
		
#line 894 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::bendblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bendblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t66 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST_in = _t;
		match(_t,BEND);
		_t = _t->getFirstChild();
#line 195 "MDTreeParser.g"
		BendStamp* currBendStamp = new BendStamp(); blockStack.push(currBendStamp);
#line 914 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_5.member(_t->getType()))) {
				bendstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop68;
			}
			
		}
		_loop68:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t66;
		_t = _t->getNextSibling();
#line 197 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBendStamp(currBendStamp); 
		
#line 942 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::torsionblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST torsionblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t79 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp26_AST_in = _t;
		match(_t,TORSION);
		_t = _t->getFirstChild();
#line 222 "MDTreeParser.g"
		TorsionStamp* currTorsionStamp = new TorsionStamp(); blockStack.push(currTorsionStamp);
#line 962 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_6.member(_t->getType()))) {
				torsionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop81;
			}
			
		}
		_loop81:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t79;
		_t = _t->getNextSibling();
#line 224 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addTorsionStamp(currTorsionStamp); 
		
#line 990 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::inversionblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST inversionblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t93 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST_in = _t;
		match(_t,INVERSION);
		_t = _t->getFirstChild();
#line 251 "MDTreeParser.g"
		InversionStamp* currInversionStamp = new InversionStamp(); blockStack.push(currInversionStamp);
#line 1010 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_7.member(_t->getType()))) {
				inversionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop95;
			}
			
		}
		_loop95:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t93;
		_t = _t->getNextSibling();
#line 253 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addInversionStamp(currInversionStamp); 
		
#line 1038 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::rigidbodyblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST rigidbodyblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 279 "MDTreeParser.g"
	
	int index;
	
#line 1054 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t105 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp30_AST_in = _t;
		match(_t,RIGIDBODY);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 283 "MDTreeParser.g"
		RigidBodyStamp* currRigidBodyStamp = new RigidBodyStamp(index); blockStack.push(currRigidBodyStamp);
#line 1065 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				rigidbodystatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop107;
			}
			
		}
		_loop107:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp31_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t105;
		_t = _t->getNextSibling();
#line 285 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addRigidBodyStamp(currRigidBodyStamp); 
		
#line 1093 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::cutoffgroupblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST cutoffgroupblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t111 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST_in = _t;
		match(_t,CUTOFFGROUP);
		_t = _t->getFirstChild();
#line 302 "MDTreeParser.g"
		CutoffGroupStamp* currCutoffGroupStamp = new CutoffGroupStamp(); blockStack.push(currCutoffGroupStamp);
#line 1113 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				cutoffgroupstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop113;
			}
			
		}
		_loop113:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t111;
		_t = _t->getNextSibling();
#line 304 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addCutoffGroupStamp(currCutoffGroupStamp); 
		
#line 1141 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::constraintblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST constraintblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t117 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp34_AST_in = _t;
		match(_t,CONSTRAINT);
		_t = _t->getFirstChild();
#line 322 "MDTreeParser.g"
		ConstraintStamp* currConstraintStamp = new ConstraintStamp(); blockStack.push(currConstraintStamp);
#line 1161 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				constraintstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop119;
			}
			
		}
		_loop119:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp35_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t117;
		_t = _t->getNextSibling();
#line 324 "MDTreeParser.g"
		
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addConstraintStamp(currConstraintStamp); 
		
#line 1188 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::sequencestring(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST sequencestring_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST id = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t39 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp36_AST_in = _t;
		match(_t,ASSIGNEQUAL);
		_t = _t->getFirstChild();
		id = _t;
		match(_t,SEQUENCE);
		_t = _t->getNextSibling();
		constant(_t,id);
		_t = _retTree;
		_t = __t39;
		_t = _t->getNextSibling();
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::fragmentstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case ATOM:
		{
			atomblock(_t);
			_t = _retTree;
			break;
		}
		case BOND:
		{
			bondblock(_t);
			_t = _retTree;
			break;
		}
		case BEND:
		{
			bendblock(_t);
			_t = _retTree;
			break;
		}
		case TORSION:
		{
			torsionblock(_t);
			_t = _retTree;
			break;
		}
		case INVERSION:
		{
			inversionblock(_t);
			_t = _retTree;
			break;
		}
		case RIGIDBODY:
		{
			rigidbodyblock(_t);
			_t = _retTree;
			break;
		}
		case CUTOFFGROUP:
		{
			cutoffgroupblock(_t);
			_t = _retTree;
			break;
		}
		case CONSTRAINT:
		{
			constraintblock(_t);
			_t = _retTree;
			break;
		}
		case NODES:
		{
			nodesblock(_t);
			_t = _retTree;
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::nodesblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST nodesblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t123 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp37_AST_in = _t;
		match(_t,NODES);
		_t = _t->getFirstChild();
#line 341 "MDTreeParser.g"
		NodesStamp* currNodesStamp = new NodesStamp(); blockStack.push(currNodesStamp);
#line 1314 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				nodesstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop125;
			}
			
		}
		_loop125:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp38_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t123;
		_t = _t->getNextSibling();
#line 343 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		FragmentStamp* currFragmentStamp = static_cast<FragmentStamp*>(blockStack.top());
		currFragmentStamp->addNodesStamp(currNodesStamp); 
		
#line 1342 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::atomstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST atomstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 153 "MDTreeParser.g"
	
	vector<RealType> dvec;
	RealType rval;
	AtomStamp* currAtomStamp =  static_cast<AtomStamp*>(blockStack.top());
	
	
#line 1361 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case POSITION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t50 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp39_AST_in = _t;
			match(_t,POSITION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t50;
			_t = _t->getNextSibling();
#line 161 "MDTreeParser.g"
			currAtomStamp->setPosition(dvec);
#line 1385 "MDTreeParser.cpp"
			break;
		}
		case ORIENTATION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t51 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp40_AST_in = _t;
			match(_t,ORIENTATION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t51;
			_t = _t->getNextSibling();
#line 162 "MDTreeParser.g"
			currAtomStamp->setOrientation(dvec);
#line 1400 "MDTreeParser.cpp"
			break;
		}
		case CHARGE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t52 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp41_AST_in = _t;
			match(_t,CHARGE);
			_t = _t->getFirstChild();
			rval=doubleNumber(_t);
			_t = _retTree;
			_t = __t52;
			_t = _t->getNextSibling();
#line 163 "MDTreeParser.g"
			currAtomStamp->overrideCharge(rval);
#line 1415 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

RealType  MDTreeParser::doubleNumber(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 382 "MDTreeParser.g"
	RealType dval;
#line 1435 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST doubleNumber_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST ic = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fc = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case NUM_INT:
		case NUM_LONG:
		{
			ic = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
			intConst(_t);
			_t = _retTree;
#line 383 "MDTreeParser.g"
			dval = lexi_cast<RealType>(ic->getText());
#line 1452 "MDTreeParser.cpp"
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			fc = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
			floatConst(_t);
			_t = _retTree;
#line 384 "MDTreeParser.g"
			dval = lexi_cast<RealType>(fc->getText());
#line 1463 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
	return dval;
}

void MDTreeParser::bondstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bondstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 178 "MDTreeParser.g"
	
	vector<int> ivec;
	RealType rval;
	vector<RealType> dvec;
	BondStamp* currBondStamp = static_cast<BondStamp*>(blockStack.top());
	
#line 1490 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t58 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp42_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t58;
			_t = _t->getNextSibling();
#line 186 "MDTreeParser.g"
			currBondStamp->setMembers(ivec);
#line 1514 "MDTreeParser.cpp"
			break;
		}
		case FIXED:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t59 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp43_AST_in = _t;
			match(_t,FIXED);
			_t = _t->getFirstChild();
			rval=doubleNumber(_t);
			_t = _retTree;
			_t = __t59;
			_t = _t->getNextSibling();
#line 187 "MDTreeParser.g"
			currBondStamp->overrideType("Fixed", rval);
#line 1529 "MDTreeParser.cpp"
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t60 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp44_AST_in = _t;
			match(_t,HARMONIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t60;
			_t = _t->getNextSibling();
#line 188 "MDTreeParser.g"
			currBondStamp->overrideType("Harmonic", dvec);
#line 1544 "MDTreeParser.cpp"
			break;
		}
		case CUBIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t61 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp45_AST_in = _t;
			match(_t,CUBIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t61;
			_t = _t->getNextSibling();
#line 189 "MDTreeParser.g"
			currBondStamp->overrideType("Cubic", dvec);
#line 1559 "MDTreeParser.cpp"
			break;
		}
		case QUARTIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t62 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp46_AST_in = _t;
			match(_t,QUARTIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t62;
			_t = _t->getNextSibling();
#line 190 "MDTreeParser.g"
			currBondStamp->overrideType("Quartic", dvec);
#line 1574 "MDTreeParser.cpp"
			break;
		}
		case POLYNOMIAL:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t63 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp47_AST_in = _t;
			match(_t,POLYNOMIAL);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t63;
			_t = _t->getNextSibling();
#line 191 "MDTreeParser.g"
			currBondStamp->overrideType("Polynomial", dvec);
#line 1589 "MDTreeParser.cpp"
			break;
		}
		case MORSE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t64 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp48_AST_in = _t;
			match(_t,MORSE);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t64;
			_t = _t->getNextSibling();
#line 192 "MDTreeParser.g"
			currBondStamp->overrideType("Morse", dvec);
#line 1604 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

vector<int>  MDTreeParser::inttuple(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 368 "MDTreeParser.g"
	vector<int> ivec;
#line 1624 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST inttuple_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 368 "MDTreeParser.g"
	
	int ival;
	
#line 1630 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt133=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == NUM_INT || _t->getType() == NUM_LONG)) {
				ival=intConst(_t);
				_t = _retTree;
#line 372 "MDTreeParser.g"
				ivec.push_back(ival);
#line 1643 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt133>=1 ) { goto _loop133; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt133++;
		}
		_loop133:;
		}  // ( ... )+
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
	return ivec;
}

void MDTreeParser::bendstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bendstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 205 "MDTreeParser.g"
	
	vector<int> ivec;
	vector<RealType> dvec;
	BendStamp* currBendStamp = static_cast<BendStamp*>(blockStack.top());
	
#line 1671 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t70 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp49_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t70;
			_t = _t->getNextSibling();
#line 212 "MDTreeParser.g"
			currBendStamp->setMembers(ivec);
#line 1695 "MDTreeParser.cpp"
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t71 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp50_AST_in = _t;
			match(_t,HARMONIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t71;
			_t = _t->getNextSibling();
#line 213 "MDTreeParser.g"
			currBendStamp->overrideType("Harmonic", dvec);
#line 1710 "MDTreeParser.cpp"
			break;
		}
		case GHOSTBEND:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t72 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp51_AST_in = _t;
			match(_t,GHOSTBEND);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t72;
			_t = _t->getNextSibling();
#line 214 "MDTreeParser.g"
			currBendStamp->overrideType("GhostBend", dvec);
#line 1725 "MDTreeParser.cpp"
			break;
		}
		case UREYBRADLEY:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t73 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp52_AST_in = _t;
			match(_t,UREYBRADLEY);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t73;
			_t = _t->getNextSibling();
#line 215 "MDTreeParser.g"
			currBendStamp->overrideType("UreyBradley", dvec);
#line 1740 "MDTreeParser.cpp"
			break;
		}
		case CUBIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t74 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp53_AST_in = _t;
			match(_t,CUBIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t74;
			_t = _t->getNextSibling();
#line 216 "MDTreeParser.g"
			currBendStamp->overrideType("Cubic", dvec);
#line 1755 "MDTreeParser.cpp"
			break;
		}
		case QUARTIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t75 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp54_AST_in = _t;
			match(_t,QUARTIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t75;
			_t = _t->getNextSibling();
#line 217 "MDTreeParser.g"
			currBendStamp->overrideType("Quartic", dvec);
#line 1770 "MDTreeParser.cpp"
			break;
		}
		case POLYNOMIAL:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t76 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp55_AST_in = _t;
			match(_t,POLYNOMIAL);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t76;
			_t = _t->getNextSibling();
#line 218 "MDTreeParser.g"
			currBendStamp->overrideType("Polynomial", dvec);
#line 1785 "MDTreeParser.cpp"
			break;
		}
		case COSINE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t77 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp56_AST_in = _t;
			match(_t,COSINE);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t77;
			_t = _t->getNextSibling();
#line 219 "MDTreeParser.g"
			currBendStamp->overrideType("Cosine", dvec);
#line 1800 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::torsionstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST torsionstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 232 "MDTreeParser.g"
	
	vector<int> ivec;
	vector<RealType> dvec;
	TorsionStamp* currTorsionStamp = static_cast<TorsionStamp*>(blockStack.top());
	
#line 1825 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t83 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp57_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t83;
			_t = _t->getNextSibling();
#line 239 "MDTreeParser.g"
			currTorsionStamp->setMembers(ivec);
#line 1849 "MDTreeParser.cpp"
			break;
		}
		case GHOSTTORSION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t84 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp58_AST_in = _t;
			match(_t,GHOSTTORSION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t84;
			_t = _t->getNextSibling();
#line 240 "MDTreeParser.g"
			currTorsionStamp->overrideType("GhostTorsion", dvec);
#line 1864 "MDTreeParser.cpp"
			break;
		}
		case CUBIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t85 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp59_AST_in = _t;
			match(_t,CUBIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t85;
			_t = _t->getNextSibling();
#line 241 "MDTreeParser.g"
			currTorsionStamp->overrideType("Cubic", dvec);
#line 1879 "MDTreeParser.cpp"
			break;
		}
		case QUARTIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t86 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp60_AST_in = _t;
			match(_t,QUARTIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t86;
			_t = _t->getNextSibling();
#line 242 "MDTreeParser.g"
			currTorsionStamp->overrideType("Quartic", dvec);
#line 1894 "MDTreeParser.cpp"
			break;
		}
		case POLYNOMIAL:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t87 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp61_AST_in = _t;
			match(_t,POLYNOMIAL);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t87;
			_t = _t->getNextSibling();
#line 243 "MDTreeParser.g"
			currTorsionStamp->overrideType("Polynomial", dvec);
#line 1909 "MDTreeParser.cpp"
			break;
		}
		case CHARMM:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t88 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp62_AST_in = _t;
			match(_t,CHARMM);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t88;
			_t = _t->getNextSibling();
#line 244 "MDTreeParser.g"
			currTorsionStamp->overrideType("Charmm", dvec);
#line 1924 "MDTreeParser.cpp"
			break;
		}
		case OPLS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t89 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp63_AST_in = _t;
			match(_t,OPLS);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t89;
			_t = _t->getNextSibling();
#line 245 "MDTreeParser.g"
			currTorsionStamp->overrideType("Opls", dvec);
#line 1939 "MDTreeParser.cpp"
			break;
		}
		case TRAPPE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t90 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp64_AST_in = _t;
			match(_t,TRAPPE);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t90;
			_t = _t->getNextSibling();
#line 246 "MDTreeParser.g"
			currTorsionStamp->overrideType("Trappe", dvec);
#line 1954 "MDTreeParser.cpp"
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t91 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp65_AST_in = _t;
			match(_t,HARMONIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t91;
			_t = _t->getNextSibling();
#line 247 "MDTreeParser.g"
			currTorsionStamp->overrideType("Harmonic", dvec);
#line 1969 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::inversionstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST inversionstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 261 "MDTreeParser.g"
	
	int icent;
	vector<int> ivec;
	vector<RealType> dvec;
	InversionStamp* currInversionStamp = static_cast<InversionStamp*>(blockStack.top());
	
#line 1995 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case CENTER:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t97 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp66_AST_in = _t;
			match(_t,CENTER);
			_t = _t->getFirstChild();
			icent=intConst(_t);
			_t = _retTree;
			_t = __t97;
			_t = _t->getNextSibling();
#line 269 "MDTreeParser.g"
			currInversionStamp->setCenter(icent);
#line 2019 "MDTreeParser.cpp"
			break;
		}
		case SATELLITES:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t98 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp67_AST_in = _t;
			match(_t,SATELLITES);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t98;
			_t = _t->getNextSibling();
#line 270 "MDTreeParser.g"
			currInversionStamp->setSatellites(ivec);
#line 2034 "MDTreeParser.cpp"
			break;
		}
		case AMBERIMPROPER:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t99 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp68_AST_in = _t;
			match(_t,AMBERIMPROPER);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t99;
			_t = _t->getNextSibling();
#line 271 "MDTreeParser.g"
			currInversionStamp->overrideType("AmberImproper", dvec);
#line 2049 "MDTreeParser.cpp"
			break;
		}
		case IMPROPERCOSINE:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t100 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp69_AST_in = _t;
			match(_t,IMPROPERCOSINE);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t100;
			_t = _t->getNextSibling();
#line 272 "MDTreeParser.g"
			currInversionStamp->overrideType("ImproperCosine", dvec);
#line 2064 "MDTreeParser.cpp"
			break;
		}
		case HARMONIC:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t101 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp70_AST_in = _t;
			match(_t,HARMONIC);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t101;
			_t = _t->getNextSibling();
#line 273 "MDTreeParser.g"
			currInversionStamp->overrideType("Harmonic", dvec);
#line 2079 "MDTreeParser.cpp"
			break;
		}
		case CENTRALATOMHEIGHT:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t102 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp71_AST_in = _t;
			match(_t,CENTRALATOMHEIGHT);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t102;
			_t = _t->getNextSibling();
#line 274 "MDTreeParser.g"
			currInversionStamp->overrideType("CentralAtomHeight", dvec);
#line 2094 "MDTreeParser.cpp"
			break;
		}
		case DREIDING:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t103 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp72_AST_in = _t;
			match(_t,DREIDING);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t103;
			_t = _t->getNextSibling();
#line 275 "MDTreeParser.g"
			currInversionStamp->overrideType("Dreiding", dvec);
#line 2109 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::rigidbodystatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST rigidbodystatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 293 "MDTreeParser.g"
	
	vector<int> ivec;
	RigidBodyStamp* currRigidBodyStamp = static_cast<RigidBodyStamp*>(blockStack.top());
	
#line 2133 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t109 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp73_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t109;
			_t = _t->getNextSibling();
#line 299 "MDTreeParser.g"
			currRigidBodyStamp->setMembers(ivec);
#line 2157 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::cutoffgroupstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST cutoffgroupstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 312 "MDTreeParser.g"
	
	vector<int> ivec;
	CutoffGroupStamp* currCutoffGroupStamp = static_cast<CutoffGroupStamp*>(blockStack.top());
	
#line 2181 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t115 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp74_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t115;
			_t = _t->getNextSibling();
#line 318 "MDTreeParser.g"
			currCutoffGroupStamp->setMembers(ivec);
#line 2205 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::constraintstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST constraintstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 331 "MDTreeParser.g"
	
	vector<int> ivec;
	ConstraintStamp* currConstraintStamp = static_cast<ConstraintStamp*>(blockStack.top());
	
#line 2229 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t121 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp75_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t121;
			_t = _t->getNextSibling();
#line 337 "MDTreeParser.g"
			currConstraintStamp->setMembers(ivec);
#line 2253 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::nodesstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST nodesstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 351 "MDTreeParser.g"
	
	vector<int> ivec;
	NodesStamp* currNodesStamp = static_cast<NodesStamp*>(blockStack.top());
	
#line 2277 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case ASSIGNEQUAL:
		{
			assignment(_t);
			_t = _retTree;
			break;
		}
		case MEMBERS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t127 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp76_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t127;
			_t = _t->getNextSibling();
#line 357 "MDTreeParser.g"
			currNodesStamp->setMembers(ivec);
#line 2301 "MDTreeParser.cpp"
			break;
		}
		default:
		{
			throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
		}
		}
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::initializeASTFactory( ANTLR_USE_NAMESPACE(antlr)ASTFactory& )
{
}
const char* MDTreeParser::tokenNames[] = {
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
	"\"sequence\"",
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
	"\"nodes\"",
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

const unsigned long MDTreeParser::_tokenSet_0_data_[] = { 117571824UL, 65536UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "fragment" "flucQ" 
// "RNEMD" "minimizer" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDTreeParser::_tokenSet_1_data_[] = { 65280UL, 65536UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "constraint" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDTreeParser::_tokenSet_2_data_[] = { 65280UL, 73728UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "constraint" "nodes" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_2(_tokenSet_2_data_,4);
const unsigned long MDTreeParser::_tokenSet_3_data_[] = { 12582912UL, 69632UL, 0UL, 0UL };
// "position" "orientation" "charge" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_3(_tokenSet_3_data_,4);
const unsigned long MDTreeParser::_tokenSet_4_data_[] = { 4161273856UL, 65537UL, 0UL, 0UL };
// "members" "Fixed" "Harmonic" "Cubic" "Quartic" "Polynomial" "Morse" 
// ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_4(_tokenSet_4_data_,4);
const unsigned long MDTreeParser::_tokenSet_5_data_[] = { 4027056128UL, 65550UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostBend" "UreyBradley" 
// "Cosine" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_5(_tokenSet_5_data_,4);
const unsigned long MDTreeParser::_tokenSet_6_data_[] = { 4027056128UL, 65776UL, 0UL, 0UL };
// "members" "Harmonic" "Cubic" "Quartic" "Polynomial" "GhostTorsion" "Charmm" 
// "Opls" "Trappe" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_6(_tokenSet_6_data_,4);
const unsigned long MDTreeParser::_tokenSet_7_data_[] = { 271581184UL, 69376UL, 0UL, 0UL };
// "center" "satellites" "Harmonic" "AmberImproper" "ImproperCosine" "CentralAtomHeight" 
// "Dreiding" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_7(_tokenSet_7_data_,4);


