/* $ANTLR 2.7.7 (20120725): "MDTreeParser.g" -> "MDTreeParser.cpp"$ */
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
#line 72 "MDTreeParser.g"
		Component* currComponet = new Component(); blockStack.push(currComponet);
#line 151 "MDTreeParser.cpp"
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
#line 74 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addComponent(currComponet);
#line 174 "MDTreeParser.cpp"
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
#line 103 "MDTreeParser.g"
		MoleculeStamp* currMoleculeStamp = new MoleculeStamp(); blockStack.push(currMoleculeStamp);
#line 194 "MDTreeParser.cpp"
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
#line 105 "MDTreeParser.g"
		blockStack.top()->validate(); blockStack.pop(); currConf->addMoleculeStamp(currMoleculeStamp);
#line 217 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp6_AST_in = _t;
		match(_t,ZCONSTRAINT);
		_t = _t->getFirstChild();
#line 77 "MDTreeParser.g"
		ZConsStamp* currZConsStamp = new ZConsStamp(); blockStack.push(currZConsStamp);
#line 237 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp7_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t14;
		_t = _t->getNextSibling();
#line 79 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addZConsStamp(currZConsStamp);
#line 260 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp8_AST_in = _t;
		match(_t,RESTRAINT);
		_t = _t->getFirstChild();
#line 82 "MDTreeParser.g"
		RestraintStamp* currRestraintStamp = new RestraintStamp(); blockStack.push(currRestraintStamp);
#line 280 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t18;
		_t = _t->getNextSibling();
#line 84 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addRestraintStamp(currRestraintStamp);
#line 303 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp10_AST_in = _t;
		match(_t,FLUCQ);
		_t = _t->getFirstChild();
#line 87 "MDTreeParser.g"
		FluctuatingChargeParameters* flucQpars = new FluctuatingChargeParameters(); blockStack.push(flucQpars);
#line 323 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp11_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t22;
		_t = _t->getNextSibling();
#line 89 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addFluctuatingChargeParameters(flucQpars);
#line 346 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp12_AST_in = _t;
		match(_t,RNEMD);
		_t = _t->getFirstChild();
#line 92 "MDTreeParser.g"
		RNEMDParameters* rnemdPars = new RNEMDParameters(); blockStack.push(rnemdPars);
#line 366 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp13_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t26;
		_t = _t->getNextSibling();
#line 94 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addRNEMDParameters(rnemdPars);
#line 389 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp14_AST_in = _t;
		match(_t,MINIMIZER);
		_t = _t->getFirstChild();
#line 97 "MDTreeParser.g"
		MinimizerParameters* minimizerPars = new MinimizerParameters(); blockStack.push(minimizerPars);
#line 409 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp15_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t30;
		_t = _t->getNextSibling();
#line 99 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addMinimizerParameters(minimizerPars);
#line 432 "MDTreeParser.cpp"
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
#line 50 "MDTreeParser.g"
	
	int ival;
	RealType dval, x, y, z;
	Vector3d dvec;
	
#line 454 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case NUM_INT:
		case NUM_LONG:
		{
			ival=intConst(_t);
			_t = _retTree;
#line 56 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), ival);
#line 467 "MDTreeParser.cpp"
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			dval=floatConst(_t);
			_t = _retTree;
#line 57 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), dval);
#line 477 "MDTreeParser.cpp"
			break;
		}
		case LPAREN:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t8 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp16_AST_in = _t;
			match(_t,LPAREN);
			_t = _t->getFirstChild();
			x=doubleNumber(_t);
			_t = _retTree;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp17_AST_in = _t;
			match(_t,COMMA);
			_t = _t->getNextSibling();
			y=doubleNumber(_t);
			_t = _retTree;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp18_AST_in = _t;
			match(_t,COMMA);
			_t = _t->getNextSibling();
			z=doubleNumber(_t);
			_t = _retTree;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp19_AST_in = _t;
			match(_t,RPAREN);
			_t = _t->getNextSibling();
			_t = __t8;
			_t = _t->getNextSibling();
#line 59 "MDTreeParser.g"
			dvec.x() = x; 
			dvec.y() = y; 
			dvec.z() = z; 
			blockStack.top()->assign(id->getText(), dvec);
			
#line 509 "MDTreeParser.cpp"
			break;
		}
		case ID:
		{
			str1 = _t;
			match(_t,ID);
			_t = _t->getNextSibling();
#line 64 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), str1->getText());
#line 519 "MDTreeParser.cpp"
			break;
		}
		case StringLiteral:
		{
			str2 = _t;
			match(_t,StringLiteral);
			_t = _t->getNextSibling();
#line 65 "MDTreeParser.g"
			std::string s =  str2->getText();
			s = s.substr(1, s.length()-2);
			blockStack.top()->assign(id->getText(),s);
			
#line 532 "MDTreeParser.cpp"
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
#line 294 "MDTreeParser.g"
	int ival;
#line 552 "MDTreeParser.cpp"
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
#line 295 "MDTreeParser.g"
			ival = lexi_cast<int>(i1->getText());
#line 568 "MDTreeParser.cpp"
			break;
		}
		case NUM_LONG:
		{
			i2 = _t;
			match(_t,NUM_LONG);
			_t = _t->getNextSibling();
#line 296 "MDTreeParser.g"
			ival = lexi_cast<int>(i2->getText());
#line 578 "MDTreeParser.cpp"
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
#line 306 "MDTreeParser.g"
	RealType dval;
#line 599 "MDTreeParser.cpp"
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
#line 307 "MDTreeParser.g"
			dval = lexi_cast<RealType>(d1->getText());
#line 615 "MDTreeParser.cpp"
			break;
		}
		case NUM_DOUBLE:
		{
			d2 = _t;
			match(_t,NUM_DOUBLE);
			_t = _t->getNextSibling();
#line 308 "MDTreeParser.g"
			dval = lexi_cast<RealType>(d2->getText());
#line 625 "MDTreeParser.cpp"
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

RealType  MDTreeParser::doubleNumber(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 300 "MDTreeParser.g"
	RealType dval;
#line 646 "MDTreeParser.cpp"
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
#line 301 "MDTreeParser.g"
			dval = lexi_cast<RealType>(ic->getText());
#line 663 "MDTreeParser.cpp"
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			fc = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
			floatConst(_t);
			_t = _retTree;
#line 302 "MDTreeParser.g"
			dval = lexi_cast<RealType>(fc->getText());
#line 674 "MDTreeParser.cpp"
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

void MDTreeParser::moleculestatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculestatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
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
		case FRAGMENT:
		{
			fragmentblock(_t);
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

void MDTreeParser::atomblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST atomblock_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 119 "MDTreeParser.g"
	
	int index;
	
#line 773 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t39 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp20_AST_in = _t;
		match(_t,ATOM);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 123 "MDTreeParser.g"
		AtomStamp* currAtomStamp = new AtomStamp(index); blockStack.push(currAtomStamp);
#line 784 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == POSITION || _t->getType() == ORIENTATION || _t->getType() == ASSIGNEQUAL)) {
				atomstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop41;
			}
			
		}
		_loop41:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp21_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t39;
		_t = _t->getNextSibling();
#line 125 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addAtomStamp(currAtomStamp); 
		
#line 812 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t46 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp22_AST_in = _t;
		match(_t,BOND);
		_t = _t->getFirstChild();
#line 145 "MDTreeParser.g"
		BondStamp* currBondStamp = new BondStamp(); blockStack.push(currBondStamp);
#line 832 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				bondstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop48;
			}
			
		}
		_loop48:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp23_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t46;
		_t = _t->getNextSibling();
#line 147 "MDTreeParser.g"
		
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBondStamp(currBondStamp); 
		
#line 859 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t52 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST_in = _t;
		match(_t,BEND);
		_t = _t->getFirstChild();
#line 163 "MDTreeParser.g"
		BendStamp* currBendStamp = new BendStamp(); blockStack.push(currBendStamp);
#line 879 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				bendstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop54;
			}
			
		}
		_loop54:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t52;
		_t = _t->getNextSibling();
#line 165 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBendStamp(currBendStamp); 
		
#line 907 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t58 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp26_AST_in = _t;
		match(_t,TORSION);
		_t = _t->getFirstChild();
#line 182 "MDTreeParser.g"
		TorsionStamp* currTorsionStamp = new TorsionStamp(); blockStack.push(currTorsionStamp);
#line 927 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				torsionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop60;
			}
			
		}
		_loop60:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t58;
		_t = _t->getNextSibling();
#line 184 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addTorsionStamp(currTorsionStamp); 
		
#line 955 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t64 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST_in = _t;
		match(_t,INVERSION);
		_t = _t->getFirstChild();
#line 201 "MDTreeParser.g"
		InversionStamp* currInversionStamp = new InversionStamp(); blockStack.push(currInversionStamp);
#line 975 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == CENTER || _t->getType() == ASSIGNEQUAL)) {
				inversionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop66;
			}
			
		}
		_loop66:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t64;
		_t = _t->getNextSibling();
#line 203 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addInversionStamp(currInversionStamp); 
		
#line 1003 "MDTreeParser.cpp"
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
#line 220 "MDTreeParser.g"
	
	int index;
	
#line 1019 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t70 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp30_AST_in = _t;
		match(_t,RIGIDBODY);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 224 "MDTreeParser.g"
		RigidBodyStamp* currRigidBodyStamp = new RigidBodyStamp(index); blockStack.push(currRigidBodyStamp);
#line 1030 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				rigidbodystatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop72;
			}
			
		}
		_loop72:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp31_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t70;
		_t = _t->getNextSibling();
#line 226 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addRigidBodyStamp(currRigidBodyStamp); 
		
#line 1058 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t76 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST_in = _t;
		match(_t,CUTOFFGROUP);
		_t = _t->getFirstChild();
#line 243 "MDTreeParser.g"
		CutoffGroupStamp* currCutoffGroupStamp = new CutoffGroupStamp(); blockStack.push(currCutoffGroupStamp);
#line 1078 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				cutoffgroupstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop78;
			}
			
		}
		_loop78:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t76;
		_t = _t->getNextSibling();
#line 245 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addCutoffGroupStamp(currCutoffGroupStamp); 
		
#line 1106 "MDTreeParser.cpp"
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
#line 262 "MDTreeParser.g"
	int ival;
#line 1120 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t82 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp34_AST_in = _t;
		match(_t,FRAGMENT);
		_t = _t->getFirstChild();
		ival=intConst(_t);
		_t = _retTree;
#line 263 "MDTreeParser.g"
		FragmentStamp* currFragmentStamp = new FragmentStamp(ival); blockStack.push(currFragmentStamp);
#line 1131 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				fragmentstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop84;
			}
			
		}
		_loop84:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp35_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t82;
		_t = _t->getNextSibling();
#line 265 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addFragmentStamp(currFragmentStamp); 
		
#line 1159 "MDTreeParser.cpp"
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
#line 133 "MDTreeParser.g"
	
	vector<RealType> dvec;
	AtomStamp* currAtomStamp =  static_cast<AtomStamp*>(blockStack.top());
	
	
#line 1177 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t43 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp36_AST_in = _t;
			match(_t,POSITION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t43;
			_t = _t->getNextSibling();
#line 140 "MDTreeParser.g"
			currAtomStamp->setPosition(dvec);
#line 1201 "MDTreeParser.cpp"
			break;
		}
		case ORIENTATION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t44 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp37_AST_in = _t;
			match(_t,ORIENTATION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t44;
			_t = _t->getNextSibling();
#line 141 "MDTreeParser.g"
			currAtomStamp->setOrientation(dvec);
#line 1216 "MDTreeParser.cpp"
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

vector<RealType>  MDTreeParser::doubleNumberTuple(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 278 "MDTreeParser.g"
	vector<RealType> dvec;
#line 1236 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST doubleNumberTuple_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 278 "MDTreeParser.g"
	
	RealType dval;
	
#line 1242 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt88=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if (((_t->getType() >= NUM_INT && _t->getType() <= NUM_DOUBLE))) {
				dval=doubleNumber(_t);
				_t = _retTree;
#line 282 "MDTreeParser.g"
				dvec.push_back(dval);
#line 1255 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt88>=1 ) { goto _loop88; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt88++;
		}
		_loop88:;
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

void MDTreeParser::bondstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bondstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 154 "MDTreeParser.g"
	
	vector<int> ivec;
	BondStamp* currBondStamp = static_cast<BondStamp*>(blockStack.top());
	
#line 1282 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t50 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp38_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t50;
			_t = _t->getNextSibling();
#line 160 "MDTreeParser.g"
			currBondStamp->setMembers(ivec);
#line 1306 "MDTreeParser.cpp"
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
#line 286 "MDTreeParser.g"
	vector<int> ivec;
#line 1326 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST inttuple_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 286 "MDTreeParser.g"
	
	int ival;
	
#line 1332 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt91=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == NUM_INT || _t->getType() == NUM_LONG)) {
				ival=intConst(_t);
				_t = _retTree;
#line 290 "MDTreeParser.g"
				ivec.push_back(ival);
#line 1345 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt91>=1 ) { goto _loop91; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt91++;
		}
		_loop91:;
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
#line 173 "MDTreeParser.g"
	
	vector<int> ivec;
	BendStamp* currBendStamp = static_cast<BendStamp*>(blockStack.top());
	
#line 1372 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t56 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp39_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t56;
			_t = _t->getNextSibling();
#line 179 "MDTreeParser.g"
			currBendStamp->setMembers(ivec);
#line 1396 "MDTreeParser.cpp"
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
#line 192 "MDTreeParser.g"
	
	vector<int> ivec;
	TorsionStamp* currTorsionStamp = static_cast<TorsionStamp*>(blockStack.top());
	
#line 1420 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t62 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp40_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t62;
			_t = _t->getNextSibling();
#line 198 "MDTreeParser.g"
			currTorsionStamp->setMembers(ivec);
#line 1444 "MDTreeParser.cpp"
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
#line 211 "MDTreeParser.g"
	
	int icent;
	InversionStamp* currInversionStamp = static_cast<InversionStamp*>(blockStack.top());
	
#line 1468 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t68 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp41_AST_in = _t;
			match(_t,CENTER);
			_t = _t->getFirstChild();
			icent=intConst(_t);
			_t = _retTree;
			_t = __t68;
			_t = _t->getNextSibling();
#line 217 "MDTreeParser.g"
			currInversionStamp->setCenter(icent);
#line 1492 "MDTreeParser.cpp"
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
#line 234 "MDTreeParser.g"
	
	vector<int> ivec;
	RigidBodyStamp* currRigidBodyStamp = static_cast<RigidBodyStamp*>(blockStack.top());
	
#line 1516 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t74 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp42_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t74;
			_t = _t->getNextSibling();
#line 240 "MDTreeParser.g"
			currRigidBodyStamp->setMembers(ivec);
#line 1540 "MDTreeParser.cpp"
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
#line 253 "MDTreeParser.g"
	
	vector<int> ivec;
	CutoffGroupStamp* currCutoffGroupStamp = static_cast<CutoffGroupStamp*>(blockStack.top());
	
#line 1564 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t80 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp43_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t80;
			_t = _t->getNextSibling();
#line 259 "MDTreeParser.g"
			currCutoffGroupStamp->setMembers(ivec);
#line 1588 "MDTreeParser.cpp"
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

void MDTreeParser::fragmentstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentstatement_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		assignment(_t);
		_t = _retTree;
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
	"\"fragment\"",
	"\"members\"",
	"\"center\"",
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

const unsigned long MDTreeParser::_tokenSet_0_data_[] = { 40894704UL, 0UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" "minimizer" 
// ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDTreeParser::_tokenSet_1_data_[] = { 33619712UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_1(_tokenSet_1_data_,4);


