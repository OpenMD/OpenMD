/* $ANTLR 2.7.7 (20110725): "MDTreeParser.g" -> "MDTreeParser.cpp"$ */
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t9 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp2_AST_in = _t;
		match(_t,COMPONENT);
		_t = _t->getFirstChild();
#line 64 "MDTreeParser.g"
		Component* currComponet = new Component(); blockStack.push(currComponet);
#line 145 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop11;
			}
			
		}
		_loop11:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp3_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t9;
		_t = _t->getNextSibling();
#line 66 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addComponent(currComponet);
#line 168 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t29 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp4_AST_in = _t;
		match(_t,MOLECULE);
		_t = _t->getFirstChild();
#line 90 "MDTreeParser.g"
		MoleculeStamp* currMoleculeStamp = new MoleculeStamp(); blockStack.push(currMoleculeStamp);
#line 188 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_1.member(_t->getType()))) {
				moleculestatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop31;
			}
			
		}
		_loop31:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp5_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t29;
		_t = _t->getNextSibling();
#line 92 "MDTreeParser.g"
		blockStack.top()->validate(); blockStack.pop(); currConf->addMoleculeStamp(currMoleculeStamp);
#line 211 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t13 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp6_AST_in = _t;
		match(_t,ZCONSTRAINT);
		_t = _t->getFirstChild();
#line 69 "MDTreeParser.g"
		ZConsStamp* currZConsStamp = new ZConsStamp(); blockStack.push(currZConsStamp);
#line 231 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop15;
			}
			
		}
		_loop15:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp7_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t13;
		_t = _t->getNextSibling();
#line 71 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addZConsStamp(currZConsStamp);
#line 254 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t17 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp8_AST_in = _t;
		match(_t,RESTRAINT);
		_t = _t->getFirstChild();
#line 74 "MDTreeParser.g"
		RestraintStamp* currRestraintStamp = new RestraintStamp(); blockStack.push(currRestraintStamp);
#line 274 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop19;
			}
			
		}
		_loop19:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t17;
		_t = _t->getNextSibling();
#line 76 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addRestraintStamp(currRestraintStamp);
#line 297 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t21 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp10_AST_in = _t;
		match(_t,FLUCQ);
		_t = _t->getFirstChild();
#line 79 "MDTreeParser.g"
		FluctuatingChargeParameters* flucQpars = new FluctuatingChargeParameters(); blockStack.push(flucQpars);
#line 317 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop23;
			}
			
		}
		_loop23:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp11_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t21;
		_t = _t->getNextSibling();
#line 81 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addFluctuatingChargeParameters(flucQpars);
#line 340 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t25 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp12_AST_in = _t;
		match(_t,RNEMD);
		_t = _t->getFirstChild();
#line 84 "MDTreeParser.g"
		RNEMDParameters* rnemdPars = new RNEMDParameters(); blockStack.push(rnemdPars);
#line 360 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				assignment(_t);
				_t = _retTree;
			}
			else {
				goto _loop27;
			}
			
		}
		_loop27:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp13_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t25;
		_t = _t->getNextSibling();
#line 86 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addRNEMDParameters(rnemdPars);
#line 383 "MDTreeParser.cpp"
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
#line 49 "MDTreeParser.g"
	
	int ival;
	RealType dval;
	
#line 404 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case NUM_INT:
		case NUM_LONG:
		{
			ival=intConst(_t);
			_t = _retTree;
#line 54 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), ival);
#line 417 "MDTreeParser.cpp"
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			dval=floatConst(_t);
			_t = _retTree;
#line 55 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), dval);
#line 427 "MDTreeParser.cpp"
			break;
		}
		case ID:
		{
			str1 = _t;
			match(_t,ID);
			_t = _t->getNextSibling();
#line 56 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), str1->getText());
#line 437 "MDTreeParser.cpp"
			break;
		}
		case StringLiteral:
		{
			str2 = _t;
			match(_t,StringLiteral);
			_t = _t->getNextSibling();
#line 57 "MDTreeParser.g"
			std::string s =  str2->getText();
			s = s.substr(1, s.length()-2);
			blockStack.top()->assign(id->getText(),s);
			
#line 450 "MDTreeParser.cpp"
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
#line 280 "MDTreeParser.g"
	int ival;
#line 470 "MDTreeParser.cpp"
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
#line 281 "MDTreeParser.g"
			ival = lexi_cast<int>(i1->getText());
#line 486 "MDTreeParser.cpp"
			break;
		}
		case NUM_LONG:
		{
			i2 = _t;
			match(_t,NUM_LONG);
			_t = _t->getNextSibling();
#line 282 "MDTreeParser.g"
			ival = lexi_cast<int>(i2->getText());
#line 496 "MDTreeParser.cpp"
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
#line 294 "MDTreeParser.g"
	RealType dval;
#line 517 "MDTreeParser.cpp"
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
#line 295 "MDTreeParser.g"
			dval = lexi_cast<RealType>(d1->getText());
#line 533 "MDTreeParser.cpp"
			break;
		}
		case NUM_DOUBLE:
		{
			d2 = _t;
			match(_t,NUM_DOUBLE);
			_t = _t->getNextSibling();
#line 296 "MDTreeParser.g"
			dval = lexi_cast<RealType>(d2->getText());
#line 543 "MDTreeParser.cpp"
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
#line 106 "MDTreeParser.g"
	
	int index;
	
#line 642 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t34 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp14_AST_in = _t;
		match(_t,ATOM);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 110 "MDTreeParser.g"
		AtomStamp* currAtomStamp = new AtomStamp(index); blockStack.push(currAtomStamp);
#line 653 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == POSITION || _t->getType() == ORIENTATION || _t->getType() == ASSIGNEQUAL)) {
				atomstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop36;
			}
			
		}
		_loop36:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp15_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t34;
		_t = _t->getNextSibling();
#line 112 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addAtomStamp(currAtomStamp); 
		
#line 681 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t41 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp16_AST_in = _t;
		match(_t,BOND);
		_t = _t->getFirstChild();
#line 132 "MDTreeParser.g"
		BondStamp* currBondStamp = new BondStamp(); blockStack.push(currBondStamp);
#line 701 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				bondstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop43;
			}
			
		}
		_loop43:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp17_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t41;
		_t = _t->getNextSibling();
#line 134 "MDTreeParser.g"
		
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBondStamp(currBondStamp); 
		
#line 728 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t47 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp18_AST_in = _t;
		match(_t,BEND);
		_t = _t->getFirstChild();
#line 150 "MDTreeParser.g"
		BendStamp* currBendStamp = new BendStamp(); blockStack.push(currBendStamp);
#line 748 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				bendstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop49;
			}
			
		}
		_loop49:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp19_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t47;
		_t = _t->getNextSibling();
#line 152 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBendStamp(currBendStamp); 
		
#line 776 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t53 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp20_AST_in = _t;
		match(_t,TORSION);
		_t = _t->getFirstChild();
#line 169 "MDTreeParser.g"
		TorsionStamp* currTorsionStamp = new TorsionStamp(); blockStack.push(currTorsionStamp);
#line 796 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				torsionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop55;
			}
			
		}
		_loop55:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp21_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t53;
		_t = _t->getNextSibling();
#line 171 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addTorsionStamp(currTorsionStamp); 
		
#line 824 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t59 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp22_AST_in = _t;
		match(_t,INVERSION);
		_t = _t->getFirstChild();
#line 188 "MDTreeParser.g"
		InversionStamp* currInversionStamp = new InversionStamp(); blockStack.push(currInversionStamp);
#line 844 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == CENTER || _t->getType() == ASSIGNEQUAL)) {
				inversionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop61;
			}
			
		}
		_loop61:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp23_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t59;
		_t = _t->getNextSibling();
#line 190 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addInversionStamp(currInversionStamp); 
		
#line 872 "MDTreeParser.cpp"
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
#line 207 "MDTreeParser.g"
	
	int index;
	
#line 888 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t65 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST_in = _t;
		match(_t,RIGIDBODY);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 211 "MDTreeParser.g"
		RigidBodyStamp* currRigidBodyStamp = new RigidBodyStamp(index); blockStack.push(currRigidBodyStamp);
#line 899 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				rigidbodystatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop67;
			}
			
		}
		_loop67:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t65;
		_t = _t->getNextSibling();
#line 213 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addRigidBodyStamp(currRigidBodyStamp); 
		
#line 927 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST __t71 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp26_AST_in = _t;
		match(_t,CUTOFFGROUP);
		_t = _t->getFirstChild();
#line 230 "MDTreeParser.g"
		CutoffGroupStamp* currCutoffGroupStamp = new CutoffGroupStamp(); blockStack.push(currCutoffGroupStamp);
#line 947 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				cutoffgroupstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop73;
			}
			
		}
		_loop73:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t71;
		_t = _t->getNextSibling();
#line 232 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addCutoffGroupStamp(currCutoffGroupStamp); 
		
#line 975 "MDTreeParser.cpp"
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
#line 249 "MDTreeParser.g"
	int ival;
#line 989 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t77 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST_in = _t;
		match(_t,FRAGMENT);
		_t = _t->getFirstChild();
		ival=intConst(_t);
		_t = _retTree;
#line 250 "MDTreeParser.g"
		FragmentStamp* currFragmentStamp = new FragmentStamp(ival); blockStack.push(currFragmentStamp);
#line 1000 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				fragmentstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop79;
			}
			
		}
		_loop79:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t77;
		_t = _t->getNextSibling();
#line 252 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addFragmentStamp(currFragmentStamp); 
		
#line 1028 "MDTreeParser.cpp"
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
#line 120 "MDTreeParser.g"
	
	vector<RealType> dvec;
	AtomStamp* currAtomStamp =  static_cast<AtomStamp*>(blockStack.top());
	
	
#line 1046 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t38 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp30_AST_in = _t;
			match(_t,POSITION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t38;
			_t = _t->getNextSibling();
#line 127 "MDTreeParser.g"
			currAtomStamp->setPosition(dvec);
#line 1070 "MDTreeParser.cpp"
			break;
		}
		case ORIENTATION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t39 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp31_AST_in = _t;
			match(_t,ORIENTATION);
			_t = _t->getFirstChild();
			dvec=doubleNumberTuple(_t);
			_t = _retTree;
			_t = __t39;
			_t = _t->getNextSibling();
#line 128 "MDTreeParser.g"
			currAtomStamp->setOrientation(dvec);
#line 1085 "MDTreeParser.cpp"
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
#line 265 "MDTreeParser.g"
	vector<RealType> dvec;
#line 1105 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST doubleNumberTuple_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 265 "MDTreeParser.g"
	
	RealType dval;
	
#line 1111 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt83=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if (((_t->getType() >= NUM_INT && _t->getType() <= NUM_DOUBLE))) {
				dval=doubleNumber(_t);
				_t = _retTree;
#line 269 "MDTreeParser.g"
				dvec.push_back(dval);
#line 1124 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt83>=1 ) { goto _loop83; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt83++;
		}
		_loop83:;
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
#line 141 "MDTreeParser.g"
	
	vector<int> ivec;
	BondStamp* currBondStamp = static_cast<BondStamp*>(blockStack.top());
	
#line 1151 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t45 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp32_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t45;
			_t = _t->getNextSibling();
#line 147 "MDTreeParser.g"
			currBondStamp->setMembers(ivec);
#line 1175 "MDTreeParser.cpp"
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
#line 272 "MDTreeParser.g"
	vector<int> ivec;
#line 1195 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST inttuple_AST_in = (_t == ANTLR_USE_NAMESPACE(antlr)RefAST(ASTNULL)) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 272 "MDTreeParser.g"
	
	int ival;
	
#line 1201 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt86=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == NUM_INT || _t->getType() == NUM_LONG)) {
				ival=intConst(_t);
				_t = _retTree;
#line 276 "MDTreeParser.g"
				ivec.push_back(ival);
#line 1214 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt86>=1 ) { goto _loop86; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt86++;
		}
		_loop86:;
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
#line 160 "MDTreeParser.g"
	
	vector<int> ivec;
	BendStamp* currBendStamp = static_cast<BendStamp*>(blockStack.top());
	
#line 1241 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t51 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp33_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t51;
			_t = _t->getNextSibling();
#line 166 "MDTreeParser.g"
			currBendStamp->setMembers(ivec);
#line 1265 "MDTreeParser.cpp"
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
#line 179 "MDTreeParser.g"
	
	vector<int> ivec;
	TorsionStamp* currTorsionStamp = static_cast<TorsionStamp*>(blockStack.top());
	
#line 1289 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t57 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp34_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t57;
			_t = _t->getNextSibling();
#line 185 "MDTreeParser.g"
			currTorsionStamp->setMembers(ivec);
#line 1313 "MDTreeParser.cpp"
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
#line 198 "MDTreeParser.g"
	
	int icent;
	InversionStamp* currInversionStamp = static_cast<InversionStamp*>(blockStack.top());
	
#line 1337 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t63 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp35_AST_in = _t;
			match(_t,CENTER);
			_t = _t->getFirstChild();
			icent=intConst(_t);
			_t = _retTree;
			_t = __t63;
			_t = _t->getNextSibling();
#line 204 "MDTreeParser.g"
			currInversionStamp->setCenter(icent);
#line 1361 "MDTreeParser.cpp"
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
#line 221 "MDTreeParser.g"
	
	vector<int> ivec;
	RigidBodyStamp* currRigidBodyStamp = static_cast<RigidBodyStamp*>(blockStack.top());
	
#line 1385 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t69 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp36_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t69;
			_t = _t->getNextSibling();
#line 227 "MDTreeParser.g"
			currRigidBodyStamp->setMembers(ivec);
#line 1409 "MDTreeParser.cpp"
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
#line 240 "MDTreeParser.g"
	
	vector<int> ivec;
	CutoffGroupStamp* currCutoffGroupStamp = static_cast<CutoffGroupStamp*>(blockStack.top());
	
#line 1433 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t75 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp37_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t75;
			_t = _t->getNextSibling();
#line 246 "MDTreeParser.g"
			currCutoffGroupStamp->setMembers(ivec);
#line 1457 "MDTreeParser.cpp"
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

RealType  MDTreeParser::doubleNumber(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 286 "MDTreeParser.g"
	RealType dval;
#line 1492 "MDTreeParser.cpp"
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
#line 288 "MDTreeParser.g"
			dval = lexi_cast<RealType>(ic->getText());
#line 1509 "MDTreeParser.cpp"
			break;
		}
		case NUM_FLOAT:
		case NUM_DOUBLE:
		{
			fc = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
			floatConst(_t);
			_t = _retTree;
#line 289 "MDTreeParser.g"
			dval = lexi_cast<RealType>(fc->getText());
#line 1520 "MDTreeParser.cpp"
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

const unsigned long MDTreeParser::_tokenSet_0_data_[] = { 19923184UL, 0UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" "restraint" "flucQ" "RNEMD" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDTreeParser::_tokenSet_1_data_[] = { 16842496UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "inversion" "rigidBody" "cutoffGroup" 
// "fragment" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_1(_tokenSet_1_data_,4);


