/* $ANTLR 2.7.5 (20050406): "MDTreeParser.g" -> "MDTreeParser.cpp"$ */
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
	ANTLR_USE_NAMESPACE(antlr)RefAST mdfile_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
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
	ANTLR_USE_NAMESPACE(antlr)RefAST statement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
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
	ANTLR_USE_NAMESPACE(antlr)RefAST assignment_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
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
	ANTLR_USE_NAMESPACE(antlr)RefAST componentblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t13 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp2_AST_in = _t;
		match(_t,COMPONENT);
		_t = _t->getFirstChild();
#line 69 "MDTreeParser.g"
		Component* currComponet = new Component(); blockStack.push(currComponet);
#line 127 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp3_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t13;
		_t = _t->getNextSibling();
#line 71 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addComponent(currComponet);
#line 150 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::moleculeblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculeblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t21 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp4_AST_in = _t;
		match(_t,MOLECULE);
		_t = _t->getFirstChild();
#line 79 "MDTreeParser.g"
		MoleculeStamp* currMoleculeStamp = new MoleculeStamp(); blockStack.push(currMoleculeStamp);
#line 170 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_1.member(_t->getType()))) {
				moleculestatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop23;
			}
			
		}
		_loop23:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp5_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t21;
		_t = _t->getNextSibling();
#line 81 "MDTreeParser.g"
		blockStack.top()->validate(); blockStack.pop(); currConf->addMoleculeStamp(currMoleculeStamp);
#line 193 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::zconstraintblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST zconstraintblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t17 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp6_AST_in = _t;
		match(_t,ZCONSTRAINT);
		_t = _t->getFirstChild();
#line 74 "MDTreeParser.g"
		ZConsStamp* currZConsStamp = new ZConsStamp(); blockStack.push(currZConsStamp);
#line 213 "MDTreeParser.cpp"
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
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp7_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t17;
		_t = _t->getNextSibling();
#line 76 "MDTreeParser.g"
		blockStack.top()->validate();blockStack.pop(); currConf->addZConsStamp(currZConsStamp);
#line 236 "MDTreeParser.cpp"
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
	ANTLR_USE_NAMESPACE(antlr)RefAST constant_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST str1 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST str2 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case OCTALINT:
		case DECIMALINT:
		case HEXADECIMALINT:
		case MINUS:
		case FLOATONE:
		case FLOATTWO:
		{
			signedIntOrFloat(_t,id);
			_t = _retTree;
			break;
		}
		case ID:
		{
			str1 = _t;
			match(_t,ID);
			_t = _t->getNextSibling();
#line 49 "MDTreeParser.g"
			blockStack.top()->assign(id->getText(), str1->getText());
#line 275 "MDTreeParser.cpp"
			break;
		}
		case StringLiteral:
		{
			str2 = _t;
			match(_t,StringLiteral);
			_t = _t->getNextSibling();
#line 50 "MDTreeParser.g"
			std::string s =  str2->getText();
			s = s.substr(1, s.length()-2);
			blockStack.top()->assign(id->getText(),s);
			
#line 288 "MDTreeParser.cpp"
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

void MDTreeParser::signedIntOrFloat(ANTLR_USE_NAMESPACE(antlr)RefAST _t,
	ANTLR_USE_NAMESPACE(antlr)RefAST id
) {
	ANTLR_USE_NAMESPACE(antlr)RefAST signedIntOrFloat_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST icMinus = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fcMinus = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST ic = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fc = ANTLR_USE_NAMESPACE(antlr)nullAST;
#line 56 "MDTreeParser.g"
	
	int ival;
	double dval;
	
#line 318 "MDTreeParser.cpp"
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case MINUS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t9 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp8_AST_in = _t;
			match(_t,MINUS);
			_t = _t->getFirstChild();
			{
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			switch ( _t->getType()) {
			case OCTALINT:
			case DECIMALINT:
			case HEXADECIMALINT:
			{
				icMinus = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				intConst(_t);
				_t = _retTree;
#line 61 "MDTreeParser.g"
				ival = lexi_cast<int>(icMinus->getText()); ival = -ival; blockStack.top()->assign(id->getText(), ival);
#line 343 "MDTreeParser.cpp"
				break;
			}
			case FLOATONE:
			case FLOATTWO:
			{
				fcMinus = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				floatConst(_t);
				_t = _retTree;
				break;
			}
			default:
			{
				throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
			}
			}
			}
#line 62 "MDTreeParser.g"
			dval = lexi_cast<double>(fcMinus->getText());dval = -dval;  blockStack.top()->assign(id->getText(), dval);
#line 362 "MDTreeParser.cpp"
			_t = __t9;
			_t = _t->getNextSibling();
			break;
		}
		case OCTALINT:
		case DECIMALINT:
		case HEXADECIMALINT:
		case FLOATONE:
		case FLOATTWO:
		{
			{
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			switch ( _t->getType()) {
			case OCTALINT:
			case DECIMALINT:
			case HEXADECIMALINT:
			{
				ic = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				intConst(_t);
				_t = _retTree;
#line 64 "MDTreeParser.g"
				ival = lexi_cast<int>(ic->getText()); blockStack.top()->assign(id->getText(), ival);
#line 386 "MDTreeParser.cpp"
				break;
			}
			case FLOATONE:
			case FLOATTWO:
			{
				fc = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				floatConst(_t);
				_t = _retTree;
#line 65 "MDTreeParser.g"
				dval = lexi_cast<double>(fc->getText());  blockStack.top()->assign(id->getText(), dval);
#line 397 "MDTreeParser.cpp"
				break;
			}
			default:
			{
				throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
			}
			}
			}
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
#line 249 "MDTreeParser.g"
	int ival;
#line 425 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST intConst_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST oival = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST dival = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST hival = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case OCTALINT:
		{
			oival = _t;
			match(_t,OCTALINT);
			_t = _t->getNextSibling();
#line 250 "MDTreeParser.g"
			ival = lexi_cast<int>(oival->getText());
#line 442 "MDTreeParser.cpp"
			break;
		}
		case DECIMALINT:
		{
			dival = _t;
			match(_t,DECIMALINT);
			_t = _t->getNextSibling();
#line 251 "MDTreeParser.g"
			ival = lexi_cast<int>(dival->getText());
#line 452 "MDTreeParser.cpp"
			break;
		}
		case HEXADECIMALINT:
		{
			hival = _t;
			match(_t,HEXADECIMALINT);
			_t = _t->getNextSibling();
#line 252 "MDTreeParser.g"
			ival = lexi_cast<int>(hival->getText());
#line 462 "MDTreeParser.cpp"
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

double  MDTreeParser::floatConst(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 266 "MDTreeParser.g"
	double dval;
#line 483 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST floatConst_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST d1 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST d2 = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case FLOATONE:
		{
			d1 = _t;
			match(_t,FLOATONE);
			_t = _t->getNextSibling();
#line 267 "MDTreeParser.g"
			dval = lexi_cast<double>(d1->getText());
#line 499 "MDTreeParser.cpp"
			break;
		}
		case FLOATTWO:
		{
			d2 = _t;
			match(_t,FLOATTWO);
			_t = _t->getNextSibling();
#line 268 "MDTreeParser.g"
			dval = lexi_cast<double>(d2->getText());
#line 509 "MDTreeParser.cpp"
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
	ANTLR_USE_NAMESPACE(antlr)RefAST moleculestatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
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
	ANTLR_USE_NAMESPACE(antlr)RefAST atomblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 94 "MDTreeParser.g"
	
	int index;
	
#line 602 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t26 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp9_AST_in = _t;
		match(_t,ATOM);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 98 "MDTreeParser.g"
		AtomStamp* currAtomStamp = new AtomStamp(index); blockStack.push(currAtomStamp);
#line 613 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == POSITION || _t->getType() == ORIENTATION || _t->getType() == ASSIGNEQUAL)) {
				atomstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop28;
			}
			
		}
		_loop28:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp10_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t26;
		_t = _t->getNextSibling();
#line 100 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addAtomStamp(currAtomStamp); 
		
#line 641 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::bondblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bondblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t33 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp11_AST_in = _t;
		match(_t,BOND);
		_t = _t->getFirstChild();
#line 120 "MDTreeParser.g"
		BondStamp* currBondStamp = new BondStamp(); blockStack.push(currBondStamp);
#line 661 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				bondstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop35;
			}
			
		}
		_loop35:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp12_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t33;
		_t = _t->getNextSibling();
#line 122 "MDTreeParser.g"
		
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBondStamp(currBondStamp); 
		
#line 688 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::bendblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST bendblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t39 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp13_AST_in = _t;
		match(_t,BEND);
		_t = _t->getFirstChild();
#line 138 "MDTreeParser.g"
		BendStamp* currBendStamp = new BendStamp(); blockStack.push(currBendStamp);
#line 708 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				bendstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop41;
			}
			
		}
		_loop41:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp14_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t39;
		_t = _t->getNextSibling();
#line 140 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addBendStamp(currBendStamp); 
		
#line 736 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::torsionblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST torsionblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t45 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp15_AST_in = _t;
		match(_t,TORSION);
		_t = _t->getFirstChild();
#line 157 "MDTreeParser.g"
		TorsionStamp* currTorsionStamp = new TorsionStamp(); blockStack.push(currTorsionStamp);
#line 756 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				torsionstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop47;
			}
			
		}
		_loop47:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp16_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t45;
		_t = _t->getNextSibling();
#line 159 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addTorsionStamp(currTorsionStamp); 
		
#line 784 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::rigidbodyblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST rigidbodyblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 176 "MDTreeParser.g"
	
	int index;
	
#line 800 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t51 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp17_AST_in = _t;
		match(_t,RIGIDBODY);
		_t = _t->getFirstChild();
		index=intConst(_t);
		_t = _retTree;
#line 180 "MDTreeParser.g"
		RigidBodyStamp* currRigidBodyStamp = new RigidBodyStamp(index); blockStack.push(currRigidBodyStamp);
#line 811 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				rigidbodystatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop53;
			}
			
		}
		_loop53:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp18_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t51;
		_t = _t->getNextSibling();
#line 182 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addRigidBodyStamp(currRigidBodyStamp); 
		
#line 839 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::cutoffgroupblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST cutoffgroupblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t57 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp19_AST_in = _t;
		match(_t,CUTOFFGROUP);
		_t = _t->getFirstChild();
#line 199 "MDTreeParser.g"
		CutoffGroupStamp* currCutoffGroupStamp = new CutoffGroupStamp(); blockStack.push(currCutoffGroupStamp);
#line 859 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == MEMBERS || _t->getType() == ASSIGNEQUAL)) {
				cutoffgroupstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop59;
			}
			
		}
		_loop59:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp20_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t57;
		_t = _t->getNextSibling();
#line 201 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addCutoffGroupStamp(currCutoffGroupStamp); 
		
#line 887 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::fragmentblock(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentblock_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 218 "MDTreeParser.g"
	int ival;
#line 901 "MDTreeParser.cpp"
	
	try {      // for error handling
		ANTLR_USE_NAMESPACE(antlr)RefAST __t63 = _t;
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp21_AST_in = _t;
		match(_t,FRAGMENT);
		_t = _t->getFirstChild();
		ival=intConst(_t);
		_t = _retTree;
#line 219 "MDTreeParser.g"
		FragmentStamp* currFragmentStamp = new FragmentStamp(ival); blockStack.push(currFragmentStamp);
#line 912 "MDTreeParser.cpp"
		{ // ( ... )*
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_t->getType() == ASSIGNEQUAL)) {
				fragmentstatement(_t);
				_t = _retTree;
			}
			else {
				goto _loop65;
			}
			
		}
		_loop65:;
		} // ( ... )*
		ANTLR_USE_NAMESPACE(antlr)RefAST tmp22_AST_in = _t;
		match(_t,ENDBLOCK);
		_t = _t->getNextSibling();
		_t = __t63;
		_t = _t->getNextSibling();
#line 221 "MDTreeParser.g"
		
		blockStack.top()->validate();
		blockStack.pop(); 
		MoleculeStamp* currMoleculeStamp = static_cast<MoleculeStamp*>(blockStack.top());
		currMoleculeStamp->addFragmentStamp(currFragmentStamp); 
		
#line 940 "MDTreeParser.cpp"
	}
	catch (ANTLR_USE_NAMESPACE(antlr)RecognitionException& ex) {
		reportError(ex);
		if ( _t != ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = _t->getNextSibling();
	}
	_retTree = _t;
}

void MDTreeParser::atomstatement(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
	ANTLR_USE_NAMESPACE(antlr)RefAST atomstatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 108 "MDTreeParser.g"
	
	vector<double> dvec;
	AtomStamp* currAtomStamp =  static_cast<AtomStamp*>(blockStack.top());
	
	
#line 958 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t30 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp23_AST_in = _t;
			match(_t,POSITION);
			_t = _t->getFirstChild();
			dvec=signedNumberTuple(_t);
			_t = _retTree;
			_t = __t30;
			_t = _t->getNextSibling();
#line 115 "MDTreeParser.g"
			currAtomStamp->setPosition(dvec);
#line 982 "MDTreeParser.cpp"
			break;
		}
		case ORIENTATION:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t31 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp24_AST_in = _t;
			match(_t,ORIENTATION);
			_t = _t->getFirstChild();
			dvec=signedNumberTuple(_t);
			_t = _retTree;
			_t = __t31;
			_t = _t->getNextSibling();
#line 116 "MDTreeParser.g"
			currAtomStamp->setOrientation(dvec);
#line 997 "MDTreeParser.cpp"
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

vector<double>  MDTreeParser::signedNumberTuple(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 234 "MDTreeParser.g"
	vector<double> dvec;
#line 1017 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST signedNumberTuple_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 234 "MDTreeParser.g"
	
	double dval;
	
#line 1023 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt69=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if ((_tokenSet_2.member(_t->getType()))) {
				dval=signedNumber(_t);
				_t = _retTree;
#line 238 "MDTreeParser.g"
				dvec.push_back(dval);
#line 1036 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt69>=1 ) { goto _loop69; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt69++;
		}
		_loop69:;
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
	ANTLR_USE_NAMESPACE(antlr)RefAST bondstatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 129 "MDTreeParser.g"
	
	vector<int> ivec;
	BondStamp* currBondStamp = static_cast<BondStamp*>(blockStack.top());
	
#line 1063 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t37 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp25_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t37;
			_t = _t->getNextSibling();
#line 135 "MDTreeParser.g"
			currBondStamp->setMembers(ivec);
#line 1087 "MDTreeParser.cpp"
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
#line 241 "MDTreeParser.g"
	vector<int> ivec;
#line 1107 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST inttuple_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 241 "MDTreeParser.g"
	
	int ival;
	
#line 1113 "MDTreeParser.cpp"
	
	try {      // for error handling
		{ // ( ... )+
		int _cnt72=0;
		for (;;) {
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			if (((_t->getType() >= OCTALINT && _t->getType() <= HEXADECIMALINT))) {
				ival=intConst(_t);
				_t = _retTree;
#line 245 "MDTreeParser.g"
				ivec.push_back(ival);
#line 1126 "MDTreeParser.cpp"
			}
			else {
				if ( _cnt72>=1 ) { goto _loop72; } else {throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);}
			}
			
			_cnt72++;
		}
		_loop72:;
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
	ANTLR_USE_NAMESPACE(antlr)RefAST bendstatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 148 "MDTreeParser.g"
	
	vector<int> ivec;
	BendStamp* currBendStamp = static_cast<BendStamp*>(blockStack.top());
	
#line 1153 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t43 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp26_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t43;
			_t = _t->getNextSibling();
#line 154 "MDTreeParser.g"
			currBendStamp->setMembers(ivec);
#line 1177 "MDTreeParser.cpp"
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
	ANTLR_USE_NAMESPACE(antlr)RefAST torsionstatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 167 "MDTreeParser.g"
	
	vector<int> ivec;
	TorsionStamp* currTorsionStamp = static_cast<TorsionStamp*>(blockStack.top());
	
#line 1201 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t49 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp27_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t49;
			_t = _t->getNextSibling();
#line 173 "MDTreeParser.g"
			currTorsionStamp->setMembers(ivec);
#line 1225 "MDTreeParser.cpp"
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
	ANTLR_USE_NAMESPACE(antlr)RefAST rigidbodystatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 190 "MDTreeParser.g"
	
	vector<int> ivec;
	RigidBodyStamp* currRigidBodyStamp = static_cast<RigidBodyStamp*>(blockStack.top());
	
#line 1249 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t55 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp28_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t55;
			_t = _t->getNextSibling();
#line 196 "MDTreeParser.g"
			currRigidBodyStamp->setMembers(ivec);
#line 1273 "MDTreeParser.cpp"
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
	ANTLR_USE_NAMESPACE(antlr)RefAST cutoffgroupstatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
#line 209 "MDTreeParser.g"
	
	vector<int> ivec;
	CutoffGroupStamp* currCutoffGroupStamp = static_cast<CutoffGroupStamp*>(blockStack.top());
	
#line 1297 "MDTreeParser.cpp"
	
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
			ANTLR_USE_NAMESPACE(antlr)RefAST __t61 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp29_AST_in = _t;
			match(_t,MEMBERS);
			_t = _t->getFirstChild();
			ivec=inttuple(_t);
			_t = _retTree;
			_t = __t61;
			_t = _t->getNextSibling();
#line 215 "MDTreeParser.g"
			currCutoffGroupStamp->setMembers(ivec);
#line 1321 "MDTreeParser.cpp"
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
	ANTLR_USE_NAMESPACE(antlr)RefAST fragmentstatement_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	
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

double  MDTreeParser::signedNumber(ANTLR_USE_NAMESPACE(antlr)RefAST _t) {
#line 256 "MDTreeParser.g"
	double dval;
#line 1356 "MDTreeParser.cpp"
	ANTLR_USE_NAMESPACE(antlr)RefAST signedNumber_AST_in = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
	ANTLR_USE_NAMESPACE(antlr)RefAST icMinus = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fcMinus = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST ic = ANTLR_USE_NAMESPACE(antlr)nullAST;
	ANTLR_USE_NAMESPACE(antlr)RefAST fc = ANTLR_USE_NAMESPACE(antlr)nullAST;
	
	try {      // for error handling
		if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
			_t = ASTNULL;
		switch ( _t->getType()) {
		case MINUS:
		{
			ANTLR_USE_NAMESPACE(antlr)RefAST __t75 = _t;
			ANTLR_USE_NAMESPACE(antlr)RefAST tmp30_AST_in = _t;
			match(_t,MINUS);
			_t = _t->getFirstChild();
			{
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			switch ( _t->getType()) {
			case OCTALINT:
			case DECIMALINT:
			case HEXADECIMALINT:
			{
				icMinus = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				intConst(_t);
				_t = _retTree;
#line 257 "MDTreeParser.g"
				dval = lexi_cast<double>(icMinus->getText()); dval = -dval;
#line 1386 "MDTreeParser.cpp"
				break;
			}
			case FLOATONE:
			case FLOATTWO:
			{
				fcMinus = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				floatConst(_t);
				_t = _retTree;
				break;
			}
			default:
			{
				throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
			}
			}
			}
#line 258 "MDTreeParser.g"
			dval = lexi_cast<double>(fcMinus->getText());dval = -dval;
#line 1405 "MDTreeParser.cpp"
			_t = __t75;
			_t = _t->getNextSibling();
			break;
		}
		case OCTALINT:
		case DECIMALINT:
		case HEXADECIMALINT:
		case FLOATONE:
		case FLOATTWO:
		{
			{
			if (_t == ANTLR_USE_NAMESPACE(antlr)nullAST )
				_t = ASTNULL;
			switch ( _t->getType()) {
			case OCTALINT:
			case DECIMALINT:
			case HEXADECIMALINT:
			{
				ic = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				intConst(_t);
				_t = _retTree;
#line 260 "MDTreeParser.g"
				dval = lexi_cast<double>(ic->getText());
#line 1429 "MDTreeParser.cpp"
				break;
			}
			case FLOATONE:
			case FLOATTWO:
			{
				fc = (_t == ASTNULL) ? ANTLR_USE_NAMESPACE(antlr)nullAST : _t;
				floatConst(_t);
				_t = _retTree;
#line 261 "MDTreeParser.g"
				dval = lexi_cast<double>(fc->getText());
#line 1440 "MDTreeParser.cpp"
				break;
			}
			default:
			{
				throw ANTLR_USE_NAMESPACE(antlr)NoViableAltException(_t);
			}
			}
			}
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
	"\"atom\"",
	"\"bond\"",
	"\"bend\"",
	"\"torsion\"",
	"\"rigidBody\"",
	"\"cutoffGroup\"",
	"\"fragment\"",
	"\"members\"",
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
	"OCTALINT",
	"DECIMALINT",
	"HEXADECIMALINT",
	"PLUS",
	"MINUS",
	"FLOATONE",
	"FLOATTWO",
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
	"Digit",
	"Decimal",
	"LongSuffix",
	"UnsignedSuffix",
	"FloatSuffix",
	"Exponent",
	"Vocabulary",
	"Number",
	0
};

const unsigned long MDTreeParser::_tokenSet_0_data_[] = { 524400UL, 0UL, 0UL, 0UL };
// "component" "molecule" "zconstraint" ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_0(_tokenSet_0_data_,4);
const unsigned long MDTreeParser::_tokenSet_1_data_[] = { 540544UL, 0UL, 0UL, 0UL };
// "atom" "bond" "bend" "torsion" "rigidBody" "cutoffGroup" "fragment" 
// ASSIGNEQUAL 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_1(_tokenSet_1_data_,4);
const unsigned long MDTreeParser::_tokenSet_2_data_[] = { 3758096384UL, 14UL, 0UL, 0UL };
// OCTALINT DECIMALINT HEXADECIMALINT MINUS FLOATONE FLOATTWO 
const ANTLR_USE_NAMESPACE(antlr)BitSet MDTreeParser::_tokenSet_2(_tokenSet_2_data_,4);


