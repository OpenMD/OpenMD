#include "utils/StringTokenizerTestCase.hpp"
#include <iostream>
#include <algorithm>
#include "utils/StringTokenizer.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( StringTokenizerTestCase );


void StringTokenizerTestCase::testStringTokenizer(){

    
    std::string str1 = "    \t  Hello \r World \n";
    StringTokenizer tokenizer1(str1);

    CPPUNIT_ASSERT(tokenizer1.getOriginal() == str1);

    CPPUNIT_ASSERT(tokenizer1.countTokens() == 2);

    std::string token1 = tokenizer1.nextToken();
    std::string token2 = tokenizer1.nextToken();

    CPPUNIT_ASSERT(token1 == "Hello");
    CPPUNIT_ASSERT(token2 == "World");

    //test reading and converting tokens to other data type
    std::string str2 = "1991.2\t129\t1e2 1 OOPSE\n";
    StringTokenizer tokenizer2(str2);

    CPPUNIT_ASSERT(tokenizer2.countTokens() == 5);

    float floatVal = tokenizer2.nextTokenAsFloat();
    int intVal = tokenizer2.nextTokenAsInt();
    double doubleVal = tokenizer2.nextTokenAsDouble();
    bool boolVal = tokenizer2.nextTokenAsBool();
    std::string stringVal = tokenizer2.nextToken();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(floatVal, 1991.2, 0.0001);
    CPPUNIT_ASSERT(intVal == 129);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(doubleVal, 100, 0.0001);
    CPPUNIT_ASSERT(boolVal);
    CPPUNIT_ASSERT(stringVal == "OOPSE");

    CPPUNIT_ASSERT(!tokenizer2.hasMoreTokens());

    //test peekNextToken and using different delimeters
    StringTokenizer tokenizer3(str2, " \n");

    CPPUNIT_ASSERT(tokenizer3.getDelimiters() == " \n");
    
    CPPUNIT_ASSERT(tokenizer3.countTokens() == 3);

    CPPUNIT_ASSERT(tokenizer3.peekNextToken() == "1991.2\t129\t1e2");

    CPPUNIT_ASSERT(tokenizer3.countTokens() == 3);

    CPPUNIT_ASSERT(tokenizer3.nextToken() == "1991.2\t129\t1e2");

    CPPUNIT_ASSERT(tokenizer3.countTokens() == 2);

    //test return tokens
    StringTokenizer tokenizer4(str2, " \n", true);
    CPPUNIT_ASSERT(tokenizer4.countTokens() == 6);
    CPPUNIT_ASSERT(tokenizer4.nextToken() == "1991.2\t129\t1e2");
    CPPUNIT_ASSERT(tokenizer4.nextToken() == " ");
    CPPUNIT_ASSERT(tokenizer4.nextToken() == "1");
    CPPUNIT_ASSERT(tokenizer4.nextToken() == " ");
    CPPUNIT_ASSERT(tokenizer4.nextToken() == "OOPSE");
    CPPUNIT_ASSERT(tokenizer4.nextToken() == "\n");
        
}

