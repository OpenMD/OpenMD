
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"




class  OMDLexer : public antlr4::Lexer {
public:
  enum {
    COMPONENT = 1, MOLECULE = 2, ZCONSTRAINT = 3, RESTRAINT = 4, ATOM = 5, 
    BOND = 6, BEND = 7, TORSION = 8, INVERSION = 9, RIGIDBODY = 10, CUTOFFGROUP = 11, 
    CONSTRAINT = 12, DISTANCE = 13, FRAGMENT = 14, SEQUENCE = 15, MEMBERS = 16, 
    CENTER = 17, SATELLITES = 18, POSITION = 19, ORIENTATION = 20, FLUCQ = 21, 
    RNEMD = 22, LIGHT = 23, MINIMIZER = 24, FIXED = 25, HARMONIC = 26, CUBIC = 27, 
    QUARTIC = 28, POLYNOMIAL = 29, MORSE = 30, GHOSTBEND = 31, UREYBRADLEY = 32, 
    COSINE = 33, GHOSTTORSION = 34, CHARMM = 35, OPLS = 36, TRAPPE = 37, 
    AMBERIMPROPER = 38, IMPROPERCOSINE = 39, CENTRALATOMHEIGHT = 40, DREIDING = 41, 
    CHARGE = 42, NODES = 43, ASSIGNEQUAL = 44, COLON = 45, COMMA = 46, QUESTIONMARK = 47, 
    SEMICOLON = 48, DOT = 49, LPAREN = 50, RPAREN = 51, LBRACKET = 52, RBRACKET = 53, 
    LCURLY = 54, RCURLY = 55, ID = 56, NUM_LONG = 57, NUM_INT = 58, NUM_FLOAT = 59, 
    NUM_DOUBLE = 60, CharLiteral = 61, StringLiteral = 62, Whitespace = 63, 
    Newline = 64, LineContinuation = 65, Comment = 66, CPPComment = 67, 
    PREPROC_DIRECTIVE = 68
  };

  explicit OMDLexer(antlr4::CharStream *input);

  ~OMDLexer() override;


      int deferredLineCount = 0;
      FilenameObserver* observer;
      
      void setObserver(FilenameObserver* osv) {
          observer = osv;
      }
      
      void initDeferredLineCount() {
          deferredLineCount = 0;
      }
      
      void deferredNewline() {
          deferredLineCount++;
      }


  std::string getGrammarFileName() const override;

  const std::vector<std::string>& getRuleNames() const override;

  const std::vector<std::string>& getChannelNames() const override;

  const std::vector<std::string>& getModeNames() const override;

  const antlr4::dfa::Vocabulary& getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;

  const antlr4::atn::ATN& getATN() const override;

  // By default the static state used to implement the lexer is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

private:

  // Individual action functions triggered by action() above.

  // Individual semantic predicate functions triggered by sempred() above.

};

