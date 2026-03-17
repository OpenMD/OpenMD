
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"




class  OMDLexer : public antlr4::Lexer {
public:
  enum {
    TRUE = 1, FALSE = 2, COMPONENT = 3, MOLECULE = 4, ZCONSTRAINT = 5, RESTRAINT = 6, 
    ATOM = 7, BOND = 8, BEND = 9, TORSION = 10, INVERSION = 11, RIGIDBODY = 12, 
    CUTOFFGROUP = 13, CONSTRAINT = 14, DISTANCE = 15, FRAGMENT = 16, SEQUENCE = 17, 
    MEMBERS = 18, CENTER = 19, SATELLITES = 20, POSITION = 21, ORIENTATION = 22, 
    FLUCQ = 23, RNEMD = 24, LIGHT = 25, MINIMIZER = 26, FIXED = 27, HARMONIC = 28, 
    CUBIC = 29, QUARTIC = 30, POLYNOMIAL = 31, MORSE = 32, GHOSTBEND = 33, 
    UREYBRADLEY = 34, COSINE = 35, GHOSTTORSION = 36, CHARMM = 37, OPLS = 38, 
    TRAPPE = 39, AMBERIMPROPER = 40, IMPROPERCOSINE = 41, CENTRALATOMHEIGHT = 42, 
    DREIDING = 43, CHARGE = 44, NODES = 45, ASSIGNEQUAL = 46, COLON = 47, 
    COMMA = 48, QUESTIONMARK = 49, SEMICOLON = 50, DOT = 51, LPAREN = 52, 
    RPAREN = 53, LBRACKET = 54, RBRACKET = 55, LCURLY = 56, RCURLY = 57, 
    NUM_LONG = 58, NUM_INT = 59, NUM_FLOAT = 60, NUM_DOUBLE = 61, CharLiteral = 62, 
    StringLiteral = 63, ID = 64, Whitespace = 65, Newline = 66, LineContinuation = 67, 
    Comment = 68, CPPComment = 69, PREPROC_DIRECTIVE = 70
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

