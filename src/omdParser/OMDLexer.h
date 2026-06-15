
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
    FLUCQ = 23, RNEMD = 24, LIGHT = 25, VELOCITYFIELD = 26, MINIMIZER = 27, 
    FIXED = 28, HARMONIC = 29, CUBIC = 30, QUARTIC = 31, POLYNOMIAL = 32, 
    MORSE = 33, GHOSTBEND = 34, UREYBRADLEY = 35, COSINE = 36, GHOSTTORSION = 37, 
    CHARMM = 38, OPLS = 39, TRAPPE = 40, AMBERIMPROPER = 41, IMPROPERCOSINE = 42, 
    CENTRALATOMHEIGHT = 43, DREIDING = 44, CHARGE = 45, NODES = 46, ASSIGNEQUAL = 47, 
    COLON = 48, COMMA = 49, QUESTIONMARK = 50, SEMICOLON = 51, DOT = 52, 
    LPAREN = 53, RPAREN = 54, LBRACKET = 55, RBRACKET = 56, LCURLY = 57, 
    RCURLY = 58, NUM_LONG = 59, NUM_INT = 60, NUM_FLOAT = 61, NUM_DOUBLE = 62, 
    CharLiteral = 63, StringLiteral = 64, ID = 65, Whitespace = 66, Newline = 67, 
    LineContinuation = 68, Comment = 69, CPPComment = 70, PREPROC_DIRECTIVE = 71
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

