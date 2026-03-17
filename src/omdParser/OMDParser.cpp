
#include "antlr4-runtime.h"
#include "utils/StringUtils.hpp"
#include "omdParser/FilenameObserver.hpp"


// Generated from OMD.g4 by ANTLR 4.13.2


#include "OMDListener.h"
#include "OMDVisitor.h"

#include "OMDParser.h"


using namespace antlrcpp;

using namespace antlr4;

namespace {

struct OMDParserStaticData final {
  OMDParserStaticData(std::vector<std::string> ruleNames,
                        std::vector<std::string> literalNames,
                        std::vector<std::string> symbolicNames)
      : ruleNames(std::move(ruleNames)), literalNames(std::move(literalNames)),
        symbolicNames(std::move(symbolicNames)),
        vocabulary(this->literalNames, this->symbolicNames) {}

  OMDParserStaticData(const OMDParserStaticData&) = delete;
  OMDParserStaticData(OMDParserStaticData&&) = delete;
  OMDParserStaticData& operator=(const OMDParserStaticData&) = delete;
  OMDParserStaticData& operator=(OMDParserStaticData&&) = delete;

  std::vector<antlr4::dfa::DFA> decisionToDFA;
  antlr4::atn::PredictionContextCache sharedContextCache;
  const std::vector<std::string> ruleNames;
  const std::vector<std::string> literalNames;
  const std::vector<std::string> symbolicNames;
  const antlr4::dfa::Vocabulary vocabulary;
  antlr4::atn::SerializedATNView serializedATN;
  std::unique_ptr<antlr4::atn::ATN> atn;
};

::antlr4::internal::OnceFlag omdParserOnceFlag;
#if ANTLR4_USE_THREAD_LOCAL_CACHE
static thread_local
#endif
std::unique_ptr<OMDParserStaticData> omdParserStaticData = nullptr;

void omdParserInitialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  if (omdParserStaticData != nullptr) {
    return;
  }
#else
  assert(omdParserStaticData == nullptr);
#endif
  auto staticData = std::make_unique<OMDParserStaticData>(
    std::vector<std::string>{
      "omdfile", "statement", "assignment", "constant", "componentblock", 
      "zconstraintblock", "restraintblock", "flucqblock", "rnemdblock", 
      "lightblock", "minimizerblock", "moleculeblock", "moleculestatement", 
      "atomblock", "atomstatement", "bondblock", "bondstatement", "bendblock", 
      "bendstatement", "torsionblock", "torsionstatement", "inversionblock", 
      "inversionstatement", "rigidbodyblock", "rigidbodystatement", "cutoffgroupblock", 
      "cutoffgroupstatement", "nodesblock", "nodesstatement", "fragmentblock", 
      "fragmentstatement", "constraintblock", "constraintstatement", "sequencestring", 
      "doubleNumberTuple", "inttuple", "intConst", "doubleNumber", "floatConst", 
      "vectorConst"
    },
    std::vector<std::string>{
      "", "'true'", "'false'", "'component'", "'molecule'", "'zconstraint'", 
      "'restraint'", "'atom'", "'bond'", "'bend'", "'torsion'", "'inversion'", 
      "'rigidBody'", "'cutoffGroup'", "'constraint'", "'distance'", "'fragment'", 
      "'sequence'", "'members'", "'center'", "'satellites'", "'position'", 
      "'orientation'", "'flucQ'", "'RNEMD'", "'light'", "'minimizer'", "'Fixed'", 
      "'Harmonic'", "'Cubic'", "'Quartic'", "'Polynomial'", "'Morse'", "'GhostBend'", 
      "'UreyBradley'", "'Cosine'", "'GhostTorsion'", "'Charmm'", "'Opls'", 
      "'Trappe'", "'AmberImproper'", "'ImproperCosine'", "'CentralAtomHeight'", 
      "'Dreiding'", "'charge'", "'nodes'", "'='", "':'", "','", "'\\u003F'", 
      "';'", "'.'", "'('", "')'", "'['", "']'", "'{'", "'}'"
    },
    std::vector<std::string>{
      "", "TRUE", "FALSE", "COMPONENT", "MOLECULE", "ZCONSTRAINT", "RESTRAINT", 
      "ATOM", "BOND", "BEND", "TORSION", "INVERSION", "RIGIDBODY", "CUTOFFGROUP", 
      "CONSTRAINT", "DISTANCE", "FRAGMENT", "SEQUENCE", "MEMBERS", "CENTER", 
      "SATELLITES", "POSITION", "ORIENTATION", "FLUCQ", "RNEMD", "LIGHT", 
      "MINIMIZER", "FIXED", "HARMONIC", "CUBIC", "QUARTIC", "POLYNOMIAL", 
      "MORSE", "GHOSTBEND", "UREYBRADLEY", "COSINE", "GHOSTTORSION", "CHARMM", 
      "OPLS", "TRAPPE", "AMBERIMPROPER", "IMPROPERCOSINE", "CENTRALATOMHEIGHT", 
      "DREIDING", "CHARGE", "NODES", "ASSIGNEQUAL", "COLON", "COMMA", "QUESTIONMARK", 
      "SEMICOLON", "DOT", "LPAREN", "RPAREN", "LBRACKET", "RBRACKET", "LCURLY", 
      "RCURLY", "NUM_LONG", "NUM_INT", "NUM_FLOAT", "NUM_DOUBLE", "CharLiteral", 
      "StringLiteral", "ID", "Whitespace", "Newline", "LineContinuation", 
      "Comment", "CPPComment", "PREPROC_DIRECTIVE"
    }
  );
  static const int32_t serializedATNSegment[] = {
  	4,1,70,655,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,6,2,
  	7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,13,7,13,2,14,7,
  	14,2,15,7,15,2,16,7,16,2,17,7,17,2,18,7,18,2,19,7,19,2,20,7,20,2,21,7,
  	21,2,22,7,22,2,23,7,23,2,24,7,24,2,25,7,25,2,26,7,26,2,27,7,27,2,28,7,
  	28,2,29,7,29,2,30,7,30,2,31,7,31,2,32,7,32,2,33,7,33,2,34,7,34,2,35,7,
  	35,2,36,7,36,2,37,7,37,2,38,7,38,2,39,7,39,1,0,5,0,82,8,0,10,0,12,0,85,
  	9,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,99,8,1,1,2,1,
  	2,1,2,1,2,1,2,1,3,1,3,1,3,1,3,1,3,1,3,1,3,3,3,113,8,3,1,4,1,4,1,4,5,4,
  	118,8,4,10,4,12,4,121,9,4,1,4,1,4,1,5,1,5,1,5,5,5,128,8,5,10,5,12,5,131,
  	9,5,1,5,1,5,1,6,1,6,1,6,5,6,138,8,6,10,6,12,6,141,9,6,1,6,1,6,1,7,1,7,
  	1,7,5,7,148,8,7,10,7,12,7,151,9,7,1,7,1,7,1,8,1,8,1,8,5,8,158,8,8,10,
  	8,12,8,161,9,8,1,8,1,8,1,9,1,9,1,9,5,9,168,8,9,10,9,12,9,171,9,9,1,9,
  	1,9,1,10,1,10,1,10,5,10,178,8,10,10,10,12,10,181,9,10,1,10,1,10,1,11,
  	1,11,1,11,5,11,188,8,11,10,11,12,11,191,9,11,1,11,1,11,1,12,1,12,1,12,
  	1,12,1,12,1,12,1,12,1,12,1,12,1,12,3,12,205,8,12,1,13,1,13,1,13,1,13,
  	1,13,1,13,5,13,213,8,13,10,13,12,13,216,9,13,1,13,1,13,1,14,1,14,1,14,
  	1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,
  	1,14,1,14,3,14,239,8,14,1,15,1,15,1,15,1,15,1,15,3,15,246,8,15,1,15,1,
  	15,5,15,250,8,15,10,15,12,15,253,9,15,1,15,1,15,1,16,1,16,1,16,1,16,1,
  	16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,
  	16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,
  	16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,3,16,300,8,16,1,
  	17,1,17,1,17,1,17,1,17,3,17,307,8,17,1,17,1,17,5,17,311,8,17,10,17,12,
  	17,314,9,17,1,17,1,17,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,
  	18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,
  	18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,
  	18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,3,18,367,8,
  	18,1,19,1,19,1,19,1,19,1,19,3,19,374,8,19,1,19,1,19,5,19,378,8,19,10,
  	19,12,19,381,9,19,1,19,1,19,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,
  	20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,
  	20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,
  	20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,
  	20,1,20,1,20,1,20,1,20,3,20,440,8,20,1,21,1,21,1,21,1,21,1,21,3,21,447,
  	8,21,1,21,1,21,5,21,451,8,21,10,21,12,21,454,9,21,1,21,1,21,1,22,1,22,
  	1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,
  	1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,
  	1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,3,22,
  	501,8,22,1,23,1,23,1,23,1,23,1,23,1,23,5,23,509,8,23,10,23,12,23,512,
  	9,23,1,23,1,23,1,24,1,24,1,24,1,24,1,24,1,24,1,24,3,24,523,8,24,1,25,
  	1,25,1,25,1,25,1,25,3,25,530,8,25,1,25,1,25,5,25,534,8,25,10,25,12,25,
  	537,9,25,1,25,1,25,1,26,1,26,1,26,1,26,1,26,1,26,1,26,3,26,548,8,26,1,
  	27,1,27,1,27,1,27,1,27,3,27,555,8,27,1,27,1,27,5,27,559,8,27,10,27,12,
  	27,562,9,27,1,27,1,27,1,28,1,28,1,28,1,28,1,28,1,28,1,28,3,28,573,8,28,
  	1,29,1,29,1,29,5,29,578,8,29,10,29,12,29,581,9,29,1,29,1,29,1,30,1,30,
  	1,30,1,30,1,30,1,30,1,30,1,30,1,30,1,30,3,30,595,8,30,1,31,1,31,1,31,
  	1,31,1,31,3,31,602,8,31,1,31,1,31,5,31,606,8,31,10,31,12,31,609,9,31,
  	1,31,1,31,1,32,1,32,1,32,1,32,1,32,1,32,1,32,3,32,620,8,32,1,33,1,33,
  	1,33,1,33,1,33,1,34,1,34,1,34,5,34,630,8,34,10,34,12,34,633,9,34,1,35,
  	1,35,1,35,5,35,638,8,35,10,35,12,35,641,9,35,1,36,1,36,1,37,1,37,3,37,
  	647,8,37,1,38,1,38,1,39,1,39,1,39,1,39,1,39,0,0,40,0,2,4,6,8,10,12,14,
  	16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,
  	62,64,66,68,70,72,74,76,78,0,2,1,0,58,59,1,0,60,61,714,0,83,1,0,0,0,2,
  	98,1,0,0,0,4,100,1,0,0,0,6,112,1,0,0,0,8,114,1,0,0,0,10,124,1,0,0,0,12,
  	134,1,0,0,0,14,144,1,0,0,0,16,154,1,0,0,0,18,164,1,0,0,0,20,174,1,0,0,
  	0,22,184,1,0,0,0,24,204,1,0,0,0,26,206,1,0,0,0,28,238,1,0,0,0,30,240,
  	1,0,0,0,32,299,1,0,0,0,34,301,1,0,0,0,36,366,1,0,0,0,38,368,1,0,0,0,40,
  	439,1,0,0,0,42,441,1,0,0,0,44,500,1,0,0,0,46,502,1,0,0,0,48,522,1,0,0,
  	0,50,524,1,0,0,0,52,547,1,0,0,0,54,549,1,0,0,0,56,572,1,0,0,0,58,574,
  	1,0,0,0,60,594,1,0,0,0,62,596,1,0,0,0,64,619,1,0,0,0,66,621,1,0,0,0,68,
  	626,1,0,0,0,70,634,1,0,0,0,72,642,1,0,0,0,74,646,1,0,0,0,76,648,1,0,0,
  	0,78,650,1,0,0,0,80,82,3,2,1,0,81,80,1,0,0,0,82,85,1,0,0,0,83,81,1,0,
  	0,0,83,84,1,0,0,0,84,86,1,0,0,0,85,83,1,0,0,0,86,87,5,0,0,1,87,1,1,0,
  	0,0,88,99,3,4,2,0,89,99,3,8,4,0,90,99,3,22,11,0,91,99,3,58,29,0,92,99,
  	3,10,5,0,93,99,3,12,6,0,94,99,3,14,7,0,95,99,3,16,8,0,96,99,3,18,9,0,
  	97,99,3,20,10,0,98,88,1,0,0,0,98,89,1,0,0,0,98,90,1,0,0,0,98,91,1,0,0,
  	0,98,92,1,0,0,0,98,93,1,0,0,0,98,94,1,0,0,0,98,95,1,0,0,0,98,96,1,0,0,
  	0,98,97,1,0,0,0,99,3,1,0,0,0,100,101,5,64,0,0,101,102,5,46,0,0,102,103,
  	3,6,3,0,103,104,5,50,0,0,104,5,1,0,0,0,105,113,3,72,36,0,106,113,3,76,
  	38,0,107,113,3,78,39,0,108,113,5,1,0,0,109,113,5,2,0,0,110,113,5,64,0,
  	0,111,113,5,63,0,0,112,105,1,0,0,0,112,106,1,0,0,0,112,107,1,0,0,0,112,
  	108,1,0,0,0,112,109,1,0,0,0,112,110,1,0,0,0,112,111,1,0,0,0,113,7,1,0,
  	0,0,114,115,5,3,0,0,115,119,5,56,0,0,116,118,3,4,2,0,117,116,1,0,0,0,
  	118,121,1,0,0,0,119,117,1,0,0,0,119,120,1,0,0,0,120,122,1,0,0,0,121,119,
  	1,0,0,0,122,123,5,57,0,0,123,9,1,0,0,0,124,125,5,5,0,0,125,129,5,56,0,
  	0,126,128,3,4,2,0,127,126,1,0,0,0,128,131,1,0,0,0,129,127,1,0,0,0,129,
  	130,1,0,0,0,130,132,1,0,0,0,131,129,1,0,0,0,132,133,5,57,0,0,133,11,1,
  	0,0,0,134,135,5,6,0,0,135,139,5,56,0,0,136,138,3,4,2,0,137,136,1,0,0,
  	0,138,141,1,0,0,0,139,137,1,0,0,0,139,140,1,0,0,0,140,142,1,0,0,0,141,
  	139,1,0,0,0,142,143,5,57,0,0,143,13,1,0,0,0,144,145,5,23,0,0,145,149,
  	5,56,0,0,146,148,3,4,2,0,147,146,1,0,0,0,148,151,1,0,0,0,149,147,1,0,
  	0,0,149,150,1,0,0,0,150,152,1,0,0,0,151,149,1,0,0,0,152,153,5,57,0,0,
  	153,15,1,0,0,0,154,155,5,24,0,0,155,159,5,56,0,0,156,158,3,4,2,0,157,
  	156,1,0,0,0,158,161,1,0,0,0,159,157,1,0,0,0,159,160,1,0,0,0,160,162,1,
  	0,0,0,161,159,1,0,0,0,162,163,5,57,0,0,163,17,1,0,0,0,164,165,5,25,0,
  	0,165,169,5,56,0,0,166,168,3,4,2,0,167,166,1,0,0,0,168,171,1,0,0,0,169,
  	167,1,0,0,0,169,170,1,0,0,0,170,172,1,0,0,0,171,169,1,0,0,0,172,173,5,
  	57,0,0,173,19,1,0,0,0,174,175,5,26,0,0,175,179,5,56,0,0,176,178,3,4,2,
  	0,177,176,1,0,0,0,178,181,1,0,0,0,179,177,1,0,0,0,179,180,1,0,0,0,180,
  	182,1,0,0,0,181,179,1,0,0,0,182,183,5,57,0,0,183,21,1,0,0,0,184,185,5,
  	4,0,0,185,189,5,56,0,0,186,188,3,24,12,0,187,186,1,0,0,0,188,191,1,0,
  	0,0,189,187,1,0,0,0,189,190,1,0,0,0,190,192,1,0,0,0,191,189,1,0,0,0,192,
  	193,5,57,0,0,193,23,1,0,0,0,194,205,3,4,2,0,195,205,3,26,13,0,196,205,
  	3,30,15,0,197,205,3,34,17,0,198,205,3,38,19,0,199,205,3,42,21,0,200,205,
  	3,46,23,0,201,205,3,50,25,0,202,205,3,62,31,0,203,205,3,66,33,0,204,194,
  	1,0,0,0,204,195,1,0,0,0,204,196,1,0,0,0,204,197,1,0,0,0,204,198,1,0,0,
  	0,204,199,1,0,0,0,204,200,1,0,0,0,204,201,1,0,0,0,204,202,1,0,0,0,204,
  	203,1,0,0,0,205,25,1,0,0,0,206,207,5,7,0,0,207,208,5,54,0,0,208,209,3,
  	72,36,0,209,210,5,55,0,0,210,214,5,56,0,0,211,213,3,28,14,0,212,211,1,
  	0,0,0,213,216,1,0,0,0,214,212,1,0,0,0,214,215,1,0,0,0,215,217,1,0,0,0,
  	216,214,1,0,0,0,217,218,5,57,0,0,218,27,1,0,0,0,219,239,3,4,2,0,220,221,
  	5,21,0,0,221,222,5,52,0,0,222,223,3,68,34,0,223,224,5,53,0,0,224,225,
  	5,50,0,0,225,239,1,0,0,0,226,227,5,22,0,0,227,228,5,52,0,0,228,229,3,
  	68,34,0,229,230,5,53,0,0,230,231,5,50,0,0,231,239,1,0,0,0,232,233,5,44,
  	0,0,233,234,5,52,0,0,234,235,3,76,38,0,235,236,5,53,0,0,236,237,5,50,
  	0,0,237,239,1,0,0,0,238,219,1,0,0,0,238,220,1,0,0,0,238,226,1,0,0,0,238,
  	232,1,0,0,0,239,29,1,0,0,0,240,245,5,8,0,0,241,242,5,54,0,0,242,243,3,
  	72,36,0,243,244,5,55,0,0,244,246,1,0,0,0,245,241,1,0,0,0,245,246,1,0,
  	0,0,246,247,1,0,0,0,247,251,5,56,0,0,248,250,3,32,16,0,249,248,1,0,0,
  	0,250,253,1,0,0,0,251,249,1,0,0,0,251,252,1,0,0,0,252,254,1,0,0,0,253,
  	251,1,0,0,0,254,255,5,57,0,0,255,31,1,0,0,0,256,300,3,4,2,0,257,258,5,
  	18,0,0,258,259,5,52,0,0,259,260,3,70,35,0,260,261,5,53,0,0,261,262,5,
  	50,0,0,262,300,1,0,0,0,263,264,5,27,0,0,264,265,5,52,0,0,265,266,3,76,
  	38,0,266,267,5,53,0,0,267,268,5,50,0,0,268,300,1,0,0,0,269,270,5,28,0,
  	0,270,271,5,52,0,0,271,272,3,68,34,0,272,273,5,53,0,0,273,274,5,50,0,
  	0,274,300,1,0,0,0,275,276,5,29,0,0,276,277,5,52,0,0,277,278,3,68,34,0,
  	278,279,5,53,0,0,279,280,5,50,0,0,280,300,1,0,0,0,281,282,5,30,0,0,282,
  	283,5,52,0,0,283,284,3,68,34,0,284,285,5,53,0,0,285,286,5,50,0,0,286,
  	300,1,0,0,0,287,288,5,31,0,0,288,289,5,52,0,0,289,290,3,68,34,0,290,291,
  	5,53,0,0,291,292,5,50,0,0,292,300,1,0,0,0,293,294,5,32,0,0,294,295,5,
  	52,0,0,295,296,3,68,34,0,296,297,5,53,0,0,297,298,5,50,0,0,298,300,1,
  	0,0,0,299,256,1,0,0,0,299,257,1,0,0,0,299,263,1,0,0,0,299,269,1,0,0,0,
  	299,275,1,0,0,0,299,281,1,0,0,0,299,287,1,0,0,0,299,293,1,0,0,0,300,33,
  	1,0,0,0,301,306,5,9,0,0,302,303,5,54,0,0,303,304,3,72,36,0,304,305,5,
  	55,0,0,305,307,1,0,0,0,306,302,1,0,0,0,306,307,1,0,0,0,307,308,1,0,0,
  	0,308,312,5,56,0,0,309,311,3,36,18,0,310,309,1,0,0,0,311,314,1,0,0,0,
  	312,310,1,0,0,0,312,313,1,0,0,0,313,315,1,0,0,0,314,312,1,0,0,0,315,316,
  	5,57,0,0,316,35,1,0,0,0,317,367,3,4,2,0,318,319,5,18,0,0,319,320,5,52,
  	0,0,320,321,3,70,35,0,321,322,5,53,0,0,322,323,5,50,0,0,323,367,1,0,0,
  	0,324,325,5,28,0,0,325,326,5,52,0,0,326,327,3,68,34,0,327,328,5,53,0,
  	0,328,329,5,50,0,0,329,367,1,0,0,0,330,331,5,33,0,0,331,332,5,52,0,0,
  	332,333,3,68,34,0,333,334,5,53,0,0,334,335,5,50,0,0,335,367,1,0,0,0,336,
  	337,5,34,0,0,337,338,5,52,0,0,338,339,3,68,34,0,339,340,5,53,0,0,340,
  	341,5,50,0,0,341,367,1,0,0,0,342,343,5,29,0,0,343,344,5,52,0,0,344,345,
  	3,68,34,0,345,346,5,53,0,0,346,347,5,50,0,0,347,367,1,0,0,0,348,349,5,
  	30,0,0,349,350,5,52,0,0,350,351,3,68,34,0,351,352,5,53,0,0,352,353,5,
  	50,0,0,353,367,1,0,0,0,354,355,5,31,0,0,355,356,5,52,0,0,356,357,3,68,
  	34,0,357,358,5,53,0,0,358,359,5,50,0,0,359,367,1,0,0,0,360,361,5,35,0,
  	0,361,362,5,52,0,0,362,363,3,68,34,0,363,364,5,53,0,0,364,365,5,50,0,
  	0,365,367,1,0,0,0,366,317,1,0,0,0,366,318,1,0,0,0,366,324,1,0,0,0,366,
  	330,1,0,0,0,366,336,1,0,0,0,366,342,1,0,0,0,366,348,1,0,0,0,366,354,1,
  	0,0,0,366,360,1,0,0,0,367,37,1,0,0,0,368,373,5,10,0,0,369,370,5,54,0,
  	0,370,371,3,72,36,0,371,372,5,55,0,0,372,374,1,0,0,0,373,369,1,0,0,0,
  	373,374,1,0,0,0,374,375,1,0,0,0,375,379,5,56,0,0,376,378,3,40,20,0,377,
  	376,1,0,0,0,378,381,1,0,0,0,379,377,1,0,0,0,379,380,1,0,0,0,380,382,1,
  	0,0,0,381,379,1,0,0,0,382,383,5,57,0,0,383,39,1,0,0,0,384,440,3,4,2,0,
  	385,386,5,18,0,0,386,387,5,52,0,0,387,388,3,70,35,0,388,389,5,53,0,0,
  	389,390,5,50,0,0,390,440,1,0,0,0,391,392,5,36,0,0,392,393,5,52,0,0,393,
  	394,3,68,34,0,394,395,5,53,0,0,395,396,5,50,0,0,396,440,1,0,0,0,397,398,
  	5,29,0,0,398,399,5,52,0,0,399,400,3,68,34,0,400,401,5,53,0,0,401,402,
  	5,50,0,0,402,440,1,0,0,0,403,404,5,30,0,0,404,405,5,52,0,0,405,406,3,
  	68,34,0,406,407,5,53,0,0,407,408,5,50,0,0,408,440,1,0,0,0,409,410,5,31,
  	0,0,410,411,5,52,0,0,411,412,3,68,34,0,412,413,5,53,0,0,413,414,5,50,
  	0,0,414,440,1,0,0,0,415,416,5,37,0,0,416,417,5,52,0,0,417,418,3,68,34,
  	0,418,419,5,53,0,0,419,420,5,50,0,0,420,440,1,0,0,0,421,422,5,38,0,0,
  	422,423,5,52,0,0,423,424,3,68,34,0,424,425,5,53,0,0,425,426,5,50,0,0,
  	426,440,1,0,0,0,427,428,5,39,0,0,428,429,5,52,0,0,429,430,3,68,34,0,430,
  	431,5,53,0,0,431,432,5,50,0,0,432,440,1,0,0,0,433,434,5,28,0,0,434,435,
  	5,52,0,0,435,436,3,68,34,0,436,437,5,53,0,0,437,438,5,50,0,0,438,440,
  	1,0,0,0,439,384,1,0,0,0,439,385,1,0,0,0,439,391,1,0,0,0,439,397,1,0,0,
  	0,439,403,1,0,0,0,439,409,1,0,0,0,439,415,1,0,0,0,439,421,1,0,0,0,439,
  	427,1,0,0,0,439,433,1,0,0,0,440,41,1,0,0,0,441,446,5,11,0,0,442,443,5,
  	54,0,0,443,444,3,72,36,0,444,445,5,55,0,0,445,447,1,0,0,0,446,442,1,0,
  	0,0,446,447,1,0,0,0,447,448,1,0,0,0,448,452,5,56,0,0,449,451,3,44,22,
  	0,450,449,1,0,0,0,451,454,1,0,0,0,452,450,1,0,0,0,452,453,1,0,0,0,453,
  	455,1,0,0,0,454,452,1,0,0,0,455,456,5,57,0,0,456,43,1,0,0,0,457,501,3,
  	4,2,0,458,459,5,19,0,0,459,460,5,52,0,0,460,461,3,72,36,0,461,462,5,53,
  	0,0,462,463,5,50,0,0,463,501,1,0,0,0,464,465,5,20,0,0,465,466,5,52,0,
  	0,466,467,3,70,35,0,467,468,5,53,0,0,468,469,5,50,0,0,469,501,1,0,0,0,
  	470,471,5,40,0,0,471,472,5,52,0,0,472,473,3,68,34,0,473,474,5,53,0,0,
  	474,475,5,50,0,0,475,501,1,0,0,0,476,477,5,41,0,0,477,478,5,52,0,0,478,
  	479,3,68,34,0,479,480,5,53,0,0,480,481,5,50,0,0,481,501,1,0,0,0,482,483,
  	5,28,0,0,483,484,5,52,0,0,484,485,3,68,34,0,485,486,5,53,0,0,486,487,
  	5,50,0,0,487,501,1,0,0,0,488,489,5,42,0,0,489,490,5,52,0,0,490,491,3,
  	68,34,0,491,492,5,53,0,0,492,493,5,50,0,0,493,501,1,0,0,0,494,495,5,43,
  	0,0,495,496,5,52,0,0,496,497,3,68,34,0,497,498,5,53,0,0,498,499,5,50,
  	0,0,499,501,1,0,0,0,500,457,1,0,0,0,500,458,1,0,0,0,500,464,1,0,0,0,500,
  	470,1,0,0,0,500,476,1,0,0,0,500,482,1,0,0,0,500,488,1,0,0,0,500,494,1,
  	0,0,0,501,45,1,0,0,0,502,503,5,12,0,0,503,504,5,54,0,0,504,505,3,72,36,
  	0,505,506,5,55,0,0,506,510,5,56,0,0,507,509,3,48,24,0,508,507,1,0,0,0,
  	509,512,1,0,0,0,510,508,1,0,0,0,510,511,1,0,0,0,511,513,1,0,0,0,512,510,
  	1,0,0,0,513,514,5,57,0,0,514,47,1,0,0,0,515,523,3,4,2,0,516,517,5,18,
  	0,0,517,518,5,52,0,0,518,519,3,70,35,0,519,520,5,53,0,0,520,521,5,50,
  	0,0,521,523,1,0,0,0,522,515,1,0,0,0,522,516,1,0,0,0,523,49,1,0,0,0,524,
  	529,5,13,0,0,525,526,5,54,0,0,526,527,3,72,36,0,527,528,5,55,0,0,528,
  	530,1,0,0,0,529,525,1,0,0,0,529,530,1,0,0,0,530,531,1,0,0,0,531,535,5,
  	56,0,0,532,534,3,52,26,0,533,532,1,0,0,0,534,537,1,0,0,0,535,533,1,0,
  	0,0,535,536,1,0,0,0,536,538,1,0,0,0,537,535,1,0,0,0,538,539,5,57,0,0,
  	539,51,1,0,0,0,540,548,3,4,2,0,541,542,5,18,0,0,542,543,5,52,0,0,543,
  	544,3,70,35,0,544,545,5,53,0,0,545,546,5,50,0,0,546,548,1,0,0,0,547,540,
  	1,0,0,0,547,541,1,0,0,0,548,53,1,0,0,0,549,554,5,45,0,0,550,551,5,54,
  	0,0,551,552,3,72,36,0,552,553,5,55,0,0,553,555,1,0,0,0,554,550,1,0,0,
  	0,554,555,1,0,0,0,555,556,1,0,0,0,556,560,5,56,0,0,557,559,3,56,28,0,
  	558,557,1,0,0,0,559,562,1,0,0,0,560,558,1,0,0,0,560,561,1,0,0,0,561,563,
  	1,0,0,0,562,560,1,0,0,0,563,564,5,57,0,0,564,55,1,0,0,0,565,573,3,4,2,
  	0,566,567,5,18,0,0,567,568,5,52,0,0,568,569,3,70,35,0,569,570,5,53,0,
  	0,570,571,5,50,0,0,571,573,1,0,0,0,572,565,1,0,0,0,572,566,1,0,0,0,573,
  	57,1,0,0,0,574,575,5,16,0,0,575,579,5,56,0,0,576,578,3,60,30,0,577,576,
  	1,0,0,0,578,581,1,0,0,0,579,577,1,0,0,0,579,580,1,0,0,0,580,582,1,0,0,
  	0,581,579,1,0,0,0,582,583,5,57,0,0,583,59,1,0,0,0,584,595,3,4,2,0,585,
  	595,3,26,13,0,586,595,3,30,15,0,587,595,3,34,17,0,588,595,3,38,19,0,589,
  	595,3,42,21,0,590,595,3,46,23,0,591,595,3,50,25,0,592,595,3,62,31,0,593,
  	595,3,54,27,0,594,584,1,0,0,0,594,585,1,0,0,0,594,586,1,0,0,0,594,587,
  	1,0,0,0,594,588,1,0,0,0,594,589,1,0,0,0,594,590,1,0,0,0,594,591,1,0,0,
  	0,594,592,1,0,0,0,594,593,1,0,0,0,595,61,1,0,0,0,596,601,5,14,0,0,597,
  	598,5,54,0,0,598,599,3,72,36,0,599,600,5,55,0,0,600,602,1,0,0,0,601,597,
  	1,0,0,0,601,602,1,0,0,0,602,603,1,0,0,0,603,607,5,56,0,0,604,606,3,64,
  	32,0,605,604,1,0,0,0,606,609,1,0,0,0,607,605,1,0,0,0,607,608,1,0,0,0,
  	608,610,1,0,0,0,609,607,1,0,0,0,610,611,5,57,0,0,611,63,1,0,0,0,612,620,
  	3,4,2,0,613,614,5,18,0,0,614,615,5,52,0,0,615,616,3,70,35,0,616,617,5,
  	53,0,0,617,618,5,50,0,0,618,620,1,0,0,0,619,612,1,0,0,0,619,613,1,0,0,
  	0,620,65,1,0,0,0,621,622,5,17,0,0,622,623,5,46,0,0,623,624,3,6,3,0,624,
  	625,5,50,0,0,625,67,1,0,0,0,626,631,3,74,37,0,627,628,5,48,0,0,628,630,
  	3,74,37,0,629,627,1,0,0,0,630,633,1,0,0,0,631,629,1,0,0,0,631,632,1,0,
  	0,0,632,69,1,0,0,0,633,631,1,0,0,0,634,639,3,72,36,0,635,636,5,48,0,0,
  	636,638,3,72,36,0,637,635,1,0,0,0,638,641,1,0,0,0,639,637,1,0,0,0,639,
  	640,1,0,0,0,640,71,1,0,0,0,641,639,1,0,0,0,642,643,7,0,0,0,643,73,1,0,
  	0,0,644,647,3,72,36,0,645,647,3,76,38,0,646,644,1,0,0,0,646,645,1,0,0,
  	0,647,75,1,0,0,0,648,649,7,1,0,0,649,77,1,0,0,0,650,651,5,52,0,0,651,
  	652,3,68,34,0,652,653,5,53,0,0,653,79,1,0,0,0,42,83,98,112,119,129,139,
  	149,159,169,179,189,204,214,238,245,251,299,306,312,366,373,379,439,446,
  	452,500,510,522,529,535,547,554,560,572,579,594,601,607,619,631,639,646
  };
  staticData->serializedATN = antlr4::atn::SerializedATNView(serializedATNSegment, sizeof(serializedATNSegment) / sizeof(serializedATNSegment[0]));

  antlr4::atn::ATNDeserializer deserializer;
  staticData->atn = deserializer.deserialize(staticData->serializedATN);

  const size_t count = staticData->atn->getNumberOfDecisions();
  staticData->decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) { 
    staticData->decisionToDFA.emplace_back(staticData->atn->getDecisionState(i), i);
  }
  omdParserStaticData = std::move(staticData);
}

}

OMDParser::OMDParser(TokenStream *input) : OMDParser(input, antlr4::atn::ParserATNSimulatorOptions()) {}

OMDParser::OMDParser(TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options) : Parser(input) {
  OMDParser::initialize();
  _interpreter = new atn::ParserATNSimulator(this, *omdParserStaticData->atn, omdParserStaticData->decisionToDFA, omdParserStaticData->sharedContextCache, options);
}

OMDParser::~OMDParser() {
  delete _interpreter;
}

const atn::ATN& OMDParser::getATN() const {
  return *omdParserStaticData->atn;
}

std::string OMDParser::getGrammarFileName() const {
  return "OMD.g4";
}

const std::vector<std::string>& OMDParser::getRuleNames() const {
  return omdParserStaticData->ruleNames;
}

const dfa::Vocabulary& OMDParser::getVocabulary() const {
  return omdParserStaticData->vocabulary;
}

antlr4::atn::SerializedATNView OMDParser::getSerializedATN() const {
  return omdParserStaticData->serializedATN;
}


//----------------- OmdfileContext ------------------------------------------------------------------

OMDParser::OmdfileContext::OmdfileContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::OmdfileContext::EOF() {
  return getToken(OMDParser::EOF, 0);
}

std::vector<OMDParser::StatementContext *> OMDParser::OmdfileContext::statement() {
  return getRuleContexts<OMDParser::StatementContext>();
}

OMDParser::StatementContext* OMDParser::OmdfileContext::statement(size_t i) {
  return getRuleContext<OMDParser::StatementContext>(i);
}


size_t OMDParser::OmdfileContext::getRuleIndex() const {
  return OMDParser::RuleOmdfile;
}

void OMDParser::OmdfileContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterOmdfile(this);
}

void OMDParser::OmdfileContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitOmdfile(this);
}


std::any OMDParser::OmdfileContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitOmdfile(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::OmdfileContext* OMDParser::omdfile() {
  OmdfileContext *_localctx = _tracker.createInstance<OmdfileContext>(_ctx, getState());
  enterRule(_localctx, 0, OMDParser::RuleOmdfile);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(83);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 3) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 3)) & 2305843009229430799) != 0)) {
      setState(80);
      statement();
      setState(85);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(86);
    match(OMDParser::EOF);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- StatementContext ------------------------------------------------------------------

OMDParser::StatementContext::StatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::StatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

OMDParser::ComponentblockContext* OMDParser::StatementContext::componentblock() {
  return getRuleContext<OMDParser::ComponentblockContext>(0);
}

OMDParser::MoleculeblockContext* OMDParser::StatementContext::moleculeblock() {
  return getRuleContext<OMDParser::MoleculeblockContext>(0);
}

OMDParser::FragmentblockContext* OMDParser::StatementContext::fragmentblock() {
  return getRuleContext<OMDParser::FragmentblockContext>(0);
}

OMDParser::ZconstraintblockContext* OMDParser::StatementContext::zconstraintblock() {
  return getRuleContext<OMDParser::ZconstraintblockContext>(0);
}

OMDParser::RestraintblockContext* OMDParser::StatementContext::restraintblock() {
  return getRuleContext<OMDParser::RestraintblockContext>(0);
}

OMDParser::FlucqblockContext* OMDParser::StatementContext::flucqblock() {
  return getRuleContext<OMDParser::FlucqblockContext>(0);
}

OMDParser::RnemdblockContext* OMDParser::StatementContext::rnemdblock() {
  return getRuleContext<OMDParser::RnemdblockContext>(0);
}

OMDParser::LightblockContext* OMDParser::StatementContext::lightblock() {
  return getRuleContext<OMDParser::LightblockContext>(0);
}

OMDParser::MinimizerblockContext* OMDParser::StatementContext::minimizerblock() {
  return getRuleContext<OMDParser::MinimizerblockContext>(0);
}


size_t OMDParser::StatementContext::getRuleIndex() const {
  return OMDParser::RuleStatement;
}

void OMDParser::StatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterStatement(this);
}

void OMDParser::StatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitStatement(this);
}


std::any OMDParser::StatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitStatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::StatementContext* OMDParser::statement() {
  StatementContext *_localctx = _tracker.createInstance<StatementContext>(_ctx, getState());
  enterRule(_localctx, 2, OMDParser::RuleStatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(98);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(88);
        assignment();
        break;
      }

      case OMDParser::COMPONENT: {
        enterOuterAlt(_localctx, 2);
        setState(89);
        componentblock();
        break;
      }

      case OMDParser::MOLECULE: {
        enterOuterAlt(_localctx, 3);
        setState(90);
        moleculeblock();
        break;
      }

      case OMDParser::FRAGMENT: {
        enterOuterAlt(_localctx, 4);
        setState(91);
        fragmentblock();
        break;
      }

      case OMDParser::ZCONSTRAINT: {
        enterOuterAlt(_localctx, 5);
        setState(92);
        zconstraintblock();
        break;
      }

      case OMDParser::RESTRAINT: {
        enterOuterAlt(_localctx, 6);
        setState(93);
        restraintblock();
        break;
      }

      case OMDParser::FLUCQ: {
        enterOuterAlt(_localctx, 7);
        setState(94);
        flucqblock();
        break;
      }

      case OMDParser::RNEMD: {
        enterOuterAlt(_localctx, 8);
        setState(95);
        rnemdblock();
        break;
      }

      case OMDParser::LIGHT: {
        enterOuterAlt(_localctx, 9);
        setState(96);
        lightblock();
        break;
      }

      case OMDParser::MINIMIZER: {
        enterOuterAlt(_localctx, 10);
        setState(97);
        minimizerblock();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- AssignmentContext ------------------------------------------------------------------

OMDParser::AssignmentContext::AssignmentContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::AssignmentContext::ID() {
  return getToken(OMDParser::ID, 0);
}

tree::TerminalNode* OMDParser::AssignmentContext::ASSIGNEQUAL() {
  return getToken(OMDParser::ASSIGNEQUAL, 0);
}

OMDParser::ConstantContext* OMDParser::AssignmentContext::constant() {
  return getRuleContext<OMDParser::ConstantContext>(0);
}

tree::TerminalNode* OMDParser::AssignmentContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}


size_t OMDParser::AssignmentContext::getRuleIndex() const {
  return OMDParser::RuleAssignment;
}

void OMDParser::AssignmentContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterAssignment(this);
}

void OMDParser::AssignmentContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitAssignment(this);
}


std::any OMDParser::AssignmentContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitAssignment(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::AssignmentContext* OMDParser::assignment() {
  AssignmentContext *_localctx = _tracker.createInstance<AssignmentContext>(_ctx, getState());
  enterRule(_localctx, 4, OMDParser::RuleAssignment);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(100);
    match(OMDParser::ID);
    setState(101);
    match(OMDParser::ASSIGNEQUAL);
    setState(102);
    constant();
    setState(103);
    match(OMDParser::SEMICOLON);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ConstantContext ------------------------------------------------------------------

OMDParser::ConstantContext::ConstantContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::IntConstContext* OMDParser::ConstantContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

OMDParser::FloatConstContext* OMDParser::ConstantContext::floatConst() {
  return getRuleContext<OMDParser::FloatConstContext>(0);
}

OMDParser::VectorConstContext* OMDParser::ConstantContext::vectorConst() {
  return getRuleContext<OMDParser::VectorConstContext>(0);
}

tree::TerminalNode* OMDParser::ConstantContext::TRUE() {
  return getToken(OMDParser::TRUE, 0);
}

tree::TerminalNode* OMDParser::ConstantContext::FALSE() {
  return getToken(OMDParser::FALSE, 0);
}

tree::TerminalNode* OMDParser::ConstantContext::ID() {
  return getToken(OMDParser::ID, 0);
}

tree::TerminalNode* OMDParser::ConstantContext::StringLiteral() {
  return getToken(OMDParser::StringLiteral, 0);
}


size_t OMDParser::ConstantContext::getRuleIndex() const {
  return OMDParser::RuleConstant;
}

void OMDParser::ConstantContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterConstant(this);
}

void OMDParser::ConstantContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitConstant(this);
}


std::any OMDParser::ConstantContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitConstant(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::ConstantContext* OMDParser::constant() {
  ConstantContext *_localctx = _tracker.createInstance<ConstantContext>(_ctx, getState());
  enterRule(_localctx, 6, OMDParser::RuleConstant);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(112);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::NUM_LONG:
      case OMDParser::NUM_INT: {
        enterOuterAlt(_localctx, 1);
        setState(105);
        intConst();
        break;
      }

      case OMDParser::NUM_FLOAT:
      case OMDParser::NUM_DOUBLE: {
        enterOuterAlt(_localctx, 2);
        setState(106);
        floatConst();
        break;
      }

      case OMDParser::LPAREN: {
        enterOuterAlt(_localctx, 3);
        setState(107);
        vectorConst();
        break;
      }

      case OMDParser::TRUE: {
        enterOuterAlt(_localctx, 4);
        setState(108);
        match(OMDParser::TRUE);
        break;
      }

      case OMDParser::FALSE: {
        enterOuterAlt(_localctx, 5);
        setState(109);
        match(OMDParser::FALSE);
        break;
      }

      case OMDParser::ID: {
        enterOuterAlt(_localctx, 6);
        setState(110);
        match(OMDParser::ID);
        break;
      }

      case OMDParser::StringLiteral: {
        enterOuterAlt(_localctx, 7);
        setState(111);
        match(OMDParser::StringLiteral);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ComponentblockContext ------------------------------------------------------------------

OMDParser::ComponentblockContext::ComponentblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::ComponentblockContext::COMPONENT() {
  return getToken(OMDParser::COMPONENT, 0);
}

tree::TerminalNode* OMDParser::ComponentblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::ComponentblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::ComponentblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::ComponentblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::ComponentblockContext::getRuleIndex() const {
  return OMDParser::RuleComponentblock;
}

void OMDParser::ComponentblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterComponentblock(this);
}

void OMDParser::ComponentblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitComponentblock(this);
}


std::any OMDParser::ComponentblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitComponentblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::ComponentblockContext* OMDParser::componentblock() {
  ComponentblockContext *_localctx = _tracker.createInstance<ComponentblockContext>(_ctx, getState());
  enterRule(_localctx, 8, OMDParser::RuleComponentblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(114);
    match(OMDParser::COMPONENT);
    setState(115);
    match(OMDParser::LCURLY);
    setState(119);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(116);
      assignment();
      setState(121);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(122);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ZconstraintblockContext ------------------------------------------------------------------

OMDParser::ZconstraintblockContext::ZconstraintblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::ZconstraintblockContext::ZCONSTRAINT() {
  return getToken(OMDParser::ZCONSTRAINT, 0);
}

tree::TerminalNode* OMDParser::ZconstraintblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::ZconstraintblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::ZconstraintblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::ZconstraintblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::ZconstraintblockContext::getRuleIndex() const {
  return OMDParser::RuleZconstraintblock;
}

void OMDParser::ZconstraintblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterZconstraintblock(this);
}

void OMDParser::ZconstraintblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitZconstraintblock(this);
}


std::any OMDParser::ZconstraintblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitZconstraintblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::ZconstraintblockContext* OMDParser::zconstraintblock() {
  ZconstraintblockContext *_localctx = _tracker.createInstance<ZconstraintblockContext>(_ctx, getState());
  enterRule(_localctx, 10, OMDParser::RuleZconstraintblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(124);
    match(OMDParser::ZCONSTRAINT);
    setState(125);
    match(OMDParser::LCURLY);
    setState(129);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(126);
      assignment();
      setState(131);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(132);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- RestraintblockContext ------------------------------------------------------------------

OMDParser::RestraintblockContext::RestraintblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::RestraintblockContext::RESTRAINT() {
  return getToken(OMDParser::RESTRAINT, 0);
}

tree::TerminalNode* OMDParser::RestraintblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::RestraintblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::RestraintblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::RestraintblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::RestraintblockContext::getRuleIndex() const {
  return OMDParser::RuleRestraintblock;
}

void OMDParser::RestraintblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterRestraintblock(this);
}

void OMDParser::RestraintblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitRestraintblock(this);
}


std::any OMDParser::RestraintblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitRestraintblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::RestraintblockContext* OMDParser::restraintblock() {
  RestraintblockContext *_localctx = _tracker.createInstance<RestraintblockContext>(_ctx, getState());
  enterRule(_localctx, 12, OMDParser::RuleRestraintblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(134);
    match(OMDParser::RESTRAINT);
    setState(135);
    match(OMDParser::LCURLY);
    setState(139);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(136);
      assignment();
      setState(141);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(142);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- FlucqblockContext ------------------------------------------------------------------

OMDParser::FlucqblockContext::FlucqblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::FlucqblockContext::FLUCQ() {
  return getToken(OMDParser::FLUCQ, 0);
}

tree::TerminalNode* OMDParser::FlucqblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::FlucqblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::FlucqblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::FlucqblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::FlucqblockContext::getRuleIndex() const {
  return OMDParser::RuleFlucqblock;
}

void OMDParser::FlucqblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterFlucqblock(this);
}

void OMDParser::FlucqblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitFlucqblock(this);
}


std::any OMDParser::FlucqblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitFlucqblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::FlucqblockContext* OMDParser::flucqblock() {
  FlucqblockContext *_localctx = _tracker.createInstance<FlucqblockContext>(_ctx, getState());
  enterRule(_localctx, 14, OMDParser::RuleFlucqblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(144);
    match(OMDParser::FLUCQ);
    setState(145);
    match(OMDParser::LCURLY);
    setState(149);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(146);
      assignment();
      setState(151);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(152);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- RnemdblockContext ------------------------------------------------------------------

OMDParser::RnemdblockContext::RnemdblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::RnemdblockContext::RNEMD() {
  return getToken(OMDParser::RNEMD, 0);
}

tree::TerminalNode* OMDParser::RnemdblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::RnemdblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::RnemdblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::RnemdblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::RnemdblockContext::getRuleIndex() const {
  return OMDParser::RuleRnemdblock;
}

void OMDParser::RnemdblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterRnemdblock(this);
}

void OMDParser::RnemdblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitRnemdblock(this);
}


std::any OMDParser::RnemdblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitRnemdblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::RnemdblockContext* OMDParser::rnemdblock() {
  RnemdblockContext *_localctx = _tracker.createInstance<RnemdblockContext>(_ctx, getState());
  enterRule(_localctx, 16, OMDParser::RuleRnemdblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(154);
    match(OMDParser::RNEMD);
    setState(155);
    match(OMDParser::LCURLY);
    setState(159);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(156);
      assignment();
      setState(161);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(162);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- LightblockContext ------------------------------------------------------------------

OMDParser::LightblockContext::LightblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::LightblockContext::LIGHT() {
  return getToken(OMDParser::LIGHT, 0);
}

tree::TerminalNode* OMDParser::LightblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::LightblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::LightblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::LightblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::LightblockContext::getRuleIndex() const {
  return OMDParser::RuleLightblock;
}

void OMDParser::LightblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterLightblock(this);
}

void OMDParser::LightblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitLightblock(this);
}


std::any OMDParser::LightblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitLightblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::LightblockContext* OMDParser::lightblock() {
  LightblockContext *_localctx = _tracker.createInstance<LightblockContext>(_ctx, getState());
  enterRule(_localctx, 18, OMDParser::RuleLightblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(164);
    match(OMDParser::LIGHT);
    setState(165);
    match(OMDParser::LCURLY);
    setState(169);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(166);
      assignment();
      setState(171);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(172);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- MinimizerblockContext ------------------------------------------------------------------

OMDParser::MinimizerblockContext::MinimizerblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::MinimizerblockContext::MINIMIZER() {
  return getToken(OMDParser::MINIMIZER, 0);
}

tree::TerminalNode* OMDParser::MinimizerblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::MinimizerblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::MinimizerblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::MinimizerblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::MinimizerblockContext::getRuleIndex() const {
  return OMDParser::RuleMinimizerblock;
}

void OMDParser::MinimizerblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterMinimizerblock(this);
}

void OMDParser::MinimizerblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitMinimizerblock(this);
}


std::any OMDParser::MinimizerblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitMinimizerblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::MinimizerblockContext* OMDParser::minimizerblock() {
  MinimizerblockContext *_localctx = _tracker.createInstance<MinimizerblockContext>(_ctx, getState());
  enterRule(_localctx, 20, OMDParser::RuleMinimizerblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(174);
    match(OMDParser::MINIMIZER);
    setState(175);
    match(OMDParser::LCURLY);
    setState(179);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(176);
      assignment();
      setState(181);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(182);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- MoleculeblockContext ------------------------------------------------------------------

OMDParser::MoleculeblockContext::MoleculeblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::MoleculeblockContext::MOLECULE() {
  return getToken(OMDParser::MOLECULE, 0);
}

tree::TerminalNode* OMDParser::MoleculeblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::MoleculeblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::MoleculestatementContext *> OMDParser::MoleculeblockContext::moleculestatement() {
  return getRuleContexts<OMDParser::MoleculestatementContext>();
}

OMDParser::MoleculestatementContext* OMDParser::MoleculeblockContext::moleculestatement(size_t i) {
  return getRuleContext<OMDParser::MoleculestatementContext>(i);
}


size_t OMDParser::MoleculeblockContext::getRuleIndex() const {
  return OMDParser::RuleMoleculeblock;
}

void OMDParser::MoleculeblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterMoleculeblock(this);
}

void OMDParser::MoleculeblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitMoleculeblock(this);
}


std::any OMDParser::MoleculeblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitMoleculeblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::MoleculeblockContext* OMDParser::moleculeblock() {
  MoleculeblockContext *_localctx = _tracker.createInstance<MoleculeblockContext>(_ctx, getState());
  enterRule(_localctx, 22, OMDParser::RuleMoleculeblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(184);
    match(OMDParser::MOLECULE);
    setState(185);
    match(OMDParser::LCURLY);
    setState(189);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 7) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 7)) & 144115188075857151) != 0)) {
      setState(186);
      moleculestatement();
      setState(191);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(192);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- MoleculestatementContext ------------------------------------------------------------------

OMDParser::MoleculestatementContext::MoleculestatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::MoleculestatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

OMDParser::AtomblockContext* OMDParser::MoleculestatementContext::atomblock() {
  return getRuleContext<OMDParser::AtomblockContext>(0);
}

OMDParser::BondblockContext* OMDParser::MoleculestatementContext::bondblock() {
  return getRuleContext<OMDParser::BondblockContext>(0);
}

OMDParser::BendblockContext* OMDParser::MoleculestatementContext::bendblock() {
  return getRuleContext<OMDParser::BendblockContext>(0);
}

OMDParser::TorsionblockContext* OMDParser::MoleculestatementContext::torsionblock() {
  return getRuleContext<OMDParser::TorsionblockContext>(0);
}

OMDParser::InversionblockContext* OMDParser::MoleculestatementContext::inversionblock() {
  return getRuleContext<OMDParser::InversionblockContext>(0);
}

OMDParser::RigidbodyblockContext* OMDParser::MoleculestatementContext::rigidbodyblock() {
  return getRuleContext<OMDParser::RigidbodyblockContext>(0);
}

OMDParser::CutoffgroupblockContext* OMDParser::MoleculestatementContext::cutoffgroupblock() {
  return getRuleContext<OMDParser::CutoffgroupblockContext>(0);
}

OMDParser::ConstraintblockContext* OMDParser::MoleculestatementContext::constraintblock() {
  return getRuleContext<OMDParser::ConstraintblockContext>(0);
}

OMDParser::SequencestringContext* OMDParser::MoleculestatementContext::sequencestring() {
  return getRuleContext<OMDParser::SequencestringContext>(0);
}


size_t OMDParser::MoleculestatementContext::getRuleIndex() const {
  return OMDParser::RuleMoleculestatement;
}

void OMDParser::MoleculestatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterMoleculestatement(this);
}

void OMDParser::MoleculestatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitMoleculestatement(this);
}


std::any OMDParser::MoleculestatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitMoleculestatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::MoleculestatementContext* OMDParser::moleculestatement() {
  MoleculestatementContext *_localctx = _tracker.createInstance<MoleculestatementContext>(_ctx, getState());
  enterRule(_localctx, 24, OMDParser::RuleMoleculestatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(204);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(194);
        assignment();
        break;
      }

      case OMDParser::ATOM: {
        enterOuterAlt(_localctx, 2);
        setState(195);
        atomblock();
        break;
      }

      case OMDParser::BOND: {
        enterOuterAlt(_localctx, 3);
        setState(196);
        bondblock();
        break;
      }

      case OMDParser::BEND: {
        enterOuterAlt(_localctx, 4);
        setState(197);
        bendblock();
        break;
      }

      case OMDParser::TORSION: {
        enterOuterAlt(_localctx, 5);
        setState(198);
        torsionblock();
        break;
      }

      case OMDParser::INVERSION: {
        enterOuterAlt(_localctx, 6);
        setState(199);
        inversionblock();
        break;
      }

      case OMDParser::RIGIDBODY: {
        enterOuterAlt(_localctx, 7);
        setState(200);
        rigidbodyblock();
        break;
      }

      case OMDParser::CUTOFFGROUP: {
        enterOuterAlt(_localctx, 8);
        setState(201);
        cutoffgroupblock();
        break;
      }

      case OMDParser::CONSTRAINT: {
        enterOuterAlt(_localctx, 9);
        setState(202);
        constraintblock();
        break;
      }

      case OMDParser::SEQUENCE: {
        enterOuterAlt(_localctx, 10);
        setState(203);
        sequencestring();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- AtomblockContext ------------------------------------------------------------------

OMDParser::AtomblockContext::AtomblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::AtomblockContext::ATOM() {
  return getToken(OMDParser::ATOM, 0);
}

tree::TerminalNode* OMDParser::AtomblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::AtomblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::AtomblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

tree::TerminalNode* OMDParser::AtomblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::AtomblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AtomstatementContext *> OMDParser::AtomblockContext::atomstatement() {
  return getRuleContexts<OMDParser::AtomstatementContext>();
}

OMDParser::AtomstatementContext* OMDParser::AtomblockContext::atomstatement(size_t i) {
  return getRuleContext<OMDParser::AtomstatementContext>(i);
}


size_t OMDParser::AtomblockContext::getRuleIndex() const {
  return OMDParser::RuleAtomblock;
}

void OMDParser::AtomblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterAtomblock(this);
}

void OMDParser::AtomblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitAtomblock(this);
}


std::any OMDParser::AtomblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitAtomblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::AtomblockContext* OMDParser::atomblock() {
  AtomblockContext *_localctx = _tracker.createInstance<AtomblockContext>(_ctx, getState());
  enterRule(_localctx, 26, OMDParser::RuleAtomblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(206);
    match(OMDParser::ATOM);
    setState(207);
    match(OMDParser::LBRACKET);
    setState(208);
    intConst();
    setState(209);
    match(OMDParser::RBRACKET);
    setState(210);
    match(OMDParser::LCURLY);
    setState(214);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 21) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 21)) & 8796101410819) != 0)) {
      setState(211);
      atomstatement();
      setState(216);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(217);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- AtomstatementContext ------------------------------------------------------------------

OMDParser::AtomstatementContext::AtomstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::AtomstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::AtomstatementContext::POSITION() {
  return getToken(OMDParser::POSITION, 0);
}

tree::TerminalNode* OMDParser::AtomstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::DoubleNumberTupleContext* OMDParser::AtomstatementContext::doubleNumberTuple() {
  return getRuleContext<OMDParser::DoubleNumberTupleContext>(0);
}

tree::TerminalNode* OMDParser::AtomstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::AtomstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}

tree::TerminalNode* OMDParser::AtomstatementContext::ORIENTATION() {
  return getToken(OMDParser::ORIENTATION, 0);
}

tree::TerminalNode* OMDParser::AtomstatementContext::CHARGE() {
  return getToken(OMDParser::CHARGE, 0);
}

OMDParser::FloatConstContext* OMDParser::AtomstatementContext::floatConst() {
  return getRuleContext<OMDParser::FloatConstContext>(0);
}


size_t OMDParser::AtomstatementContext::getRuleIndex() const {
  return OMDParser::RuleAtomstatement;
}

void OMDParser::AtomstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterAtomstatement(this);
}

void OMDParser::AtomstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitAtomstatement(this);
}


std::any OMDParser::AtomstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitAtomstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::AtomstatementContext* OMDParser::atomstatement() {
  AtomstatementContext *_localctx = _tracker.createInstance<AtomstatementContext>(_ctx, getState());
  enterRule(_localctx, 28, OMDParser::RuleAtomstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(238);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(219);
        assignment();
        break;
      }

      case OMDParser::POSITION: {
        enterOuterAlt(_localctx, 2);
        setState(220);
        match(OMDParser::POSITION);
        setState(221);
        match(OMDParser::LPAREN);
        setState(222);
        doubleNumberTuple();
        setState(223);
        match(OMDParser::RPAREN);
        setState(224);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::ORIENTATION: {
        enterOuterAlt(_localctx, 3);
        setState(226);
        match(OMDParser::ORIENTATION);
        setState(227);
        match(OMDParser::LPAREN);
        setState(228);
        doubleNumberTuple();
        setState(229);
        match(OMDParser::RPAREN);
        setState(230);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CHARGE: {
        enterOuterAlt(_localctx, 4);
        setState(232);
        match(OMDParser::CHARGE);
        setState(233);
        match(OMDParser::LPAREN);
        setState(234);
        floatConst();
        setState(235);
        match(OMDParser::RPAREN);
        setState(236);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- BondblockContext ------------------------------------------------------------------

OMDParser::BondblockContext::BondblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::BondblockContext::BOND() {
  return getToken(OMDParser::BOND, 0);
}

tree::TerminalNode* OMDParser::BondblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::BondblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::BondblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::BondblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::BondblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::BondstatementContext *> OMDParser::BondblockContext::bondstatement() {
  return getRuleContexts<OMDParser::BondstatementContext>();
}

OMDParser::BondstatementContext* OMDParser::BondblockContext::bondstatement(size_t i) {
  return getRuleContext<OMDParser::BondstatementContext>(i);
}


size_t OMDParser::BondblockContext::getRuleIndex() const {
  return OMDParser::RuleBondblock;
}

void OMDParser::BondblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterBondblock(this);
}

void OMDParser::BondblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitBondblock(this);
}


std::any OMDParser::BondblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitBondblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::BondblockContext* OMDParser::bondblock() {
  BondblockContext *_localctx = _tracker.createInstance<BondblockContext>(_ctx, getState());
  enterRule(_localctx, 30, OMDParser::RuleBondblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(240);
    match(OMDParser::BOND);
    setState(245);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(241);
      match(OMDParser::LBRACKET);
      setState(242);
      intConst();
      setState(243);
      match(OMDParser::RBRACKET);
    }
    setState(247);
    match(OMDParser::LCURLY);
    setState(251);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 18) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 18)) & 70368744209921) != 0)) {
      setState(248);
      bondstatement();
      setState(253);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(254);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- BondstatementContext ------------------------------------------------------------------

OMDParser::BondstatementContext::BondstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::BondstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::BondstatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::BondstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::BondstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::BondstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::BondstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}

tree::TerminalNode* OMDParser::BondstatementContext::FIXED() {
  return getToken(OMDParser::FIXED, 0);
}

OMDParser::FloatConstContext* OMDParser::BondstatementContext::floatConst() {
  return getRuleContext<OMDParser::FloatConstContext>(0);
}

tree::TerminalNode* OMDParser::BondstatementContext::HARMONIC() {
  return getToken(OMDParser::HARMONIC, 0);
}

OMDParser::DoubleNumberTupleContext* OMDParser::BondstatementContext::doubleNumberTuple() {
  return getRuleContext<OMDParser::DoubleNumberTupleContext>(0);
}

tree::TerminalNode* OMDParser::BondstatementContext::CUBIC() {
  return getToken(OMDParser::CUBIC, 0);
}

tree::TerminalNode* OMDParser::BondstatementContext::QUARTIC() {
  return getToken(OMDParser::QUARTIC, 0);
}

tree::TerminalNode* OMDParser::BondstatementContext::POLYNOMIAL() {
  return getToken(OMDParser::POLYNOMIAL, 0);
}

tree::TerminalNode* OMDParser::BondstatementContext::MORSE() {
  return getToken(OMDParser::MORSE, 0);
}


size_t OMDParser::BondstatementContext::getRuleIndex() const {
  return OMDParser::RuleBondstatement;
}

void OMDParser::BondstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterBondstatement(this);
}

void OMDParser::BondstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitBondstatement(this);
}


std::any OMDParser::BondstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitBondstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::BondstatementContext* OMDParser::bondstatement() {
  BondstatementContext *_localctx = _tracker.createInstance<BondstatementContext>(_ctx, getState());
  enterRule(_localctx, 32, OMDParser::RuleBondstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(299);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(256);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(257);
        match(OMDParser::MEMBERS);
        setState(258);
        match(OMDParser::LPAREN);
        setState(259);
        inttuple();
        setState(260);
        match(OMDParser::RPAREN);
        setState(261);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::FIXED: {
        enterOuterAlt(_localctx, 3);
        setState(263);
        match(OMDParser::FIXED);
        setState(264);
        match(OMDParser::LPAREN);
        setState(265);
        floatConst();
        setState(266);
        match(OMDParser::RPAREN);
        setState(267);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 4);
        setState(269);
        match(OMDParser::HARMONIC);
        setState(270);
        match(OMDParser::LPAREN);
        setState(271);
        doubleNumberTuple();
        setState(272);
        match(OMDParser::RPAREN);
        setState(273);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 5);
        setState(275);
        match(OMDParser::CUBIC);
        setState(276);
        match(OMDParser::LPAREN);
        setState(277);
        doubleNumberTuple();
        setState(278);
        match(OMDParser::RPAREN);
        setState(279);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 6);
        setState(281);
        match(OMDParser::QUARTIC);
        setState(282);
        match(OMDParser::LPAREN);
        setState(283);
        doubleNumberTuple();
        setState(284);
        match(OMDParser::RPAREN);
        setState(285);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 7);
        setState(287);
        match(OMDParser::POLYNOMIAL);
        setState(288);
        match(OMDParser::LPAREN);
        setState(289);
        doubleNumberTuple();
        setState(290);
        match(OMDParser::RPAREN);
        setState(291);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::MORSE: {
        enterOuterAlt(_localctx, 8);
        setState(293);
        match(OMDParser::MORSE);
        setState(294);
        match(OMDParser::LPAREN);
        setState(295);
        doubleNumberTuple();
        setState(296);
        match(OMDParser::RPAREN);
        setState(297);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- BendblockContext ------------------------------------------------------------------

OMDParser::BendblockContext::BendblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::BendblockContext::BEND() {
  return getToken(OMDParser::BEND, 0);
}

tree::TerminalNode* OMDParser::BendblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::BendblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::BendblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::BendblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::BendblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::BendstatementContext *> OMDParser::BendblockContext::bendstatement() {
  return getRuleContexts<OMDParser::BendstatementContext>();
}

OMDParser::BendstatementContext* OMDParser::BendblockContext::bendstatement(size_t i) {
  return getRuleContext<OMDParser::BendstatementContext>(i);
}


size_t OMDParser::BendblockContext::getRuleIndex() const {
  return OMDParser::RuleBendblock;
}

void OMDParser::BendblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterBendblock(this);
}

void OMDParser::BendblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitBendblock(this);
}


std::any OMDParser::BendblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitBendblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::BendblockContext* OMDParser::bendblock() {
  BendblockContext *_localctx = _tracker.createInstance<BendblockContext>(_ctx, getState());
  enterRule(_localctx, 34, OMDParser::RuleBendblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(301);
    match(OMDParser::BEND);
    setState(306);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(302);
      match(OMDParser::LBRACKET);
      setState(303);
      intConst();
      setState(304);
      match(OMDParser::RBRACKET);
    }
    setState(308);
    match(OMDParser::LCURLY);
    setState(312);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 18) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 18)) & 70368744422401) != 0)) {
      setState(309);
      bendstatement();
      setState(314);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(315);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- BendstatementContext ------------------------------------------------------------------

OMDParser::BendstatementContext::BendstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::BendstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::BendstatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::BendstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::BendstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::HARMONIC() {
  return getToken(OMDParser::HARMONIC, 0);
}

OMDParser::DoubleNumberTupleContext* OMDParser::BendstatementContext::doubleNumberTuple() {
  return getRuleContext<OMDParser::DoubleNumberTupleContext>(0);
}

tree::TerminalNode* OMDParser::BendstatementContext::GHOSTBEND() {
  return getToken(OMDParser::GHOSTBEND, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::UREYBRADLEY() {
  return getToken(OMDParser::UREYBRADLEY, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::CUBIC() {
  return getToken(OMDParser::CUBIC, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::QUARTIC() {
  return getToken(OMDParser::QUARTIC, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::POLYNOMIAL() {
  return getToken(OMDParser::POLYNOMIAL, 0);
}

tree::TerminalNode* OMDParser::BendstatementContext::COSINE() {
  return getToken(OMDParser::COSINE, 0);
}


size_t OMDParser::BendstatementContext::getRuleIndex() const {
  return OMDParser::RuleBendstatement;
}

void OMDParser::BendstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterBendstatement(this);
}

void OMDParser::BendstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitBendstatement(this);
}


std::any OMDParser::BendstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitBendstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::BendstatementContext* OMDParser::bendstatement() {
  BendstatementContext *_localctx = _tracker.createInstance<BendstatementContext>(_ctx, getState());
  enterRule(_localctx, 36, OMDParser::RuleBendstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(366);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(317);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(318);
        match(OMDParser::MEMBERS);
        setState(319);
        match(OMDParser::LPAREN);
        setState(320);
        inttuple();
        setState(321);
        match(OMDParser::RPAREN);
        setState(322);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 3);
        setState(324);
        match(OMDParser::HARMONIC);
        setState(325);
        match(OMDParser::LPAREN);
        setState(326);
        doubleNumberTuple();
        setState(327);
        match(OMDParser::RPAREN);
        setState(328);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::GHOSTBEND: {
        enterOuterAlt(_localctx, 4);
        setState(330);
        match(OMDParser::GHOSTBEND);
        setState(331);
        match(OMDParser::LPAREN);
        setState(332);
        doubleNumberTuple();
        setState(333);
        match(OMDParser::RPAREN);
        setState(334);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::UREYBRADLEY: {
        enterOuterAlt(_localctx, 5);
        setState(336);
        match(OMDParser::UREYBRADLEY);
        setState(337);
        match(OMDParser::LPAREN);
        setState(338);
        doubleNumberTuple();
        setState(339);
        match(OMDParser::RPAREN);
        setState(340);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 6);
        setState(342);
        match(OMDParser::CUBIC);
        setState(343);
        match(OMDParser::LPAREN);
        setState(344);
        doubleNumberTuple();
        setState(345);
        match(OMDParser::RPAREN);
        setState(346);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 7);
        setState(348);
        match(OMDParser::QUARTIC);
        setState(349);
        match(OMDParser::LPAREN);
        setState(350);
        doubleNumberTuple();
        setState(351);
        match(OMDParser::RPAREN);
        setState(352);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 8);
        setState(354);
        match(OMDParser::POLYNOMIAL);
        setState(355);
        match(OMDParser::LPAREN);
        setState(356);
        doubleNumberTuple();
        setState(357);
        match(OMDParser::RPAREN);
        setState(358);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::COSINE: {
        enterOuterAlt(_localctx, 9);
        setState(360);
        match(OMDParser::COSINE);
        setState(361);
        match(OMDParser::LPAREN);
        setState(362);
        doubleNumberTuple();
        setState(363);
        match(OMDParser::RPAREN);
        setState(364);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- TorsionblockContext ------------------------------------------------------------------

OMDParser::TorsionblockContext::TorsionblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::TorsionblockContext::TORSION() {
  return getToken(OMDParser::TORSION, 0);
}

tree::TerminalNode* OMDParser::TorsionblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::TorsionblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::TorsionblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::TorsionblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::TorsionblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::TorsionstatementContext *> OMDParser::TorsionblockContext::torsionstatement() {
  return getRuleContexts<OMDParser::TorsionstatementContext>();
}

OMDParser::TorsionstatementContext* OMDParser::TorsionblockContext::torsionstatement(size_t i) {
  return getRuleContext<OMDParser::TorsionstatementContext>(i);
}


size_t OMDParser::TorsionblockContext::getRuleIndex() const {
  return OMDParser::RuleTorsionblock;
}

void OMDParser::TorsionblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterTorsionblock(this);
}

void OMDParser::TorsionblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitTorsionblock(this);
}


std::any OMDParser::TorsionblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitTorsionblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::TorsionblockContext* OMDParser::torsionblock() {
  TorsionblockContext *_localctx = _tracker.createInstance<TorsionblockContext>(_ctx, getState());
  enterRule(_localctx, 38, OMDParser::RuleTorsionblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(368);
    match(OMDParser::TORSION);
    setState(373);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(369);
      match(OMDParser::LBRACKET);
      setState(370);
      intConst();
      setState(371);
      match(OMDParser::RBRACKET);
    }
    setState(375);
    match(OMDParser::LCURLY);
    setState(379);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 18) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 18)) & 70368748125185) != 0)) {
      setState(376);
      torsionstatement();
      setState(381);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(382);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- TorsionstatementContext ------------------------------------------------------------------

OMDParser::TorsionstatementContext::TorsionstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::TorsionstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::TorsionstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::GHOSTTORSION() {
  return getToken(OMDParser::GHOSTTORSION, 0);
}

OMDParser::DoubleNumberTupleContext* OMDParser::TorsionstatementContext::doubleNumberTuple() {
  return getRuleContext<OMDParser::DoubleNumberTupleContext>(0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::CUBIC() {
  return getToken(OMDParser::CUBIC, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::QUARTIC() {
  return getToken(OMDParser::QUARTIC, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::POLYNOMIAL() {
  return getToken(OMDParser::POLYNOMIAL, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::CHARMM() {
  return getToken(OMDParser::CHARMM, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::OPLS() {
  return getToken(OMDParser::OPLS, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::TRAPPE() {
  return getToken(OMDParser::TRAPPE, 0);
}

tree::TerminalNode* OMDParser::TorsionstatementContext::HARMONIC() {
  return getToken(OMDParser::HARMONIC, 0);
}


size_t OMDParser::TorsionstatementContext::getRuleIndex() const {
  return OMDParser::RuleTorsionstatement;
}

void OMDParser::TorsionstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterTorsionstatement(this);
}

void OMDParser::TorsionstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitTorsionstatement(this);
}


std::any OMDParser::TorsionstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitTorsionstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::TorsionstatementContext* OMDParser::torsionstatement() {
  TorsionstatementContext *_localctx = _tracker.createInstance<TorsionstatementContext>(_ctx, getState());
  enterRule(_localctx, 40, OMDParser::RuleTorsionstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(439);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(384);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(385);
        match(OMDParser::MEMBERS);
        setState(386);
        match(OMDParser::LPAREN);
        setState(387);
        inttuple();
        setState(388);
        match(OMDParser::RPAREN);
        setState(389);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::GHOSTTORSION: {
        enterOuterAlt(_localctx, 3);
        setState(391);
        match(OMDParser::GHOSTTORSION);
        setState(392);
        match(OMDParser::LPAREN);
        setState(393);
        doubleNumberTuple();
        setState(394);
        match(OMDParser::RPAREN);
        setState(395);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 4);
        setState(397);
        match(OMDParser::CUBIC);
        setState(398);
        match(OMDParser::LPAREN);
        setState(399);
        doubleNumberTuple();
        setState(400);
        match(OMDParser::RPAREN);
        setState(401);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 5);
        setState(403);
        match(OMDParser::QUARTIC);
        setState(404);
        match(OMDParser::LPAREN);
        setState(405);
        doubleNumberTuple();
        setState(406);
        match(OMDParser::RPAREN);
        setState(407);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 6);
        setState(409);
        match(OMDParser::POLYNOMIAL);
        setState(410);
        match(OMDParser::LPAREN);
        setState(411);
        doubleNumberTuple();
        setState(412);
        match(OMDParser::RPAREN);
        setState(413);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CHARMM: {
        enterOuterAlt(_localctx, 7);
        setState(415);
        match(OMDParser::CHARMM);
        setState(416);
        match(OMDParser::LPAREN);
        setState(417);
        doubleNumberTuple();
        setState(418);
        match(OMDParser::RPAREN);
        setState(419);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::OPLS: {
        enterOuterAlt(_localctx, 8);
        setState(421);
        match(OMDParser::OPLS);
        setState(422);
        match(OMDParser::LPAREN);
        setState(423);
        doubleNumberTuple();
        setState(424);
        match(OMDParser::RPAREN);
        setState(425);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::TRAPPE: {
        enterOuterAlt(_localctx, 9);
        setState(427);
        match(OMDParser::TRAPPE);
        setState(428);
        match(OMDParser::LPAREN);
        setState(429);
        doubleNumberTuple();
        setState(430);
        match(OMDParser::RPAREN);
        setState(431);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 10);
        setState(433);
        match(OMDParser::HARMONIC);
        setState(434);
        match(OMDParser::LPAREN);
        setState(435);
        doubleNumberTuple();
        setState(436);
        match(OMDParser::RPAREN);
        setState(437);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- InversionblockContext ------------------------------------------------------------------

OMDParser::InversionblockContext::InversionblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::InversionblockContext::INVERSION() {
  return getToken(OMDParser::INVERSION, 0);
}

tree::TerminalNode* OMDParser::InversionblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::InversionblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::InversionblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::InversionblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::InversionblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::InversionstatementContext *> OMDParser::InversionblockContext::inversionstatement() {
  return getRuleContexts<OMDParser::InversionstatementContext>();
}

OMDParser::InversionstatementContext* OMDParser::InversionblockContext::inversionstatement(size_t i) {
  return getRuleContext<OMDParser::InversionstatementContext>(i);
}


size_t OMDParser::InversionblockContext::getRuleIndex() const {
  return OMDParser::RuleInversionblock;
}

void OMDParser::InversionblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterInversionblock(this);
}

void OMDParser::InversionblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitInversionblock(this);
}


std::any OMDParser::InversionblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitInversionblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::InversionblockContext* OMDParser::inversionblock() {
  InversionblockContext *_localctx = _tracker.createInstance<InversionblockContext>(_ctx, getState());
  enterRule(_localctx, 42, OMDParser::RuleInversionblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(441);
    match(OMDParser::INVERSION);
    setState(446);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(442);
      match(OMDParser::LBRACKET);
      setState(443);
      intConst();
      setState(444);
      match(OMDParser::RBRACKET);
    }
    setState(448);
    match(OMDParser::LCURLY);
    setState(452);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 19) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 19)) & 35184403546627) != 0)) {
      setState(449);
      inversionstatement();
      setState(454);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(455);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- InversionstatementContext ------------------------------------------------------------------

OMDParser::InversionstatementContext::InversionstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::InversionstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::CENTER() {
  return getToken(OMDParser::CENTER, 0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::IntConstContext* OMDParser::InversionstatementContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::SATELLITES() {
  return getToken(OMDParser::SATELLITES, 0);
}

OMDParser::InttupleContext* OMDParser::InversionstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::AMBERIMPROPER() {
  return getToken(OMDParser::AMBERIMPROPER, 0);
}

OMDParser::DoubleNumberTupleContext* OMDParser::InversionstatementContext::doubleNumberTuple() {
  return getRuleContext<OMDParser::DoubleNumberTupleContext>(0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::IMPROPERCOSINE() {
  return getToken(OMDParser::IMPROPERCOSINE, 0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::HARMONIC() {
  return getToken(OMDParser::HARMONIC, 0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::CENTRALATOMHEIGHT() {
  return getToken(OMDParser::CENTRALATOMHEIGHT, 0);
}

tree::TerminalNode* OMDParser::InversionstatementContext::DREIDING() {
  return getToken(OMDParser::DREIDING, 0);
}


size_t OMDParser::InversionstatementContext::getRuleIndex() const {
  return OMDParser::RuleInversionstatement;
}

void OMDParser::InversionstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterInversionstatement(this);
}

void OMDParser::InversionstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitInversionstatement(this);
}


std::any OMDParser::InversionstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitInversionstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::InversionstatementContext* OMDParser::inversionstatement() {
  InversionstatementContext *_localctx = _tracker.createInstance<InversionstatementContext>(_ctx, getState());
  enterRule(_localctx, 44, OMDParser::RuleInversionstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(500);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(457);
        assignment();
        break;
      }

      case OMDParser::CENTER: {
        enterOuterAlt(_localctx, 2);
        setState(458);
        match(OMDParser::CENTER);
        setState(459);
        match(OMDParser::LPAREN);
        setState(460);
        intConst();
        setState(461);
        match(OMDParser::RPAREN);
        setState(462);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::SATELLITES: {
        enterOuterAlt(_localctx, 3);
        setState(464);
        match(OMDParser::SATELLITES);
        setState(465);
        match(OMDParser::LPAREN);
        setState(466);
        inttuple();
        setState(467);
        match(OMDParser::RPAREN);
        setState(468);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::AMBERIMPROPER: {
        enterOuterAlt(_localctx, 4);
        setState(470);
        match(OMDParser::AMBERIMPROPER);
        setState(471);
        match(OMDParser::LPAREN);
        setState(472);
        doubleNumberTuple();
        setState(473);
        match(OMDParser::RPAREN);
        setState(474);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::IMPROPERCOSINE: {
        enterOuterAlt(_localctx, 5);
        setState(476);
        match(OMDParser::IMPROPERCOSINE);
        setState(477);
        match(OMDParser::LPAREN);
        setState(478);
        doubleNumberTuple();
        setState(479);
        match(OMDParser::RPAREN);
        setState(480);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 6);
        setState(482);
        match(OMDParser::HARMONIC);
        setState(483);
        match(OMDParser::LPAREN);
        setState(484);
        doubleNumberTuple();
        setState(485);
        match(OMDParser::RPAREN);
        setState(486);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CENTRALATOMHEIGHT: {
        enterOuterAlt(_localctx, 7);
        setState(488);
        match(OMDParser::CENTRALATOMHEIGHT);
        setState(489);
        match(OMDParser::LPAREN);
        setState(490);
        doubleNumberTuple();
        setState(491);
        match(OMDParser::RPAREN);
        setState(492);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::DREIDING: {
        enterOuterAlt(_localctx, 8);
        setState(494);
        match(OMDParser::DREIDING);
        setState(495);
        match(OMDParser::LPAREN);
        setState(496);
        doubleNumberTuple();
        setState(497);
        match(OMDParser::RPAREN);
        setState(498);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- RigidbodyblockContext ------------------------------------------------------------------

OMDParser::RigidbodyblockContext::RigidbodyblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::RigidbodyblockContext::RIGIDBODY() {
  return getToken(OMDParser::RIGIDBODY, 0);
}

tree::TerminalNode* OMDParser::RigidbodyblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::RigidbodyblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::RigidbodyblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

tree::TerminalNode* OMDParser::RigidbodyblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::RigidbodyblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::RigidbodystatementContext *> OMDParser::RigidbodyblockContext::rigidbodystatement() {
  return getRuleContexts<OMDParser::RigidbodystatementContext>();
}

OMDParser::RigidbodystatementContext* OMDParser::RigidbodyblockContext::rigidbodystatement(size_t i) {
  return getRuleContext<OMDParser::RigidbodystatementContext>(i);
}


size_t OMDParser::RigidbodyblockContext::getRuleIndex() const {
  return OMDParser::RuleRigidbodyblock;
}

void OMDParser::RigidbodyblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterRigidbodyblock(this);
}

void OMDParser::RigidbodyblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitRigidbodyblock(this);
}


std::any OMDParser::RigidbodyblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitRigidbodyblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::RigidbodyblockContext* OMDParser::rigidbodyblock() {
  RigidbodyblockContext *_localctx = _tracker.createInstance<RigidbodyblockContext>(_ctx, getState());
  enterRule(_localctx, 46, OMDParser::RuleRigidbodyblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(502);
    match(OMDParser::RIGIDBODY);
    setState(503);
    match(OMDParser::LBRACKET);
    setState(504);
    intConst();
    setState(505);
    match(OMDParser::RBRACKET);
    setState(506);
    match(OMDParser::LCURLY);
    setState(510);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(507);
      rigidbodystatement();
      setState(512);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(513);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- RigidbodystatementContext ------------------------------------------------------------------

OMDParser::RigidbodystatementContext::RigidbodystatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::RigidbodystatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::RigidbodystatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::RigidbodystatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::RigidbodystatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::RigidbodystatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::RigidbodystatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}


size_t OMDParser::RigidbodystatementContext::getRuleIndex() const {
  return OMDParser::RuleRigidbodystatement;
}

void OMDParser::RigidbodystatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterRigidbodystatement(this);
}

void OMDParser::RigidbodystatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitRigidbodystatement(this);
}


std::any OMDParser::RigidbodystatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitRigidbodystatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::RigidbodystatementContext* OMDParser::rigidbodystatement() {
  RigidbodystatementContext *_localctx = _tracker.createInstance<RigidbodystatementContext>(_ctx, getState());
  enterRule(_localctx, 48, OMDParser::RuleRigidbodystatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(522);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(515);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(516);
        match(OMDParser::MEMBERS);
        setState(517);
        match(OMDParser::LPAREN);
        setState(518);
        inttuple();
        setState(519);
        match(OMDParser::RPAREN);
        setState(520);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- CutoffgroupblockContext ------------------------------------------------------------------

OMDParser::CutoffgroupblockContext::CutoffgroupblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::CutoffgroupblockContext::CUTOFFGROUP() {
  return getToken(OMDParser::CUTOFFGROUP, 0);
}

tree::TerminalNode* OMDParser::CutoffgroupblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::CutoffgroupblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::CutoffgroupblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::CutoffgroupblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::CutoffgroupblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::CutoffgroupstatementContext *> OMDParser::CutoffgroupblockContext::cutoffgroupstatement() {
  return getRuleContexts<OMDParser::CutoffgroupstatementContext>();
}

OMDParser::CutoffgroupstatementContext* OMDParser::CutoffgroupblockContext::cutoffgroupstatement(size_t i) {
  return getRuleContext<OMDParser::CutoffgroupstatementContext>(i);
}


size_t OMDParser::CutoffgroupblockContext::getRuleIndex() const {
  return OMDParser::RuleCutoffgroupblock;
}

void OMDParser::CutoffgroupblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterCutoffgroupblock(this);
}

void OMDParser::CutoffgroupblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitCutoffgroupblock(this);
}


std::any OMDParser::CutoffgroupblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitCutoffgroupblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::CutoffgroupblockContext* OMDParser::cutoffgroupblock() {
  CutoffgroupblockContext *_localctx = _tracker.createInstance<CutoffgroupblockContext>(_ctx, getState());
  enterRule(_localctx, 50, OMDParser::RuleCutoffgroupblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(524);
    match(OMDParser::CUTOFFGROUP);
    setState(529);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(525);
      match(OMDParser::LBRACKET);
      setState(526);
      intConst();
      setState(527);
      match(OMDParser::RBRACKET);
    }
    setState(531);
    match(OMDParser::LCURLY);
    setState(535);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(532);
      cutoffgroupstatement();
      setState(537);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(538);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- CutoffgroupstatementContext ------------------------------------------------------------------

OMDParser::CutoffgroupstatementContext::CutoffgroupstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::CutoffgroupstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::CutoffgroupstatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::CutoffgroupstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::CutoffgroupstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::CutoffgroupstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::CutoffgroupstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}


size_t OMDParser::CutoffgroupstatementContext::getRuleIndex() const {
  return OMDParser::RuleCutoffgroupstatement;
}

void OMDParser::CutoffgroupstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterCutoffgroupstatement(this);
}

void OMDParser::CutoffgroupstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitCutoffgroupstatement(this);
}


std::any OMDParser::CutoffgroupstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitCutoffgroupstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::CutoffgroupstatementContext* OMDParser::cutoffgroupstatement() {
  CutoffgroupstatementContext *_localctx = _tracker.createInstance<CutoffgroupstatementContext>(_ctx, getState());
  enterRule(_localctx, 52, OMDParser::RuleCutoffgroupstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(547);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(540);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(541);
        match(OMDParser::MEMBERS);
        setState(542);
        match(OMDParser::LPAREN);
        setState(543);
        inttuple();
        setState(544);
        match(OMDParser::RPAREN);
        setState(545);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- NodesblockContext ------------------------------------------------------------------

OMDParser::NodesblockContext::NodesblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::NodesblockContext::NODES() {
  return getToken(OMDParser::NODES, 0);
}

tree::TerminalNode* OMDParser::NodesblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::NodesblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::NodesblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::NodesblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::NodesblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::NodesstatementContext *> OMDParser::NodesblockContext::nodesstatement() {
  return getRuleContexts<OMDParser::NodesstatementContext>();
}

OMDParser::NodesstatementContext* OMDParser::NodesblockContext::nodesstatement(size_t i) {
  return getRuleContext<OMDParser::NodesstatementContext>(i);
}


size_t OMDParser::NodesblockContext::getRuleIndex() const {
  return OMDParser::RuleNodesblock;
}

void OMDParser::NodesblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterNodesblock(this);
}

void OMDParser::NodesblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitNodesblock(this);
}


std::any OMDParser::NodesblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitNodesblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::NodesblockContext* OMDParser::nodesblock() {
  NodesblockContext *_localctx = _tracker.createInstance<NodesblockContext>(_ctx, getState());
  enterRule(_localctx, 54, OMDParser::RuleNodesblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(549);
    match(OMDParser::NODES);
    setState(554);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(550);
      match(OMDParser::LBRACKET);
      setState(551);
      intConst();
      setState(552);
      match(OMDParser::RBRACKET);
    }
    setState(556);
    match(OMDParser::LCURLY);
    setState(560);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(557);
      nodesstatement();
      setState(562);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(563);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- NodesstatementContext ------------------------------------------------------------------

OMDParser::NodesstatementContext::NodesstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::NodesstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::NodesstatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::NodesstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::NodesstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::NodesstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::NodesstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}


size_t OMDParser::NodesstatementContext::getRuleIndex() const {
  return OMDParser::RuleNodesstatement;
}

void OMDParser::NodesstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterNodesstatement(this);
}

void OMDParser::NodesstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitNodesstatement(this);
}


std::any OMDParser::NodesstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitNodesstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::NodesstatementContext* OMDParser::nodesstatement() {
  NodesstatementContext *_localctx = _tracker.createInstance<NodesstatementContext>(_ctx, getState());
  enterRule(_localctx, 56, OMDParser::RuleNodesstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(572);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(565);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(566);
        match(OMDParser::MEMBERS);
        setState(567);
        match(OMDParser::LPAREN);
        setState(568);
        inttuple();
        setState(569);
        match(OMDParser::RPAREN);
        setState(570);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- FragmentblockContext ------------------------------------------------------------------

OMDParser::FragmentblockContext::FragmentblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::FragmentblockContext::FRAGMENT() {
  return getToken(OMDParser::FRAGMENT, 0);
}

tree::TerminalNode* OMDParser::FragmentblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::FragmentblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::FragmentstatementContext *> OMDParser::FragmentblockContext::fragmentstatement() {
  return getRuleContexts<OMDParser::FragmentstatementContext>();
}

OMDParser::FragmentstatementContext* OMDParser::FragmentblockContext::fragmentstatement(size_t i) {
  return getRuleContext<OMDParser::FragmentstatementContext>(i);
}


size_t OMDParser::FragmentblockContext::getRuleIndex() const {
  return OMDParser::RuleFragmentblock;
}

void OMDParser::FragmentblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterFragmentblock(this);
}

void OMDParser::FragmentblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitFragmentblock(this);
}


std::any OMDParser::FragmentblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitFragmentblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::FragmentblockContext* OMDParser::fragmentblock() {
  FragmentblockContext *_localctx = _tracker.createInstance<FragmentblockContext>(_ctx, getState());
  enterRule(_localctx, 58, OMDParser::RuleFragmentblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(574);
    match(OMDParser::FRAGMENT);
    setState(575);
    match(OMDParser::LCURLY);
    setState(579);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 7) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 7)) & 144115462953763071) != 0)) {
      setState(576);
      fragmentstatement();
      setState(581);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(582);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- FragmentstatementContext ------------------------------------------------------------------

OMDParser::FragmentstatementContext::FragmentstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::FragmentstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

OMDParser::AtomblockContext* OMDParser::FragmentstatementContext::atomblock() {
  return getRuleContext<OMDParser::AtomblockContext>(0);
}

OMDParser::BondblockContext* OMDParser::FragmentstatementContext::bondblock() {
  return getRuleContext<OMDParser::BondblockContext>(0);
}

OMDParser::BendblockContext* OMDParser::FragmentstatementContext::bendblock() {
  return getRuleContext<OMDParser::BendblockContext>(0);
}

OMDParser::TorsionblockContext* OMDParser::FragmentstatementContext::torsionblock() {
  return getRuleContext<OMDParser::TorsionblockContext>(0);
}

OMDParser::InversionblockContext* OMDParser::FragmentstatementContext::inversionblock() {
  return getRuleContext<OMDParser::InversionblockContext>(0);
}

OMDParser::RigidbodyblockContext* OMDParser::FragmentstatementContext::rigidbodyblock() {
  return getRuleContext<OMDParser::RigidbodyblockContext>(0);
}

OMDParser::CutoffgroupblockContext* OMDParser::FragmentstatementContext::cutoffgroupblock() {
  return getRuleContext<OMDParser::CutoffgroupblockContext>(0);
}

OMDParser::ConstraintblockContext* OMDParser::FragmentstatementContext::constraintblock() {
  return getRuleContext<OMDParser::ConstraintblockContext>(0);
}

OMDParser::NodesblockContext* OMDParser::FragmentstatementContext::nodesblock() {
  return getRuleContext<OMDParser::NodesblockContext>(0);
}


size_t OMDParser::FragmentstatementContext::getRuleIndex() const {
  return OMDParser::RuleFragmentstatement;
}

void OMDParser::FragmentstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterFragmentstatement(this);
}

void OMDParser::FragmentstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitFragmentstatement(this);
}


std::any OMDParser::FragmentstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitFragmentstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::FragmentstatementContext* OMDParser::fragmentstatement() {
  FragmentstatementContext *_localctx = _tracker.createInstance<FragmentstatementContext>(_ctx, getState());
  enterRule(_localctx, 60, OMDParser::RuleFragmentstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(594);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(584);
        assignment();
        break;
      }

      case OMDParser::ATOM: {
        enterOuterAlt(_localctx, 2);
        setState(585);
        atomblock();
        break;
      }

      case OMDParser::BOND: {
        enterOuterAlt(_localctx, 3);
        setState(586);
        bondblock();
        break;
      }

      case OMDParser::BEND: {
        enterOuterAlt(_localctx, 4);
        setState(587);
        bendblock();
        break;
      }

      case OMDParser::TORSION: {
        enterOuterAlt(_localctx, 5);
        setState(588);
        torsionblock();
        break;
      }

      case OMDParser::INVERSION: {
        enterOuterAlt(_localctx, 6);
        setState(589);
        inversionblock();
        break;
      }

      case OMDParser::RIGIDBODY: {
        enterOuterAlt(_localctx, 7);
        setState(590);
        rigidbodyblock();
        break;
      }

      case OMDParser::CUTOFFGROUP: {
        enterOuterAlt(_localctx, 8);
        setState(591);
        cutoffgroupblock();
        break;
      }

      case OMDParser::CONSTRAINT: {
        enterOuterAlt(_localctx, 9);
        setState(592);
        constraintblock();
        break;
      }

      case OMDParser::NODES: {
        enterOuterAlt(_localctx, 10);
        setState(593);
        nodesblock();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ConstraintblockContext ------------------------------------------------------------------

OMDParser::ConstraintblockContext::ConstraintblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::ConstraintblockContext::CONSTRAINT() {
  return getToken(OMDParser::CONSTRAINT, 0);
}

tree::TerminalNode* OMDParser::ConstraintblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::ConstraintblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

tree::TerminalNode* OMDParser::ConstraintblockContext::LBRACKET() {
  return getToken(OMDParser::LBRACKET, 0);
}

OMDParser::IntConstContext* OMDParser::ConstraintblockContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

tree::TerminalNode* OMDParser::ConstraintblockContext::RBRACKET() {
  return getToken(OMDParser::RBRACKET, 0);
}

std::vector<OMDParser::ConstraintstatementContext *> OMDParser::ConstraintblockContext::constraintstatement() {
  return getRuleContexts<OMDParser::ConstraintstatementContext>();
}

OMDParser::ConstraintstatementContext* OMDParser::ConstraintblockContext::constraintstatement(size_t i) {
  return getRuleContext<OMDParser::ConstraintstatementContext>(i);
}


size_t OMDParser::ConstraintblockContext::getRuleIndex() const {
  return OMDParser::RuleConstraintblock;
}

void OMDParser::ConstraintblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterConstraintblock(this);
}

void OMDParser::ConstraintblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitConstraintblock(this);
}


std::any OMDParser::ConstraintblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitConstraintblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::ConstraintblockContext* OMDParser::constraintblock() {
  ConstraintblockContext *_localctx = _tracker.createInstance<ConstraintblockContext>(_ctx, getState());
  enterRule(_localctx, 62, OMDParser::RuleConstraintblock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(596);
    match(OMDParser::CONSTRAINT);
    setState(601);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(597);
      match(OMDParser::LBRACKET);
      setState(598);
      intConst();
      setState(599);
      match(OMDParser::RBRACKET);
    }
    setState(603);
    match(OMDParser::LCURLY);
    setState(607);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(604);
      constraintstatement();
      setState(609);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(610);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ConstraintstatementContext ------------------------------------------------------------------

OMDParser::ConstraintstatementContext::ConstraintstatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::AssignmentContext* OMDParser::ConstraintstatementContext::assignment() {
  return getRuleContext<OMDParser::AssignmentContext>(0);
}

tree::TerminalNode* OMDParser::ConstraintstatementContext::MEMBERS() {
  return getToken(OMDParser::MEMBERS, 0);
}

tree::TerminalNode* OMDParser::ConstraintstatementContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::InttupleContext* OMDParser::ConstraintstatementContext::inttuple() {
  return getRuleContext<OMDParser::InttupleContext>(0);
}

tree::TerminalNode* OMDParser::ConstraintstatementContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}

tree::TerminalNode* OMDParser::ConstraintstatementContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}


size_t OMDParser::ConstraintstatementContext::getRuleIndex() const {
  return OMDParser::RuleConstraintstatement;
}

void OMDParser::ConstraintstatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterConstraintstatement(this);
}

void OMDParser::ConstraintstatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitConstraintstatement(this);
}


std::any OMDParser::ConstraintstatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitConstraintstatement(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::ConstraintstatementContext* OMDParser::constraintstatement() {
  ConstraintstatementContext *_localctx = _tracker.createInstance<ConstraintstatementContext>(_ctx, getState());
  enterRule(_localctx, 64, OMDParser::RuleConstraintstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(619);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(612);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(613);
        match(OMDParser::MEMBERS);
        setState(614);
        match(OMDParser::LPAREN);
        setState(615);
        inttuple();
        setState(616);
        match(OMDParser::RPAREN);
        setState(617);
        match(OMDParser::SEMICOLON);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- SequencestringContext ------------------------------------------------------------------

OMDParser::SequencestringContext::SequencestringContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::SequencestringContext::SEQUENCE() {
  return getToken(OMDParser::SEQUENCE, 0);
}

tree::TerminalNode* OMDParser::SequencestringContext::ASSIGNEQUAL() {
  return getToken(OMDParser::ASSIGNEQUAL, 0);
}

OMDParser::ConstantContext* OMDParser::SequencestringContext::constant() {
  return getRuleContext<OMDParser::ConstantContext>(0);
}

tree::TerminalNode* OMDParser::SequencestringContext::SEMICOLON() {
  return getToken(OMDParser::SEMICOLON, 0);
}


size_t OMDParser::SequencestringContext::getRuleIndex() const {
  return OMDParser::RuleSequencestring;
}

void OMDParser::SequencestringContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterSequencestring(this);
}

void OMDParser::SequencestringContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitSequencestring(this);
}


std::any OMDParser::SequencestringContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitSequencestring(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::SequencestringContext* OMDParser::sequencestring() {
  SequencestringContext *_localctx = _tracker.createInstance<SequencestringContext>(_ctx, getState());
  enterRule(_localctx, 66, OMDParser::RuleSequencestring);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(621);
    match(OMDParser::SEQUENCE);
    setState(622);
    match(OMDParser::ASSIGNEQUAL);
    setState(623);
    constant();
    setState(624);
    match(OMDParser::SEMICOLON);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DoubleNumberTupleContext ------------------------------------------------------------------

OMDParser::DoubleNumberTupleContext::DoubleNumberTupleContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<OMDParser::DoubleNumberContext *> OMDParser::DoubleNumberTupleContext::doubleNumber() {
  return getRuleContexts<OMDParser::DoubleNumberContext>();
}

OMDParser::DoubleNumberContext* OMDParser::DoubleNumberTupleContext::doubleNumber(size_t i) {
  return getRuleContext<OMDParser::DoubleNumberContext>(i);
}

std::vector<tree::TerminalNode *> OMDParser::DoubleNumberTupleContext::COMMA() {
  return getTokens(OMDParser::COMMA);
}

tree::TerminalNode* OMDParser::DoubleNumberTupleContext::COMMA(size_t i) {
  return getToken(OMDParser::COMMA, i);
}


size_t OMDParser::DoubleNumberTupleContext::getRuleIndex() const {
  return OMDParser::RuleDoubleNumberTuple;
}

void OMDParser::DoubleNumberTupleContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterDoubleNumberTuple(this);
}

void OMDParser::DoubleNumberTupleContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitDoubleNumberTuple(this);
}


std::any OMDParser::DoubleNumberTupleContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitDoubleNumberTuple(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::DoubleNumberTupleContext* OMDParser::doubleNumberTuple() {
  DoubleNumberTupleContext *_localctx = _tracker.createInstance<DoubleNumberTupleContext>(_ctx, getState());
  enterRule(_localctx, 68, OMDParser::RuleDoubleNumberTuple);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(626);
    doubleNumber();
    setState(631);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::COMMA) {
      setState(627);
      match(OMDParser::COMMA);
      setState(628);
      doubleNumber();
      setState(633);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- InttupleContext ------------------------------------------------------------------

OMDParser::InttupleContext::InttupleContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<OMDParser::IntConstContext *> OMDParser::InttupleContext::intConst() {
  return getRuleContexts<OMDParser::IntConstContext>();
}

OMDParser::IntConstContext* OMDParser::InttupleContext::intConst(size_t i) {
  return getRuleContext<OMDParser::IntConstContext>(i);
}

std::vector<tree::TerminalNode *> OMDParser::InttupleContext::COMMA() {
  return getTokens(OMDParser::COMMA);
}

tree::TerminalNode* OMDParser::InttupleContext::COMMA(size_t i) {
  return getToken(OMDParser::COMMA, i);
}


size_t OMDParser::InttupleContext::getRuleIndex() const {
  return OMDParser::RuleInttuple;
}

void OMDParser::InttupleContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterInttuple(this);
}

void OMDParser::InttupleContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitInttuple(this);
}


std::any OMDParser::InttupleContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitInttuple(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::InttupleContext* OMDParser::inttuple() {
  InttupleContext *_localctx = _tracker.createInstance<InttupleContext>(_ctx, getState());
  enterRule(_localctx, 70, OMDParser::RuleInttuple);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(634);
    intConst();
    setState(639);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::COMMA) {
      setState(635);
      match(OMDParser::COMMA);
      setState(636);
      intConst();
      setState(641);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- IntConstContext ------------------------------------------------------------------

OMDParser::IntConstContext::IntConstContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::IntConstContext::NUM_INT() {
  return getToken(OMDParser::NUM_INT, 0);
}

tree::TerminalNode* OMDParser::IntConstContext::NUM_LONG() {
  return getToken(OMDParser::NUM_LONG, 0);
}


size_t OMDParser::IntConstContext::getRuleIndex() const {
  return OMDParser::RuleIntConst;
}

void OMDParser::IntConstContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterIntConst(this);
}

void OMDParser::IntConstContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitIntConst(this);
}


std::any OMDParser::IntConstContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitIntConst(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::IntConstContext* OMDParser::intConst() {
  IntConstContext *_localctx = _tracker.createInstance<IntConstContext>(_ctx, getState());
  enterRule(_localctx, 72, OMDParser::RuleIntConst);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(642);
    _la = _input->LA(1);
    if (!(_la == OMDParser::NUM_LONG

    || _la == OMDParser::NUM_INT)) {
    _errHandler->recoverInline(this);
    }
    else {
      _errHandler->reportMatch(this);
      consume();
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DoubleNumberContext ------------------------------------------------------------------

OMDParser::DoubleNumberContext::DoubleNumberContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

OMDParser::IntConstContext* OMDParser::DoubleNumberContext::intConst() {
  return getRuleContext<OMDParser::IntConstContext>(0);
}

OMDParser::FloatConstContext* OMDParser::DoubleNumberContext::floatConst() {
  return getRuleContext<OMDParser::FloatConstContext>(0);
}


size_t OMDParser::DoubleNumberContext::getRuleIndex() const {
  return OMDParser::RuleDoubleNumber;
}

void OMDParser::DoubleNumberContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterDoubleNumber(this);
}

void OMDParser::DoubleNumberContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitDoubleNumber(this);
}


std::any OMDParser::DoubleNumberContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitDoubleNumber(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::DoubleNumberContext* OMDParser::doubleNumber() {
  DoubleNumberContext *_localctx = _tracker.createInstance<DoubleNumberContext>(_ctx, getState());
  enterRule(_localctx, 74, OMDParser::RuleDoubleNumber);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(646);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::NUM_LONG:
      case OMDParser::NUM_INT: {
        enterOuterAlt(_localctx, 1);
        setState(644);
        intConst();
        break;
      }

      case OMDParser::NUM_FLOAT:
      case OMDParser::NUM_DOUBLE: {
        enterOuterAlt(_localctx, 2);
        setState(645);
        floatConst();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- FloatConstContext ------------------------------------------------------------------

OMDParser::FloatConstContext::FloatConstContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::FloatConstContext::NUM_FLOAT() {
  return getToken(OMDParser::NUM_FLOAT, 0);
}

tree::TerminalNode* OMDParser::FloatConstContext::NUM_DOUBLE() {
  return getToken(OMDParser::NUM_DOUBLE, 0);
}


size_t OMDParser::FloatConstContext::getRuleIndex() const {
  return OMDParser::RuleFloatConst;
}

void OMDParser::FloatConstContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterFloatConst(this);
}

void OMDParser::FloatConstContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitFloatConst(this);
}


std::any OMDParser::FloatConstContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitFloatConst(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::FloatConstContext* OMDParser::floatConst() {
  FloatConstContext *_localctx = _tracker.createInstance<FloatConstContext>(_ctx, getState());
  enterRule(_localctx, 76, OMDParser::RuleFloatConst);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(648);
    _la = _input->LA(1);
    if (!(_la == OMDParser::NUM_FLOAT

    || _la == OMDParser::NUM_DOUBLE)) {
    _errHandler->recoverInline(this);
    }
    else {
      _errHandler->reportMatch(this);
      consume();
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- VectorConstContext ------------------------------------------------------------------

OMDParser::VectorConstContext::VectorConstContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::VectorConstContext::LPAREN() {
  return getToken(OMDParser::LPAREN, 0);
}

OMDParser::DoubleNumberTupleContext* OMDParser::VectorConstContext::doubleNumberTuple() {
  return getRuleContext<OMDParser::DoubleNumberTupleContext>(0);
}

tree::TerminalNode* OMDParser::VectorConstContext::RPAREN() {
  return getToken(OMDParser::RPAREN, 0);
}


size_t OMDParser::VectorConstContext::getRuleIndex() const {
  return OMDParser::RuleVectorConst;
}

void OMDParser::VectorConstContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterVectorConst(this);
}

void OMDParser::VectorConstContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitVectorConst(this);
}


std::any OMDParser::VectorConstContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitVectorConst(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::VectorConstContext* OMDParser::vectorConst() {
  VectorConstContext *_localctx = _tracker.createInstance<VectorConstContext>(_ctx, getState());
  enterRule(_localctx, 78, OMDParser::RuleVectorConst);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(650);
    match(OMDParser::LPAREN);
    setState(651);
    doubleNumberTuple();
    setState(652);
    match(OMDParser::RPAREN);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

void OMDParser::initialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  omdParserInitialize();
#else
  ::antlr4::internal::call_once(omdParserOnceFlag, omdParserInitialize);
#endif
}
