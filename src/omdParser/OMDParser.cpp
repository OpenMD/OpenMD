
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
      "", "'component'", "'molecule'", "'zconstraint'", "'restraint'", "'atom'", 
      "'bond'", "'bend'", "'torsion'", "'inversion'", "'rigidBody'", "'cutoffGroup'", 
      "'constraint'", "'distance'", "'fragment'", "'sequence'", "'members'", 
      "'center'", "'satellites'", "'position'", "'orientation'", "'flucQ'", 
      "'RNEMD'", "'light'", "'minimizer'", "'Fixed'", "'Harmonic'", "'Cubic'", 
      "'Quartic'", "'Polynomial'", "'Morse'", "'GhostBend'", "'UreyBradley'", 
      "'Cosine'", "'GhostTorsion'", "'Charmm'", "'Opls'", "'Trappe'", "'AmberImproper'", 
      "'ImproperCosine'", "'CentralAtomHeight'", "'Dreiding'", "'charge'", 
      "'nodes'", "'='", "':'", "','", "'\\u003F'", "';'", "'.'", "'('", 
      "')'", "'['", "']'", "'{'", "'}'"
    },
    std::vector<std::string>{
      "", "COMPONENT", "MOLECULE", "ZCONSTRAINT", "RESTRAINT", "ATOM", "BOND", 
      "BEND", "TORSION", "INVERSION", "RIGIDBODY", "CUTOFFGROUP", "CONSTRAINT", 
      "DISTANCE", "FRAGMENT", "SEQUENCE", "MEMBERS", "CENTER", "SATELLITES", 
      "POSITION", "ORIENTATION", "FLUCQ", "RNEMD", "LIGHT", "MINIMIZER", 
      "FIXED", "HARMONIC", "CUBIC", "QUARTIC", "POLYNOMIAL", "MORSE", "GHOSTBEND", 
      "UREYBRADLEY", "COSINE", "GHOSTTORSION", "CHARMM", "OPLS", "TRAPPE", 
      "AMBERIMPROPER", "IMPROPERCOSINE", "CENTRALATOMHEIGHT", "DREIDING", 
      "CHARGE", "NODES", "ASSIGNEQUAL", "COLON", "COMMA", "QUESTIONMARK", 
      "SEMICOLON", "DOT", "LPAREN", "RPAREN", "LBRACKET", "RBRACKET", "LCURLY", 
      "RCURLY", "ID", "NUM_LONG", "NUM_INT", "NUM_FLOAT", "NUM_DOUBLE", 
      "CharLiteral", "StringLiteral", "Whitespace", "Newline", "LineContinuation", 
      "Comment", "CPPComment", "PREPROC_DIRECTIVE"
    }
  );
  static const int32_t serializedATNSegment[] = {
  	4,1,68,653,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,6,2,
  	7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,13,7,13,2,14,7,
  	14,2,15,7,15,2,16,7,16,2,17,7,17,2,18,7,18,2,19,7,19,2,20,7,20,2,21,7,
  	21,2,22,7,22,2,23,7,23,2,24,7,24,2,25,7,25,2,26,7,26,2,27,7,27,2,28,7,
  	28,2,29,7,29,2,30,7,30,2,31,7,31,2,32,7,32,2,33,7,33,2,34,7,34,2,35,7,
  	35,2,36,7,36,2,37,7,37,2,38,7,38,2,39,7,39,1,0,5,0,82,8,0,10,0,12,0,85,
  	9,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,99,8,1,1,2,1,
  	2,1,2,1,2,1,2,1,3,1,3,1,3,1,3,1,3,3,3,111,8,3,1,4,1,4,1,4,5,4,116,8,4,
  	10,4,12,4,119,9,4,1,4,1,4,1,5,1,5,1,5,5,5,126,8,5,10,5,12,5,129,9,5,1,
  	5,1,5,1,6,1,6,1,6,5,6,136,8,6,10,6,12,6,139,9,6,1,6,1,6,1,7,1,7,1,7,5,
  	7,146,8,7,10,7,12,7,149,9,7,1,7,1,7,1,8,1,8,1,8,5,8,156,8,8,10,8,12,8,
  	159,9,8,1,8,1,8,1,9,1,9,1,9,5,9,166,8,9,10,9,12,9,169,9,9,1,9,1,9,1,10,
  	1,10,1,10,5,10,176,8,10,10,10,12,10,179,9,10,1,10,1,10,1,11,1,11,1,11,
  	5,11,186,8,11,10,11,12,11,189,9,11,1,11,1,11,1,12,1,12,1,12,1,12,1,12,
  	1,12,1,12,1,12,1,12,1,12,3,12,203,8,12,1,13,1,13,1,13,1,13,1,13,1,13,
  	5,13,211,8,13,10,13,12,13,214,9,13,1,13,1,13,1,14,1,14,1,14,1,14,1,14,
  	1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,1,14,
  	3,14,237,8,14,1,15,1,15,1,15,1,15,1,15,3,15,244,8,15,1,15,1,15,5,15,248,
  	8,15,10,15,12,15,251,9,15,1,15,1,15,1,16,1,16,1,16,1,16,1,16,1,16,1,16,
  	1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,
  	1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,
  	1,16,1,16,1,16,1,16,1,16,1,16,1,16,1,16,3,16,298,8,16,1,17,1,17,1,17,
  	1,17,1,17,3,17,305,8,17,1,17,1,17,5,17,309,8,17,10,17,12,17,312,9,17,
  	1,17,1,17,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,
  	1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,
  	1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,
  	1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,1,18,3,18,365,8,18,1,19,1,19,
  	1,19,1,19,1,19,3,19,372,8,19,1,19,1,19,5,19,376,8,19,10,19,12,19,379,
  	9,19,1,19,1,19,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,
  	1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,
  	1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,
  	1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,1,20,
  	1,20,1,20,3,20,438,8,20,1,21,1,21,1,21,1,21,1,21,3,21,445,8,21,1,21,1,
  	21,5,21,449,8,21,10,21,12,21,452,9,21,1,21,1,21,1,22,1,22,1,22,1,22,1,
  	22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,
  	22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,
  	22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,3,22,499,8,22,1,
  	23,1,23,1,23,1,23,1,23,1,23,5,23,507,8,23,10,23,12,23,510,9,23,1,23,1,
  	23,1,24,1,24,1,24,1,24,1,24,1,24,1,24,3,24,521,8,24,1,25,1,25,1,25,1,
  	25,1,25,3,25,528,8,25,1,25,1,25,5,25,532,8,25,10,25,12,25,535,9,25,1,
  	25,1,25,1,26,1,26,1,26,1,26,1,26,1,26,1,26,3,26,546,8,26,1,27,1,27,1,
  	27,1,27,1,27,3,27,553,8,27,1,27,1,27,5,27,557,8,27,10,27,12,27,560,9,
  	27,1,27,1,27,1,28,1,28,1,28,1,28,1,28,1,28,1,28,3,28,571,8,28,1,29,1,
  	29,1,29,5,29,576,8,29,10,29,12,29,579,9,29,1,29,1,29,1,30,1,30,1,30,1,
  	30,1,30,1,30,1,30,1,30,1,30,1,30,3,30,593,8,30,1,31,1,31,1,31,1,31,1,
  	31,3,31,600,8,31,1,31,1,31,5,31,604,8,31,10,31,12,31,607,9,31,1,31,1,
  	31,1,32,1,32,1,32,1,32,1,32,1,32,1,32,3,32,618,8,32,1,33,1,33,1,33,1,
  	33,1,33,1,34,1,34,1,34,5,34,628,8,34,10,34,12,34,631,9,34,1,35,1,35,1,
  	35,5,35,636,8,35,10,35,12,35,639,9,35,1,36,1,36,1,37,1,37,3,37,645,8,
  	37,1,38,1,38,1,39,1,39,1,39,1,39,1,39,0,0,40,0,2,4,6,8,10,12,14,16,18,
  	20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,
  	66,68,70,72,74,76,78,0,2,1,0,57,58,1,0,59,60,710,0,83,1,0,0,0,2,98,1,
  	0,0,0,4,100,1,0,0,0,6,110,1,0,0,0,8,112,1,0,0,0,10,122,1,0,0,0,12,132,
  	1,0,0,0,14,142,1,0,0,0,16,152,1,0,0,0,18,162,1,0,0,0,20,172,1,0,0,0,22,
  	182,1,0,0,0,24,202,1,0,0,0,26,204,1,0,0,0,28,236,1,0,0,0,30,238,1,0,0,
  	0,32,297,1,0,0,0,34,299,1,0,0,0,36,364,1,0,0,0,38,366,1,0,0,0,40,437,
  	1,0,0,0,42,439,1,0,0,0,44,498,1,0,0,0,46,500,1,0,0,0,48,520,1,0,0,0,50,
  	522,1,0,0,0,52,545,1,0,0,0,54,547,1,0,0,0,56,570,1,0,0,0,58,572,1,0,0,
  	0,60,592,1,0,0,0,62,594,1,0,0,0,64,617,1,0,0,0,66,619,1,0,0,0,68,624,
  	1,0,0,0,70,632,1,0,0,0,72,640,1,0,0,0,74,644,1,0,0,0,76,646,1,0,0,0,78,
  	648,1,0,0,0,80,82,3,2,1,0,81,80,1,0,0,0,82,85,1,0,0,0,83,81,1,0,0,0,83,
  	84,1,0,0,0,84,86,1,0,0,0,85,83,1,0,0,0,86,87,5,0,0,1,87,1,1,0,0,0,88,
  	99,3,4,2,0,89,99,3,8,4,0,90,99,3,22,11,0,91,99,3,58,29,0,92,99,3,10,5,
  	0,93,99,3,12,6,0,94,99,3,14,7,0,95,99,3,16,8,0,96,99,3,18,9,0,97,99,3,
  	20,10,0,98,88,1,0,0,0,98,89,1,0,0,0,98,90,1,0,0,0,98,91,1,0,0,0,98,92,
  	1,0,0,0,98,93,1,0,0,0,98,94,1,0,0,0,98,95,1,0,0,0,98,96,1,0,0,0,98,97,
  	1,0,0,0,99,3,1,0,0,0,100,101,5,56,0,0,101,102,5,44,0,0,102,103,3,6,3,
  	0,103,104,5,48,0,0,104,5,1,0,0,0,105,111,3,72,36,0,106,111,3,76,38,0,
  	107,111,3,78,39,0,108,111,5,56,0,0,109,111,5,62,0,0,110,105,1,0,0,0,110,
  	106,1,0,0,0,110,107,1,0,0,0,110,108,1,0,0,0,110,109,1,0,0,0,111,7,1,0,
  	0,0,112,113,5,1,0,0,113,117,5,54,0,0,114,116,3,4,2,0,115,114,1,0,0,0,
  	116,119,1,0,0,0,117,115,1,0,0,0,117,118,1,0,0,0,118,120,1,0,0,0,119,117,
  	1,0,0,0,120,121,5,55,0,0,121,9,1,0,0,0,122,123,5,3,0,0,123,127,5,54,0,
  	0,124,126,3,4,2,0,125,124,1,0,0,0,126,129,1,0,0,0,127,125,1,0,0,0,127,
  	128,1,0,0,0,128,130,1,0,0,0,129,127,1,0,0,0,130,131,5,55,0,0,131,11,1,
  	0,0,0,132,133,5,4,0,0,133,137,5,54,0,0,134,136,3,4,2,0,135,134,1,0,0,
  	0,136,139,1,0,0,0,137,135,1,0,0,0,137,138,1,0,0,0,138,140,1,0,0,0,139,
  	137,1,0,0,0,140,141,5,55,0,0,141,13,1,0,0,0,142,143,5,21,0,0,143,147,
  	5,54,0,0,144,146,3,4,2,0,145,144,1,0,0,0,146,149,1,0,0,0,147,145,1,0,
  	0,0,147,148,1,0,0,0,148,150,1,0,0,0,149,147,1,0,0,0,150,151,5,55,0,0,
  	151,15,1,0,0,0,152,153,5,22,0,0,153,157,5,54,0,0,154,156,3,4,2,0,155,
  	154,1,0,0,0,156,159,1,0,0,0,157,155,1,0,0,0,157,158,1,0,0,0,158,160,1,
  	0,0,0,159,157,1,0,0,0,160,161,5,55,0,0,161,17,1,0,0,0,162,163,5,23,0,
  	0,163,167,5,54,0,0,164,166,3,4,2,0,165,164,1,0,0,0,166,169,1,0,0,0,167,
  	165,1,0,0,0,167,168,1,0,0,0,168,170,1,0,0,0,169,167,1,0,0,0,170,171,5,
  	55,0,0,171,19,1,0,0,0,172,173,5,24,0,0,173,177,5,54,0,0,174,176,3,4,2,
  	0,175,174,1,0,0,0,176,179,1,0,0,0,177,175,1,0,0,0,177,178,1,0,0,0,178,
  	180,1,0,0,0,179,177,1,0,0,0,180,181,5,55,0,0,181,21,1,0,0,0,182,183,5,
  	2,0,0,183,187,5,54,0,0,184,186,3,24,12,0,185,184,1,0,0,0,186,189,1,0,
  	0,0,187,185,1,0,0,0,187,188,1,0,0,0,188,190,1,0,0,0,189,187,1,0,0,0,190,
  	191,5,55,0,0,191,23,1,0,0,0,192,203,3,4,2,0,193,203,3,26,13,0,194,203,
  	3,30,15,0,195,203,3,34,17,0,196,203,3,38,19,0,197,203,3,42,21,0,198,203,
  	3,46,23,0,199,203,3,50,25,0,200,203,3,62,31,0,201,203,3,66,33,0,202,192,
  	1,0,0,0,202,193,1,0,0,0,202,194,1,0,0,0,202,195,1,0,0,0,202,196,1,0,0,
  	0,202,197,1,0,0,0,202,198,1,0,0,0,202,199,1,0,0,0,202,200,1,0,0,0,202,
  	201,1,0,0,0,203,25,1,0,0,0,204,205,5,5,0,0,205,206,5,52,0,0,206,207,3,
  	72,36,0,207,208,5,53,0,0,208,212,5,54,0,0,209,211,3,28,14,0,210,209,1,
  	0,0,0,211,214,1,0,0,0,212,210,1,0,0,0,212,213,1,0,0,0,213,215,1,0,0,0,
  	214,212,1,0,0,0,215,216,5,55,0,0,216,27,1,0,0,0,217,237,3,4,2,0,218,219,
  	5,19,0,0,219,220,5,50,0,0,220,221,3,68,34,0,221,222,5,51,0,0,222,223,
  	5,48,0,0,223,237,1,0,0,0,224,225,5,20,0,0,225,226,5,50,0,0,226,227,3,
  	68,34,0,227,228,5,51,0,0,228,229,5,48,0,0,229,237,1,0,0,0,230,231,5,42,
  	0,0,231,232,5,50,0,0,232,233,3,76,38,0,233,234,5,51,0,0,234,235,5,48,
  	0,0,235,237,1,0,0,0,236,217,1,0,0,0,236,218,1,0,0,0,236,224,1,0,0,0,236,
  	230,1,0,0,0,237,29,1,0,0,0,238,243,5,6,0,0,239,240,5,52,0,0,240,241,3,
  	72,36,0,241,242,5,53,0,0,242,244,1,0,0,0,243,239,1,0,0,0,243,244,1,0,
  	0,0,244,245,1,0,0,0,245,249,5,54,0,0,246,248,3,32,16,0,247,246,1,0,0,
  	0,248,251,1,0,0,0,249,247,1,0,0,0,249,250,1,0,0,0,250,252,1,0,0,0,251,
  	249,1,0,0,0,252,253,5,55,0,0,253,31,1,0,0,0,254,298,3,4,2,0,255,256,5,
  	16,0,0,256,257,5,50,0,0,257,258,3,70,35,0,258,259,5,51,0,0,259,260,5,
  	48,0,0,260,298,1,0,0,0,261,262,5,25,0,0,262,263,5,50,0,0,263,264,3,76,
  	38,0,264,265,5,51,0,0,265,266,5,48,0,0,266,298,1,0,0,0,267,268,5,26,0,
  	0,268,269,5,50,0,0,269,270,3,68,34,0,270,271,5,51,0,0,271,272,5,48,0,
  	0,272,298,1,0,0,0,273,274,5,27,0,0,274,275,5,50,0,0,275,276,3,68,34,0,
  	276,277,5,51,0,0,277,278,5,48,0,0,278,298,1,0,0,0,279,280,5,28,0,0,280,
  	281,5,50,0,0,281,282,3,68,34,0,282,283,5,51,0,0,283,284,5,48,0,0,284,
  	298,1,0,0,0,285,286,5,29,0,0,286,287,5,50,0,0,287,288,3,68,34,0,288,289,
  	5,51,0,0,289,290,5,48,0,0,290,298,1,0,0,0,291,292,5,30,0,0,292,293,5,
  	50,0,0,293,294,3,68,34,0,294,295,5,51,0,0,295,296,5,48,0,0,296,298,1,
  	0,0,0,297,254,1,0,0,0,297,255,1,0,0,0,297,261,1,0,0,0,297,267,1,0,0,0,
  	297,273,1,0,0,0,297,279,1,0,0,0,297,285,1,0,0,0,297,291,1,0,0,0,298,33,
  	1,0,0,0,299,304,5,7,0,0,300,301,5,52,0,0,301,302,3,72,36,0,302,303,5,
  	53,0,0,303,305,1,0,0,0,304,300,1,0,0,0,304,305,1,0,0,0,305,306,1,0,0,
  	0,306,310,5,54,0,0,307,309,3,36,18,0,308,307,1,0,0,0,309,312,1,0,0,0,
  	310,308,1,0,0,0,310,311,1,0,0,0,311,313,1,0,0,0,312,310,1,0,0,0,313,314,
  	5,55,0,0,314,35,1,0,0,0,315,365,3,4,2,0,316,317,5,16,0,0,317,318,5,50,
  	0,0,318,319,3,70,35,0,319,320,5,51,0,0,320,321,5,48,0,0,321,365,1,0,0,
  	0,322,323,5,26,0,0,323,324,5,50,0,0,324,325,3,68,34,0,325,326,5,51,0,
  	0,326,327,5,48,0,0,327,365,1,0,0,0,328,329,5,31,0,0,329,330,5,50,0,0,
  	330,331,3,68,34,0,331,332,5,51,0,0,332,333,5,48,0,0,333,365,1,0,0,0,334,
  	335,5,32,0,0,335,336,5,50,0,0,336,337,3,68,34,0,337,338,5,51,0,0,338,
  	339,5,48,0,0,339,365,1,0,0,0,340,341,5,27,0,0,341,342,5,50,0,0,342,343,
  	3,68,34,0,343,344,5,51,0,0,344,345,5,48,0,0,345,365,1,0,0,0,346,347,5,
  	28,0,0,347,348,5,50,0,0,348,349,3,68,34,0,349,350,5,51,0,0,350,351,5,
  	48,0,0,351,365,1,0,0,0,352,353,5,29,0,0,353,354,5,50,0,0,354,355,3,68,
  	34,0,355,356,5,51,0,0,356,357,5,48,0,0,357,365,1,0,0,0,358,359,5,33,0,
  	0,359,360,5,50,0,0,360,361,3,68,34,0,361,362,5,51,0,0,362,363,5,48,0,
  	0,363,365,1,0,0,0,364,315,1,0,0,0,364,316,1,0,0,0,364,322,1,0,0,0,364,
  	328,1,0,0,0,364,334,1,0,0,0,364,340,1,0,0,0,364,346,1,0,0,0,364,352,1,
  	0,0,0,364,358,1,0,0,0,365,37,1,0,0,0,366,371,5,8,0,0,367,368,5,52,0,0,
  	368,369,3,72,36,0,369,370,5,53,0,0,370,372,1,0,0,0,371,367,1,0,0,0,371,
  	372,1,0,0,0,372,373,1,0,0,0,373,377,5,54,0,0,374,376,3,40,20,0,375,374,
  	1,0,0,0,376,379,1,0,0,0,377,375,1,0,0,0,377,378,1,0,0,0,378,380,1,0,0,
  	0,379,377,1,0,0,0,380,381,5,55,0,0,381,39,1,0,0,0,382,438,3,4,2,0,383,
  	384,5,16,0,0,384,385,5,50,0,0,385,386,3,70,35,0,386,387,5,51,0,0,387,
  	388,5,48,0,0,388,438,1,0,0,0,389,390,5,34,0,0,390,391,5,50,0,0,391,392,
  	3,68,34,0,392,393,5,51,0,0,393,394,5,48,0,0,394,438,1,0,0,0,395,396,5,
  	27,0,0,396,397,5,50,0,0,397,398,3,68,34,0,398,399,5,51,0,0,399,400,5,
  	48,0,0,400,438,1,0,0,0,401,402,5,28,0,0,402,403,5,50,0,0,403,404,3,68,
  	34,0,404,405,5,51,0,0,405,406,5,48,0,0,406,438,1,0,0,0,407,408,5,29,0,
  	0,408,409,5,50,0,0,409,410,3,68,34,0,410,411,5,51,0,0,411,412,5,48,0,
  	0,412,438,1,0,0,0,413,414,5,35,0,0,414,415,5,50,0,0,415,416,3,68,34,0,
  	416,417,5,51,0,0,417,418,5,48,0,0,418,438,1,0,0,0,419,420,5,36,0,0,420,
  	421,5,50,0,0,421,422,3,68,34,0,422,423,5,51,0,0,423,424,5,48,0,0,424,
  	438,1,0,0,0,425,426,5,37,0,0,426,427,5,50,0,0,427,428,3,68,34,0,428,429,
  	5,51,0,0,429,430,5,48,0,0,430,438,1,0,0,0,431,432,5,26,0,0,432,433,5,
  	50,0,0,433,434,3,68,34,0,434,435,5,51,0,0,435,436,5,48,0,0,436,438,1,
  	0,0,0,437,382,1,0,0,0,437,383,1,0,0,0,437,389,1,0,0,0,437,395,1,0,0,0,
  	437,401,1,0,0,0,437,407,1,0,0,0,437,413,1,0,0,0,437,419,1,0,0,0,437,425,
  	1,0,0,0,437,431,1,0,0,0,438,41,1,0,0,0,439,444,5,9,0,0,440,441,5,52,0,
  	0,441,442,3,72,36,0,442,443,5,53,0,0,443,445,1,0,0,0,444,440,1,0,0,0,
  	444,445,1,0,0,0,445,446,1,0,0,0,446,450,5,54,0,0,447,449,3,44,22,0,448,
  	447,1,0,0,0,449,452,1,0,0,0,450,448,1,0,0,0,450,451,1,0,0,0,451,453,1,
  	0,0,0,452,450,1,0,0,0,453,454,5,55,0,0,454,43,1,0,0,0,455,499,3,4,2,0,
  	456,457,5,17,0,0,457,458,5,50,0,0,458,459,3,72,36,0,459,460,5,51,0,0,
  	460,461,5,48,0,0,461,499,1,0,0,0,462,463,5,18,0,0,463,464,5,50,0,0,464,
  	465,3,70,35,0,465,466,5,51,0,0,466,467,5,48,0,0,467,499,1,0,0,0,468,469,
  	5,38,0,0,469,470,5,50,0,0,470,471,3,68,34,0,471,472,5,51,0,0,472,473,
  	5,48,0,0,473,499,1,0,0,0,474,475,5,39,0,0,475,476,5,50,0,0,476,477,3,
  	68,34,0,477,478,5,51,0,0,478,479,5,48,0,0,479,499,1,0,0,0,480,481,5,26,
  	0,0,481,482,5,50,0,0,482,483,3,68,34,0,483,484,5,51,0,0,484,485,5,48,
  	0,0,485,499,1,0,0,0,486,487,5,40,0,0,487,488,5,50,0,0,488,489,3,68,34,
  	0,489,490,5,51,0,0,490,491,5,48,0,0,491,499,1,0,0,0,492,493,5,41,0,0,
  	493,494,5,50,0,0,494,495,3,68,34,0,495,496,5,51,0,0,496,497,5,48,0,0,
  	497,499,1,0,0,0,498,455,1,0,0,0,498,456,1,0,0,0,498,462,1,0,0,0,498,468,
  	1,0,0,0,498,474,1,0,0,0,498,480,1,0,0,0,498,486,1,0,0,0,498,492,1,0,0,
  	0,499,45,1,0,0,0,500,501,5,10,0,0,501,502,5,52,0,0,502,503,3,72,36,0,
  	503,504,5,53,0,0,504,508,5,54,0,0,505,507,3,48,24,0,506,505,1,0,0,0,507,
  	510,1,0,0,0,508,506,1,0,0,0,508,509,1,0,0,0,509,511,1,0,0,0,510,508,1,
  	0,0,0,511,512,5,55,0,0,512,47,1,0,0,0,513,521,3,4,2,0,514,515,5,16,0,
  	0,515,516,5,50,0,0,516,517,3,70,35,0,517,518,5,51,0,0,518,519,5,48,0,
  	0,519,521,1,0,0,0,520,513,1,0,0,0,520,514,1,0,0,0,521,49,1,0,0,0,522,
  	527,5,11,0,0,523,524,5,52,0,0,524,525,3,72,36,0,525,526,5,53,0,0,526,
  	528,1,0,0,0,527,523,1,0,0,0,527,528,1,0,0,0,528,529,1,0,0,0,529,533,5,
  	54,0,0,530,532,3,52,26,0,531,530,1,0,0,0,532,535,1,0,0,0,533,531,1,0,
  	0,0,533,534,1,0,0,0,534,536,1,0,0,0,535,533,1,0,0,0,536,537,5,55,0,0,
  	537,51,1,0,0,0,538,546,3,4,2,0,539,540,5,16,0,0,540,541,5,50,0,0,541,
  	542,3,70,35,0,542,543,5,51,0,0,543,544,5,48,0,0,544,546,1,0,0,0,545,538,
  	1,0,0,0,545,539,1,0,0,0,546,53,1,0,0,0,547,552,5,43,0,0,548,549,5,52,
  	0,0,549,550,3,72,36,0,550,551,5,53,0,0,551,553,1,0,0,0,552,548,1,0,0,
  	0,552,553,1,0,0,0,553,554,1,0,0,0,554,558,5,54,0,0,555,557,3,56,28,0,
  	556,555,1,0,0,0,557,560,1,0,0,0,558,556,1,0,0,0,558,559,1,0,0,0,559,561,
  	1,0,0,0,560,558,1,0,0,0,561,562,5,55,0,0,562,55,1,0,0,0,563,571,3,4,2,
  	0,564,565,5,16,0,0,565,566,5,50,0,0,566,567,3,70,35,0,567,568,5,51,0,
  	0,568,569,5,48,0,0,569,571,1,0,0,0,570,563,1,0,0,0,570,564,1,0,0,0,571,
  	57,1,0,0,0,572,573,5,14,0,0,573,577,5,54,0,0,574,576,3,60,30,0,575,574,
  	1,0,0,0,576,579,1,0,0,0,577,575,1,0,0,0,577,578,1,0,0,0,578,580,1,0,0,
  	0,579,577,1,0,0,0,580,581,5,55,0,0,581,59,1,0,0,0,582,593,3,4,2,0,583,
  	593,3,26,13,0,584,593,3,30,15,0,585,593,3,34,17,0,586,593,3,38,19,0,587,
  	593,3,42,21,0,588,593,3,46,23,0,589,593,3,50,25,0,590,593,3,62,31,0,591,
  	593,3,54,27,0,592,582,1,0,0,0,592,583,1,0,0,0,592,584,1,0,0,0,592,585,
  	1,0,0,0,592,586,1,0,0,0,592,587,1,0,0,0,592,588,1,0,0,0,592,589,1,0,0,
  	0,592,590,1,0,0,0,592,591,1,0,0,0,593,61,1,0,0,0,594,599,5,12,0,0,595,
  	596,5,52,0,0,596,597,3,72,36,0,597,598,5,53,0,0,598,600,1,0,0,0,599,595,
  	1,0,0,0,599,600,1,0,0,0,600,601,1,0,0,0,601,605,5,54,0,0,602,604,3,64,
  	32,0,603,602,1,0,0,0,604,607,1,0,0,0,605,603,1,0,0,0,605,606,1,0,0,0,
  	606,608,1,0,0,0,607,605,1,0,0,0,608,609,5,55,0,0,609,63,1,0,0,0,610,618,
  	3,4,2,0,611,612,5,16,0,0,612,613,5,50,0,0,613,614,3,70,35,0,614,615,5,
  	51,0,0,615,616,5,48,0,0,616,618,1,0,0,0,617,610,1,0,0,0,617,611,1,0,0,
  	0,618,65,1,0,0,0,619,620,5,15,0,0,620,621,5,44,0,0,621,622,3,6,3,0,622,
  	623,5,48,0,0,623,67,1,0,0,0,624,629,3,74,37,0,625,626,5,46,0,0,626,628,
  	3,74,37,0,627,625,1,0,0,0,628,631,1,0,0,0,629,627,1,0,0,0,629,630,1,0,
  	0,0,630,69,1,0,0,0,631,629,1,0,0,0,632,637,3,72,36,0,633,634,5,46,0,0,
  	634,636,3,72,36,0,635,633,1,0,0,0,636,639,1,0,0,0,637,635,1,0,0,0,637,
  	638,1,0,0,0,638,71,1,0,0,0,639,637,1,0,0,0,640,641,7,0,0,0,641,73,1,0,
  	0,0,642,645,3,72,36,0,643,645,3,76,38,0,644,642,1,0,0,0,644,643,1,0,0,
  	0,645,75,1,0,0,0,646,647,7,1,0,0,647,77,1,0,0,0,648,649,5,50,0,0,649,
  	650,3,68,34,0,650,651,5,51,0,0,651,79,1,0,0,0,42,83,98,110,117,127,137,
  	147,157,167,177,187,202,212,236,243,249,297,304,310,364,371,377,437,444,
  	450,498,508,520,527,533,545,552,558,570,577,592,599,605,617,629,637,644
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
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72057594069401630) != 0)) {
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
    setState(110);
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

      case OMDParser::ID: {
        enterOuterAlt(_localctx, 4);
        setState(108);
        match(OMDParser::ID);
        break;
      }

      case OMDParser::StringLiteral: {
        enterOuterAlt(_localctx, 5);
        setState(109);
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
    setState(112);
    match(OMDParser::COMPONENT);
    setState(113);
    match(OMDParser::LCURLY);
    setState(117);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(114);
      assignment();
      setState(119);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(120);
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
    setState(122);
    match(OMDParser::ZCONSTRAINT);
    setState(123);
    match(OMDParser::LCURLY);
    setState(127);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(124);
      assignment();
      setState(129);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(130);
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
    setState(132);
    match(OMDParser::RESTRAINT);
    setState(133);
    match(OMDParser::LCURLY);
    setState(137);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(134);
      assignment();
      setState(139);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(140);
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
    setState(142);
    match(OMDParser::FLUCQ);
    setState(143);
    match(OMDParser::LCURLY);
    setState(147);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(144);
      assignment();
      setState(149);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(150);
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
    setState(152);
    match(OMDParser::RNEMD);
    setState(153);
    match(OMDParser::LCURLY);
    setState(157);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(154);
      assignment();
      setState(159);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(160);
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
    setState(162);
    match(OMDParser::LIGHT);
    setState(163);
    match(OMDParser::LCURLY);
    setState(167);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(164);
      assignment();
      setState(169);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(170);
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
    setState(172);
    match(OMDParser::MINIMIZER);
    setState(173);
    match(OMDParser::LCURLY);
    setState(177);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(174);
      assignment();
      setState(179);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(180);
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
    setState(182);
    match(OMDParser::MOLECULE);
    setState(183);
    match(OMDParser::LCURLY);
    setState(187);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72057594037968864) != 0)) {
      setState(184);
      moleculestatement();
      setState(189);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(190);
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
    setState(202);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(192);
        assignment();
        break;
      }

      case OMDParser::ATOM: {
        enterOuterAlt(_localctx, 2);
        setState(193);
        atomblock();
        break;
      }

      case OMDParser::BOND: {
        enterOuterAlt(_localctx, 3);
        setState(194);
        bondblock();
        break;
      }

      case OMDParser::BEND: {
        enterOuterAlt(_localctx, 4);
        setState(195);
        bendblock();
        break;
      }

      case OMDParser::TORSION: {
        enterOuterAlt(_localctx, 5);
        setState(196);
        torsionblock();
        break;
      }

      case OMDParser::INVERSION: {
        enterOuterAlt(_localctx, 6);
        setState(197);
        inversionblock();
        break;
      }

      case OMDParser::RIGIDBODY: {
        enterOuterAlt(_localctx, 7);
        setState(198);
        rigidbodyblock();
        break;
      }

      case OMDParser::CUTOFFGROUP: {
        enterOuterAlt(_localctx, 8);
        setState(199);
        cutoffgroupblock();
        break;
      }

      case OMDParser::CONSTRAINT: {
        enterOuterAlt(_localctx, 9);
        setState(200);
        constraintblock();
        break;
      }

      case OMDParser::SEQUENCE: {
        enterOuterAlt(_localctx, 10);
        setState(201);
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
    setState(204);
    match(OMDParser::ATOM);
    setState(205);
    match(OMDParser::LBRACKET);
    setState(206);
    intConst();
    setState(207);
    match(OMDParser::RBRACKET);
    setState(208);
    match(OMDParser::LCURLY);
    setState(212);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72061992086011904) != 0)) {
      setState(209);
      atomstatement();
      setState(214);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(215);
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
    setState(236);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(217);
        assignment();
        break;
      }

      case OMDParser::POSITION: {
        enterOuterAlt(_localctx, 2);
        setState(218);
        match(OMDParser::POSITION);
        setState(219);
        match(OMDParser::LPAREN);
        setState(220);
        doubleNumberTuple();
        setState(221);
        match(OMDParser::RPAREN);
        setState(222);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::ORIENTATION: {
        enterOuterAlt(_localctx, 3);
        setState(224);
        match(OMDParser::ORIENTATION);
        setState(225);
        match(OMDParser::LPAREN);
        setState(226);
        doubleNumberTuple();
        setState(227);
        match(OMDParser::RPAREN);
        setState(228);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CHARGE: {
        enterOuterAlt(_localctx, 4);
        setState(230);
        match(OMDParser::CHARGE);
        setState(231);
        match(OMDParser::LPAREN);
        setState(232);
        floatConst();
        setState(233);
        match(OMDParser::RPAREN);
        setState(234);
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
    setState(238);
    match(OMDParser::BOND);
    setState(243);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(239);
      match(OMDParser::LBRACKET);
      setState(240);
      intConst();
      setState(241);
      match(OMDParser::RBRACKET);
    }
    setState(245);
    match(OMDParser::LCURLY);
    setState(249);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72057596151922688) != 0)) {
      setState(246);
      bondstatement();
      setState(251);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(252);
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
    setState(297);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(254);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(255);
        match(OMDParser::MEMBERS);
        setState(256);
        match(OMDParser::LPAREN);
        setState(257);
        inttuple();
        setState(258);
        match(OMDParser::RPAREN);
        setState(259);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::FIXED: {
        enterOuterAlt(_localctx, 3);
        setState(261);
        match(OMDParser::FIXED);
        setState(262);
        match(OMDParser::LPAREN);
        setState(263);
        floatConst();
        setState(264);
        match(OMDParser::RPAREN);
        setState(265);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 4);
        setState(267);
        match(OMDParser::HARMONIC);
        setState(268);
        match(OMDParser::LPAREN);
        setState(269);
        doubleNumberTuple();
        setState(270);
        match(OMDParser::RPAREN);
        setState(271);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 5);
        setState(273);
        match(OMDParser::CUBIC);
        setState(274);
        match(OMDParser::LPAREN);
        setState(275);
        doubleNumberTuple();
        setState(276);
        match(OMDParser::RPAREN);
        setState(277);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 6);
        setState(279);
        match(OMDParser::QUARTIC);
        setState(280);
        match(OMDParser::LPAREN);
        setState(281);
        doubleNumberTuple();
        setState(282);
        match(OMDParser::RPAREN);
        setState(283);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 7);
        setState(285);
        match(OMDParser::POLYNOMIAL);
        setState(286);
        match(OMDParser::LPAREN);
        setState(287);
        doubleNumberTuple();
        setState(288);
        match(OMDParser::RPAREN);
        setState(289);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::MORSE: {
        enterOuterAlt(_localctx, 8);
        setState(291);
        match(OMDParser::MORSE);
        setState(292);
        match(OMDParser::LPAREN);
        setState(293);
        doubleNumberTuple();
        setState(294);
        match(OMDParser::RPAREN);
        setState(295);
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
    setState(299);
    match(OMDParser::BEND);
    setState(304);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(300);
      match(OMDParser::LBRACKET);
      setState(301);
      intConst();
      setState(302);
      match(OMDParser::RBRACKET);
    }
    setState(306);
    match(OMDParser::LCURLY);
    setState(310);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72057610077011968) != 0)) {
      setState(307);
      bendstatement();
      setState(312);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(313);
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
    setState(364);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(315);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(316);
        match(OMDParser::MEMBERS);
        setState(317);
        match(OMDParser::LPAREN);
        setState(318);
        inttuple();
        setState(319);
        match(OMDParser::RPAREN);
        setState(320);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 3);
        setState(322);
        match(OMDParser::HARMONIC);
        setState(323);
        match(OMDParser::LPAREN);
        setState(324);
        doubleNumberTuple();
        setState(325);
        match(OMDParser::RPAREN);
        setState(326);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::GHOSTBEND: {
        enterOuterAlt(_localctx, 4);
        setState(328);
        match(OMDParser::GHOSTBEND);
        setState(329);
        match(OMDParser::LPAREN);
        setState(330);
        doubleNumberTuple();
        setState(331);
        match(OMDParser::RPAREN);
        setState(332);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::UREYBRADLEY: {
        enterOuterAlt(_localctx, 5);
        setState(334);
        match(OMDParser::UREYBRADLEY);
        setState(335);
        match(OMDParser::LPAREN);
        setState(336);
        doubleNumberTuple();
        setState(337);
        match(OMDParser::RPAREN);
        setState(338);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 6);
        setState(340);
        match(OMDParser::CUBIC);
        setState(341);
        match(OMDParser::LPAREN);
        setState(342);
        doubleNumberTuple();
        setState(343);
        match(OMDParser::RPAREN);
        setState(344);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 7);
        setState(346);
        match(OMDParser::QUARTIC);
        setState(347);
        match(OMDParser::LPAREN);
        setState(348);
        doubleNumberTuple();
        setState(349);
        match(OMDParser::RPAREN);
        setState(350);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 8);
        setState(352);
        match(OMDParser::POLYNOMIAL);
        setState(353);
        match(OMDParser::LPAREN);
        setState(354);
        doubleNumberTuple();
        setState(355);
        match(OMDParser::RPAREN);
        setState(356);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::COSINE: {
        enterOuterAlt(_localctx, 9);
        setState(358);
        match(OMDParser::COSINE);
        setState(359);
        match(OMDParser::LPAREN);
        setState(360);
        doubleNumberTuple();
        setState(361);
        match(OMDParser::RPAREN);
        setState(362);
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
    setState(366);
    match(OMDParser::TORSION);
    setState(371);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(367);
      match(OMDParser::LBRACKET);
      setState(368);
      intConst();
      setState(369);
      match(OMDParser::RBRACKET);
    }
    setState(373);
    match(OMDParser::LCURLY);
    setState(377);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72057852742664192) != 0)) {
      setState(374);
      torsionstatement();
      setState(379);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(380);
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
    setState(437);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(382);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(383);
        match(OMDParser::MEMBERS);
        setState(384);
        match(OMDParser::LPAREN);
        setState(385);
        inttuple();
        setState(386);
        match(OMDParser::RPAREN);
        setState(387);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::GHOSTTORSION: {
        enterOuterAlt(_localctx, 3);
        setState(389);
        match(OMDParser::GHOSTTORSION);
        setState(390);
        match(OMDParser::LPAREN);
        setState(391);
        doubleNumberTuple();
        setState(392);
        match(OMDParser::RPAREN);
        setState(393);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 4);
        setState(395);
        match(OMDParser::CUBIC);
        setState(396);
        match(OMDParser::LPAREN);
        setState(397);
        doubleNumberTuple();
        setState(398);
        match(OMDParser::RPAREN);
        setState(399);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 5);
        setState(401);
        match(OMDParser::QUARTIC);
        setState(402);
        match(OMDParser::LPAREN);
        setState(403);
        doubleNumberTuple();
        setState(404);
        match(OMDParser::RPAREN);
        setState(405);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 6);
        setState(407);
        match(OMDParser::POLYNOMIAL);
        setState(408);
        match(OMDParser::LPAREN);
        setState(409);
        doubleNumberTuple();
        setState(410);
        match(OMDParser::RPAREN);
        setState(411);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CHARMM: {
        enterOuterAlt(_localctx, 7);
        setState(413);
        match(OMDParser::CHARMM);
        setState(414);
        match(OMDParser::LPAREN);
        setState(415);
        doubleNumberTuple();
        setState(416);
        match(OMDParser::RPAREN);
        setState(417);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::OPLS: {
        enterOuterAlt(_localctx, 8);
        setState(419);
        match(OMDParser::OPLS);
        setState(420);
        match(OMDParser::LPAREN);
        setState(421);
        doubleNumberTuple();
        setState(422);
        match(OMDParser::RPAREN);
        setState(423);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::TRAPPE: {
        enterOuterAlt(_localctx, 9);
        setState(425);
        match(OMDParser::TRAPPE);
        setState(426);
        match(OMDParser::LPAREN);
        setState(427);
        doubleNumberTuple();
        setState(428);
        match(OMDParser::RPAREN);
        setState(429);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 10);
        setState(431);
        match(OMDParser::HARMONIC);
        setState(432);
        match(OMDParser::LPAREN);
        setState(433);
        doubleNumberTuple();
        setState(434);
        match(OMDParser::RPAREN);
        setState(435);
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
    setState(439);
    match(OMDParser::INVERSION);
    setState(444);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(440);
      match(OMDParser::LBRACKET);
      setState(441);
      intConst();
      setState(442);
      match(OMDParser::RBRACKET);
    }
    setState(446);
    match(OMDParser::LCURLY);
    setState(450);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72061717274034176) != 0)) {
      setState(447);
      inversionstatement();
      setState(452);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(453);
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
    setState(498);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(455);
        assignment();
        break;
      }

      case OMDParser::CENTER: {
        enterOuterAlt(_localctx, 2);
        setState(456);
        match(OMDParser::CENTER);
        setState(457);
        match(OMDParser::LPAREN);
        setState(458);
        intConst();
        setState(459);
        match(OMDParser::RPAREN);
        setState(460);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::SATELLITES: {
        enterOuterAlt(_localctx, 3);
        setState(462);
        match(OMDParser::SATELLITES);
        setState(463);
        match(OMDParser::LPAREN);
        setState(464);
        inttuple();
        setState(465);
        match(OMDParser::RPAREN);
        setState(466);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::AMBERIMPROPER: {
        enterOuterAlt(_localctx, 4);
        setState(468);
        match(OMDParser::AMBERIMPROPER);
        setState(469);
        match(OMDParser::LPAREN);
        setState(470);
        doubleNumberTuple();
        setState(471);
        match(OMDParser::RPAREN);
        setState(472);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::IMPROPERCOSINE: {
        enterOuterAlt(_localctx, 5);
        setState(474);
        match(OMDParser::IMPROPERCOSINE);
        setState(475);
        match(OMDParser::LPAREN);
        setState(476);
        doubleNumberTuple();
        setState(477);
        match(OMDParser::RPAREN);
        setState(478);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 6);
        setState(480);
        match(OMDParser::HARMONIC);
        setState(481);
        match(OMDParser::LPAREN);
        setState(482);
        doubleNumberTuple();
        setState(483);
        match(OMDParser::RPAREN);
        setState(484);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CENTRALATOMHEIGHT: {
        enterOuterAlt(_localctx, 7);
        setState(486);
        match(OMDParser::CENTRALATOMHEIGHT);
        setState(487);
        match(OMDParser::LPAREN);
        setState(488);
        doubleNumberTuple();
        setState(489);
        match(OMDParser::RPAREN);
        setState(490);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::DREIDING: {
        enterOuterAlt(_localctx, 8);
        setState(492);
        match(OMDParser::DREIDING);
        setState(493);
        match(OMDParser::LPAREN);
        setState(494);
        doubleNumberTuple();
        setState(495);
        match(OMDParser::RPAREN);
        setState(496);
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
    setState(500);
    match(OMDParser::RIGIDBODY);
    setState(501);
    match(OMDParser::LBRACKET);
    setState(502);
    intConst();
    setState(503);
    match(OMDParser::RBRACKET);
    setState(504);
    match(OMDParser::LCURLY);
    setState(508);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(505);
      rigidbodystatement();
      setState(510);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(511);
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
    setState(520);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(513);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(514);
        match(OMDParser::MEMBERS);
        setState(515);
        match(OMDParser::LPAREN);
        setState(516);
        inttuple();
        setState(517);
        match(OMDParser::RPAREN);
        setState(518);
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
    setState(522);
    match(OMDParser::CUTOFFGROUP);
    setState(527);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(523);
      match(OMDParser::LBRACKET);
      setState(524);
      intConst();
      setState(525);
      match(OMDParser::RBRACKET);
    }
    setState(529);
    match(OMDParser::LCURLY);
    setState(533);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(530);
      cutoffgroupstatement();
      setState(535);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(536);
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
    setState(545);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(538);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(539);
        match(OMDParser::MEMBERS);
        setState(540);
        match(OMDParser::LPAREN);
        setState(541);
        inttuple();
        setState(542);
        match(OMDParser::RPAREN);
        setState(543);
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
    setState(547);
    match(OMDParser::NODES);
    setState(552);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(548);
      match(OMDParser::LBRACKET);
      setState(549);
      intConst();
      setState(550);
      match(OMDParser::RBRACKET);
    }
    setState(554);
    match(OMDParser::LCURLY);
    setState(558);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(555);
      nodesstatement();
      setState(560);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(561);
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
    setState(570);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(563);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(564);
        match(OMDParser::MEMBERS);
        setState(565);
        match(OMDParser::LPAREN);
        setState(566);
        inttuple();
        setState(567);
        match(OMDParser::RPAREN);
        setState(568);
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
    setState(572);
    match(OMDParser::FRAGMENT);
    setState(573);
    match(OMDParser::LCURLY);
    setState(577);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & 72066390130958304) != 0)) {
      setState(574);
      fragmentstatement();
      setState(579);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(580);
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
    setState(592);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(582);
        assignment();
        break;
      }

      case OMDParser::ATOM: {
        enterOuterAlt(_localctx, 2);
        setState(583);
        atomblock();
        break;
      }

      case OMDParser::BOND: {
        enterOuterAlt(_localctx, 3);
        setState(584);
        bondblock();
        break;
      }

      case OMDParser::BEND: {
        enterOuterAlt(_localctx, 4);
        setState(585);
        bendblock();
        break;
      }

      case OMDParser::TORSION: {
        enterOuterAlt(_localctx, 5);
        setState(586);
        torsionblock();
        break;
      }

      case OMDParser::INVERSION: {
        enterOuterAlt(_localctx, 6);
        setState(587);
        inversionblock();
        break;
      }

      case OMDParser::RIGIDBODY: {
        enterOuterAlt(_localctx, 7);
        setState(588);
        rigidbodyblock();
        break;
      }

      case OMDParser::CUTOFFGROUP: {
        enterOuterAlt(_localctx, 8);
        setState(589);
        cutoffgroupblock();
        break;
      }

      case OMDParser::CONSTRAINT: {
        enterOuterAlt(_localctx, 9);
        setState(590);
        constraintblock();
        break;
      }

      case OMDParser::NODES: {
        enterOuterAlt(_localctx, 10);
        setState(591);
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
    setState(594);
    match(OMDParser::CONSTRAINT);
    setState(599);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(595);
      match(OMDParser::LBRACKET);
      setState(596);
      intConst();
      setState(597);
      match(OMDParser::RBRACKET);
    }
    setState(601);
    match(OMDParser::LCURLY);
    setState(605);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(602);
      constraintstatement();
      setState(607);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(608);
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
    setState(617);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(610);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(611);
        match(OMDParser::MEMBERS);
        setState(612);
        match(OMDParser::LPAREN);
        setState(613);
        inttuple();
        setState(614);
        match(OMDParser::RPAREN);
        setState(615);
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
    setState(619);
    match(OMDParser::SEQUENCE);
    setState(620);
    match(OMDParser::ASSIGNEQUAL);
    setState(621);
    constant();
    setState(622);
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
    setState(624);
    doubleNumber();
    setState(629);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::COMMA) {
      setState(625);
      match(OMDParser::COMMA);
      setState(626);
      doubleNumber();
      setState(631);
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
    setState(632);
    intConst();
    setState(637);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::COMMA) {
      setState(633);
      match(OMDParser::COMMA);
      setState(634);
      intConst();
      setState(639);
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
    setState(640);
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
    setState(644);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::NUM_LONG:
      case OMDParser::NUM_INT: {
        enterOuterAlt(_localctx, 1);
        setState(642);
        intConst();
        break;
      }

      case OMDParser::NUM_FLOAT:
      case OMDParser::NUM_DOUBLE: {
        enterOuterAlt(_localctx, 2);
        setState(643);
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
    setState(646);
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
    setState(648);
    match(OMDParser::LPAREN);
    setState(649);
    doubleNumberTuple();
    setState(650);
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
