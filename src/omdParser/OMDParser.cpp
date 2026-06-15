
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
      "lightblock", "velocityfieldblock", "minimizerblock", "moleculeblock", 
      "moleculestatement", "atomblock", "atomstatement", "bondblock", "bondstatement", 
      "bendblock", "bendstatement", "torsionblock", "torsionstatement", 
      "inversionblock", "inversionstatement", "rigidbodyblock", "rigidbodystatement", 
      "cutoffgroupblock", "cutoffgroupstatement", "nodesblock", "nodesstatement", 
      "fragmentblock", "fragmentstatement", "constraintblock", "constraintstatement", 
      "sequencestring", "doubleNumberTuple", "inttuple", "intConst", "doubleNumber", 
      "floatConst", "vectorConst"
    },
    std::vector<std::string>{
      "", "'true'", "'false'", "'component'", "'molecule'", "'zconstraint'", 
      "'restraint'", "'atom'", "'bond'", "'bend'", "'torsion'", "'inversion'", 
      "'rigidBody'", "'cutoffGroup'", "'constraint'", "'distance'", "'fragment'", 
      "'sequence'", "'members'", "'center'", "'satellites'", "'position'", 
      "'orientation'", "'flucQ'", "'RNEMD'", "'light'", "'velocityField'", 
      "'minimizer'", "'Fixed'", "'Harmonic'", "'Cubic'", "'Quartic'", "'Polynomial'", 
      "'Morse'", "'GhostBend'", "'UreyBradley'", "'Cosine'", "'GhostTorsion'", 
      "'Charmm'", "'Opls'", "'Trappe'", "'AmberImproper'", "'ImproperCosine'", 
      "'CentralAtomHeight'", "'Dreiding'", "'charge'", "'nodes'", "'='", 
      "':'", "','", "'\\u003F'", "';'", "'.'", "'('", "')'", "'['", "']'", 
      "'{'", "'}'"
    },
    std::vector<std::string>{
      "", "TRUE", "FALSE", "COMPONENT", "MOLECULE", "ZCONSTRAINT", "RESTRAINT", 
      "ATOM", "BOND", "BEND", "TORSION", "INVERSION", "RIGIDBODY", "CUTOFFGROUP", 
      "CONSTRAINT", "DISTANCE", "FRAGMENT", "SEQUENCE", "MEMBERS", "CENTER", 
      "SATELLITES", "POSITION", "ORIENTATION", "FLUCQ", "RNEMD", "LIGHT", 
      "VELOCITYFIELD", "MINIMIZER", "FIXED", "HARMONIC", "CUBIC", "QUARTIC", 
      "POLYNOMIAL", "MORSE", "GHOSTBEND", "UREYBRADLEY", "COSINE", "GHOSTTORSION", 
      "CHARMM", "OPLS", "TRAPPE", "AMBERIMPROPER", "IMPROPERCOSINE", "CENTRALATOMHEIGHT", 
      "DREIDING", "CHARGE", "NODES", "ASSIGNEQUAL", "COLON", "COMMA", "QUESTIONMARK", 
      "SEMICOLON", "DOT", "LPAREN", "RPAREN", "LBRACKET", "RBRACKET", "LCURLY", 
      "RCURLY", "NUM_LONG", "NUM_INT", "NUM_FLOAT", "NUM_DOUBLE", "CharLiteral", 
      "StringLiteral", "ID", "Whitespace", "Newline", "LineContinuation", 
      "Comment", "CPPComment", "PREPROC_DIRECTIVE"
    }
  );
  static const int32_t serializedATNSegment[] = {
  	4,1,71,668,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,6,2,
  	7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,13,7,13,2,14,7,
  	14,2,15,7,15,2,16,7,16,2,17,7,17,2,18,7,18,2,19,7,19,2,20,7,20,2,21,7,
  	21,2,22,7,22,2,23,7,23,2,24,7,24,2,25,7,25,2,26,7,26,2,27,7,27,2,28,7,
  	28,2,29,7,29,2,30,7,30,2,31,7,31,2,32,7,32,2,33,7,33,2,34,7,34,2,35,7,
  	35,2,36,7,36,2,37,7,37,2,38,7,38,2,39,7,39,2,40,7,40,1,0,5,0,84,8,0,10,
  	0,12,0,87,9,0,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,
  	102,8,1,1,2,1,2,1,2,1,2,1,2,1,3,1,3,1,3,1,3,1,3,1,3,1,3,3,3,116,8,3,1,
  	4,1,4,1,4,5,4,121,8,4,10,4,12,4,124,9,4,1,4,1,4,1,5,1,5,1,5,5,5,131,8,
  	5,10,5,12,5,134,9,5,1,5,1,5,1,6,1,6,1,6,5,6,141,8,6,10,6,12,6,144,9,6,
  	1,6,1,6,1,7,1,7,1,7,5,7,151,8,7,10,7,12,7,154,9,7,1,7,1,7,1,8,1,8,1,8,
  	5,8,161,8,8,10,8,12,8,164,9,8,1,8,1,8,1,9,1,9,1,9,5,9,171,8,9,10,9,12,
  	9,174,9,9,1,9,1,9,1,10,1,10,1,10,5,10,181,8,10,10,10,12,10,184,9,10,1,
  	10,1,10,1,11,1,11,1,11,5,11,191,8,11,10,11,12,11,194,9,11,1,11,1,11,1,
  	12,1,12,1,12,5,12,201,8,12,10,12,12,12,204,9,12,1,12,1,12,1,13,1,13,1,
  	13,1,13,1,13,1,13,1,13,1,13,1,13,1,13,3,13,218,8,13,1,14,1,14,1,14,1,
  	14,1,14,1,14,5,14,226,8,14,10,14,12,14,229,9,14,1,14,1,14,1,15,1,15,1,
  	15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,15,1,
  	15,1,15,1,15,3,15,252,8,15,1,16,1,16,1,16,1,16,1,16,3,16,259,8,16,1,16,
  	1,16,5,16,263,8,16,10,16,12,16,266,9,16,1,16,1,16,1,17,1,17,1,17,1,17,
  	1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,
  	1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,
  	1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,3,17,313,8,17,
  	1,18,1,18,1,18,1,18,1,18,3,18,320,8,18,1,18,1,18,5,18,324,8,18,10,18,
  	12,18,327,9,18,1,18,1,18,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,
  	1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,
  	1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,
  	1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,1,19,3,19,380,
  	8,19,1,20,1,20,1,20,1,20,1,20,3,20,387,8,20,1,20,1,20,5,20,391,8,20,10,
  	20,12,20,394,9,20,1,20,1,20,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,
  	21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,
  	21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,
  	21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,21,1,
  	21,1,21,1,21,1,21,1,21,3,21,453,8,21,1,22,1,22,1,22,1,22,1,22,3,22,460,
  	8,22,1,22,1,22,5,22,464,8,22,10,22,12,22,467,9,22,1,22,1,22,1,23,1,23,
  	1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,
  	1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,
  	1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,1,23,3,23,
  	514,8,23,1,24,1,24,1,24,1,24,1,24,1,24,5,24,522,8,24,10,24,12,24,525,
  	9,24,1,24,1,24,1,25,1,25,1,25,1,25,1,25,1,25,1,25,3,25,536,8,25,1,26,
  	1,26,1,26,1,26,1,26,3,26,543,8,26,1,26,1,26,5,26,547,8,26,10,26,12,26,
  	550,9,26,1,26,1,26,1,27,1,27,1,27,1,27,1,27,1,27,1,27,3,27,561,8,27,1,
  	28,1,28,1,28,1,28,1,28,3,28,568,8,28,1,28,1,28,5,28,572,8,28,10,28,12,
  	28,575,9,28,1,28,1,28,1,29,1,29,1,29,1,29,1,29,1,29,1,29,3,29,586,8,29,
  	1,30,1,30,1,30,5,30,591,8,30,10,30,12,30,594,9,30,1,30,1,30,1,31,1,31,
  	1,31,1,31,1,31,1,31,1,31,1,31,1,31,1,31,3,31,608,8,31,1,32,1,32,1,32,
  	1,32,1,32,3,32,615,8,32,1,32,1,32,5,32,619,8,32,10,32,12,32,622,9,32,
  	1,32,1,32,1,33,1,33,1,33,1,33,1,33,1,33,1,33,3,33,633,8,33,1,34,1,34,
  	1,34,1,34,1,34,1,35,1,35,1,35,5,35,643,8,35,10,35,12,35,646,9,35,1,36,
  	1,36,1,36,5,36,651,8,36,10,36,12,36,654,9,36,1,37,1,37,1,38,1,38,3,38,
  	660,8,38,1,39,1,39,1,40,1,40,1,40,1,40,1,40,0,0,41,0,2,4,6,8,10,12,14,
  	16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,
  	62,64,66,68,70,72,74,76,78,80,0,2,1,0,59,60,1,0,61,62,728,0,85,1,0,0,
  	0,2,101,1,0,0,0,4,103,1,0,0,0,6,115,1,0,0,0,8,117,1,0,0,0,10,127,1,0,
  	0,0,12,137,1,0,0,0,14,147,1,0,0,0,16,157,1,0,0,0,18,167,1,0,0,0,20,177,
  	1,0,0,0,22,187,1,0,0,0,24,197,1,0,0,0,26,217,1,0,0,0,28,219,1,0,0,0,30,
  	251,1,0,0,0,32,253,1,0,0,0,34,312,1,0,0,0,36,314,1,0,0,0,38,379,1,0,0,
  	0,40,381,1,0,0,0,42,452,1,0,0,0,44,454,1,0,0,0,46,513,1,0,0,0,48,515,
  	1,0,0,0,50,535,1,0,0,0,52,537,1,0,0,0,54,560,1,0,0,0,56,562,1,0,0,0,58,
  	585,1,0,0,0,60,587,1,0,0,0,62,607,1,0,0,0,64,609,1,0,0,0,66,632,1,0,0,
  	0,68,634,1,0,0,0,70,639,1,0,0,0,72,647,1,0,0,0,74,655,1,0,0,0,76,659,
  	1,0,0,0,78,661,1,0,0,0,80,663,1,0,0,0,82,84,3,2,1,0,83,82,1,0,0,0,84,
  	87,1,0,0,0,85,83,1,0,0,0,85,86,1,0,0,0,86,88,1,0,0,0,87,85,1,0,0,0,88,
  	89,5,0,0,1,89,1,1,0,0,0,90,102,3,4,2,0,91,102,3,8,4,0,92,102,3,24,12,
  	0,93,102,3,60,30,0,94,102,3,10,5,0,95,102,3,12,6,0,96,102,3,14,7,0,97,
  	102,3,16,8,0,98,102,3,18,9,0,99,102,3,20,10,0,100,102,3,22,11,0,101,90,
  	1,0,0,0,101,91,1,0,0,0,101,92,1,0,0,0,101,93,1,0,0,0,101,94,1,0,0,0,101,
  	95,1,0,0,0,101,96,1,0,0,0,101,97,1,0,0,0,101,98,1,0,0,0,101,99,1,0,0,
  	0,101,100,1,0,0,0,102,3,1,0,0,0,103,104,5,65,0,0,104,105,5,47,0,0,105,
  	106,3,6,3,0,106,107,5,51,0,0,107,5,1,0,0,0,108,116,3,74,37,0,109,116,
  	3,78,39,0,110,116,3,80,40,0,111,116,5,1,0,0,112,116,5,2,0,0,113,116,5,
  	65,0,0,114,116,5,64,0,0,115,108,1,0,0,0,115,109,1,0,0,0,115,110,1,0,0,
  	0,115,111,1,0,0,0,115,112,1,0,0,0,115,113,1,0,0,0,115,114,1,0,0,0,116,
  	7,1,0,0,0,117,118,5,3,0,0,118,122,5,57,0,0,119,121,3,4,2,0,120,119,1,
  	0,0,0,121,124,1,0,0,0,122,120,1,0,0,0,122,123,1,0,0,0,123,125,1,0,0,0,
  	124,122,1,0,0,0,125,126,5,58,0,0,126,9,1,0,0,0,127,128,5,5,0,0,128,132,
  	5,57,0,0,129,131,3,4,2,0,130,129,1,0,0,0,131,134,1,0,0,0,132,130,1,0,
  	0,0,132,133,1,0,0,0,133,135,1,0,0,0,134,132,1,0,0,0,135,136,5,58,0,0,
  	136,11,1,0,0,0,137,138,5,6,0,0,138,142,5,57,0,0,139,141,3,4,2,0,140,139,
  	1,0,0,0,141,144,1,0,0,0,142,140,1,0,0,0,142,143,1,0,0,0,143,145,1,0,0,
  	0,144,142,1,0,0,0,145,146,5,58,0,0,146,13,1,0,0,0,147,148,5,23,0,0,148,
  	152,5,57,0,0,149,151,3,4,2,0,150,149,1,0,0,0,151,154,1,0,0,0,152,150,
  	1,0,0,0,152,153,1,0,0,0,153,155,1,0,0,0,154,152,1,0,0,0,155,156,5,58,
  	0,0,156,15,1,0,0,0,157,158,5,24,0,0,158,162,5,57,0,0,159,161,3,4,2,0,
  	160,159,1,0,0,0,161,164,1,0,0,0,162,160,1,0,0,0,162,163,1,0,0,0,163,165,
  	1,0,0,0,164,162,1,0,0,0,165,166,5,58,0,0,166,17,1,0,0,0,167,168,5,25,
  	0,0,168,172,5,57,0,0,169,171,3,4,2,0,170,169,1,0,0,0,171,174,1,0,0,0,
  	172,170,1,0,0,0,172,173,1,0,0,0,173,175,1,0,0,0,174,172,1,0,0,0,175,176,
  	5,58,0,0,176,19,1,0,0,0,177,178,5,26,0,0,178,182,5,57,0,0,179,181,3,4,
  	2,0,180,179,1,0,0,0,181,184,1,0,0,0,182,180,1,0,0,0,182,183,1,0,0,0,183,
  	185,1,0,0,0,184,182,1,0,0,0,185,186,5,58,0,0,186,21,1,0,0,0,187,188,5,
  	27,0,0,188,192,5,57,0,0,189,191,3,4,2,0,190,189,1,0,0,0,191,194,1,0,0,
  	0,192,190,1,0,0,0,192,193,1,0,0,0,193,195,1,0,0,0,194,192,1,0,0,0,195,
  	196,5,58,0,0,196,23,1,0,0,0,197,198,5,4,0,0,198,202,5,57,0,0,199,201,
  	3,26,13,0,200,199,1,0,0,0,201,204,1,0,0,0,202,200,1,0,0,0,202,203,1,0,
  	0,0,203,205,1,0,0,0,204,202,1,0,0,0,205,206,5,58,0,0,206,25,1,0,0,0,207,
  	218,3,4,2,0,208,218,3,28,14,0,209,218,3,32,16,0,210,218,3,36,18,0,211,
  	218,3,40,20,0,212,218,3,44,22,0,213,218,3,48,24,0,214,218,3,52,26,0,215,
  	218,3,64,32,0,216,218,3,68,34,0,217,207,1,0,0,0,217,208,1,0,0,0,217,209,
  	1,0,0,0,217,210,1,0,0,0,217,211,1,0,0,0,217,212,1,0,0,0,217,213,1,0,0,
  	0,217,214,1,0,0,0,217,215,1,0,0,0,217,216,1,0,0,0,218,27,1,0,0,0,219,
  	220,5,7,0,0,220,221,5,55,0,0,221,222,3,74,37,0,222,223,5,56,0,0,223,227,
  	5,57,0,0,224,226,3,30,15,0,225,224,1,0,0,0,226,229,1,0,0,0,227,225,1,
  	0,0,0,227,228,1,0,0,0,228,230,1,0,0,0,229,227,1,0,0,0,230,231,5,58,0,
  	0,231,29,1,0,0,0,232,252,3,4,2,0,233,234,5,21,0,0,234,235,5,53,0,0,235,
  	236,3,70,35,0,236,237,5,54,0,0,237,238,5,51,0,0,238,252,1,0,0,0,239,240,
  	5,22,0,0,240,241,5,53,0,0,241,242,3,70,35,0,242,243,5,54,0,0,243,244,
  	5,51,0,0,244,252,1,0,0,0,245,246,5,45,0,0,246,247,5,53,0,0,247,248,3,
  	78,39,0,248,249,5,54,0,0,249,250,5,51,0,0,250,252,1,0,0,0,251,232,1,0,
  	0,0,251,233,1,0,0,0,251,239,1,0,0,0,251,245,1,0,0,0,252,31,1,0,0,0,253,
  	258,5,8,0,0,254,255,5,55,0,0,255,256,3,74,37,0,256,257,5,56,0,0,257,259,
  	1,0,0,0,258,254,1,0,0,0,258,259,1,0,0,0,259,260,1,0,0,0,260,264,5,57,
  	0,0,261,263,3,34,17,0,262,261,1,0,0,0,263,266,1,0,0,0,264,262,1,0,0,0,
  	264,265,1,0,0,0,265,267,1,0,0,0,266,264,1,0,0,0,267,268,5,58,0,0,268,
  	33,1,0,0,0,269,313,3,4,2,0,270,271,5,18,0,0,271,272,5,53,0,0,272,273,
  	3,72,36,0,273,274,5,54,0,0,274,275,5,51,0,0,275,313,1,0,0,0,276,277,5,
  	28,0,0,277,278,5,53,0,0,278,279,3,78,39,0,279,280,5,54,0,0,280,281,5,
  	51,0,0,281,313,1,0,0,0,282,283,5,29,0,0,283,284,5,53,0,0,284,285,3,70,
  	35,0,285,286,5,54,0,0,286,287,5,51,0,0,287,313,1,0,0,0,288,289,5,30,0,
  	0,289,290,5,53,0,0,290,291,3,70,35,0,291,292,5,54,0,0,292,293,5,51,0,
  	0,293,313,1,0,0,0,294,295,5,31,0,0,295,296,5,53,0,0,296,297,3,70,35,0,
  	297,298,5,54,0,0,298,299,5,51,0,0,299,313,1,0,0,0,300,301,5,32,0,0,301,
  	302,5,53,0,0,302,303,3,70,35,0,303,304,5,54,0,0,304,305,5,51,0,0,305,
  	313,1,0,0,0,306,307,5,33,0,0,307,308,5,53,0,0,308,309,3,70,35,0,309,310,
  	5,54,0,0,310,311,5,51,0,0,311,313,1,0,0,0,312,269,1,0,0,0,312,270,1,0,
  	0,0,312,276,1,0,0,0,312,282,1,0,0,0,312,288,1,0,0,0,312,294,1,0,0,0,312,
  	300,1,0,0,0,312,306,1,0,0,0,313,35,1,0,0,0,314,319,5,9,0,0,315,316,5,
  	55,0,0,316,317,3,74,37,0,317,318,5,56,0,0,318,320,1,0,0,0,319,315,1,0,
  	0,0,319,320,1,0,0,0,320,321,1,0,0,0,321,325,5,57,0,0,322,324,3,38,19,
  	0,323,322,1,0,0,0,324,327,1,0,0,0,325,323,1,0,0,0,325,326,1,0,0,0,326,
  	328,1,0,0,0,327,325,1,0,0,0,328,329,5,58,0,0,329,37,1,0,0,0,330,380,3,
  	4,2,0,331,332,5,18,0,0,332,333,5,53,0,0,333,334,3,72,36,0,334,335,5,54,
  	0,0,335,336,5,51,0,0,336,380,1,0,0,0,337,338,5,29,0,0,338,339,5,53,0,
  	0,339,340,3,70,35,0,340,341,5,54,0,0,341,342,5,51,0,0,342,380,1,0,0,0,
  	343,344,5,34,0,0,344,345,5,53,0,0,345,346,3,70,35,0,346,347,5,54,0,0,
  	347,348,5,51,0,0,348,380,1,0,0,0,349,350,5,35,0,0,350,351,5,53,0,0,351,
  	352,3,70,35,0,352,353,5,54,0,0,353,354,5,51,0,0,354,380,1,0,0,0,355,356,
  	5,30,0,0,356,357,5,53,0,0,357,358,3,70,35,0,358,359,5,54,0,0,359,360,
  	5,51,0,0,360,380,1,0,0,0,361,362,5,31,0,0,362,363,5,53,0,0,363,364,3,
  	70,35,0,364,365,5,54,0,0,365,366,5,51,0,0,366,380,1,0,0,0,367,368,5,32,
  	0,0,368,369,5,53,0,0,369,370,3,70,35,0,370,371,5,54,0,0,371,372,5,51,
  	0,0,372,380,1,0,0,0,373,374,5,36,0,0,374,375,5,53,0,0,375,376,3,70,35,
  	0,376,377,5,54,0,0,377,378,5,51,0,0,378,380,1,0,0,0,379,330,1,0,0,0,379,
  	331,1,0,0,0,379,337,1,0,0,0,379,343,1,0,0,0,379,349,1,0,0,0,379,355,1,
  	0,0,0,379,361,1,0,0,0,379,367,1,0,0,0,379,373,1,0,0,0,380,39,1,0,0,0,
  	381,386,5,10,0,0,382,383,5,55,0,0,383,384,3,74,37,0,384,385,5,56,0,0,
  	385,387,1,0,0,0,386,382,1,0,0,0,386,387,1,0,0,0,387,388,1,0,0,0,388,392,
  	5,57,0,0,389,391,3,42,21,0,390,389,1,0,0,0,391,394,1,0,0,0,392,390,1,
  	0,0,0,392,393,1,0,0,0,393,395,1,0,0,0,394,392,1,0,0,0,395,396,5,58,0,
  	0,396,41,1,0,0,0,397,453,3,4,2,0,398,399,5,18,0,0,399,400,5,53,0,0,400,
  	401,3,72,36,0,401,402,5,54,0,0,402,403,5,51,0,0,403,453,1,0,0,0,404,405,
  	5,37,0,0,405,406,5,53,0,0,406,407,3,70,35,0,407,408,5,54,0,0,408,409,
  	5,51,0,0,409,453,1,0,0,0,410,411,5,30,0,0,411,412,5,53,0,0,412,413,3,
  	70,35,0,413,414,5,54,0,0,414,415,5,51,0,0,415,453,1,0,0,0,416,417,5,31,
  	0,0,417,418,5,53,0,0,418,419,3,70,35,0,419,420,5,54,0,0,420,421,5,51,
  	0,0,421,453,1,0,0,0,422,423,5,32,0,0,423,424,5,53,0,0,424,425,3,70,35,
  	0,425,426,5,54,0,0,426,427,5,51,0,0,427,453,1,0,0,0,428,429,5,38,0,0,
  	429,430,5,53,0,0,430,431,3,70,35,0,431,432,5,54,0,0,432,433,5,51,0,0,
  	433,453,1,0,0,0,434,435,5,39,0,0,435,436,5,53,0,0,436,437,3,70,35,0,437,
  	438,5,54,0,0,438,439,5,51,0,0,439,453,1,0,0,0,440,441,5,40,0,0,441,442,
  	5,53,0,0,442,443,3,70,35,0,443,444,5,54,0,0,444,445,5,51,0,0,445,453,
  	1,0,0,0,446,447,5,29,0,0,447,448,5,53,0,0,448,449,3,70,35,0,449,450,5,
  	54,0,0,450,451,5,51,0,0,451,453,1,0,0,0,452,397,1,0,0,0,452,398,1,0,0,
  	0,452,404,1,0,0,0,452,410,1,0,0,0,452,416,1,0,0,0,452,422,1,0,0,0,452,
  	428,1,0,0,0,452,434,1,0,0,0,452,440,1,0,0,0,452,446,1,0,0,0,453,43,1,
  	0,0,0,454,459,5,11,0,0,455,456,5,55,0,0,456,457,3,74,37,0,457,458,5,56,
  	0,0,458,460,1,0,0,0,459,455,1,0,0,0,459,460,1,0,0,0,460,461,1,0,0,0,461,
  	465,5,57,0,0,462,464,3,46,23,0,463,462,1,0,0,0,464,467,1,0,0,0,465,463,
  	1,0,0,0,465,466,1,0,0,0,466,468,1,0,0,0,467,465,1,0,0,0,468,469,5,58,
  	0,0,469,45,1,0,0,0,470,514,3,4,2,0,471,472,5,19,0,0,472,473,5,53,0,0,
  	473,474,3,74,37,0,474,475,5,54,0,0,475,476,5,51,0,0,476,514,1,0,0,0,477,
  	478,5,20,0,0,478,479,5,53,0,0,479,480,3,72,36,0,480,481,5,54,0,0,481,
  	482,5,51,0,0,482,514,1,0,0,0,483,484,5,41,0,0,484,485,5,53,0,0,485,486,
  	3,70,35,0,486,487,5,54,0,0,487,488,5,51,0,0,488,514,1,0,0,0,489,490,5,
  	42,0,0,490,491,5,53,0,0,491,492,3,70,35,0,492,493,5,54,0,0,493,494,5,
  	51,0,0,494,514,1,0,0,0,495,496,5,29,0,0,496,497,5,53,0,0,497,498,3,70,
  	35,0,498,499,5,54,0,0,499,500,5,51,0,0,500,514,1,0,0,0,501,502,5,43,0,
  	0,502,503,5,53,0,0,503,504,3,70,35,0,504,505,5,54,0,0,505,506,5,51,0,
  	0,506,514,1,0,0,0,507,508,5,44,0,0,508,509,5,53,0,0,509,510,3,70,35,0,
  	510,511,5,54,0,0,511,512,5,51,0,0,512,514,1,0,0,0,513,470,1,0,0,0,513,
  	471,1,0,0,0,513,477,1,0,0,0,513,483,1,0,0,0,513,489,1,0,0,0,513,495,1,
  	0,0,0,513,501,1,0,0,0,513,507,1,0,0,0,514,47,1,0,0,0,515,516,5,12,0,0,
  	516,517,5,55,0,0,517,518,3,74,37,0,518,519,5,56,0,0,519,523,5,57,0,0,
  	520,522,3,50,25,0,521,520,1,0,0,0,522,525,1,0,0,0,523,521,1,0,0,0,523,
  	524,1,0,0,0,524,526,1,0,0,0,525,523,1,0,0,0,526,527,5,58,0,0,527,49,1,
  	0,0,0,528,536,3,4,2,0,529,530,5,18,0,0,530,531,5,53,0,0,531,532,3,72,
  	36,0,532,533,5,54,0,0,533,534,5,51,0,0,534,536,1,0,0,0,535,528,1,0,0,
  	0,535,529,1,0,0,0,536,51,1,0,0,0,537,542,5,13,0,0,538,539,5,55,0,0,539,
  	540,3,74,37,0,540,541,5,56,0,0,541,543,1,0,0,0,542,538,1,0,0,0,542,543,
  	1,0,0,0,543,544,1,0,0,0,544,548,5,57,0,0,545,547,3,54,27,0,546,545,1,
  	0,0,0,547,550,1,0,0,0,548,546,1,0,0,0,548,549,1,0,0,0,549,551,1,0,0,0,
  	550,548,1,0,0,0,551,552,5,58,0,0,552,53,1,0,0,0,553,561,3,4,2,0,554,555,
  	5,18,0,0,555,556,5,53,0,0,556,557,3,72,36,0,557,558,5,54,0,0,558,559,
  	5,51,0,0,559,561,1,0,0,0,560,553,1,0,0,0,560,554,1,0,0,0,561,55,1,0,0,
  	0,562,567,5,46,0,0,563,564,5,55,0,0,564,565,3,74,37,0,565,566,5,56,0,
  	0,566,568,1,0,0,0,567,563,1,0,0,0,567,568,1,0,0,0,568,569,1,0,0,0,569,
  	573,5,57,0,0,570,572,3,58,29,0,571,570,1,0,0,0,572,575,1,0,0,0,573,571,
  	1,0,0,0,573,574,1,0,0,0,574,576,1,0,0,0,575,573,1,0,0,0,576,577,5,58,
  	0,0,577,57,1,0,0,0,578,586,3,4,2,0,579,580,5,18,0,0,580,581,5,53,0,0,
  	581,582,3,72,36,0,582,583,5,54,0,0,583,584,5,51,0,0,584,586,1,0,0,0,585,
  	578,1,0,0,0,585,579,1,0,0,0,586,59,1,0,0,0,587,588,5,16,0,0,588,592,5,
  	57,0,0,589,591,3,62,31,0,590,589,1,0,0,0,591,594,1,0,0,0,592,590,1,0,
  	0,0,592,593,1,0,0,0,593,595,1,0,0,0,594,592,1,0,0,0,595,596,5,58,0,0,
  	596,61,1,0,0,0,597,608,3,4,2,0,598,608,3,28,14,0,599,608,3,32,16,0,600,
  	608,3,36,18,0,601,608,3,40,20,0,602,608,3,44,22,0,603,608,3,48,24,0,604,
  	608,3,52,26,0,605,608,3,64,32,0,606,608,3,56,28,0,607,597,1,0,0,0,607,
  	598,1,0,0,0,607,599,1,0,0,0,607,600,1,0,0,0,607,601,1,0,0,0,607,602,1,
  	0,0,0,607,603,1,0,0,0,607,604,1,0,0,0,607,605,1,0,0,0,607,606,1,0,0,0,
  	608,63,1,0,0,0,609,614,5,14,0,0,610,611,5,55,0,0,611,612,3,74,37,0,612,
  	613,5,56,0,0,613,615,1,0,0,0,614,610,1,0,0,0,614,615,1,0,0,0,615,616,
  	1,0,0,0,616,620,5,57,0,0,617,619,3,66,33,0,618,617,1,0,0,0,619,622,1,
  	0,0,0,620,618,1,0,0,0,620,621,1,0,0,0,621,623,1,0,0,0,622,620,1,0,0,0,
  	623,624,5,58,0,0,624,65,1,0,0,0,625,633,3,4,2,0,626,627,5,18,0,0,627,
  	628,5,53,0,0,628,629,3,72,36,0,629,630,5,54,0,0,630,631,5,51,0,0,631,
  	633,1,0,0,0,632,625,1,0,0,0,632,626,1,0,0,0,633,67,1,0,0,0,634,635,5,
  	17,0,0,635,636,5,47,0,0,636,637,3,6,3,0,637,638,5,51,0,0,638,69,1,0,0,
  	0,639,644,3,76,38,0,640,641,5,49,0,0,641,643,3,76,38,0,642,640,1,0,0,
  	0,643,646,1,0,0,0,644,642,1,0,0,0,644,645,1,0,0,0,645,71,1,0,0,0,646,
  	644,1,0,0,0,647,652,3,74,37,0,648,649,5,49,0,0,649,651,3,74,37,0,650,
  	648,1,0,0,0,651,654,1,0,0,0,652,650,1,0,0,0,652,653,1,0,0,0,653,73,1,
  	0,0,0,654,652,1,0,0,0,655,656,7,0,0,0,656,75,1,0,0,0,657,660,3,74,37,
  	0,658,660,3,78,39,0,659,657,1,0,0,0,659,658,1,0,0,0,660,77,1,0,0,0,661,
  	662,7,1,0,0,662,79,1,0,0,0,663,664,5,53,0,0,664,665,3,70,35,0,665,666,
  	5,54,0,0,666,81,1,0,0,0,43,85,101,115,122,132,142,152,162,172,182,192,
  	202,217,227,251,258,264,312,319,325,379,386,392,452,459,465,513,523,535,
  	542,548,560,567,573,585,592,607,614,620,632,644,652,659
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
    setState(85);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 3) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 3)) & 4611686018459901967) != 0)) {
      setState(82);
      statement();
      setState(87);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(88);
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

OMDParser::VelocityfieldblockContext* OMDParser::StatementContext::velocityfieldblock() {
  return getRuleContext<OMDParser::VelocityfieldblockContext>(0);
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
    setState(101);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(90);
        assignment();
        break;
      }

      case OMDParser::COMPONENT: {
        enterOuterAlt(_localctx, 2);
        setState(91);
        componentblock();
        break;
      }

      case OMDParser::MOLECULE: {
        enterOuterAlt(_localctx, 3);
        setState(92);
        moleculeblock();
        break;
      }

      case OMDParser::FRAGMENT: {
        enterOuterAlt(_localctx, 4);
        setState(93);
        fragmentblock();
        break;
      }

      case OMDParser::ZCONSTRAINT: {
        enterOuterAlt(_localctx, 5);
        setState(94);
        zconstraintblock();
        break;
      }

      case OMDParser::RESTRAINT: {
        enterOuterAlt(_localctx, 6);
        setState(95);
        restraintblock();
        break;
      }

      case OMDParser::FLUCQ: {
        enterOuterAlt(_localctx, 7);
        setState(96);
        flucqblock();
        break;
      }

      case OMDParser::RNEMD: {
        enterOuterAlt(_localctx, 8);
        setState(97);
        rnemdblock();
        break;
      }

      case OMDParser::LIGHT: {
        enterOuterAlt(_localctx, 9);
        setState(98);
        lightblock();
        break;
      }

      case OMDParser::VELOCITYFIELD: {
        enterOuterAlt(_localctx, 10);
        setState(99);
        velocityfieldblock();
        break;
      }

      case OMDParser::MINIMIZER: {
        enterOuterAlt(_localctx, 11);
        setState(100);
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
    setState(103);
    match(OMDParser::ID);
    setState(104);
    match(OMDParser::ASSIGNEQUAL);
    setState(105);
    constant();
    setState(106);
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
    setState(115);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::NUM_LONG:
      case OMDParser::NUM_INT: {
        enterOuterAlt(_localctx, 1);
        setState(108);
        intConst();
        break;
      }

      case OMDParser::NUM_FLOAT:
      case OMDParser::NUM_DOUBLE: {
        enterOuterAlt(_localctx, 2);
        setState(109);
        floatConst();
        break;
      }

      case OMDParser::LPAREN: {
        enterOuterAlt(_localctx, 3);
        setState(110);
        vectorConst();
        break;
      }

      case OMDParser::TRUE: {
        enterOuterAlt(_localctx, 4);
        setState(111);
        match(OMDParser::TRUE);
        break;
      }

      case OMDParser::FALSE: {
        enterOuterAlt(_localctx, 5);
        setState(112);
        match(OMDParser::FALSE);
        break;
      }

      case OMDParser::ID: {
        enterOuterAlt(_localctx, 6);
        setState(113);
        match(OMDParser::ID);
        break;
      }

      case OMDParser::StringLiteral: {
        enterOuterAlt(_localctx, 7);
        setState(114);
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
    setState(117);
    match(OMDParser::COMPONENT);
    setState(118);
    match(OMDParser::LCURLY);
    setState(122);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(119);
      assignment();
      setState(124);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(125);
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
    setState(127);
    match(OMDParser::ZCONSTRAINT);
    setState(128);
    match(OMDParser::LCURLY);
    setState(132);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(129);
      assignment();
      setState(134);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(135);
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
    setState(137);
    match(OMDParser::RESTRAINT);
    setState(138);
    match(OMDParser::LCURLY);
    setState(142);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(139);
      assignment();
      setState(144);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(145);
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
    setState(147);
    match(OMDParser::FLUCQ);
    setState(148);
    match(OMDParser::LCURLY);
    setState(152);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(149);
      assignment();
      setState(154);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(155);
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
    setState(157);
    match(OMDParser::RNEMD);
    setState(158);
    match(OMDParser::LCURLY);
    setState(162);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(159);
      assignment();
      setState(164);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(165);
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
    setState(167);
    match(OMDParser::LIGHT);
    setState(168);
    match(OMDParser::LCURLY);
    setState(172);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(169);
      assignment();
      setState(174);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(175);
    match(OMDParser::RCURLY);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- VelocityfieldblockContext ------------------------------------------------------------------

OMDParser::VelocityfieldblockContext::VelocityfieldblockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* OMDParser::VelocityfieldblockContext::VELOCITYFIELD() {
  return getToken(OMDParser::VELOCITYFIELD, 0);
}

tree::TerminalNode* OMDParser::VelocityfieldblockContext::LCURLY() {
  return getToken(OMDParser::LCURLY, 0);
}

tree::TerminalNode* OMDParser::VelocityfieldblockContext::RCURLY() {
  return getToken(OMDParser::RCURLY, 0);
}

std::vector<OMDParser::AssignmentContext *> OMDParser::VelocityfieldblockContext::assignment() {
  return getRuleContexts<OMDParser::AssignmentContext>();
}

OMDParser::AssignmentContext* OMDParser::VelocityfieldblockContext::assignment(size_t i) {
  return getRuleContext<OMDParser::AssignmentContext>(i);
}


size_t OMDParser::VelocityfieldblockContext::getRuleIndex() const {
  return OMDParser::RuleVelocityfieldblock;
}

void OMDParser::VelocityfieldblockContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterVelocityfieldblock(this);
}

void OMDParser::VelocityfieldblockContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<OMDListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitVelocityfieldblock(this);
}


std::any OMDParser::VelocityfieldblockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<OMDVisitor*>(visitor))
    return parserVisitor->visitVelocityfieldblock(this);
  else
    return visitor->visitChildren(this);
}

OMDParser::VelocityfieldblockContext* OMDParser::velocityfieldblock() {
  VelocityfieldblockContext *_localctx = _tracker.createInstance<VelocityfieldblockContext>(_ctx, getState());
  enterRule(_localctx, 20, OMDParser::RuleVelocityfieldblock);
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
    setState(177);
    match(OMDParser::VELOCITYFIELD);
    setState(178);
    match(OMDParser::LCURLY);
    setState(182);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(179);
      assignment();
      setState(184);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(185);
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
  enterRule(_localctx, 22, OMDParser::RuleMinimizerblock);
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
    setState(187);
    match(OMDParser::MINIMIZER);
    setState(188);
    match(OMDParser::LCURLY);
    setState(192);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::ID) {
      setState(189);
      assignment();
      setState(194);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(195);
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
  enterRule(_localctx, 24, OMDParser::RuleMoleculeblock);
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
    setState(197);
    match(OMDParser::MOLECULE);
    setState(198);
    match(OMDParser::LCURLY);
    setState(202);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 7) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 7)) & 288230376151713023) != 0)) {
      setState(199);
      moleculestatement();
      setState(204);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(205);
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
  enterRule(_localctx, 26, OMDParser::RuleMoleculestatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(217);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(207);
        assignment();
        break;
      }

      case OMDParser::ATOM: {
        enterOuterAlt(_localctx, 2);
        setState(208);
        atomblock();
        break;
      }

      case OMDParser::BOND: {
        enterOuterAlt(_localctx, 3);
        setState(209);
        bondblock();
        break;
      }

      case OMDParser::BEND: {
        enterOuterAlt(_localctx, 4);
        setState(210);
        bendblock();
        break;
      }

      case OMDParser::TORSION: {
        enterOuterAlt(_localctx, 5);
        setState(211);
        torsionblock();
        break;
      }

      case OMDParser::INVERSION: {
        enterOuterAlt(_localctx, 6);
        setState(212);
        inversionblock();
        break;
      }

      case OMDParser::RIGIDBODY: {
        enterOuterAlt(_localctx, 7);
        setState(213);
        rigidbodyblock();
        break;
      }

      case OMDParser::CUTOFFGROUP: {
        enterOuterAlt(_localctx, 8);
        setState(214);
        cutoffgroupblock();
        break;
      }

      case OMDParser::CONSTRAINT: {
        enterOuterAlt(_localctx, 9);
        setState(215);
        constraintblock();
        break;
      }

      case OMDParser::SEQUENCE: {
        enterOuterAlt(_localctx, 10);
        setState(216);
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
  enterRule(_localctx, 28, OMDParser::RuleAtomblock);
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
    setState(219);
    match(OMDParser::ATOM);
    setState(220);
    match(OMDParser::LBRACKET);
    setState(221);
    intConst();
    setState(222);
    match(OMDParser::RBRACKET);
    setState(223);
    match(OMDParser::LCURLY);
    setState(227);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 21) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 21)) & 17592202821635) != 0)) {
      setState(224);
      atomstatement();
      setState(229);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(230);
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
  enterRule(_localctx, 30, OMDParser::RuleAtomstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(251);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(232);
        assignment();
        break;
      }

      case OMDParser::POSITION: {
        enterOuterAlt(_localctx, 2);
        setState(233);
        match(OMDParser::POSITION);
        setState(234);
        match(OMDParser::LPAREN);
        setState(235);
        doubleNumberTuple();
        setState(236);
        match(OMDParser::RPAREN);
        setState(237);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::ORIENTATION: {
        enterOuterAlt(_localctx, 3);
        setState(239);
        match(OMDParser::ORIENTATION);
        setState(240);
        match(OMDParser::LPAREN);
        setState(241);
        doubleNumberTuple();
        setState(242);
        match(OMDParser::RPAREN);
        setState(243);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CHARGE: {
        enterOuterAlt(_localctx, 4);
        setState(245);
        match(OMDParser::CHARGE);
        setState(246);
        match(OMDParser::LPAREN);
        setState(247);
        floatConst();
        setState(248);
        match(OMDParser::RPAREN);
        setState(249);
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
  enterRule(_localctx, 32, OMDParser::RuleBondblock);
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
    setState(253);
    match(OMDParser::BOND);
    setState(258);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(254);
      match(OMDParser::LBRACKET);
      setState(255);
      intConst();
      setState(256);
      match(OMDParser::RBRACKET);
    }
    setState(260);
    match(OMDParser::LCURLY);
    setState(264);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 18) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 18)) & 140737488419841) != 0)) {
      setState(261);
      bondstatement();
      setState(266);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(267);
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
  enterRule(_localctx, 34, OMDParser::RuleBondstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(312);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(269);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(270);
        match(OMDParser::MEMBERS);
        setState(271);
        match(OMDParser::LPAREN);
        setState(272);
        inttuple();
        setState(273);
        match(OMDParser::RPAREN);
        setState(274);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::FIXED: {
        enterOuterAlt(_localctx, 3);
        setState(276);
        match(OMDParser::FIXED);
        setState(277);
        match(OMDParser::LPAREN);
        setState(278);
        floatConst();
        setState(279);
        match(OMDParser::RPAREN);
        setState(280);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 4);
        setState(282);
        match(OMDParser::HARMONIC);
        setState(283);
        match(OMDParser::LPAREN);
        setState(284);
        doubleNumberTuple();
        setState(285);
        match(OMDParser::RPAREN);
        setState(286);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 5);
        setState(288);
        match(OMDParser::CUBIC);
        setState(289);
        match(OMDParser::LPAREN);
        setState(290);
        doubleNumberTuple();
        setState(291);
        match(OMDParser::RPAREN);
        setState(292);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 6);
        setState(294);
        match(OMDParser::QUARTIC);
        setState(295);
        match(OMDParser::LPAREN);
        setState(296);
        doubleNumberTuple();
        setState(297);
        match(OMDParser::RPAREN);
        setState(298);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 7);
        setState(300);
        match(OMDParser::POLYNOMIAL);
        setState(301);
        match(OMDParser::LPAREN);
        setState(302);
        doubleNumberTuple();
        setState(303);
        match(OMDParser::RPAREN);
        setState(304);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::MORSE: {
        enterOuterAlt(_localctx, 8);
        setState(306);
        match(OMDParser::MORSE);
        setState(307);
        match(OMDParser::LPAREN);
        setState(308);
        doubleNumberTuple();
        setState(309);
        match(OMDParser::RPAREN);
        setState(310);
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
  enterRule(_localctx, 36, OMDParser::RuleBendblock);
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
    setState(314);
    match(OMDParser::BEND);
    setState(319);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(315);
      match(OMDParser::LBRACKET);
      setState(316);
      intConst();
      setState(317);
      match(OMDParser::RBRACKET);
    }
    setState(321);
    match(OMDParser::LCURLY);
    setState(325);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 18) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 18)) & 140737488844801) != 0)) {
      setState(322);
      bendstatement();
      setState(327);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(328);
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
  enterRule(_localctx, 38, OMDParser::RuleBendstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(379);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(330);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(331);
        match(OMDParser::MEMBERS);
        setState(332);
        match(OMDParser::LPAREN);
        setState(333);
        inttuple();
        setState(334);
        match(OMDParser::RPAREN);
        setState(335);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 3);
        setState(337);
        match(OMDParser::HARMONIC);
        setState(338);
        match(OMDParser::LPAREN);
        setState(339);
        doubleNumberTuple();
        setState(340);
        match(OMDParser::RPAREN);
        setState(341);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::GHOSTBEND: {
        enterOuterAlt(_localctx, 4);
        setState(343);
        match(OMDParser::GHOSTBEND);
        setState(344);
        match(OMDParser::LPAREN);
        setState(345);
        doubleNumberTuple();
        setState(346);
        match(OMDParser::RPAREN);
        setState(347);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::UREYBRADLEY: {
        enterOuterAlt(_localctx, 5);
        setState(349);
        match(OMDParser::UREYBRADLEY);
        setState(350);
        match(OMDParser::LPAREN);
        setState(351);
        doubleNumberTuple();
        setState(352);
        match(OMDParser::RPAREN);
        setState(353);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 6);
        setState(355);
        match(OMDParser::CUBIC);
        setState(356);
        match(OMDParser::LPAREN);
        setState(357);
        doubleNumberTuple();
        setState(358);
        match(OMDParser::RPAREN);
        setState(359);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 7);
        setState(361);
        match(OMDParser::QUARTIC);
        setState(362);
        match(OMDParser::LPAREN);
        setState(363);
        doubleNumberTuple();
        setState(364);
        match(OMDParser::RPAREN);
        setState(365);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 8);
        setState(367);
        match(OMDParser::POLYNOMIAL);
        setState(368);
        match(OMDParser::LPAREN);
        setState(369);
        doubleNumberTuple();
        setState(370);
        match(OMDParser::RPAREN);
        setState(371);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::COSINE: {
        enterOuterAlt(_localctx, 9);
        setState(373);
        match(OMDParser::COSINE);
        setState(374);
        match(OMDParser::LPAREN);
        setState(375);
        doubleNumberTuple();
        setState(376);
        match(OMDParser::RPAREN);
        setState(377);
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
  enterRule(_localctx, 40, OMDParser::RuleTorsionblock);
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
    setState(381);
    match(OMDParser::TORSION);
    setState(386);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(382);
      match(OMDParser::LBRACKET);
      setState(383);
      intConst();
      setState(384);
      match(OMDParser::RBRACKET);
    }
    setState(388);
    match(OMDParser::LCURLY);
    setState(392);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 18) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 18)) & 140737496250369) != 0)) {
      setState(389);
      torsionstatement();
      setState(394);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(395);
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
  enterRule(_localctx, 42, OMDParser::RuleTorsionstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(452);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(397);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(398);
        match(OMDParser::MEMBERS);
        setState(399);
        match(OMDParser::LPAREN);
        setState(400);
        inttuple();
        setState(401);
        match(OMDParser::RPAREN);
        setState(402);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::GHOSTTORSION: {
        enterOuterAlt(_localctx, 3);
        setState(404);
        match(OMDParser::GHOSTTORSION);
        setState(405);
        match(OMDParser::LPAREN);
        setState(406);
        doubleNumberTuple();
        setState(407);
        match(OMDParser::RPAREN);
        setState(408);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CUBIC: {
        enterOuterAlt(_localctx, 4);
        setState(410);
        match(OMDParser::CUBIC);
        setState(411);
        match(OMDParser::LPAREN);
        setState(412);
        doubleNumberTuple();
        setState(413);
        match(OMDParser::RPAREN);
        setState(414);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::QUARTIC: {
        enterOuterAlt(_localctx, 5);
        setState(416);
        match(OMDParser::QUARTIC);
        setState(417);
        match(OMDParser::LPAREN);
        setState(418);
        doubleNumberTuple();
        setState(419);
        match(OMDParser::RPAREN);
        setState(420);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::POLYNOMIAL: {
        enterOuterAlt(_localctx, 6);
        setState(422);
        match(OMDParser::POLYNOMIAL);
        setState(423);
        match(OMDParser::LPAREN);
        setState(424);
        doubleNumberTuple();
        setState(425);
        match(OMDParser::RPAREN);
        setState(426);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CHARMM: {
        enterOuterAlt(_localctx, 7);
        setState(428);
        match(OMDParser::CHARMM);
        setState(429);
        match(OMDParser::LPAREN);
        setState(430);
        doubleNumberTuple();
        setState(431);
        match(OMDParser::RPAREN);
        setState(432);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::OPLS: {
        enterOuterAlt(_localctx, 8);
        setState(434);
        match(OMDParser::OPLS);
        setState(435);
        match(OMDParser::LPAREN);
        setState(436);
        doubleNumberTuple();
        setState(437);
        match(OMDParser::RPAREN);
        setState(438);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::TRAPPE: {
        enterOuterAlt(_localctx, 9);
        setState(440);
        match(OMDParser::TRAPPE);
        setState(441);
        match(OMDParser::LPAREN);
        setState(442);
        doubleNumberTuple();
        setState(443);
        match(OMDParser::RPAREN);
        setState(444);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 10);
        setState(446);
        match(OMDParser::HARMONIC);
        setState(447);
        match(OMDParser::LPAREN);
        setState(448);
        doubleNumberTuple();
        setState(449);
        match(OMDParser::RPAREN);
        setState(450);
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
  enterRule(_localctx, 44, OMDParser::RuleInversionblock);
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
    setState(454);
    match(OMDParser::INVERSION);
    setState(459);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(455);
      match(OMDParser::LBRACKET);
      setState(456);
      intConst();
      setState(457);
      match(OMDParser::RBRACKET);
    }
    setState(461);
    match(OMDParser::LCURLY);
    setState(465);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 19) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 19)) & 70368807093251) != 0)) {
      setState(462);
      inversionstatement();
      setState(467);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(468);
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
  enterRule(_localctx, 46, OMDParser::RuleInversionstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(513);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(470);
        assignment();
        break;
      }

      case OMDParser::CENTER: {
        enterOuterAlt(_localctx, 2);
        setState(471);
        match(OMDParser::CENTER);
        setState(472);
        match(OMDParser::LPAREN);
        setState(473);
        intConst();
        setState(474);
        match(OMDParser::RPAREN);
        setState(475);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::SATELLITES: {
        enterOuterAlt(_localctx, 3);
        setState(477);
        match(OMDParser::SATELLITES);
        setState(478);
        match(OMDParser::LPAREN);
        setState(479);
        inttuple();
        setState(480);
        match(OMDParser::RPAREN);
        setState(481);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::AMBERIMPROPER: {
        enterOuterAlt(_localctx, 4);
        setState(483);
        match(OMDParser::AMBERIMPROPER);
        setState(484);
        match(OMDParser::LPAREN);
        setState(485);
        doubleNumberTuple();
        setState(486);
        match(OMDParser::RPAREN);
        setState(487);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::IMPROPERCOSINE: {
        enterOuterAlt(_localctx, 5);
        setState(489);
        match(OMDParser::IMPROPERCOSINE);
        setState(490);
        match(OMDParser::LPAREN);
        setState(491);
        doubleNumberTuple();
        setState(492);
        match(OMDParser::RPAREN);
        setState(493);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::HARMONIC: {
        enterOuterAlt(_localctx, 6);
        setState(495);
        match(OMDParser::HARMONIC);
        setState(496);
        match(OMDParser::LPAREN);
        setState(497);
        doubleNumberTuple();
        setState(498);
        match(OMDParser::RPAREN);
        setState(499);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::CENTRALATOMHEIGHT: {
        enterOuterAlt(_localctx, 7);
        setState(501);
        match(OMDParser::CENTRALATOMHEIGHT);
        setState(502);
        match(OMDParser::LPAREN);
        setState(503);
        doubleNumberTuple();
        setState(504);
        match(OMDParser::RPAREN);
        setState(505);
        match(OMDParser::SEMICOLON);
        break;
      }

      case OMDParser::DREIDING: {
        enterOuterAlt(_localctx, 8);
        setState(507);
        match(OMDParser::DREIDING);
        setState(508);
        match(OMDParser::LPAREN);
        setState(509);
        doubleNumberTuple();
        setState(510);
        match(OMDParser::RPAREN);
        setState(511);
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
  enterRule(_localctx, 48, OMDParser::RuleRigidbodyblock);
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
    setState(515);
    match(OMDParser::RIGIDBODY);
    setState(516);
    match(OMDParser::LBRACKET);
    setState(517);
    intConst();
    setState(518);
    match(OMDParser::RBRACKET);
    setState(519);
    match(OMDParser::LCURLY);
    setState(523);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(520);
      rigidbodystatement();
      setState(525);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(526);
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
  enterRule(_localctx, 50, OMDParser::RuleRigidbodystatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(535);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(528);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(529);
        match(OMDParser::MEMBERS);
        setState(530);
        match(OMDParser::LPAREN);
        setState(531);
        inttuple();
        setState(532);
        match(OMDParser::RPAREN);
        setState(533);
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
  enterRule(_localctx, 52, OMDParser::RuleCutoffgroupblock);
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
    setState(537);
    match(OMDParser::CUTOFFGROUP);
    setState(542);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(538);
      match(OMDParser::LBRACKET);
      setState(539);
      intConst();
      setState(540);
      match(OMDParser::RBRACKET);
    }
    setState(544);
    match(OMDParser::LCURLY);
    setState(548);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(545);
      cutoffgroupstatement();
      setState(550);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(551);
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
  enterRule(_localctx, 54, OMDParser::RuleCutoffgroupstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(560);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(553);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(554);
        match(OMDParser::MEMBERS);
        setState(555);
        match(OMDParser::LPAREN);
        setState(556);
        inttuple();
        setState(557);
        match(OMDParser::RPAREN);
        setState(558);
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
  enterRule(_localctx, 56, OMDParser::RuleNodesblock);
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
    setState(562);
    match(OMDParser::NODES);
    setState(567);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(563);
      match(OMDParser::LBRACKET);
      setState(564);
      intConst();
      setState(565);
      match(OMDParser::RBRACKET);
    }
    setState(569);
    match(OMDParser::LCURLY);
    setState(573);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(570);
      nodesstatement();
      setState(575);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(576);
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
  enterRule(_localctx, 58, OMDParser::RuleNodesstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(585);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(578);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(579);
        match(OMDParser::MEMBERS);
        setState(580);
        match(OMDParser::LPAREN);
        setState(581);
        inttuple();
        setState(582);
        match(OMDParser::RPAREN);
        setState(583);
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
  enterRule(_localctx, 60, OMDParser::RuleFragmentblock);
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
    setState(587);
    match(OMDParser::FRAGMENT);
    setState(588);
    match(OMDParser::LCURLY);
    setState(592);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (((((_la - 7) & ~ 0x3fULL) == 0) &&
      ((1ULL << (_la - 7)) & 288230925907525887) != 0)) {
      setState(589);
      fragmentstatement();
      setState(594);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(595);
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
  enterRule(_localctx, 62, OMDParser::RuleFragmentstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(607);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(597);
        assignment();
        break;
      }

      case OMDParser::ATOM: {
        enterOuterAlt(_localctx, 2);
        setState(598);
        atomblock();
        break;
      }

      case OMDParser::BOND: {
        enterOuterAlt(_localctx, 3);
        setState(599);
        bondblock();
        break;
      }

      case OMDParser::BEND: {
        enterOuterAlt(_localctx, 4);
        setState(600);
        bendblock();
        break;
      }

      case OMDParser::TORSION: {
        enterOuterAlt(_localctx, 5);
        setState(601);
        torsionblock();
        break;
      }

      case OMDParser::INVERSION: {
        enterOuterAlt(_localctx, 6);
        setState(602);
        inversionblock();
        break;
      }

      case OMDParser::RIGIDBODY: {
        enterOuterAlt(_localctx, 7);
        setState(603);
        rigidbodyblock();
        break;
      }

      case OMDParser::CUTOFFGROUP: {
        enterOuterAlt(_localctx, 8);
        setState(604);
        cutoffgroupblock();
        break;
      }

      case OMDParser::CONSTRAINT: {
        enterOuterAlt(_localctx, 9);
        setState(605);
        constraintblock();
        break;
      }

      case OMDParser::NODES: {
        enterOuterAlt(_localctx, 10);
        setState(606);
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
  enterRule(_localctx, 64, OMDParser::RuleConstraintblock);
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
    setState(609);
    match(OMDParser::CONSTRAINT);
    setState(614);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == OMDParser::LBRACKET) {
      setState(610);
      match(OMDParser::LBRACKET);
      setState(611);
      intConst();
      setState(612);
      match(OMDParser::RBRACKET);
    }
    setState(616);
    match(OMDParser::LCURLY);
    setState(620);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::MEMBERS

    || _la == OMDParser::ID) {
      setState(617);
      constraintstatement();
      setState(622);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(623);
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
  enterRule(_localctx, 66, OMDParser::RuleConstraintstatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(632);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(625);
        assignment();
        break;
      }

      case OMDParser::MEMBERS: {
        enterOuterAlt(_localctx, 2);
        setState(626);
        match(OMDParser::MEMBERS);
        setState(627);
        match(OMDParser::LPAREN);
        setState(628);
        inttuple();
        setState(629);
        match(OMDParser::RPAREN);
        setState(630);
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
  enterRule(_localctx, 68, OMDParser::RuleSequencestring);

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
    match(OMDParser::SEQUENCE);
    setState(635);
    match(OMDParser::ASSIGNEQUAL);
    setState(636);
    constant();
    setState(637);
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
  enterRule(_localctx, 70, OMDParser::RuleDoubleNumberTuple);
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
    setState(639);
    doubleNumber();
    setState(644);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::COMMA) {
      setState(640);
      match(OMDParser::COMMA);
      setState(641);
      doubleNumber();
      setState(646);
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
  enterRule(_localctx, 72, OMDParser::RuleInttuple);
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
    setState(647);
    intConst();
    setState(652);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == OMDParser::COMMA) {
      setState(648);
      match(OMDParser::COMMA);
      setState(649);
      intConst();
      setState(654);
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
  enterRule(_localctx, 74, OMDParser::RuleIntConst);
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
    setState(655);
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
  enterRule(_localctx, 76, OMDParser::RuleDoubleNumber);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(659);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case OMDParser::NUM_LONG:
      case OMDParser::NUM_INT: {
        enterOuterAlt(_localctx, 1);
        setState(657);
        intConst();
        break;
      }

      case OMDParser::NUM_FLOAT:
      case OMDParser::NUM_DOUBLE: {
        enterOuterAlt(_localctx, 2);
        setState(658);
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
  enterRule(_localctx, 78, OMDParser::RuleFloatConst);
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
    setState(661);
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
  enterRule(_localctx, 80, OMDParser::RuleVectorConst);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(663);
    match(OMDParser::LPAREN);
    setState(664);
    doubleNumberTuple();
    setState(665);
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
