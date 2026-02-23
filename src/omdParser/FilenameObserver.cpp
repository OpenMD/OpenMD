#include "FilenameObserver.hpp"
#include "OMDLexer.h"
#include "OMDParser.h"

FilenameObserver::FilenameObserver() : parser_(NULL), lexer_(NULL) {}
void FilenameObserver::setParser(OMDParser* parser) {parser_ = parser;}
void FilenameObserver::setLexer(OMDLexer* lexer) {lexer_ = lexer;}
void FilenameObserver::notify(const std::string& filename) {
  /*    if (lexer_)
        lexer_->setFilename(filename);
    if (parser_)
        parser_->setFilename(filename);
  */
}
