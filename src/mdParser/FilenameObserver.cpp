#include "FilenameObserver.hpp"
#include "MDLexer.hpp"
#include "MDParser.hpp"

FilenameObserver::FilenameObserver() : parser_(NULL), lexer_(NULL) {}
void FilenameObserver::setParser(MDParser* parser) {parser_ = parser;}
void FilenameObserver::setLexer(MDLexer* lexer) {lexer_ = lexer;}
void FilenameObserver::notify(const std::string& filename) {
    if (lexer_)
        lexer_->setFilename(filename);
    if (parser_)
        parser_->setFilename(filename);
}
