// Example usage of ANTLR v4 OMDParser
// This file shows how to integrate the new parser into OpenMD

#include "antlr4-runtime.h"
#include "OMDLexer.h"
#include "OMDParser.h"
#include "OMDTreeVisitor.hpp"
#include "io/Globals.hpp"
#include <fstream>
#include <iostream>

// Custom error listener for better error reporting
class OMDErrorListener : public antlr4::BaseErrorListener {
public:
    void syntaxError(
        antlr4::Recognizer *recognizer,
        antlr4::Token *offendingSymbol,
        size_t line,
        size_t charPositionInLine,
        const std::string &msg,
        std::exception_ptr e
    ) override {
        std::cerr << "Syntax error at line " << line 
                  << ":" << charPositionInLine 
                  << " - " << msg << std::endl;
        
        // You can throw an exception here if you want to halt parsing
        // throw std::runtime_error("Parse error");
    }
};

// Function to parse an OpenMD file and return Globals configuration
Globals* parseOmdFile(const std::string& filename) {
    try {
        // Open the file
        std::ifstream stream(filename);
        if (!stream.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return nullptr;
        }
        
        // Create ANTLR input stream from file
        antlr4::ANTLRInputStream input(stream);
        
        // Create lexer
        OMDLexer lexer(&input);
        
        // Optional: Set up filename observer for preprocessor directives
        // FilenameObserver* observer = new FilenameObserver();
        // lexer.setObserver(observer);
        
        // Create token stream
        antlr4::CommonTokenStream tokens(&lexer);
        
        // Create parser
        OMDParser parser(&tokens);
        
        // Optional: Add custom error listener
        OMDErrorListener errorListener;
        parser.removeErrorListeners();  // Remove default console error listener
        parser.addErrorListener(&errorListener);
        
        // Parse the file (top-level rule is 'omdfile')
        OMDParser::OmdfileContext* tree = parser.omdfile();
        
        // Check for parse errors
        if (parser.getNumberOfSyntaxErrors() > 0) {
            std::cerr << "Parsing failed with " 
                      << parser.getNumberOfSyntaxErrors() 
                      << " errors" << std::endl;
            return nullptr;
        }
        
        // Create visitor and walk the parse tree
        OMDTreeVisitor visitor;
        Globals* conf = visitor.walkTree(tree);
        
        return conf;
        
    } catch (const std::exception& e) {
        std::cerr << "Exception during parsing: " << e.what() << std::endl;
        return nullptr;
    }
}

// Function to parse from a string (useful for testing)
Globals* parseOmdString(const std::string& input) {
    try {
        // Create ANTLR input stream from string
        antlr4::ANTLRInputStream inputStream(input);
        
        // Create lexer
        OMDLexer lexer(&inputStream);
        
        // Create token stream
        antlr4::CommonTokenStream tokens(&lexer);
        
        // Create parser
        OMDParser parser(&tokens);
        
        // Parse
        OMDParser::OmdfileContext* tree = parser.omdfile();
        
        // Visit tree
        OMDTreeVisitor visitor;
        Globals* conf = visitor.walkTree(tree);
        
        return conf;
        
    } catch (const std::exception& e) {
        std::cerr << "Exception during parsing: " << e.what() << std::endl;
        return nullptr;
    }
}

// Example main function
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.omd>" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    
    std::cout << "Parsing file: " << filename << std::endl;
    
    Globals* config = parseOmdFile(filename);
    
    if (config) {
        std::cout << "Parsing successful!" << std::endl;
        // Use config...
        
        // Don't forget to delete when done
        delete config;
        return 0;
    } else {
        std::cerr << "Parsing failed!" << std::endl;
        return 1;
    }
}

// ============================================================================
// Advanced Usage Examples
// ============================================================================

// Example: Custom visitor for specific processing
class CustomOMDVisitor : public OMDTreeVisitor {
public:
    int moleculeCount = 0;
    int atomCount = 0;
    
    virtual antlrcpp::Any visitMoleculeblock(OMDParser::MoleculeblockContext *ctx) override {
        moleculeCount++;
        std::cout << "Found molecule #" << moleculeCount << std::endl;
        return OMDTreeVisitor::visitMoleculeblock(ctx);
    }
    
    virtual antlrcpp::Any visitAtomblock(OMDParser::AtomblockContext *ctx) override {
        atomCount++;
        std::cout << "Found atom #" << atomCount << std::endl;
        return OMDTreeVisitor::visitAtomblock(ctx);
    }
};

// Example: Print parse tree structure
void printParseTree(antlr4::tree::ParseTree* tree, const OMDParser& parser, int indent = 0) {
    std::string indentStr(indent * 2, ' ');
    
    if (auto* terminalNode = dynamic_cast<antlr4::tree::TerminalNode*>(tree)) {
        // Leaf node (token)
        std::cout << indentStr << "TOKEN: " << terminalNode->getText() << std::endl;
    } else if (auto* ruleNode = dynamic_cast<antlr4::RuleContext*>(tree)) {
        // Rule node
        std::string ruleName = parser.getRuleNames()[ruleNode->getRuleIndex()];
        std::cout << indentStr << "RULE: " << ruleName << std::endl;
        
        // Recursively print children
        for (size_t i = 0; i < ruleNode->children.size(); i++) {
            printParseTree(ruleNode->children[i], parser, indent + 1);
        }
    }
}

// Example: Validate without processing
bool validateOmdFile(const std::string& filename) {
    std::ifstream stream(filename);
    if (!stream.is_open()) {
        return false;
    }
    
    antlr4::ANTLRInputStream input(stream);
    OMDLexer lexer(&input);
    antlr4::CommonTokenStream tokens(&lexer);
    OMDParser parser(&tokens);
    
    // Just parse, don't process
    parser.omdfile();
    
    return parser.getNumberOfSyntaxErrors() == 0;
}

// Example: Get detailed error information
struct ParseError {
    size_t line;
    size_t column;
    std::string message;
};

class DetailedErrorListener : public antlr4::BaseErrorListener {
public:
    std::vector<ParseError> errors;
    
    void syntaxError(
        antlr4::Recognizer *recognizer,
        antlr4::Token *offendingSymbol,
        size_t line,
        size_t charPositionInLine,
        const std::string &msg,
        std::exception_ptr e
    ) override {
        ParseError error;
        error.line = line;
        error.column = charPositionInLine;
        error.message = msg;
        errors.push_back(error);
    }
};

std::vector<ParseError> getParseErrors(const std::string& filename) {
    std::ifstream stream(filename);
    antlr4::ANTLRInputStream input(stream);
    OMDLexer lexer(&input);
    antlr4::CommonTokenStream tokens(&lexer);
    OMDParser parser(&tokens);
    
    DetailedErrorListener errorListener;
    parser.removeErrorListeners();
    parser.addErrorListener(&errorListener);
    
    parser.omdfile();
    
    return errorListener.errors;
}
