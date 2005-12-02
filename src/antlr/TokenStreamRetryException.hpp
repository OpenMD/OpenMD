#ifndef INC_TokenStreamRetryException_hpp__
#define INC_TokenStreamRetryException_hpp__

/* ANTLR Translator Generator
 * Project led by Terence Parr at http://www.jGuru.com
 * Software rights: http://www.antlr.org/license.html
 *
 * $Id: TokenStreamRetryException.hpp,v 1.1 2005-12-02 15:38:02 tim Exp $
 */

#include <antlr/config.hpp>
#include <antlr/TokenStreamException.hpp>

#ifdef ANTLR_CXX_SUPPORTS_NAMESPACE
namespace antlr {
#endif

class TokenStreamRetryException : public TokenStreamException {
public:
	TokenStreamRetryException() {}
	~TokenStreamRetryException() throw() {}
};

#ifdef ANTLR_CXX_SUPPORTS_NAMESPACE
}
#endif

#endif //INC_TokenStreamRetryException_hpp__
