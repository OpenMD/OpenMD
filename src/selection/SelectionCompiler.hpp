/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */

#ifndef SELECTION_SELECTIONCOMPILER_HPP
#define SELECTION_SELECTIONCOMPILER_HPP

#include <vector>
namespace oopse {


/**
 * @class SelectionCompiler SelectionCompiler.hpp "selection/SelectionCompiler.hpp"
 * @brief compile a selection script to tokens
 * <pre>

    expression       :: = clauseOr

    clauseOr         ::= clauseAnd {OR clauseAnd}*

    clauseAnd        ::= clauseNot {AND clauseNot}*

    clauseNot        ::= NOT clauseNot | clausePrimitive

    clausePrimitive  ::= clauseComparator |
                         clauseWithin |
                         clauseResidueSpec |
                         none | all |
                         ( clauseOr )

    clauseComparator ::= atomproperty comparatorop integer

    clauseWithin     ::= WITHIN ( clauseDistance , expression )

    clauseDistance   ::= integer | decimal

    


    clauseResidueSpec::= { clauseResNameSpec }
                         { clauseResNumSpec }
                         { chainSpec }
                         { clauseAtomSpec }
                         { modelSpec }

    clauseResNameSpec::= * | [ resNamePattern ] | resNamePattern

    // question marks are part of identifiers
    // they get split up and dealt with as wildcards at runtime
    // and the integers which are residue number chains get bundled
    // in with the identifier and also split out at runtime
    // iff a variable of that name does not exist

    resNamePattern   ::= up to 3 alphanumeric chars with * and ?

    clauseResNumSpec ::= * | clauseSequenceRange

    clauseSequenceRange ::= clauseSequenceCode { - clauseSequenceCode }

    clauseSequenceCode ::= seqcode | {-} integer

    clauseChainSpec  ::= {:} * | identifier | integer

    clauseAtomSpec   ::= . * | . identifier {*} // note that this * is *not* a wildcard

    clauseModelSpec  ::= {:|/} * | integer
 
 * </pre>
 */
class SelectionCompiler{
    public:
        bool compile();

        std::vector<Token> getCompiledTokens();
    private:

        bool clauseOr();
        bool clauseAnd();
        bool clauseNot();
        bool clausePrimitive();
        bool clauseWithin();
        bool clauseComparator();
        
        internalCompile();

        std::vector<Token> compiledTokens_;
};

}
#endif
