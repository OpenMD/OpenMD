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
 
/**
 * @file PolynomialBondType.hpp
 * @author    teng lin
 * @date  11/16/2004
 * @version 1.0
 */ 

#ifndef TYPES_POLYNOMIALBONDTYPE_HPP
#define TYPES_POLYNOMIALBONDTYPE_HPP


#include "math/Polynomial.hpp"
#include "types/BondType.hpp"

namespace oopse {

/**
 * @class PolynomialBondType PolynomialBondType.hpp "types/PolynomialBondType.hpp"
 * @todo documentation
 */
class PolynomialBondType : public BondType{

    public:
        PolynomialBondType(double r0) : BondType(r0) {}

        void setCoefficient(int power, double coefficient) {
            polynomial_.setCoefficient(power, coefficient);
        }

        double getCoefficient(int power) {
            return polynomial_.getCoefficient(power);
        }
        
        void calcForce(double r, double & V, double & dVdr) {
            double delta = r - r0;
            V = polynomial_.evaluate(delta);
            dVdr = polynomial_.evaluateDerivative(delta);

        }

        friend std::ostream& operator <<(std::ostream& os, PolynomialBondType& pbt);
        
    private:
        
        DoublePolynomial polynomial_;
};

std::ostream& operator <<(std::ostream& os, PolynomialBondType& pbt) {
   DoublePolynomial::const_iterator i;

    i = pbt.polynomial_.begin();
    
    if (i == pbt.polynomial_.end()) {
        os << "This PolynomialBondType contains nothing" << std::endl;
        return os;
    }

    os << "This PolynomialBondType contains below terms:" << std::endl;    
    
    while(true){
        os << i->second << "*" << "(r - " << pbt.getEquilibriumBondLength() << ")" << "^" << i->first;

        if (++i == pbt.polynomial_.end()) {
            //if we reach the end of the polynomial pair, write out a newline and then escape the loop
            os << std::endl;
            break;
        } else {
            //otherwise, write out a "+"
            os << " + ";
        }
    }
    
    return os;
}


} //end namespace oopse
#endif //TYPES_POLYNOMIALBONDTYPE_HPP
