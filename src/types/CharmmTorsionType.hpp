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
 * @file TorsionType.hpp
 * @author teng lin
 * @date  11/16/2004
 * @version 1.0
 */ 

#ifndef TYPES_CHARMMTORSIONTYPE_HPP
#define TYPES_CHARMMTORSIONTYPE_HPP
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "types/TorsionType.hpp"
namespace oopse {

struct CharmmTorsionParameter {
    double kchi;
    int n;
    double delta;
};

class SamePeriodicityFunctor {
    public:
        SamePeriodicityFunctor(int n) : n_(n) {}
        bool operator()(CharmmTorsionParameter p) {
            return p.n == n_;
        }
        
   private:
        int n_;
};
/**
 * @class CharmmTorsionType CharmmTorsionType.hpp "types/CharmmTorsionType.hpp"
 */
class CharmmTorsionType : public TorsionType{
    public:

        void setCharmmTorsionParameter(CharmmTorsionParameter param) {

            assert(param.n >= 0 && param.n < 6);

            std::vector<CharmmTorsionParameter>::iterator i ;
            i = std::find_if(parameter_.begin(), parameter_.end(), SamePeriodicityFunctor(param.n));

            if (i != parameter_.end()) {
                std::cerr << "a parameter set with " << param.n <<" is already there" << std::endl;
            } else {
                parameter_.push_back(param);
            }            
        }

        void setCharmmTorsionParameter(double kchi, int n, double delta) {
            CharmmTorsionParameter param;
            param.kchi = kchi;
            param.n = n;
            param.delta = delta;
            setCharmmTorsionParameter(param);
        }
        
        virtual void calcForce(double cosPhi, double sinPhi, double& V, double& dVdPhi);

    private:
        std::vector<CharmmTorsionParameter> parameter_;
};
  
} //end namespace oopse
#endif //TYPES_CHARMMTORSIONTYPE_HPP

