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
 * @file BondType.hpp
 * @author    tlin
 * @date  11/01/2004
 * @version 1.0
 */ 

#ifndef TYPES_CUBICBONDTYPE_HPP
#define TYPES_CUBICBONDTYPE_HPP

#include "types/BondType.hpp"

namespace oopse {
/**
 * @class CubicBondType 
 * @todo document
 */
class CubicBondType : public BondType {

    public:

        CubicBondType(double r0, double k3, double k2, double k1, double k0) 
            : BondType(r0), k3_(k3), k2_(k2),  k1_(k1), k0_(k0){
        }

        void setForceConstant(double k3, double k2, double k1, double k0) {
            k3_ = k3;
            k2_ = k2;
            k1_ = k1;
            k0_ = k0;

        }

        void getForceConstant(double& k3, double& k2, double& k1, double& k0) {
            k3 = k3_;
            k2  = k2_;
            k1 = k1_;
            k0 = k0_;
        }

        virtual void calcForce(double r, double& V, double& dVdr) {
            double dr =  r- r0;
            double dr2 = dr * dr;
            double dr3 = dr2 * dr;
            
            V =k0_ + k1_ * dr + k2_*dr2 + k3_*dr3;
            dVdr = k1_ + 2.0*k2_ * dr + 3.0 * k3_*dr2;            
        }
        
    private:

        double k3_;
        double k2_;
        double k1_;
        double k0_;

};

}//end namespace oopse
#endif //TYPES_CUBICBONDTYPE_HPP
