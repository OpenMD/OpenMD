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
  * @file NVT.hpp
  * @author tlin
  * @date 11/19/2004
  * @time 13:25am
  * @version 1.0
  */

#ifndef INTEGRATOR_NVT_HPP
#define INTEGRATOR_NVT_HPP

#include "integrators/VelocityVerletIntegrator.hpp"

namespace oopse {

/**
 * @class NVT NVT.hpp "integrators/NVT.hpp"
 * Basic thermostating via Hoover, Phys.Rev.A, 1985, Vol. 31 (5) 1695-1697
 * @todo document
*/
class NVT : public VelocityVerletIntegrator{
    public:
        NVT(SimInfo* info);

        int getMaxIterationNumber() {
            return maxIterNum_;
        }
        
        void setMaxIterationNumber(int maxIter) {
            maxIterNum_ = maxIter;
        }

        double getTauThermostat() {
            return tauThermostat_;
        }

        void setTauThermostat(double tt) {
            tauThermostat_ = tt;
        }

        double getTargetTemp() {
            return targetTemp_;
        }

        void setTargetTemp(double tt) {
            targetTemp_ = tt;
        }

        double getChiTolerance() {
            return chiTolerance_;
        }

        void setChiTolerance(double tol) {
            chiTolerance_ = tol;
        }

        
    protected:
        virtual void moveA();

        virtual void moveB();

        virtual void doUpdate() ;
        
    private:
        virtual double calcConservedQuantity();              

        int maxIterNum_;
        double targetTemp_;
        double tauThermostat_;
        double chiTolerance_;

        std::vector<Vector3d> oldVel_;
        std::vector<Vector3d> oldJi_;
 };


} //end namespace oopse

#endif //INTEGRATOR_NVT_HPP
