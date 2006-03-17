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
#include "applications/hydrodynamics/AnalyticalModel.hpp"
#include "applications/hydrodynamics/Spheric.hpp"
#include "applications/hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"
#include "math/LU.hpp"
namespace oopse {
bool AnalyticalModel::calcHydroProps(Spheric* spheric, double viscosity, double temperature) {

    double radius = spheric->getRadius(); 
    HydroProps props;
    props.center =V3Zero;
    double Xitt  = 6.0 * NumericConstant::PI * viscosity * radius;
    double Xirr = 8.0 * NumericConstant::PI * viscosity * radius * radius * radius;
    props.Xi(0, 0) = Xitt;
    props.Xi(1, 1) = Xitt;
    props.Xi(2, 2) = Xitt;
    props.Xi(3, 3) = Xirr;
    props.Xi(4, 4) = Xirr;
    props.Xi(5, 5) = Xirr;
    
    const double convertConstant = 6.023; //convert poise.angstrom to amu/fs
    props.Xi *= convertConstant;
    Mat6x6d XiCopy = props.Xi;
    invertMatrix(XiCopy, props.D);
    double kt = OOPSEConstant::kB * temperature;
    props.D *= kt;
    props.Xi *= OOPSEConstant::kb * temperature;

    setCR(props);
    setCD(props);

    return true;
    
}

/**
 * calculate the ratio of friction coeffiction constant between ellipsoid and spheric 
 * with same volume.
 * @param m
 * @param n 
 * @note 
 * Reference:
 *
 * (1) Victor A. Bloomfield, On-Line Biophysics Textbook, Volume: Separations and Hydrodynamics
 * Chapter 1,Survey of Biomolecular Hydrodynamics
 * http://www.biophysics.org/education/vbloomfield.pdf 
 * (2) F. Perrin , J. Phys. Radium, [7] 5, 497-511, 1934
 * (3) F. Perrin, J. Phys. Radium, [7] 7, 1-11, 1936
 */        
bool AnalyticalModel::calcHydroProps(Ellipsoid* ellipsoid, double viscosity, double temperature) {
    double ft;
    double fra;
    double frb;
    double a = ellipsoid->getA();
    double b = ellipsoid->getB();
    double q = a/b; //?
    if (q > 1.0) {//prolate
        ft = sqrt(1-q*q)/(pow(q, 2.0/3.0)*log((1 + sqrt(1-q*q))/q));
        fra = 4*(1-q*q)/(3*(2 - 2*pow(q, 4.0/3.0)/ft)); //not sure
        frb = 4*(1-q*q*q*q) /(3*q*q*(2*pow(q, -2.0/3.0)*(2-q*q)/ft-2));
    } else {//oblate
        ft = sqrt(1-q*q)/(pow(q, 2.0/3.0)*atan(sqrt(q*q-1)));
        fra = 4*(1-q*q)/(3*(2 - 2*pow(q, 4.0/3.0)/ft)); //not sure
        frb = 4*(1-q*q*q*q) /(3*q*q*(2*pow(q, -2.0/3.0)*(2-q*q)/ft-2));
    }
                    
    double radius = pow(a*a*b, 1.0/3.0);
    HydroProps props;
    double Xitt  = 6.0 * NumericConstant::PI * viscosity * radius;
    double Xirr = 8.0 * NumericConstant::PI * viscosity * radius * radius * radius;
    props.Xi(0, 0) = Xitt;
    props.Xi(1, 1) = Xitt;
    props.Xi(2, 2) = Xitt;
    props.Xi(3, 3) = Xirr;
    props.Xi(4, 4) = Xirr;
    props.Xi(5, 5) = Xirr;
    
    const double convertConstant = 6.023; //convert poise.angstrom to amu/fs
    props.Xi *= convertConstant;    
    props.Xi(0,0) *= ft;
    props.Xi(1,1) *= ft;
    props.Xi(2,2) *= ft;
    props.Xi(3,3) *= fra;
    props.Xi(4,4) *= fra;
    props.Xi(5,5) *= frb;
    
    Mat6x6d XiCopy = props.Xi;
    XiCopy /= OOPSEConstant::kb * temperature;
    invertMatrix(XiCopy, props.D);
    double kt = OOPSEConstant::kB * temperature;
    props.D *= kt;

    setCR(props);
    setCD(props);

    return true;
}

bool AnalyticalModel::calcHydroProps(CompositeShape* compositexShape, double viscosity, double temperature) {
    return false;
}
        


}
