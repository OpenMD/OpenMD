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
#include "applications/hydrodynamics/HydrodynamicsModel.hpp"
#include "applications/hydrodynamics/Spheric.hpp"
#include "applications/hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"

namespace oopse {

bool HydrodynamicsModel::calcHydroProps(Spheric* spheric, double viscosity, double temperature) {
    return false;
}

bool HydrodynamicsModel::calcHydroProps(Ellipsoid* ellipsoid, double viscosity, double temperature) {
    return false;
}

bool HydrodynamicsModel::calcHydroProps(CompositeShape* compositexShape, double viscosity, double temperature) {
    return false;
}

void HydrodynamicsModel::writeHydroProps(std::ostream& os) {

    
    os << sd_->getType() << "\t";
    
    //center of resistance
    os << cr_.center[0] <<  "\t" << cr_.center[1] <<  "\t" << cr_.center[2] <<  "\t"; 

    //resistance tensor at center of resistance
    //translation
    os << cr_.Xi(0, 0) <<  "\t" << cr_.Xi(0, 1) <<  "\t" << cr_.Xi(0, 2) <<  "\t" 
        << cr_.Xi(1, 0) <<  "\t" << cr_.Xi(1, 1) <<  "\t" << cr_.Xi(1, 2) <<  "\t" 
        << cr_.Xi(2, 0) <<  "\t" << cr_.Xi(2, 1) <<  "\t" << cr_.Xi(2, 2) <<  "\t"; 

    //rotation-translation
    os << cr_.Xi(0, 3) <<  "\t" << cr_.Xi(0, 4) <<  "\t" << cr_.Xi(0, 5) <<  "\t" 
        << cr_.Xi(1, 3) <<  "\t" << cr_.Xi(1, 4) <<  "\t" << cr_.Xi(1, 5) <<  "\t" 
        << cr_.Xi(2, 3) <<  "\t" << cr_.Xi(2, 4) <<  "\t" << cr_.Xi(2, 5) <<  "\t";

    //translation-rotation
    os << cr_.Xi(3, 0) <<  "\t" << cr_.Xi(3, 1) <<  "\t" << cr_.Xi(3, 2) <<  "\t" 
        << cr_.Xi(4, 0) <<  "\t" << cr_.Xi(4, 1) <<  "\t" << cr_.Xi(4, 2) <<  "\t" 
        << cr_.Xi(5, 0) <<  "\t" << cr_.Xi(5, 1) <<  "\t" << cr_.Xi(5, 2) <<  "\t"; 

    //rotation
    os << cr_.Xi(3, 3) <<  "\t" << cr_.Xi(3, 4) <<  "\t" << cr_.Xi(3, 5) <<  "\t" 
        << cr_.Xi(4, 3) <<  "\t" << cr_.Xi(4, 4) <<  "\t" << cr_.Xi(4, 5) <<  "\t" 
        << cr_.Xi(5, 3) <<  "\t" << cr_.Xi(5, 4) <<  "\t" << cr_.Xi(5, 5) <<  "\t"; 


    //diffusion tensor at center of resistance
    //translation
    os << cr_.D(0, 0) <<  "\t" << cr_.D(0, 1) <<  "\t" << cr_.D(0, 2) <<  "\t" 
        << cr_.D(1, 0) <<  "\t" << cr_.D(1, 1) <<  "\t" << cr_.D(1, 2) <<  "\t" 
        << cr_.D(2, 0) <<  "\t" << cr_.D(2, 1) <<  "\t" << cr_.D(2, 2) <<  "\t"; 

    //rotation-translation
    os << cr_.D(0, 3) <<  "\t" << cr_.D(0, 4) <<  "\t" << cr_.D(0, 5) <<  "\t" 
        << cr_.D(1, 3) <<  "\t" << cr_.D(1, 4) <<  "\t" << cr_.D(1, 5) <<  "\t" 
        << cr_.D(2, 3) <<  "\t" << cr_.D(2, 4) <<  "\t" << cr_.D(2, 5) <<  "\t"; 

    //translation-rotation
    os << cr_.D(3, 0) <<  "\t" << cr_.D(3, 1) <<  "\t" << cr_.D(3, 2) <<  "\t" 
        << cr_.D(4, 0) <<  "\t" << cr_.D(4, 1) <<  "\t" << cr_.D(4, 2) <<  "\t" 
        << cr_.D(5, 0) <<  "\t" << cr_.D(5, 1) <<  "\t" << cr_.D(5, 2) <<  "\t"; 

    //rotation
    os << cr_.D(3, 3) <<  "\t" << cr_.D(3, 4) <<  "\t" << cr_.D(3, 5) <<  "\t" 
        << cr_.D(4, 3) <<  "\t" << cr_.D(4, 4) <<  "\t" << cr_.D(4, 5) <<  "\t" 
        << cr_.D(5, 3) <<  "\t" << cr_.D(5, 4) <<  "\t" << cr_.D(5, 5) <<  "\t"; 
        
    //---------------------------------------------------------------------

    //center of diffusion
    os << cd_.center[0] <<  "\t" << cd_.center[1] <<  "\t" << cd_.center[2] <<  "\t"; 

    //resistance tensor at center of diffusion
    //translation
    os << cd_.Xi(0, 0) <<  "\t" << cd_.Xi(0, 1) <<  "\t" << cd_.Xi(0, 2) <<  "\t" 
        << cd_.Xi(1, 0) <<  "\t" << cd_.Xi(1, 1) <<  "\t" << cd_.Xi(1, 2) <<  "\t" 
        << cd_.Xi(2, 0) <<  "\t" << cd_.Xi(2, 1) <<  "\t" << cd_.Xi(2, 2) <<  "\t"; 

    //rotation-translation
    os << cd_.Xi(0, 3) <<  "\t" << cd_.Xi(0, 4) <<  "\t" << cd_.Xi(0, 5) <<  "\t" 
        << cd_.Xi(1, 3) <<  "\t" << cd_.Xi(1, 4) <<  "\t" << cd_.Xi(1, 5) <<  "\t" 
        << cd_.Xi(2, 3) <<  "\t" << cd_.Xi(2, 4) <<  "\t" << cd_.Xi(2, 5) <<  "\t"; 

    //translation-rotation
    os << cd_.Xi(3, 0) <<  "\t" << cd_.Xi(3, 1) <<  "\t" << cd_.Xi(3, 2) <<  "\t" 
        << cd_.Xi(4, 0) <<  "\t" << cd_.Xi(4, 1) <<  "\t" << cd_.Xi(4, 2) <<  "\t" 
        << cd_.Xi(5, 0) <<  "\t" << cd_.Xi(5, 1) <<  "\t" << cd_.Xi(5, 2) <<  "\t"; 

    //rotation
    os << cd_.Xi(3, 3) <<  "\t" << cd_.Xi(3, 4) <<  "\t" << cd_.Xi(3, 5) <<  "\t" 
        << cd_.Xi(4, 3) <<  "\t" << cd_.Xi(4, 4) <<  "\t" << cd_.Xi(4, 5) <<  "\t" 
        << cd_.Xi(5, 3) <<  "\t" << cd_.Xi(5, 4) <<  "\t" << cd_.Xi(5, 5) <<  "\t"; 


    //diffusion tensor at center of diffusion
    //translation
    os << cd_.D(0, 0) <<  "\t" << cd_.D(0, 1) <<  "\t" << cd_.D(0, 2) <<  "\t" 
        << cd_.D(1, 0) <<  "\t" << cd_.D(1, 1) <<  "\t" << cd_.D(1, 2) <<  "\t" 
        << cd_.D(2, 0) <<  "\t" << cd_.D(2, 1) <<  "\t" << cd_.D(2, 2) <<  "\t"; 

    //rotation-translation
    os << cd_.D(0, 3) <<  "\t" << cd_.D(0, 4) <<  "\t" << cd_.D(0, 5) <<  "\t" 
        << cd_.D(1, 3) <<  "\t" << cd_.D(1, 4) <<  "\t" << cd_.D(1, 5) <<  "\t" 
        << cd_.D(2, 3) <<  "\t" << cd_.D(2, 4) <<  "\t" << cd_.D(2, 5) <<  "\t"; 

    //translation-rotation
    os << cd_.D(3, 0) <<  "\t" << cd_.D(3, 1) <<  "\t" << cd_.D(3, 2) <<  "\t" 
        << cd_.D(4, 0) <<  "\t" << cd_.D(4, 1) <<  "\t" << cd_.D(4, 2) <<  "\t" 
        << cd_.D(5, 0) <<  "\t" << cd_.D(5, 1) <<  "\t" << cd_.D(5, 2) <<  "\t"; 

    //rotation
    os << cd_.D(3, 3) <<  "\t" << cd_.D(3, 4) <<  "\t" << cd_.D(3, 5) <<  "\t" 
        << cd_.D(4, 3) <<  "\t" << cd_.D(4, 4) <<  "\t" << cd_.D(4, 5) <<  "\t" 
        << cd_.D(5, 3) <<  "\t" << cd_.D(5, 4) <<  "\t" << cd_.D(5, 5) <<  "\n"; 


}

}
