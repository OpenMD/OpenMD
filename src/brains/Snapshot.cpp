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
  * @file Snapshot.cpp
  * @author tlin
  * @date 11/11/2004
  * @time 10:56am
  * @version 1.0
  */

#include "brains/Snapshot.hpp"
#include "utils/NumericConstant.hpp"
#include "utils/simError.h"
#include "utils/Utility.hpp"
namespace oopse {

void  Snapshot::setHmat(const Mat3x3d& m) {
    const double orthoTolerance = NumericConstant::epsilon;
    hmat_ = m;
    invHmat_ = hmat_.inverse();
    
    //prepare fortran Hmat 
    double fortranHmat[9];
    double fortranInvHmat[9];
    hmat_.getArray(fortranHmat);
    invHmat_.getArray(fortranInvHmat);

    //determine whether the box is orthoTolerance or not
    int oldOrthoRhombic = orthoRhombic_;
    
    double smallDiag = fabs(hmat_(0, 0));
    if(smallDiag > fabs(hmat_(1, 1))) smallDiag = fabs(hmat_(1, 1));
    if(smallDiag > fabs(hmat_(2, 2))) smallDiag = fabs(hmat_(2, 2));    
    double tol = smallDiag * orthoTolerance;

    orthoRhombic_ = 1;

    for (int i = 0; i < 3; i++ ) {
        for (int j = 0 ; j < 3; j++) {
            if (i != j) {
                if (orthoRhombic_) {
                    if ( fabs(hmat_(i, j)) >= tol)
                        orthoRhombic_ = 0;
                }        
            }
        }
    }

    if( oldOrthoRhombic != orthoRhombic_ ){

        if( orthoRhombic_ ) {
            sprintf( painCave.errMsg,
                "OOPSE is switching from the default Non-Orthorhombic\n"
                "\tto the faster Orthorhombic periodic boundary computations.\n"
                "\tThis is usually a good thing, but if you wan't the\n"
                "\tNon-Orthorhombic computations, make the orthoBoxTolerance\n"
                "\tvariable ( currently set to %G ) smaller.\n",
                orthoTolerance);
            painCave.severity = OOPSE_INFO;
            simError();
        }
        else {
            sprintf( painCave.errMsg,
                "OOPSE is switching from the faster Orthorhombic to the more\n"
                "\tflexible Non-Orthorhombic periodic boundary computations.\n"
                "\tThis is usually because the box has deformed under\n"
                "\tNPTf integration. If you wan't to live on the edge with\n"
                "\tthe Orthorhombic computations, make the orthoBoxTolerance\n"
                "\tvariable ( currently set to %G ) larger.\n",
                orthoTolerance);
            painCave.severity = OOPSE_WARNING;
            simError();
        }
    }    

    //notify fortran simulation box has changed
    setFortranBox(fortranHmat, fortranInvHmat, &orthoRhombic_);
}


void Snapshot::wrapVector(Vector3d& pos) {

    int i;
    Vector3d scaled;

    if( !orthoRhombic_ ){

        // calc the scaled coordinates.
        scaled = invHmat_* pos;

        // wrap the scaled coordinates
        for (i = 0; i < 3; ++i) {
            scaled[i] -= roundMe(scaled[i]);
        }

        // calc the wrapped real coordinates from the wrapped scaled coordinates
        pos = hmat_ * scaled;    

    } else {//if it is orthoRhombic, we could improve efficiency by only caculating the diagonal element
    
        // calc the scaled coordinates.
        for (i=0; i<3; i++) {
            scaled[i] = pos[i] * invHmat_(i, i);
        }
        
        // wrap the scaled coordinates
        for (i = 0; i < 3; ++i) {
            scaled[i] -= roundMe(scaled[i]);
        }

        // calc the wrapped real coordinates from the wrapped scaled coordinates
        for (i=0; i<3; i++) {
            pos[i] = scaled[i] * hmat_(i, i);
        }
        
    }

}

}
  
