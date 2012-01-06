/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <cmath>

#include "minimizers/CGFamilyMinimizer.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Utility.hpp"
#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

  CGFamilyMinimizer::CGFamilyMinimizer(SimInfo *info) : Minimizer(info){
    prevG.resize(ndim);
    prevX.resize(ndim);
    direction.resize(ndim);

    stepSize = paramSet->getStepSize();
  }

  int CGFamilyMinimizer::checkConvg(){
    RealType fTol;
    RealType relativeFTol;  // relative tolerance
    RealType deltaF;
    RealType gTol;
    RealType relativeGTol;
    RealType gnorm;

    // test function tolerance test
    fTol =paramSet->getFTol();

    relativeFTol = fTol * std::max(RealType(1.0), fabs(curF));  // relative tolerance

    deltaF = prevF - curF;
  
    if (fabs(deltaF) <= relativeFTol) {

      if (bVerbose){
	std::cout << "function value tolerance test passed" << std::endl;
	std::cout << "ftol = " << fTol
		  << "\tdeltaf = " << deltaF<< std::endl;
      }
      return CONVG_FTOL;
    }
  
    //gradient tolerance test
    gTol = paramSet->getGTol();
    relativeGTol = gTol * std::max(RealType(1.0), fabs(curF));

#ifndef IS_MPI
    gnorm = sqrt(dotProduct(curG, curG));
#else
    RealType localDP;
    RealType globalDP;

    localDP = dotProduct(curG, curG);
    MPI_Allreduce(&localDP, &globalDP, 1, MPI_REALTYPE,MPI_SUM, MPI_COMM_WORLD);  
    gnorm  = sqrt(globalDP);
#endif

    if (gnorm <= relativeGTol) {
      std::cout << "gradient tolerance test" << std::endl;
      std::cout << "gnorm = " << gnorm
		<< "\trelativeGTol = " << relativeGTol<< std::endl;
      return CONVG_GTOL;
    }
  
    //absolute gradient tolerance test

    if (gnorm <= gTol) {
      std::cout << "absolute gradient tolerance test" << std::endl;
      std::cout << "gnorm = " << gnorm
		<< "\tgTol = " << gTol<< std::endl;
      return CONVG_ABSGTOL;
    }

    // did not converge yet
    return CONVG_UNCONVG;
  }

}
