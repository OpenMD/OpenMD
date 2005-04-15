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
 
#include <cmath>

#include "minimizers/CGFamilyMinimizer.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Utility.hpp"

namespace oopse {

  CGFamilyMinimizer::CGFamilyMinimizer(SimInfo *info) : Minimizer(info){
    prevG.resize(ndim);
    prevX.resize(ndim);
    direction.resize(ndim);

    stepSize = paramSet->getStepSize();
  }

  int CGFamilyMinimizer::checkConvg(){
    double fTol;
    double relativeFTol;  // relative tolerance
    double deltaF;
    double gTol;
    double relativeGTol;
    double gnorm;
  

    // test function tolerance test
    fTol =paramSet->getFTol();
    relativeFTol = fTol * std::max(1.0,fabs(curF));  // relative tolerance
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
    relativeGTol = gTol * std::max(1.0,fabs(curF));

#ifndef IS_MPI
    gnorm = sqrt(dotProduct(curG, curG));
#else
    double localDP;
    double globalDP;

    localDP = dotProduct(curG, curG);
    MPI_Allreduce(&localDP, &globalDP, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);  
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
