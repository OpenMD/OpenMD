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
 
#include "brains/Register.hpp"

#include "integrators/IntegratorFactory.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "integrators/Integrator.hpp"
#include "integrators/NVE.hpp"
#include "integrators/NVT.hpp"
#include "integrators/NPTi.hpp"
#include "integrators/NPTf.hpp"
#include "integrators/NPTxyz.hpp"
#include "integrators/NPTsz.hpp"
#include "integrators/NPAT.hpp"
#include "integrators/NPrT.hpp"
#include "integrators/NgammaT.hpp"
#include "integrators/LangevinDynamics.hpp"
#if defined(HAVE_QHULL)
#include "integrators/LangevinHullDynamics.hpp"
#endif

#include "optimization/OptimizationFactory.hpp"
#include "optimization/OptimizationCreator.hpp"
#include "optimization/Method.hpp"
#include "optimization/SteepestDescent.hpp"
#include "optimization/ConjugateGradient.hpp"
#include "optimization/BFGS.hpp"

#include "lattice/LatticeFactory.hpp"
#include "lattice/LatticeCreator.hpp"
#include "lattice/FCCLattice.hpp"

using namespace QuantLib;
namespace OpenMD {

  void registerIntegrators() {
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NVE>("NVE"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NVT>("NVT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTi>("NPTI"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTf>("NPTF"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTxyz>("NPTXYZ"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPAT>("NPAT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPrT>("NPRT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPrT>("NPGT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NgammaT>("NGT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NgammaT>("NGAMMAT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<LangevinDynamics>("LANGEVINDYNAMICS"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<LangevinDynamics>("LD"));
#if defined(HAVE_QHULL)
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<LangevinHullDynamics>("LHULL"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<LangevinHullDynamics>("LANGEVINHULL"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<LangevinHullDynamics>("SMIPD"));
#endif
  }

  void registerOptimizers() {
    OptimizationFactory::getInstance()->registerOptimization(new OptimizationBuilder<QuantLib::SteepestDescent>("SD"));
    OptimizationFactory::getInstance()->registerOptimization(new OptimizationBuilder<QuantLib::ConjugateGradient>("CG"));
    OptimizationFactory::getInstance()->registerOptimization(new OptimizationBuilder<QuantLib::BFGS>("BFGS"));
  }

  void registerLattice(){
    LatticeFactory::getInstance()->registerLattice(new LatticeBuilder<FCCLattice>("FCC"));
  }

  void registerAll() {
    registerIntegrators();
    registerOptimizers();
  }

}
