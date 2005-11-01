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
 
#include "brains/Register.hpp"

#include "integrators/IntegratorFactory.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "integrators/Integrator.hpp"
#include "integrators/NVE.hpp"
#include "integrators/NVT.hpp"
#include "integrators/NPTi.hpp"
#include "integrators/NPTf.hpp"
#include "integrators/NPTxyz.hpp"
#include "integrators/NPAT.hpp"
#include "integrators/NPrT.hpp"

#include "minimizers/MinimizerFactory.hpp"
#include "minimizers/MinimizerCreator.hpp"
#include "minimizers/PRCG.hpp"
#include "minimizers/SDMinimizer.hpp"
#include "UseTheForce/DUFF.hpp"
#include "UseTheForce/EAM_FF.hpp"
#include "UseTheForce/ForceFieldFactory.hpp"
#include "UseTheForce/ForceFieldCreator.hpp"
#include "UseTheForce/SHAPES_FF.hpp"
#include "lattice/LatticeFactory.hpp"
#include "lattice/LatticeCreator.hpp"
#include "lattice/FCCLattice.hpp"

namespace oopse {


  void registerForceFields() {
    /** @todo move to a seperate initialization module */
    //DUFF, WATER and LJ are merged into one force field
    ForceFieldFactory::getInstance()->registerForceField(new ForceFieldBuilder<DUFF>("DUFF"));
    ForceFieldFactory::getInstance()->registerForceField(new ForceFieldBuilder<DUFF>("WATER"));
    ForceFieldFactory::getInstance()->registerForceField(new ForceFieldBuilder<DUFF>("LJ"));
    //in theory, EAM can also be merged
    ForceFieldFactory::getInstance()->registerForceField(new ForceFieldBuilder<EAM_FF>("EAM"));
    //heck, that worked...  let's try merging SHAPES
    ForceFieldFactory::getInstance()->registerForceField(new ForceFieldBuilder<SHAPES_FF>("SHAPES"));

  }

  void registerIntegrators() {
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NVE>("NVE"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NVT>("NVT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTi>("NPTI"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTf>("NPTF"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTxyz>("NPTXYZ"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPAT>("NPAT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPrT>("NPRT"));

  }

  void registerMinimizers() {
    MinimizerFactory::getInstance()->registerMinimizer(new MinimizerBuilder<SDMinimizer>("SD"));
    MinimizerFactory::getInstance()->registerMinimizer(new MinimizerBuilder<PRCGMinimizer>("CG"));
  }

  void registerLattice(){
    LatticeFactory::getInstance()->registerLattice(new LatticeBuilder<FCCLattice>("FCC"));
  }

  void registerAll() {
    registerForceFields();
    registerIntegrators();
    registerMinimizers();
  }

}
