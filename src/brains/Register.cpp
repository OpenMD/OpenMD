/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "brains/Register.hpp"

#include "integrators/Integrator.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "integrators/IntegratorFactory.hpp"
#include "integrators/LangevinDynamics.hpp"
#include "integrators/LangevinPiston.hpp"
#include "integrators/NPA.hpp"
#include "integrators/NPAT.hpp"
#include "integrators/NPTf.hpp"
#include "integrators/NPTi.hpp"
#include "integrators/NPTsz.hpp"
#include "integrators/NPTxyz.hpp"
#include "integrators/NPrT.hpp"
#include "integrators/NVE.hpp"
#include "integrators/NVT.hpp"
#include "integrators/NgammaT.hpp"
#if defined(HAVE_QHULL)
#include "integrators/LangevinHullDynamics.hpp"
#endif
#include "integrators/SPFDynamics.hpp"
#include "lattice/BCCLattice.hpp"
#include "lattice/FCCLattice.hpp"
#include "lattice/LatticeCreator.hpp"
#include "lattice/LatticeFactory.hpp"
#include "lattice/SCLattice.hpp"
#include "optimization/BFGS.hpp"
#include "optimization/ConjugateGradient.hpp"
#include "optimization/Method.hpp"
#include "optimization/OptimizationCreator.hpp"
#include "optimization/OptimizationFactory.hpp"
#include "optimization/SteepestDescent.hpp"

using namespace QuantLib;
namespace OpenMD {

  void registerIntegrators() {
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NVE>("NVE"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NVT>("NVT"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPTi>("NPTI"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPTf>("NPTF"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPTxyz>("NPTXYZ"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPAT>("NPAT"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPA>("NPA"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPrT>("NPRT"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NPrT>("NPGT"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NgammaT>("NGT"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<NgammaT>("NGAMMAT"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<LangevinDynamics>("LANGEVINDYNAMICS"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<LangevinDynamics>("LD"));
#if defined(HAVE_QHULL)
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<LangevinHullDynamics>("LHULL"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<LangevinHullDynamics>("LANGEVINHULL"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<LangevinHullDynamics>("SMIPD"));
#endif
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<LangevinPiston>("LANGEVINPISTON"));
    IntegratorFactory::getInstance().registerIntegrator(
        new IntegratorBuilder<SPFDynamics>("SPF"));
  }

  void registerOptimizers() {
    OptimizationFactory::getInstance().registerOptimization(
        new OptimizationBuilder<QuantLib::SteepestDescent>("SD"));
    OptimizationFactory::getInstance().registerOptimization(
        new OptimizationBuilder<QuantLib::ConjugateGradient>("CG"));
    OptimizationFactory::getInstance().registerOptimization(
        new OptimizationBuilder<QuantLib::BFGS>("BFGS"));
  }

  void registerLattice() {
    LatticeFactory::getInstance().registerLattice(
        new LatticeBuilder<FCCLattice>("FCC"));
    LatticeFactory::getInstance().registerLattice(
        new LatticeBuilder<SCLattice>("SC"));
    LatticeFactory::getInstance().registerLattice(
        new LatticeBuilder<BCCLattice>("BCC"));
  }

  void registerAll() {
    registerIntegrators();
    registerOptimizers();
  }

}  // namespace OpenMD
