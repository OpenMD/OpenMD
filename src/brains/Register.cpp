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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
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
#include "integrators/NPA.hpp"
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
#include "lattice/SCLattice.hpp"
#include "lattice/BCCLattice.hpp"

#include "analysis/AnalyzerFactory.hpp"
#include "analysis/AnalyzerCreator.hpp"
#include "analysis/GofR.hpp"
#include "analysis/GofZ.hpp"
#include "analysis/GofRZ.hpp"
#include "analysis/GofRAngle.hpp"
#include "analysis/GofAngle2.hpp"
#include "analysis/GofRAngle2.hpp"
#include "analysis/GofXyz.hpp"
#include "analysis/TwoDGofR.hpp"
#include "analysis/P2OrderParameter.hpp"
#include "analysis/BondOrderParameter.hpp"
#include "analysis/BOPofR.hpp"
#include "analysis/RippleOP.hpp"
#include "analysis/SCDOrderParameter.hpp"
#include "analysis/DensityPlot.hpp"
#include "analysis/ObjectCount.hpp"
#include "analysis/RhoZ.hpp"
#include "analysis/PipeDensity.hpp"
#include "analysis/pAngle.hpp"
#include "analysis/BondAngleDistribution.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "analysis/Hxy.hpp"
#endif
#include "analysis/RhoR.hpp"
#include "analysis/AngleR.hpp"
#include "analysis/TetrahedralityParam.hpp"
#include "analysis/TetrahedralityParamZ.hpp"
#include "analysis/TetrahedralityParamXYZ.hpp"
#include "analysis/TetrahedralityParamDens.hpp"
#include "analysis/RNEMDStats.hpp"
#include "analysis/NitrileFrequencyMap.hpp"
#include "analysis/MultipoleSum.hpp"
#include "analysis/CoordinationNumber.hpp"
#include "analysis/HBondGeometric.hpp"
#include "analysis/PotDiff.hpp"
#include "analysis/TetrahedralityHBMatrix.hpp"
#include "analysis/Kirkwood.hpp"
#include "analysis/Field.hpp"
#include "analysis/VelocityZ.hpp"
#include "analysis/CenterOfMass.hpp"
#include "analysis/ContactAngle1.hpp"
#include "analysis/ContactAngle2.hpp"
#include "analysis/GCNSeq.hpp"
#include "analysis/NanoLength.hpp"
#include "analysis/NanoVolume.hpp"

using namespace QuantLib;
namespace OpenMD {

  void registerIntegrators() {
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NVE>("NVE"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NVT>("NVT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTi>("NPTI"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTf>("NPTF"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPTxyz>("NPTXYZ"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPAT>("NPAT"));
    IntegratorFactory::getInstance()->registerIntegrator(new IntegratorBuilder<NPA>("NPA"));
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

  void registerAnalyzers() {
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<BondOrderParameter>("bo"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<IcosahedralOfR>("ior"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<FCCOfR>("for"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<BondAngleDistribution>("bad"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<ObjectCount>("count"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofR>("gofr"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofZ>("gofz"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofRTheta>("r_theta"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofROmega>("r_omega"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofRZ>("r_z"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofAngle2>("theta_omega"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofRAngle2>("r_theta_omega"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GofXyz>("gxyz"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<TwoDGofR>("twodgofr"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<P2OrderParameter>("p2"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<RippleOP>("rp2"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<SCDOrderParameter>("scd"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<DensityPlot>("density"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<RhoZ>("slab_density"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<PipeDensity>("pipe_density"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<pAngle>("p_angle"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<Hxy>("hxy"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<RhoR>("rho_r"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<AngleR>("angle_r"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<TetrahedralityParam>("tet_param"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<TetrahedralityParamZ>("tet_param_z"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<TetrahedralityParamDens>("tet_param_dens"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<TetrahedralityParamXYZ>("tet_param_xyz"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<RNEMDZ>("rnemdz"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<RNEMDR>("rnemdr"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<RNEMDRTheta>("rnemdrt"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<NitrileFrequencyMap>("nitrile"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<MultipoleSum>("multipole"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<CoordinationNumber>("cn"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<SCN>("scn"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GCN>("gcn"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<HBondGeometric>("hbond"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<PotDiff>("potDiff"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<TetrahedralityHBMatrix>("tet_hb"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<Kirkwood>("kirkwood"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<KirkwoodQuadrupoles>("kirkwoodQ"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<DensityField>("densityfield"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<VelocityField>("velocityfield"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<VelocityZ>("velocityZ"));
    // Sequential Analyzers:
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<CenterOfMass>("com"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<GCNSeq>("gcnseq"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<ContactAngle1>("ca1"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<ContactAngle2>("ca2"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<NanoLength>("nanolength"));
    AnalyzerFactory::getInstance()->registerAnalyzer(new AnalyzerBuilder<NanoVolume>("nanovolume"));
  }
  
  void registerLattice(){
    LatticeFactory::getInstance()->registerLattice(new LatticeBuilder<FCCLattice>("FCC"));
    LatticeFactory::getInstance()->registerLattice(new LatticeBuilder<SCLattice>("SC"));
    LatticeFactory::getInstance()->registerLattice(new LatticeBuilder<BCCLattice>("BCC"));
  }
  
  void registerAll() {
    registerIntegrators();
    registerOptimizers();
    registerAnalyzers();
  }

}
