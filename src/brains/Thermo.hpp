/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef BRAINS_THERMO_HPP
#define BRAINS_THERMO_HPP

#include "brains/SimInfo.hpp"
#include "primitives/Atom.hpp"

namespace OpenMD {

  class Thermo {
  public:
    Thermo(SimInfo* info) : info_(info) {}

    // note: all the following energies are in kcal/mol

    RealType getTranslationalKinetic();  // the translational kinetic energy
    RealType getRotationalKinetic();     // the rotational kinetic energy
    RealType getElectronicKinetic();     // the electronic kinetic energy
    RealType getKinetic();               // the total kinetic energy
    RealType getPotential();             // the total potential energy
    potVec getSelectionPotentials();     // the potential energy of a selection

    RealType getTotalEnergy();  // gets the total energy

    RealType getTemperature();            // Gives the instant temp. in K
    RealType getElectronicTemperature();  // gives the instant electronic
                                          // temperature in K
    RealType getNetCharge();       // gives the total net charge on the system
    RealType getChargeMomentum();  // gives the instantaneous charge momentum in
                                   // kcal fs / e / mol
    std::vector<Vector3d> getCurrentDensity();

    RealType getPressure();  // gives the instant pressure in atm;
    RealType getPressure(Snapshot* snap); // gives the instant pressure in atm for a given snapshot;

    /** \brief gives the pressure tensor in amu*fs^-2*Ang^-1 */
    Mat3x3d getPressureTensor();
    /** \brief gives the pressure tensor in amu*fs^-2*Ang^-1 for a given snapshot */
    Mat3x3d getPressureTensor(Snapshot* snap);

    RealType getVolume();  // gives the volume in Ang^3

    /** \brief accumulate and return the simulation box dipole moment in C*m */
    Vector3d getSystemDipole();

    /** \brief accumulate and return the simulation box dipole moment in debye
     * Angstroms */
    Mat3x3d getSystemQuadrupole();

    Vector3d getHeatFlux();

    /** \brief Returns the center of the mass of the whole system.*/
    Vector3d getCom();

    /** \brief Returns the velocity of center of mass of the whole system.*/
    Vector3d getComVel();

    /** \brief Returns the center of the mass and Center of Mass velocity of
        the whole system.*/
    void getComAll(Vector3d& com, Vector3d& comVel);

    /** \brief Returns the inertia tensor and the total angular
        momentum for for the entire system
     * \param[out] inertiaTensor the inertia tensor
     * \param[out] angularMomentum the angular momentum vector
     * \ingroup surface
     */
    void getInertiaTensor(Mat3x3d& inertiaTensor, Vector3d& angularMomentum);

    /** \brief Returns the Axis-aligned bounding box for the current system.
     */
    Mat3x3d getBoundingBox();

    /** \brief Returns system angular momentum */
    Vector3d getAngularMomentum();

    /** \brief Returns volume of system as estimated by an ellipsoid defined
        by the radii of gyration */
    RealType getGyrationalVolume();

    /** \brief Overloaded version of gyrational volume that also returns
        det(I) so dV/dr can be calculated */
    void getGyrationalVolume(RealType& vol, RealType& detI);

    RealType getHullVolume();

    RealType getTaggedAtomPairDistance();

  private:
    SimInfo* info_ {nullptr};
  };
}  // namespace OpenMD

#endif
