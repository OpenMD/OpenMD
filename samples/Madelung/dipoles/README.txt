This directory contains a set of dipolar crystals that can be used to
test electrostatic energies for dipole-dipole interactions.  The
dipolar analogues to the structural Madelung constants for ionic
crystals were first worked out by Sauer who computed the energies of
certain selected dipole arrays (ordered arrays of zero magnetization)
and obtained a number of these constants.[1]

This theory was developed more completely by Luttinger & Tisza [2] and
they tabulated energy constants for the Sauer arrays (and other
periodic structures).  We have repeated the Luttinger & Tisza series
summations to much higher order and obtained the following energy
constants:

Array Type   Lattice    Dipole Direction     Energy constants
----------   -------    ----------------     ----------------
   A           SC             001               -2.67678868438
   A           BCC            001                0
   A           BCC            111               -1.770
   A           FCC            001                2.16693283503
   A           FCC            011               -1.08346641751

   B           SC             001               -2.67678868438
   B           BCC            001               -1.33839434219
   B           BCC            111               -1.770
   B           FCC            001               -1.08346641751
   B           FCC            011               -1.80757363405
 
Type "A" arrays have nearest neighbor strings of antiparallel dipoles.

Type "B" arrays have nearest neighbor strings of antiparallel dipoles
if the dipoles are contained in a plane perpendicular to the dipole
direction that passes through the dipole.

Note that these arrays are not necessarily the minimum energy
structures, and those interested in this problem should consult the
Luttinger & Tisza paper for more details.
            
The electrostatic energy for one of these dipolar arrays is

  E = C N^2 mu^2

where C is the energy constant above, N is the number of dipoles per
unit volume, and mu is the strength of the dipole.

In the units used by OpenMD, with dipoles of 1 Debye, lengths in
angstroms, and energies reported in kcal / mol, the electrostatic
energies are:

  E = 14.39325 C N^2 mu^2

E here is an energy density, so it must be multiplied by the total
volume of the box.

For example, the A_sc_001.md sample has a 8000 dipoles in a (40 angstrom)^3 box
with a lattice spacing of 2 angstroms, so the resulting energy is:

  E = Total volume * Energy per unit volume
    = 40 * 40 * 40 * 14.39325 * (-2.67678868438) * N^2 * mu^2
    = 64000 * (-38.527688) * (1/8)^2 * 1^2 
    = 38527.69 KCal/Mol

Here we've used the definition of N as the number of dipoles per unit 
volume:  N = 1/a^3 = 1/(2^3) 

---------------------------------------------------------------------
[1] J. A. Sauer, "Magnetic Energy Constants of Dipolar Lattices,"
    Phys. Rev. 57, 142–146 (1940) doi: 10.1103/PhysRev.57.142

[2] J. M. Luttinger and L. Tisza, "Theory of Dipole Interaction in
    Crystals," Phys. Rev. 70, 954–964 (1946) doi: 10.1103/PhysRev.70.954
