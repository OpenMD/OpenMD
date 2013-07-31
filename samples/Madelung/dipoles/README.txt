This directory contains a set of dipolar crystals that can be used to
test electrostatic energies for dipole-dipole interactions.  The
dipolar analogues to the Madelung constants for ionic crystals were
first worked out by Sauer who computed the energies of certain
selected dipole arrays (ordered arrays of zero magnetization) and
obtained a number of these constants.[1]

This theory was developed more completely by Luttinger & Tisza [2] who
tabulated energy constants for the Sauer arrays (and other periodic
structures).  We have repeated the Luttinger & Tisza series summations
to much higher order and obtained the following energy constants:

Array Type   Lattice    Dipole Direction     Energy constants
----------   -------    ----------------     ----------------
   A           SC             001               -2.676788684
   A           BCC            001                0
   A           BCC            111               -1.770078733
   A           FCC            001                2.166932835
   A           FCC            011               -1.083466417

   *           BCC          minimum             -1.985920929

   B           SC             001               -2.676788684
   B           BCC            001               -1.338394342
   B           BCC            111               -1.770078733
   B           FCC            001               -1.083466417
   B           FCC            011               -1.807573634
 
Type "A" arrays have nearest neighbor strings of antiparallel dipoles.

Type "B" arrays have nearest neighbor strings of antiparallel dipoles
if the dipoles are contained in a plane perpendicular to the dipole
direction that passes through the dipole.

There's also an additional minimum energy structure for the BCC
lattice that was found by Luttinger & Tisza.  All of the energy
constants can be recomputed to a very high degree of accuracy using
the "fieldValues" python script contained in this directory.  Dipolar
arrays matching these configurations which are suitable for use with
OpenMD can be generated with the "buildDipolarArray" program.
            
The electrostatic energy for one of these dipolar arrays is

  E = C N^2 mu^2

where C is the energy constant above, N is the number of dipoles per
unit volume, and mu is the strength of the dipole.

In the units used by OpenMD, dipoles are measured in Debye, lengths in
angstroms, and energies reported in kcal / mol, the electrostatic
energies are:

  E = 14.39325 C N^2 mu^2

E here is an energy density, so it must be multiplied by the total
volume of the box.

For example, the A_sc_001.md sample has a 8000 dipoles (1 Debye each)
in a (40 angstrom)^3 box with a lattice spacing of 2 angstroms, so the
resulting energy is:

  E = Total volume * Energy per unit volume
    = 40 * 40 * 40 * 14.39325 * (-2.676788684) * N^2 * mu^2
    = 64000 * (-38.52769) * (1/8)^2 * 1^2 
    = -38527.69 kcal/mol

Here we've used the definition of N as the number of dipoles per unit 
volume:  N = 1/a^3 = 1/(2^3) 

---------------------------------------------------------------------
[1] J. A. Sauer, "Magnetic Energy Constants of Dipolar Lattices,"
    Phys. Rev. 57, 142–146 (1940) doi: 10.1103/PhysRev.57.142

[2] J. M. Luttinger and L. Tisza, "Theory of Dipole Interaction in
    Crystals," Phys. Rev. 70, 954–964 (1946) doi: 10.1103/PhysRev.70.954
    Also note the errata contained in: Phys. Rev. 72, 257 (1947) 
    doi: 10.1103/PhysRev.72.257

