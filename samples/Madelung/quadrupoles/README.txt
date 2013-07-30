This directory contains a set of quadrupolar crystals that can be used
to test electrostatic energies for quadrupole-quadrupole interactions.
The quadrupolar analogues to the structural Madelung constants for
ionic crystals were first worked out by Nagai and Nakamura who
computed the energies of certain selected quadrupole arrays and
obtained a number of these constants.[1] Their work is based on
earlier work on dipolar crystals that was done by Luttinger & Tisza
[2].  We have generated only the lowest energy configurations for the
linear quadrupoles 

 Lattice    Quadrupole Direction     Energy constants   Notes
 -------    --------------------     ----------------   ------
   SC             111                     -8.30         Fig 6a
   BCC            011                    -21.7          Fig 7b
   FCC            111                    -80.5          Fig 8
             
The electrostatic energy for one of these dipolar arrays is

  E = C (3/4) N^2 Q^2

where C is the energy constant above, N is the number of quadrupoles
per unit volume, and Q is the strength of the dipole.

In the units used by OpenMD, with quadrupole momentss of 1
Debye-angstrom, lengths in angstroms, and energies reported in
kcal / mol, the electrostatic energies are:

  E = 14.39325 (3/4) C N^2 Q^2

E here is an energy density, so it must be multiplied by the total
volume of the box.


---------------------------------------------------------------------
[1] O. Nagai and T. Nakamura, "Quadrupole Interaction in Crystals,"
    Progress of Theoretical Physics 24 (2), 432-454 (1960) 
    doi: 10.1143/PTP.24.432

[2] J. M. Luttinger and L. Tisza, "Theory of Dipole Interaction in
    Crystals," Phys. Rev. 70, 954–964 (1946) doi: 10.1103/PhysRev.70.954
