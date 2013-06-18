This directory contains a set of dipolar crystals that can be used to
test electrostatic computation for dipole-dipole interactions.  The
dipolar analogues to the structural Madelung constants for ionic
crystals were first worked out by Sauer who computed the energies of
certain selected dipole arrays (ordered arrays of zero magnetization)
and obtained a number of these constants.[1]

This theory was developed more completely by Luttinger & Tisza [2] and
they tabulated these constants as follows:

Array Type   Lattice    Dipole Direction     Energy constants
----------   -------    ----------------     ----------------
   A           SC             001               -2.676
   A           BCC            001                0
   A           BCC            111               -1.770
   A           FCC            001                2.167
   A           FCC            011               -1.084

   B           SC             001               -2.676
   B           BCC            001               -1.338
   B           BCC            111               -1.770
   B           FCC            001               -1.084
   B           FCC            011               -1.808
 
Type "A" arrays have nearest neighbor strings of antiparallel dipoles.

Type "B" arrays have nearest neighbor strings of antiparallel dipoles
if the dipoles are contained in a plane perpendicular to the dipole
direction that passes through the dipole.

Note that these arrays are not necessarily the minimum energy
structures, and those interested in this problem should consult the
Luttinger & Tisza paper for more details.
            
The electrostatic energy for one of these dipolar arrays is

  E = C N^2 mu^2

where C is the energy constant above, N is the number of dipoles, and
mu is the strength of the dipole.
