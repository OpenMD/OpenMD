# Quadrupolar Crystals

## Background Information
This directory contains a set of quadrupolar crystals that can be used
to test electrostatic energies for quadrupole-quadrupole interactions.
The quadrupolar analogues to the structural Madelung constants for
ionic crystals were first worked out by Nagai and Nakamura who
computed the energies of certain selected quadrupole arrays and
obtained a number of these constants.

> O. Nagai and T. Nakamura, "Quadrupole Interaction in Crystals,"
> Progress of Theoretical Physics 24 (2), 432-454 (1960) 
> doi: 10.1143/PTP.24.432,   also see Errata: Progress of
> Theoretical Physics 30 (3), 412 (1963).  doi: 10.1143/PTP.30.412a

Their work is based on earlier work on dipolar crystals that was done 
by Luttinger & Tisza:

> J. M. Luttinger and L. Tisza, "Theory of Dipole Interaction in
> Crystals," Phys. Rev. 70, 954–964 (1946) doi: 10.1103/PhysRev.70.954

We have generated only the lowest energy configurations for the linear
quadrupoles

| Lattice | Quadrupole Direction | Energy constants | Notes  |
|---------|----------------------|------------------|--------|
| SC      |      111             |       -8.30      | Fig 6a |
| BCC     |      011             |      -21.7       | Fig 7b |
| FCC     |      111             |      -80.5       | Fig 8  |
             
The electrostatic energy for one of these dipolar arrays is

  $$ E = C \left( \frac{3}{4} \right) N^2 Q^2 $$

where $C$ is the energy constant above, N is the number of quadrupoles
per unit volume, and $Q$ is the strength of the dipole.

In the units used by OpenMD, with quadrupole moments of 1
Debye-angstrom, lengths in angstroms, and energies reported in kcal /
mol, the electrostatic energies are:

  $$ E = 14.39325 \left( \frac{3}{4} \right) C N^2 Q^2$$

$E$ here is an energy density, so it must be multiplied by the total
volume of the box.

For example, the `fcc.omd` sample has a 8000 quadrupoles (1
Debye-angstrom each) in a (40 angstrom)^3 box with a lattice spacing
of 2 angstroms, so the resulting energy is:

$$\begin{aligned}
  E &= \text{Total volume} * \text{Energy per unit volume} \\
    &= 40 * 40 * 40 * 14.39325 * \left( \frac{3}{4} \right) * (-8.30) * N^2 * Q^2 \\
    &= 64000 * (-89.598) * (1/8)^2 * 1^2  \\
    &= -89598 \text{ kcal mol}^{-1} \\
\end{aligned}$$

Here we've used the definition of $N$ as the number of quadrupoles per unit 
volume:

$$N = 1/a^3 = 1/(2^3)$$

## Instructions

To test the shifted force (DSF) (or shfited gradient) method of
carrying out electrostatic interactions with a large cutoff, and
*without damping*, we would add these lines to a sample file:

```
cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "undamped";
cutoffRadius = 20.0;
statFileFormat = "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGYTEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITYELECTROSTATIC_POTENTIAL";
```

Other options for the `cutoffMethod` include `shifted_potential`,
`reaction_field`, `ewald_full`, `taylor_shifted`, and `hard`, although
`shifted_force` provides remarkable accuracy at minimal computational
cost.

Once those lines have been set we run a single step simulation:

```
mpirun -np 4 openmd_MPI sc.omd
```

## Expected Output

The report generated from this simulation should look like:
```
###############################################################################
# Status Report:                                                              #
#              Total Time:           1 fs                                     #
#       Number of Samples:           2                                        #
#            Total Energy:      -92904  ±  0            kcal/mol              #
#        Potential Energy:      -92904  ±  0.00152901   kcal/mol              #
#          Kinetic Energy: 7.41398e-11  ±  1.02752e-10  kcal/mol              #
#             Temperature: 1.55461e-12  ±  2.15458e-12  K                     #
#                Pressure:     -166440  ±  0.00277351   atm                   #
#                  Volume:       64000  ±  0            A^3                   #
#      Conserved Quantity:      -92904  ±  0            kcal/mol              #
# Electrostatic Potential:    -89472.4  ±  0            kcal/mol              #
###############################################################################
```

One interesting test is to modify the parameters to include:
```
electrostaticScreeningMethod = "damped";
dampingAlpha = 0.18;
```
In this case, the electrostatic potential gets a bit closer to the series approximation 
at  -89490.607 kcal/mol.

