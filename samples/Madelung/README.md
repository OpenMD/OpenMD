# Madelung Energy

## Background Information
The sample in this directory provides a way of checking the value of
the Madelung Energy for a perfect crystal of NaCl using different methods of computing electrostatic interactions.  The relevant quantities are:

$$
\begin{aligned}
M \text{ (Madelung constant)} &= 1.74756 \\
a \text{ (lattice constant)} &= 5.65 \AA \\

\frac{q^2}{4 \pi \epsilon_0 a} &= 58.77233 \text{ kcal mol}^{-1} \\
M \frac{q^2}{4 \pi \epsilon_0 a} & = 102.708173 \text{ kcal mol}^{-1}\\
\end{aligned}
$$ 

The file `NaCl.omd` contains 8000 ions, so the total electrostatic energy
of the perfect crystal in this file should be:

$$
V_\text{electrostatic} = -821665.38 \text{ kcal mol}^{-1}
$$

## Instructions
Using different electrostatic calculation methods, we can get quite
close to this value.

For example, to test the Damped Shifted Force (DSF) model with a very large cutoff, we set these lines in `NaCl.omd` :

```
cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "damped";
cutoffRadius = 28;
dampingAlpha = 0.14159292;
```

Other options for the `cutoffMethod` include `shifted_potential`, `reaction_field`, `ewald_full`, and `hard`, although `shifted_force` provides remarkable accuracy at minimal computational cost.

To extract values for the electrostatic potential in OpenMD, we add the
`ELECTROSTATIC_POTENTIAL` keyword to the end of the `statFileFormat`:

```
statFileFormat = "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY|ELECTROSTATIC_POTENTIAL";
```

Once those lines have been set we run a single step simulation:

```
mpirun -np 4 openmd_MPI NaCl.omd
```

## Expected Output
The report generated from this simulation should look like:
```
###############################################################################
# Status Report:                                                              #
#              Total Time:           2 fs                                     #
#       Number of Samples:           2                                        #
#            Total Energy:     -709366  ±  0            kcal/mol              #
#        Potential Energy:     -709366  ±  0            kcal/mol              #
#          Kinetic Energy: 2.34108e-25  ±  3.24457e-25  kcal/mol              #
#             Temperature: 9.81846e-27  ±  1.36077e-26  K                     #
#                Pressure:      116088  ±  0.00118777   atm                   #
#                  Volume:      180362  ±  0            A^3                   #
#      Conserved Quantity:     -709366  ±  0            kcal/mol              #
# Electrostatic Potential:     -821668  ±  0.00589253   kcal/mol              #
###############################################################################
```

Looking into the `NaCl.stat` file shows that the resultant electrostatic potential is:
$$
V_\text{electrostatic} = -821667.61 \text{ kcal mol}^{-1}
$$
which is only a 0.0002714 % error from the true Madelung energy.