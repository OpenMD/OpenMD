Madelung Energy Sample Calculation

The sample in this directory provides a way of checking the value of
the Madelung Energy for a perfect crystal of NaCl.  The relevant
quantities are:

M (Madelung constant) = 1.74756
a (lattice constant)  = 5.65 Angstroms
q^2 / (4 pi e0 a)     = 58.77233   kcal / mol
M q^2 / (4 pi e0 a)   = 102.708173 kcal / mol

The file NaCl.md contains 8000 ions, so the total electrostatic energy
of the perfect crystal in this file should be:

V_electrostatic = -821665.38  kcal / mol

Using different electrostatic calculation methods, we can get quite
close to this value.

For example, with :
cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "damped";
cutoffRadius = 28;
dampingAlpha = 0.14159292;

The resultant electrostatic potential is:  -821667.79 kcal / mol

To obtain values for the electrostatic potential in OpenMD, we add the
ELECTROSTATIC_POTENTIAL keyword to the end of the statFileFormat:

statFileFormat = "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY|ELECTROSTATIC_POTENTIAL";
