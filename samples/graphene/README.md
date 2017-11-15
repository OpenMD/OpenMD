This sample illustrates a simulation of united-atom (UA) propylene
monomers confined between and all-atom (AA) representation of graphene
sheets.  The force field is based partially on OPLS-AA, and partially
on TraPPE-UA.

Propylene parameters from: C.D. Wick, M.G. Martin, and J.I. Siepmann,
"Transferable potentials for phase equilibria. 4. United-atom
description of linear and branched alkenes and of alkylbenzenes,"
*J. Phys. Chem. B*, **104**, pp. 8008-8016 (2000).
DOI: [10.1021/jp001044x](https://doi.org/10.1021/jp001044x)

Note that the molecule definition in `graphene.inc` includes bonds that span the
box boundaries.  If you want to extend to larger boxes, start with the molecule
definition in `graphene.raw.inc` and then add bonds that span the boundaries of
the new box.

Files included:

1. `graphene.frc`: Force field file
2. `graphene.inc`: Molecule definition include file for a graphene
   sheet (includes cross-box bonding, so assumes a particular box
   geometry)
3. `graphene.raw.inc`: Molecule definition without the additional
   cross-box bonds. Contains unterminated sp<sup>2</sup> carbon atoms.
4. `graphene.omd`: initial OpenMD file for starting a simulation
5. `propylene.omd`: Molecule definition include file for the propylene monomers
6. `propylene.xyz`: A skeletal propylene xyz file for use with Packmol
7. `system.pack`: A [packmol](http://www.ime.unicamp.br/~martinez/packmol/)
   input script that places propylene molecules inside bounds of the
   graphene sheets.

To create the initial configurations, we typically run:

`packmol < system.pack`

which will create a `system.xyz` file containing all of the propylene moleucles.
To get this into a reasonable OpenMD starting configuration, we would run:

`atom2omd -ixyz system.xyz`

This creates `system.omd`, which has periodic box guessed from the bounding box
of the solvent molecules, and this box is generally not the same size as the
box containing the graphene sheets.  To successfully combine the monomers with
the graphene sheets, the `system.omd` file must be edited to modify the Hmat
line to read:

~~~~
        Hmat: {{ 22.23, 0, 0 }, { 0, 42.78, 0 }, { 0, 0, 50 }}
~~~~

Then to combine the propylene monomers with the graphene sheets, eliminating
molecules that overlap, we would run:

`omd-solvator -u graphene.omd -v system.omd -r 3.5 -o combined.omd -n 360 -p 3`

Following this, we typically edit the `combined.omd` file to include the
molecule definitions, set the force field, and set various simulation
parameters:

~~~~
#include "graphene.inc"
#include "propylene.omd"

component{
  type = graphene;
  nMol = 2;
}
component{
  type = propylene;
nMol = 23;
}

forceField = "graphene";
ensemble = NVT;
cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "damped";
cutoffRadius = 9;
dampingAlpha = 0.2;
targetTemp = 300;
tauThermostat = 1000;
dt = 1.0;
runTime = 1e3;
sampleTime = 100;
statusTime = 10;
~~~~
