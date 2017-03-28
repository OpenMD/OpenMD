This sample illustrates a simulation of united-atom (UA) propylene monomers
confined between and all-atom (AA) representation of graphene sheets.  The force
field is based partially on OPLS-AA, and partially on TraPPE-UA.

Propylene parameters from: C.D. Wick, M.G. Martin, and J.I. Siepmann,
"Transferable potentials for phase equilibria. 4. United-atom description of
linear and branched alkenes and of alkylbenzenes," *J. Phys. Chem. B*, **104**,
8008-8016 (2000).

Note that the molecule definition in `graphene.inc` includes bonds that span the
box boundaries.  If you want to extend to larger boxes, start with the molecule
definition in `graphene.raw.inc` and then add bonds that span the boundaries of
the new box.

Files included:

1. `graphene.frc`: Force field file
2. `graphene.inc`: Molecule definition include file for a graphene sheet (includes cross-box bonding, so assumes a particular box geometry)
3. `graphene.raw.inc`: Molecule definition without the additional cross-box bonds. Contains unterminated sp<sup>2</sup> carbon atoms.
4. `graphene.omd`: initial OpenMD file for starting a simulation
