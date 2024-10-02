
# Samples included with OpenMD

- [DR-EAM](DR-EAM/README.md) - various samples using the Density-Readjusting Embedded Atom Model for metals. Subdirectories include:

    + [fieldScreening](DR-EAM/fieldScreening/) - internal screening of a metal slab in an external field
    + [imageCharge](DR-EAM/imageCharge) - response of fluctuating charge atoms to approach of an ion
    + [orderedStructures](DR-EAM/orderedStructures/) - charge distributions in  $L1_{2}$ and $L1_{0}$ alloy structures

- [LangevinHull](LangevinHull/README.md) - samples using constant Pressure and Temperature for non-periodic systems
- [Madelung](Madelung/README.md) - computes Madelung constants using different approximate electrostatic methods.  Also includes subdirectories for:

    + [dipoles](Madelung/dipoles/README.md) - electrostatic energies for dipolar crysals
    + [quadrupoles](Madelung/quadrupoles/README.md) - electrostatic energies for quadrupolar crystals

- [RNEMD](RNEMD/README.md) - Reverse Non-Equilibrium Molecular Dynamics, in bulk and at [interfaces](RNEMD/interfaces)
- [air](air/README.md) - a simple mixture of diatomic nitrogen, diatomic oxygen, and argon atoms
- [aqueousIons](aqueousIons/README.md) - some salt water boxes showing how to run with Li-Song-Merz 12-6 and 12-6-4 models for ions
- [argon](argon/README.md) - simple liquid argon boxes using Lennard-Jones potentials
- [bond-order](bond-order/README.md) - demonstration files for testing bond orientational order parameters in different structures
- [builders](builders/README.md) - demonstrations of how to create initial structures using the builder codes.
- [graphene](graphene/README.md) - models for graphene sheets and nanoporous graphene membranes
- [metals](metals/README.md) - metal structures (bulk and surface) using different interaction models.  Subdirectories include:

    + [EAM](metals/EAM/README.md) - metal simulations using the Embedded Atom Method
    + [Sutton-Chen](metals/Sutton-Chen/README.md) - metal simulations using the Sutton-Chen potential
    + [surfaces](metals/surfaces/README.md) - building and running a Au(111) / water interface and a Pt(557) surface with adsorbed carbon monoxide

- [liquid](liquid/README.md) - build and equilibrate a liquid sample using methanol as a demonstration molecule
- [minimizer](minimizer/README.md) - samples showing how to invoke the minimizer
- [nanoparticles](nanoparticles/README.md) - build and equilibrate spherical and icosahedral nanoparticles
- [protein](protein/README.md) - build and equilibrate a small pentapeptide starting from a PDB file
- [water](water/README.md) - multiple water models, and initial structures including ice, ice surfaces, liquids, and dimers. Subdirectories include:

    + [dimer](water/dimer/README.md) - simple annealing of dimer structures
    + [iceCrystals](water/iceCrystals/README.md) - build ice crystals using a basis set of orthorhombic building blocks
    + [ice_surfaces](water/ice_surfaces/README.md) - pre-equilibrated ice/water interfaces
    + [liquid](water/liquid/README.md) - liquid simulations using 5 different water models
    + [tip4p-fq](water/tip4p-fq/README.md) - the TIP4P-FQ (fluctuating charge) model for water

- [zcons](zcons/README.md) - tests for Z-constraints to keep molecules fixed at specific values of one coordinate
- [zeolite](zeolite/README.md) - diffusion of water in a ZSM5 structure simulated using the CLAYFF force field
