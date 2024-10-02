# Getting Started with OpenMD

OpenMD is a simulation engine that is started from the command line, so we're going to assume that you are reasonably familiar with unix command line environments. This tutorial will focus on using OpenMD assuming it has already been installed on your machine. For directions on installing, checkout [INSTALL.md](docs/INSTALL.md) before continuing.

Once OpenMD has been installed completely, and the `openmd` program is already in your path, you don't need to do anything special to set up your environment.  However, if you are just test-driving, and built it locally without installing it, you'll need to specify a few variables before continuing.  For example, if the code was unpacked into the directory `~/OpenMD-3.1` and you built in the subdirectory `~/OpenMD-3.1/build`, then you might need to do this:

```bash
export FORCE_PARAM_PATH="~/OpenMD-3.1/forceFields"
export PATH="${PATH}:~/OpenMD-3.1/build/bin"
```
This will set your path to include the directory where openmd was built, and will reference the default force field files.  If you use `csh` or `tcsh` as your shell, the equivalent of those lines are:

```bash
setenv FORCE_PARAM_PATH ~/OpenMD-3.1/forceFields
setenv PATH ${PATH}:~/OpenMD-3.1/build/bin
```

## Sample simulations

In the `samples` subdirectory, there are a wide range of example simulations that have already been set up as demonstrations of the various capabilities of `OpenMD`.  We're going to start with a relatively simple case, a simulation of liquid water that uses the SPC/E water model for inter-atomic interactions.  Let's start by changing into the samples directory:

```bash
cd ~/OpenMD-3.1/samples
```

Now, let's go to the `water/liquid` sample directory and run a simulation: 

```bash
cd water/liquid
openmd spce.omd                     # run the sample
```

As the simulations is progressing, you should see a progress bar indicating the time remaining in the simulation.
If you compiled the parallel-aware version of OpenMD linked to MPI libraries, you can speed up the simulations by using multiple processors:

```bash
mpirun -np 4 openmd_MPI spce.omd   # run the sample, but with 4 processors
```

At the end of your simulation, you'll see a print out of a number of relevant thermodynamic quantities. You'll also notice a few new files in your sample directory, namely: 

- `spce.dump` - the simulation trajectory file, containing multiple snapshots of atomic positions and velocities
- `spce.eor` - the 'end-of-run' for the simulation, or the final configuration.
- `spce.stat` - a 'statistics' file which records thermodynamically-interesting quantities like temperature, pressure, volume, potential, kinetic, and total energy
- `spce.report` - statistical means and confidence intervals for the same quantities recorded in the `.stat` file

In the next section, we'll explore our recommended ways of parsing and analyzing these file types.

## Visualization Tools

Suppose we want to see how the total, potential, and kinetic energies, the temperature, pressure, and volume of the box evolved as the simulation progressed. We can visualize the `.stat` file to look at these quantities:
```bash
xmgrace -nxy spce.stat         # visualize simulation properties as they evolve during the simulation
```
Here, we're  using the `xmgrace` package to plot our data, but `spce.stat` is just a tab-separated text file, so you can look at it in many different ways.  The first few lines of the `.stat` file give information about the version of the code that was used as well as the thermodynamic quantities that were recorded (and the units used).

To visualize the trajectory or the contents of the `.dump` file after the simulation has completed we use a utility to generate an xyz file:
```bash
Dump2XYZ -i spce.dump -m -b   # create an xyz file from the trajectory, mapping back to the periodic box
jmol spce.xyz                 # visualize the atomic coordinates in the dump file using Jmol
vmd spce.xyz                  # visualize the atomic coordinates in the dump file using VMD
```
The `Dump2XYZ` program is quite flexible. Note that some simulations will use atom types that are derived from "base" atom types, so we often use the `-b` flag to map atoms back to their `Base` types (usually element names). There are many trajectory visualization tools, but `Jmol` and `VMD` are relatively common.

## Statistical Analysis

At this point we have a completed simulation (256 SPC/E water molecules, simulated in NVE conditions for 10 ps) and we've tracked its progress, what's next?

Suppose you want to dive a bit deeper into your data analysis and want to compute either static properties or time correlation functions.  Once the trajectory has been generated, you can use two other utilities to study many derived properties.  

For example, this will let you compute and visualize the Oxygen-Oxygen pair distribution function, g(r), from the trajectory data, where we specify or 'select' which atoms we want to look at:
```bash
StaticProps -i spce.dump --gofr --sele1="select O_SPCE" --sele2="select O_SPCE" 
xmgrace spce.gofr              # visualize the g(r) calculated by StaticProps
```

And this will let you compute a time correlation function, in this case the mean squared displacement (for all atoms) and visualize it:
```bash
DynamicProps -i spce.dump --rcorr  # Calculate the mean squared displacement <|r(t) - r(0)|^2>
xmgrace spce.rcorr                 # visualize the MSD calculated by DynamicProps
```

Both `StaticProps` and `DynamicProps` can compute many properties. You can get a list of them using the `-h` flag:

```
StaticProps -h
DynamicProps -h
```

## Next Steps

So what's next? By now you've familiarized yourself with a simple OpenMD simulation and analyzed the result both visually and statistically. OpenMD provides a number of other [samples](#list-of-samples-included-with-openmd) which can help guide you on what kind of projects are possible in our engine. From there you're likely ready to start creating your own simulations and run experiments on them. The [OpenMD manual](https://github.com/OpenMD/OpenMD/blob/main/docs/OpenMDmanual.pdf) is a good place to learn the in-depth theory, design decisions, and configuration options available to you. If you're also looking to contribute to the codebase, check out our [contributing guidelines](.github/CONTRIBUTING.md) and [API documentation](https://openmd.org/wp-content/docs/code/index.html).

### List of Samples included with OpenMD

- [DR-EAM](samples/DR-EAM/README.md) - various samples using the Density-Readjusting Embedded Atom Model for metals
- [LangevinHull](samples/LangevinHull/README.md) - samples using constant Pressure and Temperature for non-periodic systems
- [Madelung](samples/Madelung/README.md) - computes Madelung constants using different approximate electrostatic methods
- [RNEMD](samples/RNEMD/README.md) - Reverse Non-Equilibrium Molecular Dynamics, in bulk and at interfaces
- [air](samples/air/README.md) - a simple mixture of diatomic nitrogen, diatomic oxygen, and argon atoms
- [aqueousIons](samples/aqueousIons/README.md) - some salt water boxes showing how to run with Li-Song-Merz 12-6 and 12-6-4 models for ions
- [argon](samples/argon/README.md) - simple liquid argon boxes using Lennard-Jones potentials
- [bond-order](samples/bond-order/README.md) - demonstration files for testing bond orientational order parameters in different structures
- [builders](samples/builders/README.md) - demonstrations of how to create initial structures using the builder codes.
- [graphene](samples/graphene/README.md) - models for graphene sheets and nanoporous graphene membranes
- [metals](samples/metals/README.md) - metal structures (bulk and surface) using different interaction models.
- [liquid](samples/liquid/README.md) - build and equilibrate a liquid sample using methanol as a demonstration molecule
- [minimizer](samples/minimizer/README.md) - small, simple samples showing how to invoke the minimizer
- [nanoparticles](samples/nanoparticles/README.md) - build and equilibrate spherical and icosahedral nanoparticles
- [protein](samples/protein/README.md) - build and equilibrate a small pentapeptide starting from a PDB file
- [water](samples/water/README.md) - multiple water models, and initial structures including ice, ice surfaces, liquids, and dimers
- [zcons](samples/zcons/README.md) - tests for Z-constraints to keep molecules fixed at specific values of one coordinate
- [zeolite](samples/zeolite/README.md) - diffusion in a ZSM5 structure simulated using the CLAYFF force field
