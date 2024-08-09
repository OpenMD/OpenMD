# Builders

## Background Information

This is a collection of sample commands that can be used to build
OpenMD start files.  In OpenMD, the start files have a `<MetaData>`
block to give information about the kind of simulation being performed.
The start files also contain at least one `<Snapshot>` block which contains
information about the instantaneous configuration.

One of the difficult tasks in using any simulation program is figuring
out how to format the start file correctly.  OpenMD includes a set of
"builder" programs to make that process a bit less painful.
 
Here we outline a variety of types of systems able to be constructed
using some of the builder tools included with `OpenMD`

+ fcc lattices of a given unit cell (using **simpleBuilder**)
+ nanospheres (using **nanoparticleBuilder**)
+ spherically-capped nanorods (using **nanorodBuilder**)
+ pentagonal nanorods (using **nanorod_pentBuilder**)
+ icosohedra (using **icosahedralBuilder**)
+ cuboctahedra (using **icosahedralBuilder**)
+ truncated cube particles (using **icosahedralBuilder**)

Note: to use the builders, you must have the `<MetaData>` section of
an OpenMD (.omd) file as input. Examples of such `<MetaData>` sections
are in `one_component.omd`, `three_component.omd`, `gold.omd`, and
`bimetallic.omd`.

Many of the following examples also include a *thermalizer* command which resamples the velocities from a Maxwell-Boltzmann distribution set to a desired temperature.

## Instructions

### Example 1

Builds an FCC lattice from the `<MetaData>` block in `one_component.omd`
Uses 5 unit cells in each direction, a density of 1.0 g / cm^3, and
places the output (which can be used to start an OpenMD job) in 
FCC.omd

Note that builders will rewrite the number of molecules in each component
to match the number of lattice sites.

The thermalizer command takes the FCC.omd file and resamples the velocities
from a Maxwell-Boltzmann distribution set to 100K:

```bash
simpleBuilder -o FCC.omd --nx=5 --ny=5 --nz=5 --density=1.0 one_component.omd
thermalizer -o FCC-100K.omd -t 100 FCC.omd
Dump2XYZ -i FCC-100K.omd
```

### Example 2

Builds an FCC lattice from the `<MetaData>` block in `three_component.omd`
uses 4 unit cells in each direction, a density of 1.0 g / cm^3, and
molFractions of 0.4, 0.4, and 0.2 for the three components.  Places
the output (which can be used to start an OpenMD job) in `random_FCC.omd`

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
randomBuilder -o random_FCC.omd --nx=4 --ny=4 --nz=4 --density=1.0 --molFraction=0.4 --molFraction=0.4 three_component.omd
thermalizer -o random_FCC-100K.omd -t 100 random_FCC.omd
Dump2XYZ -i random_FCC-100K.omd
```

 ### Example 3
 
 Builds a spherical nanoparticle (FCC) from the `<MetaData>` block in `gold.omd`
 using a particle radius of 30 Angstroms, and a lattice constant of 4.09
 angstroms. Places the output (which can be used to start an OpenMD job) in 
 `gold_sphere.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.
 
```bash
nanoparticleBuilder -o gold_sphere.omd --radius=30.0 --latticeConstant=4.09 gold.omd
thermalizer -o gold_sphere-500K.omd -t 500.0 gold_sphere.omd
Dump2XYZ -i gold_sphere-500K.omd
```

 ### Example 4

 Builds a random alloy spherical nanoparticle (FCC) from the `<MetaData>` 
 block in `bimetallic.omd` using a particle radius of 30 Angstroms, a 
 lattice constant of 4.09 angstroms, and a mole fraction for the gold of 0.4.
 Places the output (which can be used to start an OpenMD job) in 
 `Au_Ag_alloy.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
nanoparticleBuilder -o Au_Ag_alloy.omd --radius=30.0 --latticeConstant=4.09 --molFraction=0.4 bimetallic.omd
thermalizer -o Au_Ag_alloy-600K.omd -t 600 Au_Ag_alloy.omd
Dump2XYZ -i Au_Ag_alloy-600K.omd
```

### Example 5

 Builds a Au(core)-Ag(shell) spherical nanoparticle (FCC) from the `<MetaData>` 
 block in `bimetallic.omd` using a particle radius of 25 Angstroms, a 
 lattice constant of 4.09 angstroms, and a core radius for the gold atoms 
 of 12.5 angstroms. Places the output (which can be used to start an 
 OpenMD job) in `Au-core-Ag-shell.omd `

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
nanoparticleBuilder -o Au-core-Ag-shell.omd --radius=30.0 --latticeConstant=4.09 --shellRadius=12.5 bimetallic.omd
thermalizer -o Au-core-Ag-shell-800K.omd -t 800.0 Au-core-Ag-shell.omd
Dump2XYZ -i Au-core-Ag-shell-800K.omd
```

### Example 6

 Reverses example 5 by building a Ag(core)-Au(shell) spherical nanoparticle.
 Uses the same `<MetaData>` block from `bimetallic.omd`, 
 a particle radius of 25 Angstroms, a lattice constant of 4.09 angstroms, 
 and a core radius for the silver atoms of 12.5 angstroms.  
 Places the output (which can be used to start an OpenMD job) in 
 `Ag-core-Au-shell.omd` 

 Note that the last radius in Example 5 was taken as the particle radius,
 but since the components are reversed in this example, both are specified:
 
```bash
nanoparticleBuilder -o Ag-core-Au-shell.omd --radius=30.0 --latticeConstant=4.09 --shellRadius=30.0,12.5 bimetallic.omd
thermalizer -o Ag-core-Au-shell-800K.omd -t 800.0 Ag-core-Au-shell.omd
Dump2XYZ -i Ag-core-Au-shell-800K.omd
```

### Example 7

 Builds a Au(core)-Ag(shell) spherical nanoparticle (FCC) from the `<MetaData>` 
 block in `bimetallic.omd` using a particle radius of 25 Angstroms, a 
 lattice constant of 4.09 angstroms, and a core radius for the gold atoms 
 of 12.5 angstroms. Places the output (which can be used to start an 
 OpenMD job) in `Au-core-Ag-shell.omd` 

 This example also introduces 70% vacancies in a 6 angstrom radial band
 around the bimetallic interface:

```bash
nanoparticleBuilder -o vacancy_interface.omd --radius=20.0 --latticeConstant=4.09 --shellRadius=12.5 --vacancyPercent=70 --vacancyInnerRadius=9.5 --vacancyOuterRadius=15.5 bimetallic.omd
thermalizer -o vacancy_interface-800K.omd -t 800 vacancy_interface.omd
Dump2XYZ -i vacancy_interface-800K.omd
```

### Example 8

 Builds a random alloy spherical nanoparticle with 30% vacancies using the
 `<MetaData>` block in `bimetallic.omd`, a particle radius of 30 Angstroms, a 
 lattice constant of 4.09 angstroms, and a mole fraction for the gold of 0.4.
 Places the output (which can be used to start an OpenMD job) in 
 `vacancy_alloy.omd`

```bash
nanoparticleBuilder -o vacancy_alloy.omd --radius=30.0 --latticeConstant=4.09 --molFraction=0.4 --vacancyPercent=80 bimetallic.omd
thermalizer -o vacancy_alloy-900K.omd -t 900 vacancy_alloy.omd
Dump2XYZ -i vacancy_alloy-900K.omd
```

### Example 9

 Builds a spherically-capped nanorod (FCC) from the `<MetaData>` block in `gold.omd`
 using a nanorod radius of 20 Angstroms, a length of 50 Angstroms and a lattice constant of 4.08
 angstroms. Places the output (which can be used to start an OpenMD job) in 
 `gold_fccrod.omd `

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
nanorodBuilder -o gold_fccrod.omd --radius=20.0 --length=50.0 --latticeConstant=4.08 gold.omd
Dump2XYZ -i gold_fccrod.omd
```

### Example 10

 Builds a pentagonal nanorod from the `<MetaData>` block in gold.omd
 using a nanorod radius of 15 Angstroms, a length of 64 Angstroms and a lattice constant of 4.08
 angstroms. Places the output (which can be used to start an OpenMD job) in 
 `gold_pentrod.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
nanorod_pentBuilder -o gold_pentrod.omd --radius=15.0 --length=64.0 --latticeConstant=4.08 gold.omd
Dump2XYZ -i gold_pentrod.omd
```

### Example 11

 Builds a Mackay icosahedral nanoparticle from the `<MetaData>` block in `gold.omd`
 using a 8 shells, and a lattice constant of 4.08 angstroms.
 Places the output (which can be used to start an OpenMD job) in 
 `gold_ico.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
icosahedralBuilder --ico -o gold_ico.omd --shells=8 --latticeConstant=4.08 gold.omd
thermalizer -o gold_ico_300K.omd -t 300 gold_ico.omd
Dump2XYZ -i gold_ico_300K.omd
```

### Example 12

 Builds a regular decahedral nanoparticle from the `<MetaData>` block in `gold.omd`
 using a 10 shells, and a lattice constant of 4.08 angstroms.
 Places the output (which can be used to start an OpenMD job) in 
 `gold_deca.omd `

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
icosahedralBuilder --deca -o gold_deca.omd --shells=10 --latticeConstant=4.08 gold.omd
thermalizer -o gold_deca_300.omd -t 300 gold_deca.omd
Dump2XYZ -i gold_deca_300.omd
```

### Example 13

 Builds a ino-decahedral nanorod from the `<MetaData>` block in `gold.omd`
 using a 10 shells, 5 atoms along the twin boundary, 100 atoms along the 
 column axis, and a lattice constant of 4.08 angstroms.
 Places the output (which can be used to start an OpenMD job) in 
 `penta_rod.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
icosahedralBuilder --ino --columnAtoms=100 --twinAtoms=5 --shells=10 -d 4.08 -o penta_rod.omd gold.omd
thermalizer -o gold_penta_rod_300.omd -t 300 penta_rod.omd
Dump2XYZ -i gold_penta_rod_300.omd
```

### Example 14
 Builds a cuboctahedral particle from the `<MetaData>` block in `gold.omd`
 using a 10 unit cells, and a lattice constant of 4.08 angstroms.
 Places the output (which can be used to start an OpenMD job) in 
 `cuboctahedron.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```bash
icosahedralBuilder --cuboctahedron --lattice="FCC" --unitCells=10 -d 4.08 -o cuboctahedron.omd gold.omd
thermalizer -o cuboctahedron_300.omd -t 300 cuboctahedron.omd
Dump2XYZ -i cuboctahedron_300.omd
```

### Example 15

 Builds a truncated cube particle from the `<MetaData>` block in
 `gold.omd` using a 10 unit cells, 4 truncated {111} planes, and a lattice
 constant of 4.08 angstroms.
 Places the output (which can be used to start an OpenMD job) in 
 `truncatedCube.omd` 

 Note that builders will rewrite the number of molecules in each component
 to match the number of lattice sites.

```
icosahedralBuilder --truncatedCube --lattice="FCC" --unitCells=10 -d 4.08 -o truncatedCube.omd gold.omd
thermalizer -o truncatedCube_300.omd -t 300 truncatedCube.omd
Dump2XYZ -i truncatedCube_300.omd
```

## Expected Output

All of the examples listed above generate both an `omd` files (and `xyz` files).  The `xyz` files are easy to visualize using **Jmol** or **VMD**:
```bash
jmol truncatedCube_300.xyz
```
Likewise all of the `omd` files can be used to start OpenMD simulations:
```bash
openmd gold_sphere-500K.omd
```
