# Utility scripts to aide in construction of OpenMD (.omd) files

Much of the magic outlined below is described in greater detail in the
file **runMe.in** in this directory.

Here we outline a variety of types of systems able to be constructed
using **simpleBuilder**.

+ fcc lattices of a given unit cell
+ nanospheres
+ spherically-capped nanorods
+ pentagonal nanorods
+ icosohedra
+ cuboctahedra
+ truncated cube particles

Note: to use simpleBuilder, you must have the *<MetaData>* section of
an OpenMD (.omd) file as input. Examples of such *<MetaData>* sections
are **one_component.omd**, **three_component.omd**, **gold.omd**, and
**bimetallic.omd**.

## Examples of builder commands shown in **runMe.in** are:

1. Builds an FCC lattice from the <MetaData> block in one_component.omd
2. Builds an FCC lattice from the <MetaData> block in three_component.omd
3. Builds a spherical nanoparticle (FCC) from the <MetaData> block in gold.omd
4. Builds a random alloy spherical nanoparticle (FCC) from the <MetaData>
5. Builds a Au(core)-Ag(shell) spherical nanoparticle (FCC) from the
   <MetaData> block in bimetallic.omd
6. Reverses example 5 by building a Ag(core)-Au(shell) spherical
   nanoparticle. Uses the same <MetaData> block from bimetallic.omd
7. Builds a Au(core)-Ag(shell) spherical nanoparticle (FCC) from the
   <MetaData> block in bimetallic.omd
8. Builds a random alloy spherical nanoparticle with 30% vacancies
   using the <MetaData> block in bimetallic.omd
9. Builds a spherically-capped nanorod (FCC) from the <MetaData> block in gold.omd
10. Builds a pentagonal nanorod from the <MetaData> block in gold.omd
11. Builds a Mackay icosahedral nanoparticle from the <MetaData> block in gold.omd
12. Builds a regular decahedral nanoparticle from the <MetaData> block in gold.omd
13. Builds a ino-decahedral nanorod from the <MetaData> block in gold.omd
14. Builds a cuboctahedral particle from the <MetaData> block in gold.omd
15. Builds a truncated cube particle from the <MetaData> block in gold.omd
