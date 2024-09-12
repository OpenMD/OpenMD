# Simulating metals using OpenMD

`OpenMD` implements a number of related potentials that describe
bonding in transition metals. These potentials have an attractive
interaction which models the embedding of a positively charged
pseudo-atom core in the electron density due to the free valance sea
of electrons created by the surrounding atoms in the system.  A
pairwise part of the potential describes the interaction of the
positively charged metal core ions with one another.  This family of
potentials generally has the form:

$$
V  =  \sum_{i} F_{i}\left[\rho_{i}\right] + \sum_{i} \sum_{j \neq i}
\phi_{ij}(\mathbf{r}_{ij})
$$

where $F_{i}$ is an embedding functional that approximates the energy
required to embed a positively-charged core ion $i$ into a local
electron density given by $\rho_{i}$,

$$
\rho_{i}   =  \sum_{j \neq i} f_{j}(\mathbf{r}_{ij}),
$$

Since the density at site $i$ ($\rho_i$) must be computed before the
embedding functional can be evaluated, transition metal potentials
generally require two loops through the atom pairs to compute the
inter-atomic forces.

The pairwise portion of the potential, $\phi_{ij}$, is usually a
repulsive interaction between atoms $i$ and $j$.

The included subdirectories illustrate a few examples of metals simulated
using potentials in this family.

### Sample Directories
+ [EAM](EAM/README.md) - using the Embedded Atom Method
+ [Sutton-Chen](Sutton-Chen/README.md) - using the Sutton-Chen potential
+ [minimizer](minimizer/README.md) - optimizing a structure
+ [surfaces](surfaces/README.md) - creating and simulating metal surfaces
