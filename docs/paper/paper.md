---
title: 'OpenMD: A parallel molecular dynamics engine for complex systems and interfaces'
tags:
  - C++
  - Chemistry
  - Molecular Dynamics
authors:
  - name: Cody R. Drisko
    orcid: 0009-0006-3968-5088
    affiliation: 1
  - name: Hemanta Bhattarai
    orcid: 0000-0002-3573-9716
    affiliation: 2    
  - name: Christopher J. Fennell
    orcid: 0000-0001-8963-4103
    affiliation: 3
  - name: Kelsey M. Stocker
    orcid: 0000-0002-1799-393X
    affiliation: 4
  - name: Charles F. Vardeman II
    orcid: 0000-0003-4091-6059
    affiliation: 5
  - name: J. Daniel Gezelter
    orcid: 0000-0002-2935-3163
    corresponding: true
    affiliation: 1
affiliations:
 - name: Department of Chemistry and Biochemistry, University of Notre Dame, Notre Dame, IN 46556
   index: 1
 - name: Department of Physics, Goshen College, Goshen, IN 46526
   index: 2      
 - name: Department of Chemistry, Oklahoma State University, Stillwater, OK, 74078
   index: 3
 - name: Department of Biochemistry, Chemistry, Environment, and Physics, Suffolk University, Boston, MA 02108
   index: 4
 - name: Center for Research Computing, University of Notre Dame, Notre Dame, IN 46556
   index: 5   
date: 25 June 2024
bibliography: paper.bib
---

# Summary

Molecular dynamics (MD) simulations help bridge the gap between quantum mechanical calculations, which trade computational resources and system size for chemical accuracy, and continuum simulations, which utilize bulk and surface material properties to make predictions. MD provides the ability to study emergent properties and dynamics for relatively large systems ($\sim 10^6$ atoms) over a span of nano- to micro-seconds. MD has been used for simulations of complex systems such as liquids, proteins, membranes, nanoparticles, interfaces, and porous materials, and the approximations used in MD simulations allow us to reach experimentally-relevant time and length scales. A molecular dynamics engine leverages classical equations of motion to evolve atomic or molecular coordinates according to a well-defined potential energy function, known as a force field, which is a function of those atomic coordinates. There are a number of high quality molecular dynamics engines that specialize in materials [@LAMMPS] or biomolecular [@CHARMM; @AmberTools] simulations or are models of computational efficiency [@GROMACS; @OpenMM]. In this paper, we provide background on an open source molecular dynamics engine, `OpenMD`, which was just released into version 3.1.

`OpenMD` is capable of efficiently simulating a variety of complex systems using standard point-particle atoms, as well as atoms with orientational degrees of freedom (e.g. point multipoles and coarse-grained assemblies), and atoms with additional degrees of freedom (e.g. fluctuating charges). It provides a test-bed for new molecular simulation methodology, while being efficient and easy to use. Liquids, proteins, zeolites, lipids, inorganic nanomaterials, transition metals (bulk, flat interfaces, and nanoparticles), alloys, solid-liquid interfaces, and a wide array of other systems have all been simulated using this code. `OpenMD` works on parallel computers using the Message Passing Interface (MPI), and is packaged with trajectory analysis and utility programs that are easy to use, extend, and modify.

From the beginning, `OpenMD` has been an open source project, and has been maintained in accordance with the FAIR (findable, accessible, interoperable, and reusable) guidelines [@Wilkinson2016]. It uses a robust meta-data language that is tightly integrated with input (`.omd`) and trajectory (`.dump`) files, providing a standardized, human-readable way to completely describe molecular systems and simulation parameters. This allows `OpenMD` simulations to be easily reproduced, modified, and reused by others. The `<MetaData>` section of these files also serve as a form of rich, machine-actionable meta-data, clearly specifying the composition of the molecular system in a manner that promotes interoperability with other software tools. Additionally, all data files are stamped with the code revision that generated that data, greatly aiding reproducibility.

# Statement of need

`OpenMD` builds on the foundations of the Object-Oriented Parallel Simulation Engine (`OOPSE`) program [@Meineke2005]. Entirely rewritten in modern C++, it extends `OOPSE` in a number of key areas. One of the added capabilities is Reverse Non-Equilibrium Molecular Dynamics (RNEMD), a family of algorithms which impose a non-physical flux on a system and use linear response theory to compute transport properties as the system approaches steady state. The goal of RNEMD methods is to calculate the relevant transport property ($\lambda$) that connects the flux ($\textbf{J}$) and driving force ($\nabla X$) according to the generalized equation,
\begin{equation}\label{eq:linearConstitutive}
  \textbf{J} = - \lambda \nabla X.
\end{equation}

`OpenMD` is also capable of performing condensed phase simulations without the use of periodic boundary conditions. To do so, an external pressure and temperature bath is applied to atoms comprising the system's convex hull, rather than the interior region. This method, the Langevin Hull, allows for constant pressure, temperature, or isobaric-isothermal (NPT) simulations of explicitly non-periodic molecular systems. Other major developments are the inclusion of advanced real-space electrostatics for point multipoles, and polarizable force fields using fluctuating charges or fluctuating electron densities.

## Reverse Non-Equilibrium Molecular Dynamics

RNEMD methods impose a non-physical (heat, momentum, or particle) flux between different regions of the simulation. In response, the system develops a temperature, velocity, or concentration gradient between the two regions, and the linear coefficient connecting applied flux and measured gradient is a transport property of the material. Since the amount of the applied flux is known exactly in RNEMD, and the measurement of gradients is generally straightforward, imposed-flux methods typically take shorter simulation times to obtain converged results for transport properties, when compared with equilibrium MD or forward-NEMD approaches. If an interface lies between the two regions, these methods can also provide *interfacial* transport coefficients by mapping any spatial discontinuities in temperature, velocity, or concentration with the applied flux.

Non-equilibrium molecular dynamics is a well-developed area of research, and `OpenMD` supports many different RNEMD algorithms. The first is the original "swapping" approach by M&uuml;ller-Plathe [@Muller-Plathe1997; @Muller-Plathe1999]. Here, the entire momentum vectors of two particles in separate slabs may be exchanged to generate a thermal flux. Occasionally, non-ideal Maxwell-Boltzmann distributions will develop in velocity profiles using this approach [@Tenney2010]. `OpenMD` also introduces a number of new algorithms which extend the capabilities of RNEMD.

Rather than using momentum exchanges between individual particles in each region, the Non-Isotropic Velocity Scaling (NIVS) algorithm applies velocity scaling to all of the selected particles in both regions [@Kuang2010]. NIVS was shown to be very effective at computing thermal conductivities, but is not suitable for imposing a momentum flux or for computing shear viscosities. However, simultaneous velocity shearing and scaling (VSS) exchanges between the two regions remove all of these limitations [@Kuang2012]. The VSS-RNEMD method yields a simple set of equations which satisfy energy and linear momentum conservation constraints, while simultaneously imposing a desired flux between the two regions. The VSS approach is versatile in that it may be used to implement both thermal and shear transport either separately or simultaneously. `OpenMD` is also capable of leveraging the VSS method in *non-periodic simulations*, in which the regions have been generalized to concentric regions of space [@Stocker2014], allowing for simulations of heat transport away from nanostructures. In the following section, we explore the algorithm that makes non-periodic boundary simulations possible in `OpenMD`.

Another novel RNEMD algorithm allows for particle positions to be swapped between RNEMD regions resulting in an applied *particle flux* [@Drisko2024]. The scaled particle flux (SPF) method accurately calculates diffusion coefficients while maintaining energy and linear momentum constraints, and can map the temperature dependence of diffusion when used in tandem with a thermal flux in VSS-RNEMD. SPF-RNEMD has also been applied to interfacial systems of nanoporous graphene in a molecular fluid. In this case, permeabilities were computed by imposing a molecular flux between regions on opposite sides of the membrane and measuring both the hydraulic and osmotic pressure that develops as a result of this flux.

## Langevin Hull

In many molecular simulations, systems have near-uniform compressibility, and `OpenMD` implements a range of integrators to sample the isothermal-isobaric (NPT) ensemble using the Nose&#769;-Hoover-Andersen equations of motion. These integrators implement various forms of affine scaling to provide isotropic (NPTi) or fully-flexible (NPTf) scaling motions of the periodic box [@Andersen1980; @Hoover1986; @Sturgeon2000]. Additional constant pressure integrators use restricted affine scaling to enforce constant surface area (NPAT), constant surface tension (NP&gamma;T), or even orthorhombic box geometries (NPTxyz). For systems comprising materials of different compressibilities, such as a solvated nanoparticle, scaled coordinate transformations may cause numerical instability or poor volume sampling depending on the strength of the applied barostat. Users may also wish to represent systems without the explicit periodicity required by box-scaling constant pressure methods. For example, proteins may be in artificially crowded environments if periodic box simulations are required.

To address these problems, `OpenMD` implements a method called the Langevin Hull which samples the isothermal-isobaric (NPT) ensemble for non-periodic systems [@Vardeman2011]. The method, based on Langevin dynamics, couples an external bath at constant pressure, temperature, and effective solvent viscosity to the atomic sites on a dynamically-computed convex hull. This convex hull is defined as the set of facets that have no concave corners at an atomic site [@Edelsbrunner1994]. The hull is computed using Delaunay triangulation between coplanar neighbors [@Delaunay1934; @Lee1980]. These computations are performed by the external `Qhull` library [@Barber1996], and are computed each time step, allowing molecules to move freely between the inner region and outer hull. Atoms in the interior evolve according to Newtonian mechanics. The equations of motion for sites on the hull,
\begin{align} \label{eq:equationsOfMotion}
  m_i \dot{\textbf{v}}_i (t) &= - \nabla_i U +\textbf{F}_i^{\textrm{ext}} \\ \nonumber
  &= - \nabla_i U + \sum_{f} \frac{1}{3} \left( -\hat{\mathbf{n}}_f P A_f - \Xi_f (t) \left( \frac{1}{3} \sum_{i = 1} \textbf{v}_i \right) + \textbf{R}_f (t) \right)
\end{align}
include additional forces from the external bath. Each facet $f$ on the convex hull has a contribution from a pressure bath acting in proportion to the facet's surface area and in the direction of the surface normal $(- \hat{\mathbf{n}}_f P A_f)$. The facets of the hull are also in contact with an implicit solvent which provides a drag on the velocity of the facet according to an approximate resistance tensor, $(- \Xi_f \mathbf{v}_f)$ and which also kicks the facet via a Gaussian random force $(\mathbf{R}_f)$.

Computation of a convex hull is $\mathcal{O}(N \log N)$ for sequential machines and remains the performance bottleneck for parallelization. In parallel, the global convex hull is computed using the union of sites from all local (processor-specific) hulls. Testing and validation for this method were carried out on three unique systems, a gold nanoparticle and an SPC/E water droplet (both with uniform compressibilities), and a gold nanoparticle solvated in SPC/E water (non-uniform compressibility), shown in Fig. \ref{fig:LHull}. This method was shown to work well across all test systems [@Vardeman2011], and remains the preferred method of simulating nanoparticles in the isothermal-isobaric (NPT) ensemble with `OpenMD`.

![A Langevin Hull surrounds explicit water solvating a gold nanoparticle. The Langevin Hull imposes an external pressure and temperature bath and maintains isobaric-isothermal conditions without periodic boundaries. (Image created with the help of Dr. Kristina Davis from the Notre Dame Center for Research Computing.)\label{fig:LHull}](nanoHull_image.png)

## Real Space Electrostatics

Electrostatic interactions are one of the most important intramolecular forces and are present in all but the most basic of molecular simulations. These interactions are also long ranged, and are typically the most computationally expensive. As a result, significant effort has gone into to balancing the accuracy and efficiency of these calculations. `OpenMD` implements a number of techniques which are generally classified according to how solvent molecules are incorporated into the systems of interest. Implicit methods, which exclude solvent molecules from the simulation, offer computational efficiency at the cost of accuracy. One example would be the use of a reaction field [@Onsager1936] for electrostatics coupled with Langevin dynamics to include the hydrodynamic effects of the solvent. Explicit methods which include all solvent molecules directly are the most widely used with `OpenMD`. Explicit electrostatic methods can further be classified as either *Real Space* or *Ewald*-based methods.

The default electrostatics summation method used in `OpenMD` is a real space, damped-shifted force (DSF) model [@Fennell2006] which extends and combines the standard shifted potentials of @Wolf1999 and the damped potentials of @Zahn2002. The potential due to the damped-shifted force has the form:
\begin{equation} \label{eq:coulombic}
  U_\mathrm{Coulomb} = \frac{1}{4 \pi \epsilon_0} \left[ \sum_i \sum_{j>i} U_\mathrm{DSF}(q_i, q_j, r_{ij}) + \sum_i U_i^\mathrm{self}(q_i) \right]
\end{equation}
where the damped shifted force potential,
\begin{equation} \label{eq:DSF}
  U_{\textrm{DSF}}(q_i, q_j, r_{ij}) = \begin{cases}
  q_i q_j \left[ f(r_{ij}) - f(R_c) - f^\prime(R_c) (r_{ij} - R_c) \right], & r_{ij} \le R_c \\
  0, & r_{ij} > R_c
  \end{cases}
\end{equation}
cuts off smoothly as $r_{ij} \rightarrow R_c$, and the Coulombic kernel is damped using a complementary error function, $f(r) = \frac{\textrm{erfc}(\alpha r)}{r}$.  The shifted potential term can be thought of as the interactions of charges with neutralizing charges on the surface of the cutoff sphere. The damping parameter $(\alpha)$ can be specified directly in `OpenMD`, or is set by default from the cutoff radius, $\alpha = 0.425 - 0.02 R_c$, where the code switches to an undamped kernel for $R_c > 21.25$ &#8491;. The self potential represents the interaction of each charge with its own neutralizing charge on the cutoff sphere, and this term resembles the self-interaction in the Ewald sum,
\begin{equation}
  U_i^\mathrm{self}(q_i) =  - \left(\frac{\textrm{erfc}(\alpha R_c)}{R_c} + \frac{\alpha}{\pi^{1/2}}\right) q_i^2.
\end{equation}

DSF offers an attractive compromise between the computational efficiency of *Real Space* methods ($\mathcal{O}(N)$) and the accuracy of the full Ewald sum [@Fennell2006]. The DSF method has also been extended for use with point dipoles and quadrupoles as *Gradient Shifted* and *Taylor Shifted* potentials [@Lamichhane2014] and has been validated for potentials and atomic forces [@Lamichhane2014a], as well as dielectric properties [@Lamichhane2016], against the full Ewald sum. Note that the Ewald method was extended to point multipoles up to quadrupolar order [@Smith82; @Smith98], and this has also been implemented in `OpenMD`. The *Shifted Force* potential generalizes most directly as the *Gradient Shifted* potential for multipoles, and these are the default electrostatic summation methods in `OpenMD`.

## Fluctuating Charges and Densities

One way to include the effects of polarizability in molecular simulations is to use electronegativity equalization [@Rappe1991] or fluctuating charges on atomic sites [@Rick1994]. `OpenMD` makes it relatively easy to add extended variables (e.g. charges) on sites to support these methods. In general, the equations of motion are derived from an extended Lagrangian,
\begin{equation} \label{eq:extendedLagrangian}
  \mathcal{L} = \sum_{i = 1}^N \left[ \frac{1}{2} m_i \dot{\textbf{r}}_i^2 + \frac{1}{2} M_q \dot{q}_i^2 \right] - U(\{\textbf{r}\}, \{q\}) - \lambda \left( \sum_{i = 1}^N q_i - Q \right)
\end{equation}
where the potential depends on both atomic coordinates and the dynamic fluctuating charges on each site. The final term in Eq. \ref{eq:extendedLagrangian} constrains the total charge on the system to a fixed value, and $M_q$ is a fictitious charge mass that governs the speed of propagation of the extended variables.

A relatively new model for simulating bulk metals, the density readjusting embedded atom method (DR-EAM), has also been implemented [@Bhattarai2019]. DR-EAM allows fluctuating densities within the framework of the Embedded Atom Method (EAM) [@Daw1984], by adding an additional degree of freedom, the charge for each atomic site. The total configurational potential energy, $U$, as a function of both instantaneous positions, $\{\textbf{r}\}$, and partial charges $\{q\}$:
\begin{align} \label{eq:fluctuatingChargePotential}
  U_\mathrm{DR-EAM}(\{\textbf{r}\}, \{q\}) = &\sum_i F_i[\bar{\rho}_i] + \frac{1}{2} \sum_i \sum_{j \ne i} \phi_{ij} (r_{ij}, q_i, q_j) \nonumber \\
  &+ \frac{1}{2} \sum_i \sum_{j \ne i} q_i q_j J(r_{ij}) + \sum_i U_i^{\textrm{self}}(q_i)
\end{align}
where the cost of embedding atom $i$ in a total valence density of $\bar{\rho}_i$ is computed using the embedding functional, $F_i[\bar{\rho}_i]$. $\phi_{ij}$ is the pair potential between atoms $i$ and $j$, and $J(r_{ij})$ is the Coulomb integral (which can be computed using the DSF approximation above). Lastly, $U_{\textrm{self}}$ is an additional self potential, modeled as a sixth-order polynomial parameterized by electron affinities and ionization potentials for a wide range of metals,
\begin{equation} \label{eq:selfPotential}
  U_i^{\textrm{self}}(q_i) = \sum_{n = 1}^6 a_n q_i^n.
\end{equation}
The contribution to the local density at site $i$ depends on the instantaneous partial charges on all other atoms,
\begin{equation} \label{eq:rhoi}
  \bar{\rho}_i = \sum_{j \neq i} \left(1 - \frac{q_j}{N_j} \right) f_j(r_{ij})
\end{equation}
with $N_j$ as the valency count for atom $j$. Modifications to the pair potential used in DR-EAM are also supported.

DR-EAM was shown to perform well for bulk metals, metal surfaces, and alloys; and most importantly, it retains the strengths of the unmodified EAM in modeling bulk elastic constants and surface energies [@Bhattarai2019]. DR-EAM has similar performance to the unmodified EAM in that both approaches require a double-pass through the force loop, once to compute local densities and again to compute forces and energies.

We note that the infrastructure required to implement DR-EAM is a superset of what is required for common *fluc-q* potentials like the TIP4P-FQ model for water [@Rick1994].

# Acknowledgements

We would like to acknowledge the contributions of Matthew A. Meineke and Teng Lin to the original `OOPSE` code. Contributions to the `OpenMD` codebase have also come from: Patrick B. Louden, Joseph R. Michalka, James M. Marr, Anderson D.S. Duraes, Suzanne M. Neidhart, Shenyu Kuang, Madan Lamichhane, Xiuquan Sun, Sydney A. Shavalier, Benjamin M. Harless, Veronica Freund, Minh Nhat Pham, Chunlei Li, Kyle Daily, Alexander Mazanek, and Yang Zheng.

Development of `OpenMD` was carried out with support from the National Science Foundation under grants CHE-0848243, CHE-1362211, CHE-1663773, and CHE-191954648.

# References
