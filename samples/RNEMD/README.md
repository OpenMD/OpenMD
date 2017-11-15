# Sample (.omd) files using Reverse Non-Equilibrium Molecular Dynamics (RNEMD)

## Momentum transport in bulk fluids

The file **shearWater.omd** is a box of 1500 SPC/E water molecules
which has the momentum flux functionality of RNEMD turned on. Notice
in the RNEMD block of the (.omd) file that

```
	fluxType = "Px";
	momentumFlux = 6.0e-7;
```

With these parameters, the applied flux will be an x-axis momentum
flux transferred across the z-dimension of the box. Using this
functionality, one is able to measure the *shear viscosity* of the
liquid at the simulated temperature by relating the imposed momentum
flux to the system's gradient response of the velocity.

## Thermal transport in bulk materials

It is also possible to measure the *thermal conductivity* of a
material using the RNEMD functionality in OpenMD. As an example,
**graphene.omd** is an (.omd) file where two sheets of graphene have a
thermal flux applied accross the long axis of the sheet. In the (.omd)
file the fluxType has been set to a kinetic energy flux, and also that
the kineticFlux is defined.

```
	fluxType = "KE";
	kineticFlux = -6.55e-11;
```
	
The system responds to the thermal flux by developing a temperature
gradient across the z-axis of the system. The thermal conductivity can
be computed by relating the resulting thermal gradient to the imposed
kinetic energy flux.

## Thermal transport across an interface

While computing the thermal conductivity of a bulk material is
certainly of interest, one may also want to investigate interfacial
thermal conductance across an interface. The RNEMD functionality of
OpenMD easily allows for this, and the RNEMD definition block in the
(.omd) file only needs a minor tweaks.  Before our selection, was of
only one component, now we need to allow for more than one atom or
molecule type to be selected.

**gold_water_interface.omd** is an example system where we can compute
the thermal conductivity across an interface, here, a gold / water
interface. Notice in the RNEMD declaration block in the (.omd) file
that the objectSelection is now,

```
	objectSelection = "select SPCE_RB_0 or Au";
```

It is important to make sure your simulation cell is constructed
properly for these kinds of simulations. Since the two RNEMD exchange
regions are defined along the z-dimension at the middle of the
simulation cell and at the far edges (wrapping about the periodic
box), the gold and water need to be properly distributed throughout
the box or else your computation will not give you what you want.


## Simultaneous shearing and thermal transport of bulk materials

The file **2744_shear.omd** is a box of 2744 Argon atoms which has a
simultaneous momentum and kinetic energy flux through the box. Notice
in the RNEMD block of the (.omd) file that

```
	fluxType = "KE+Pvector";
```

and both *kineticFlux* and *momentumFluxVector* are defined.

```
	kineticFlux = -5.0e-6;
	momentumFluxVector = (-2e-7, 0, 0);
```

Application of simultaneous momentum and kinetic energy fluxes result
in both a velocity and thermal gradient response of the system,
allowing for measurement of the shear viscosity of the fluid at a
large number of temperature domains with one simulation.


## Thermal transport in non-periodic systems

OpenMD can perform non-periodic simulations using the Langevin Hull
along with the RNEMD functionality, allowing for computation of
thermal conductance across solvated nanoparticle interfaces. OpenMD
hosts a large number of builder utility scripts which aide in the
construction of these nanoparticles (nanospheres, icosohedra,
cuboctahedra), which can be found in *samples/builders/*.

**NP20_hex_KEflux.omd** is an example of a Gold nanosphere solvated in
  hexane, with a kinetic energy flux that moves thermal energy from
  the solvent into the particle. The RNEMD declaration block in the
  (.omd) file is only slightly different than above. The syntax is
  described in detail in the OpenMD manual.

```
	useRNEMD = "true";
	objectSelection = "select Au or Hexane";
	sphereAradius = 10;
	sphereBradius = 41;
	method = "VSS";
	fluxType = "KE";
	kineticFlux = 1E-5;
	exchangeTime = 10;
	outputBins = 60;

```
The only notable change to the RNEMD declaration block is the addition
of *sphereAradius* and *sphereBradius*, which define the two exchange
regions for the RNEMD moves.
