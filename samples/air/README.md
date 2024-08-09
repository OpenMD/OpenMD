# Air

## Background Information

The files here can help set up some simple gas-phase simulations of
the common components of air. Three site rigid body models for
N<sub>2</sub>, O<sub>2</sub>, and CO<sub>2</sub> are based on the
TraPPE force field, while the SPC/E water model (also rigid) is
included for simulating various humidity levels. Lennard-Jones
parameters for various noble gases (Ar, He, Ne, Kr) are also included.

Note that if you want heat capacities at high temperatures, you would
need to include vibrational motion for these molecules as well (not
included in these parameters).

Air has a density at sea level and at 15C of roughly 1.225 kg / m<sup>3</sup>
(0.001225 g / cm<sup>3</sup>).  The components of dry air are

| Gas            |  Fraction by Volume |
|----------------|---------------------|
| N<sub>2</sub>  | 0.7809              |
| O<sub>2</sub>  | 0.2095              |
| Ar             | 0.00933             |
| CO<sub>2</sub> | 0.0003              |
| Ne             | 0.000018            |
| He             | 0.000005            |
| Kr             | 0.000001            |

In addition to the monatomic models for Ar, Ne, He, and Kr, this sample also requires a few models for diatomic and triatomic molecules.  In this sample, these molecules are represented as rigid linear (N<sub>2</sub>, O<sub>2</sub>, CO<sub>2</sub>) or three site models (SPC/E water).

## Instructions

The `air.inc` and `Air.frc` files contain parameters for most of the
simple components of air.  To build a small box, one might start with
the `mix.omd` which declares the three most prevalent components
(N<sub>2</sub>, O<sub>2</sub>, and Ar).

With all 3 of these files in this sample directory, you can make an air
mixture using randomBuilder:
```
randomBuilder mix.omd -o test.omd --density=0.001225 --nx=5 --ny=5 --nz=5 --molFraction=0.78084 --molFraction=0.20946
```
This creates a `test.omd` structure containing 500 molecules in
roughly the correct proportions (and density).

To warm this mixture up to 25C (and assign initial velocities):
```
thermalizer -t 298 -i test.omd -o warm.omd
```
The simulation is relatively short:
```
openmd warm.omd
```

## Expected Output

Some missing interaction warnings are expected here (notably between Ar atoms and the non-atomic 'X' sites on the N<sub>2</sub> and O<sub>2</sub> models).  OpenMD is just warning the user that there are atom pairs present in the simulation with no corresponding interaction in the `frc` file.

This is a relatively short simulation, but the expected report would be:

```
###############################################################################
# Status Report:                                                              #
#              Total Time:       10000 fs                                     #
#       Number of Samples:         101                                        #
#            Total Energy:     735.802  ±  0.586997     kcal/mol              #
#        Potential Energy:   -0.279875  ±  0.0654789    kcal/mol              #
#          Kinetic Energy:     736.082  ±  0.617257     kcal/mol              #
#             Temperature:     298.116  ±  0.249991     K                     #
#                Pressure:     1.03131  ±  0.00126752   atm                   #
#                  Volume: 1.96509e+07  ±  0            A^3                   #
#      Conserved Quantity:     740.366  ±  0.0023359    kcal/mol              #
###############################################################################
```
This approximates a box of air at standard temperature and Pressure (298 K and 1 atm).
Then `Dump2XYZ` can output the base atom types and map back to the simulation box for visualization with `Jmol`:
```
Dump2XYZ -i warm.dump -b -m
jmol warm.xyz
```

## References

| Molecular Model| Number of Sites | DOI  |
| ----------:|:---------------:|-----:|
| N<sub>2</sub>  | 3 |[10.1002/aic.690470719](https://doi.org/10.1002/aic.690470719) |
| O<sub>2</sub>  | 3 |[10.1007/s00214-005-0073-1](https://doi.org/10.1007/s00214-005-0073-1) |
| CO<sub>2</sub> | 3 |[10.1002/aic.690470719](https://doi.org/10.1002/aic.690470719)     |
| SPC/E | 3 |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| CH<sub>4</sub> | 1 |[10.1021/jp972543+](https://doi.org/10.1021/jp972543+) |
| He, Ne, Ar, Kr | 1 | Maitland, G.C., Rigby, M., Smith, E.B., and Wakeham, W.A. (1981) *Intermolecular forces: their origin and determination*. Clarendon Press, Oxford. |

