# Water models in OpenMD

OpenMD supports a large number of water models, only some of which
have examples shown here.

### Single site models
+ SSD, SSD-E, SSD-RF
+ SSDQ, SSDQO

### Three site models
+ SPC, SPC-HW
+ SPC/E
+ TIP3P

### Four site models
+ TIP4P, TIP4P-Ew, TIP4P/2005, TIP4P/Ice

### Five site models
+ TIP5P, TIP5P-E

### Six site models
+ NE6



One of the powerful functionalities of OpenMD is the separation
between declaration of *components* and *stuntDoubles*. In every
(.omd) file, a component block is used to specifiy what types of atoms
or molecules the stuntDoubles correspond to.

```
	component{
		type = "SPCE";
		nMol = 256;
	}

```

```
    <StuntDoubles>
         0    pvqj         -13.882269          -6.595441         -10.074898  8.730000e-04  7.358000e-03  2.248000e-03  7.861640e-01 -3.843370e-01 -3.681500e-01  3.141600e-01 -1.691600e-02  1.558000e-03  1.901000e-02
         1    pvqj         -10.800233          -3.862809          -9.917968 -1.088300e-02  6.540000e-04  4.554000e-03  4.137720e-01 -7.080340e-01 -4.812670e-01 -3.096160e-01 -7.502000e-03  1.201400e-02  9.878000e-03
         2    pvqj         -13.359152         -13.264782          -5.753301  3.626000e-03  3.221000e-03  1.990000e-04  5.155220e-01 -1.470710e-01 -5.520980e-01  6.385880e-01 -3.590000e-03 -1.001200e-02 -2.305800e-02
 ```

Due to this separation, a simulation of water with one model can be
easily transformed into a simulation with a different water model
(regardless of the number of potential sites) by simply changing the
component block.

```
	component{
		type = "TIP4P";
		nMol = 256;
	}

```

Now, the same stuntDoubles will be treated as *TIP4P* water instead of
*SPC/E* water. Of course, due to the change in potential you may have
to re-equilibrate the system. However, no further changes need to be
made.

In addition to the ease of transferring between potentials, OpenMD
also has a water system builder which generates OpenMD (.omd) files
for you, **waterBoxer**.

## waterBoxer

**waterBoxer** is a perl script which generates FCC lattices of water
molecules at user-defined densities and system sizes. The user can
specifiy which water model to generate a system of, or as described
above, change the component block definition from the default once
generated.  **waterBoxer** prints a helpful discription of how its use
and functionalities when passed *-h*.


# References

| Water Model| Number of Sites | DOI  |
| ----------:|:---------------:|-----:|
| SSD        | 1 |[10.1016/S0009-2614(03)01044-3](https://doi.org/10.1016/S0009-2614(03)01044-3) |
| SSD/RF     | 1 |[10.1063/1.1697381](https://doi.org/10.1063/1.1697381)     |
| SSD/E      | 1 |[10.1063/1.1697381](https://doi.org/10.1063/1.1697381)     |
| SPC        | 3 |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| SPC-HW     | 3 |[10.1063/1.1359183](https://doi.org/10.1063/1.1359183)     |
| SPC/E      | 3 |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| TIP3P      | 3 |[10.1063/1.445869](https://doi.org/10.1063/1.445869)       |
| TIP4P      | 4 |[10.1063/1.445869](https://doi.org/10.1063/1.445869)       |
| TIP4P-Ew   | 4 |[10.1063/1.1683075](https://doi.org/10.1063/1.1683075)     |
| TIP4P/2005 | 4 |[10.1063/1.2121687](https://doi.org/10.1063/1.2121687)     |
| TIP4P/Ice  | 4 |[10.1063/1.1931662](https://doi.org/10.1063/1.1931662)     |
| TIP5P      | 5 |[10.1063/1.481505](https://doi.org/10.1063/1.481505)       |
| TIP5P-E    | 5 |[10.1063/1.1652434](https://doi.org/10.1063/1.1652434)     |
| NE6        | 6 |[10.1063/1.1562610](https://doi.org/10.1063/1.1562610)     |
