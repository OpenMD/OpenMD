# Water models in OpenMD

OpenMD supports a large number of water models, only some of which
have examples shown here. Currently OpenMD does not support
fluctuating charge water models, however, this feature is in
development.

| Water Model| Number of Sites | DOI  |
| ---------:|:-------------:|----:|
| SSD       | 1             |[10.1016/S0009-2614(03)01044-3](https://doi.org/10.1016/S0009-2614(03)01044-3) |
| SSD/RF    | 1             |[10.1063/1.1697381](https://doi.org/10.1063/1.1697381) |
| SSD/E     | 1             |[10.1063/1.1697381](https://doi.org/10.1063/1.1697381) |
| SSDQ      | 1             |
| SSDQ0     | 1             |
| SPC       | 3
| |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| SPC-HW    | 3 | |[10.1063/1.1359183](https://doi.org/10.1063/1.1359183) |
| SPC/E     | 3 | |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| TIP3P     | 3 | |[10.1063/1.445869](https://doi.org/10.1063/1.445869) |
| TIP4P     | 4 | |[10.1063/1.445869](https://doi.org/10.1063/1.445869) |
| TIP4P-Ew  | 4 | |[10.1063/1.1683075](https://doi.org/10.1063/1.1683075) |
| TIP4P/2005 | 4 | |[10.1063/1.2121687](https://doi.org/10.1063/1.2121687) |
| TIP4P/Ice | 4 | |[10.1063/1.1931662](https://doi.org/10.1063/1.1931662) |
| TIP5P     | 5 | |[10.1063/1.481505](https://doi.org/10.1063/1.481505) |
| TIP5P-E   | 5             |[10.1063/1.1652434](https://doi.org/10.1063/1.1652434) |
| NE6       | 6 | |[10.1063/1.1562610](https://doi.org/10.1063/1.1562610) |


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

## waterBoxer

OpenMD also has a utility script to aide in the generation of water
simulation cells. **waterBoxer** is a perl script which generates FCC
lattices of water molecules at user-defined densities and system
sizes. **waterBoxer** prints a helpful discription of how its use and
functionalities when passed *-h*.
