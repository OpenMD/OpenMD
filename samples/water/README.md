# Water models in OpenMD

OpenMD supports a large number of water models, only some of which
have examples shown here. Currently OpenMD does not support
fluctuating charge water models, however, this feature is in
development. 

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
