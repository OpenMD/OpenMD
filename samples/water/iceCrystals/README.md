# This directory contains unit cell and larger structures for simulations of ice crystals.

## Proton-disordered ice crystals
The following proton-disordered ice-Ih crystals,

+ 3x3x2-C2.omd
+ 3x3x2-e.omd
+ 3x3x2.omd
+ 5x3x3.omd
+ 6x4x4.omd
+ 3x3x2-CH.omd
+ 3x3x2-h.omd
+ 4x3x2.omd
+ 6x3x3.omd
+ 9x5x1.omd

were taken from the supporting information of "Unit cells for hexagonal ice" by J. A. Hayward and J. R. Reimers, JCP 106 (1997) 1518.

## Proton-ordered ice unit cells
The following structures are unit cells for proton-ordered ice-Ih crystals.

+ HO-struct1.omd
+ HO-struct2.omd
+ HO-struct3.omd
+ HO-struct4.omd
+ HO-struct5.omd
+ HO-struct6.omd
+ HO-struct7.omd
+ HO-struct8.omd
+ HO-struct9.omd
+ HO-struct10.omd
+ HO-struct11.omd
+ HO-struct12.omd
+ HO-struct3.omd
+ HO-struct14.omd
+ HO-struct15.omd
+ HO-struct16.omd

These structures come from Table 1. of "Quantum-Chemical and Force-Field Investigations of Ice Ih" by Thomas K. Hirsch and Lars Ojamae, J. Phys. Chem. B (2004) 108, 15856-15864.

NOTE: HO-struct1.omd	is actually ice XI.

When replicated, HO-struct6.omd and HO-struct7.omd create proton stripes on the basal surfaces. 


## Creating large ice crystals from these structures
In order to generate larger crystals from these structures, use omd2omd with the -x -y -z flags to replicate these unit cells in the x, y, and z dimensions. 

```
omd2omd -i HO-struct1.omd -o bigCrystal.omd -x 5 -y 3 -z 5
```

Also, while you are unable to cleave the crystals with the current OpenMD software, you are able to rotate these crystals exposing the prismatic and secondary prismatic facets using the -p -q -r functionality of omd2omd.

```
omd2omd -i bigCrystal.omd -p 90 -q 90 -r 0
```

## Sample equilibration scheme
NOTE:
These structures are ideal ice crystals, and should be gently equilibrated with whichever water model you choose. Depending on the model, these starting structures may be more or less favorable. A sample equilibration scheme might be as follows.

1. Short NPTxyz at a low temperature, approximately 10 to 50 K, with resetTime set to a small time, approximately 10 to 50 fs.
2. Once the pressure tensor elements are about zero and the volume is fluctuating around some average, affineScale the simulation cell to the average volume. Turn resetTime off.
3. Perform an NVT simulation with the targetTemperature at whatever your desired temperature is.
4. When the temperature is fluctuating about the average, and the total energy is fluctuation about some average, thermalScale the simulation cell to this average energy.
5. You can now perform NVE simulations.
