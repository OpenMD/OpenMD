These are sample files for carrying out isobaric-isothermal (NPT)
simulations for non-periodic systems using the Langevin Hull.

Input File              Composition             Bath Properties
---------------         ------------------      ---------------------------
SPC/E Water clusters:
spce_1atm.md            1372 SPC/E cluster      1 atm, 300 K, 0.0089 Poise 
spce_100atm.md          1372 SPC/E cluster      100 atm, 300 K, 0.0089 Poise

Gold nanospheres:
Au_300K.md              20 A radius (1985 Au)   4 GPa, 300K, 0.0089 Poise
Au_lowvisc.md           40 A radius (15707 Au)  4 GPa, 400K, 0.0089 Poise 
Au_highvisc.md          40 A radius (15707 Au)  4 GPa, 400K, 0.07288 Poise

18 Angstrom gold nanospheres in SPC/E Water:
spce_Au_1atm.md         1433 Au + 5000 SPC/E    1 atm, 300K, 0.0089 Poise  
spce_Au_100atm.md       1433 Au + 5000 SPC/E    100 atm, 300K, 0.0089 Poise

In general, to use the LangevinHull integrator, you'll want include
these lines in your .md file:
	     		   	
ensemble = "LangevinHull";
targetTemp = 300;
targetPressure = 100;
viscosity = 0.0089;
usePeriodicBoundaryConditions = "false";

The bath is characterized by a pressure, temperature, and viscosity,
so these keywords are required, but the values depend on what you are
trying to simulate.

Note that the time required for thermal equilibration depends on
exposed surface area and bath viscosity.

