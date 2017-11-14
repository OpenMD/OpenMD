The following proton-disordered ice-Ih crystals,

3x3x2-C2.omd	3x3x2-e.omd	3x3x2.omd	5x3x3.omd	6x4x4.omd
3x3x2-CH.omd	3x3x2-h.omd	4x3x2.omd	6x3x3.omd	9x5x1.omd

were taken from the supporting information of "Unit cells for hexagonal ice" by J. A. Hayward and J. R. Reimers, JCP 106 (1997) 1518.

NOTE:
These structures are ideal ice crystals, and should be gently equilibrated with whichever water model you choose. Depending on the model, these starting structures may be more or less favorable. A sample equilibration scheme might be as follows.

1. Short NPTxyz at a low temperature, approximately 10 to 50 K, with resetTime set to a small time, approximately 10 to 50 fs.
2. Once the pressure tensor elements are about zero and the volume is fluctuating around some average, affineScale the simulation cell to the average volume. Turn resetTime off.
3. Perform an NVT simulation with the targetTemperature at whatever your desired temperature is.
4. When the temperature is fluctuating about the average, and the total energy is fluctuation about some average, thermalScale the simulation cell to this average energy.
5. You can now perform NVE simulations.
