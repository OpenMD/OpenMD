<OOPSE version=4>
  <MetaData>
#include "metals.md"


component{
  type = "Cu";
	nMol = 1;
}

component{
  type = "Ag";
	nMol = 12;
}



forceField = "SC";
forceFieldFileName = "SuttonChen.QSC.frc";
targetTemp = 300.0;


ensemble = "LANGEVINDYNAMICS";
langevinBufferRadius = 18.0;
frozenBufferRadius = 30.0;
viscosity=0.00890;

dt = 4.0;
runTime = 1.5e6;

usePeriodicBoundaryConditions = "false";

sampleTime = 10000.0;
statusTime = 1000.0;
thermalTime = 20.0;
tempSet = "false";
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 1000000
        Hmat: {{ 30.85, 0, 0 }, { 0, 30.85, 0 }, { 0, 0, 30.85 }}
  Thermostat: 0 , 0
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0      pv                  0                  0                  0  -7.100000e-04 -1.563000e-03 -1.861000e-03
         1      pv                  0               4.09        6.617759014  -1.002000e-03  2.540000e-04 -2.680000e-04
         2      pv                  0              -4.09        6.617759014  -1.188000e-03 -4.970000e-04  7.600000e-05
         3      pv                  0               4.09       -6.617759014  -2.096000e-03  5.370000e-04  2.171000e-03
         4      pv                  0              -4.09       -6.617759014  -1.400000e-04  2.340000e-03 -8.580000e-04
         5      pv               4.09        6.617759014                  0   0.000000e+00 -3.282000e-03  4.970000e-04
         6      pv              -4.09        6.617759014                  0   0.000000e+00  1.126000e-03 -6.470000e-04
         7      pv               4.09       -6.617759014                  0   0.000000e+00 -7.150000e-04  2.150000e-04
         8      pv              -4.09       -6.617759014                  0   1.470000e-03  4.730000e-04  2.846000e-03
         9      pv        6.617759014                  0               4.09   3.213000e-03  5.880000e-04  6.100000e-04
        10      pv       -6.617759014                  0               4.09   1.312000e-03 -1.007000e-03  6.900000e-05
        11      pv        6.617759014                  0              -4.09  -6.430000e-04 -1.496000e-03  4.024000e-03
        12      pv       -6.617759014                  0              -4.09   2.067000e-03 -9.920000e-04 -1.012000e-03
    </StuntDoubles>
  </Snapshot>
</OOPSE>
