<OpenMD version=2>
  <MetaData>
#include "metals.md"


component{
  type = "Cu";
	nMol = 1;
}

component{
  type = "Au";
	nMol = 4;
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
    </FrameData>
    <StuntDoubles>
         0      pv                  0                  0                  0   0.000000e+00  0.000000e+00  0.000000e+00
         1      pv                  0                  0                  4.08   0.000000e+00  0.000000e+00  0.000000e+00
         2      pv                  0                  3.84666           -1.36   0.000000e+00  0.000000e+00  0.000000e+00
         3      pv                  -3.33131           -1.92333          -1.36   0.000000e+00  0.000000e+00  0.000000e+00
         4      pv                  3.33131            -1.92333          -1.36   0.000000e+00  0.000000e+00  0.000000e+00
    </StuntDoubles>
  </Snapshot>
