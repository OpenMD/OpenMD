<OpenMD version=2>
  <MetaData>
#include "metals.md"


component{
  type = "Cu";
	nMol = 1;
}

component{
  type = "Au";
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
    </FrameData>
    <StuntDoubles>
	0	pv	 0.00000000   0.00000000   0.00000000   	0.0	0.0	0.0
	1	pv	 1.00000000   0.00000000   0.00000000   	0.0	0.0	0.0
	2	pv	 0.50000000   0.86602540   0.00000000   	0.0	0.0	0.0
	3	pv	-0.50000000   0.86602540   0.00000000   	0.0	0.0	0.0
	4	pv	-1.00000000   0.00000000   0.00000000   	0.0	0.0	0.0
	5	pv	-0.50000000  -0.86602540   0.00000000   	0.0	0.0	0.0
	6	pv	 0.50000000  -0.86602540   0.00000000   	0.0	0.0	0.0
	7	pv	 0.50000000   0.28867513   0.81649658   	0.0	0.0	0.0
	8	pv	-0.50000000   0.28867513   0.81649658   	0.0	0.0	0.0
	9	pv	 0.00000000  -0.57735027   0.81649658   	0.0	0.0	0.0
	10	pv	 0.00000000   0.57735027  -0.81649658   	0.0	0.0	0.0
	11	pv	 0.50000000  -0.28867513  -0.81649658   	0.0	0.0	0.0
	12	pv	-0.50000000  -0.28867513  -0.81649658   	0.0	0.0	0.0
    </StuntDoubles>
  </Snapshot>
