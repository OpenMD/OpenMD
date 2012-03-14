<OpenMD version=1>
  <MetaData>
#include "gb.md"
#include "lj.md"


component{
  type = "C";
  nMol =1;
}
component{
  type = "GBlinear";
  nMol = 1;
}



ensemble = NVE;
targetTemp = 1;
tauThermostat = 1000;
forceField = "DUFF";

cutoffMethod = "switched";
cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1;
runTime = 10000;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

thermalTime = 10;
sampleTime = 100;
statusTime = 1;
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 0
        Hmat: {{ 69.7166, 0, 0 }, { 0, 69.7166, 0 }, { 0, 0, 69.7166 }}
    </FrameData>
    <StuntDoubles>
         0      pv                  0                  0                  0   0.000000e+00  0.000000e+00  0.000000e+00
         1    pvqj                  8                  0                  3  0.000000e+00  0.000000e+00  0.000000e+00  5.000000e-01  5.000000e-01  5.000000e-01  5.000000e-01  0.000000e+00  0.000000e+00  0.000000e+00
    </StuntDoubles>
  </Snapshot>
</OpenMD>
