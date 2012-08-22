<OpenMD version=1>
  <MetaData>
#include "gb.md"
#include "lj.md"


component{
  type = "GBlinear";
  nMol = 2;
}



ensemble = NVT;
targetTemp = 1;
tauThermostat = 30;
forceField = "DUFF";
forceFieldFileName = "DUFF2.frc";

cutoffRadius = 20.0;
switchingRadius = 18.0;
dt = 1.0;
runTime = 10;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

thermalTime = 10;
sampleTime = 100;
statusTime = 1;
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 0
        Hmat: {{ 63.7166, 0, 0 }, { 0, 63.7166, 0 }, { 0, 0, 63.7166 }}
    </FrameData>
    <StuntDoubles>
         0    pvqj                  5                  0                  1  0.000000e+00  0.000000e+00  0.000000e+00  7.071070e-01  7.071070e-01  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
         1    pvqj                 -5                  0                 -3  0.000000e+00  0.000000e+00  0.000000e+00  5.000000e-01  5.000000e-01  5.000000e-01  5.000000e-01  0.000000e+00  0.000000e+00  0.000000e+00
    </StuntDoubles>
  </Snapshot>
</OpenMD>
