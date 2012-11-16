<OpenMD version=2>
  <MetaData>
#include "water.md"

component{
  type = "TIP4P";
  nMol = 2;
}

ensemble = NVT;
forceField = "DUFF";
cutoffMethod = "shifted_force";
cutoffRadius = 9.0;
dampingAlpha = 0.2;

targetTemp = 0.001;
targetPressure = 1.0;

tauThermostat = 1e4;
tauBarostat = 1e4;

dt = 1.0;
runTime = 1e4;

sampleTime = 1e2;
statusTime = 1;
useInitialTime = "false";
useInitialExtendedSystemState = "false";
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 99991
        Hmat: {{ 59.7166, 0, 0 }, { 0, 59.7166, 0 }, { 0, 0, 59.7166 }}
  Thermostat: -0.00400019 , 11.42
    </FrameData>
    <StuntDoubles>
         0    pvqj           1.139005            0.02947           0.017771  8.100000e-05 -2.510000e-04  2.200000e-05  6.728890e-01 -5.173170e-01 -5.134360e-01  1.264420e-01  1.472000e-03  3.632000e-03 -9.030000e-04
         1    pvqj           3.860995           -0.02947          -0.017771 -8.100000e-05  2.510000e-04 -2.200000e-05  7.059430e-01  1.780690e-01 -5.268640e-01  4.385770e-01 -1.064400e-02  2.690000e-03  1.290000e-04
    </StuntDoubles>
  </Snapshot>
</OpenMD>
