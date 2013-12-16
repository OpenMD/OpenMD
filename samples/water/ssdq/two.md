<OpenMD version=1>
  <MetaData>
#include "water.md"

component{
  type = "SSDQ";
  nMol = 2;
}

ensemble = NVT;
forceField = "SSDQ";
cutoffMethod = "shifted_force";
cutoffRadius = 12.0;
dampingAlpha = 0.18;

targetTemp = 0.1;
targetPressure = 1.0;

tauThermostat = 1e4;
tauBarostat = 1e4;

dt = 1.0;
runTime = 1e4;

sampleTime = 1e2;
statusTime = 10;
useInitialTime = "false";
useInitialExtendedSystemState = "false";
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 100000
        Hmat: {{ 59.7166, 0, 0 }, { 0, 59.7166, 0 }, { 0, 0, 59.7166 }}
  Thermostat: 0.0413948 , 4187.83
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0    pvqj           1.192654           0.024173           0.011944 -0.000000e+00  0.000000e+00 -0.000000e+00  6.514170e-01 -6.266110e-01 -3.051170e-01  2.998650e-01  0.000000e+00  0.000000e+00 -0.000000e+00
         1    pvqj           3.837346          -0.024173          -0.011944  0.000000e+00 -0.000000e+00  0.000000e+00 -6.613910e-01  2.582810e-01 -6.644630e-01  2.331120e-01  0.000000e+00  0.000000e+00  0.000000e+00
    </StuntDoubles>
  </Snapshot>
</OpenMD>
