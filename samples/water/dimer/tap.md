<OpenMD version=1>
  <MetaData>
#include "water.md"


component{
  type = "TAP";
  nMol = 2;
}



ensemble = NVT;
forceField = "DUFF";
electrostaticSummationMethod = "none";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.8;



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
        Time: 10000
        Hmat: {{ 59.7166, 0, 0 }, { 0, 59.7166, 0 }, { 0, 0, 59.7166 }}
  Thermostat: 0.0208838 , 207.625
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0    pvqj           1.118261           0.030441           0.007695  0.000000e+00 -0.000000e+00 -0.000000e+00  6.682090e-01 -6.132790e-01 -3.089100e-01  2.862860e-01 -0.000000e+00 -0.000000e+00  0.000000e+00
         1    pvqj           3.881739          -0.030441          -0.007695 -0.000000e+00  0.000000e+00 -0.000000e+00  6.207330e-01  3.413810e-01 -6.447840e-01  2.870600e-01 -0.000000e+00 -0.000000e+00 -0.000000e+00
    </StuntDoubles>
  </Snapshot>
</OpenMD>
