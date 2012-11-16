<OpenMD version=2>
  <MetaData>
#include "water.md"

component{
  type = "Cl-";
  nMol = 1;
}
component{
  type = "SSD_E";
  nMol = 1;
}

ensemble = NVT;
forceField = "DUFF";
cutoffMethod = "shifted_force";
dampingAlpha = 0.2;
cutoffRadius = 9.0;

targetTemp = 0.001;
targetPressure = 1.0;

tauThermostat = 1e4;
tauBarostat = 1e4;

dt = 1.0;
runTime = 1e3;

sampleTime = 10;
statusTime = 1;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 100;
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 1000
        Hmat: {{ 59.7166, 0, 0 }, { 0, 59.7166, 0 }, { 0, 0, 59.7166 }}
  Thermostat: 0.0385672 , 36.2096
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0      pv           4.596242           0.717716          -0.489407   0.000000e+00  0.000000e+00 -0.000000e+00
         1    pvqj            1.81745          -1.440934           0.945913 -0.000000e+00 -0.000000e+00  0.000000e+00  2.831510e-01 -6.511590e-01 -5.879580e-01 -3.874570e-01 -0.000000e+00 -0.000000e+00  0.000000e+00
    </StuntDoubles>
  </Snapshot>
</OpenMD>
