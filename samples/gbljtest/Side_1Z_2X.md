<OpenMD version=1>
  <MetaData>
 
#include "TestAtom.md"
 
component{
  type = "TEST";
  nMol = 2;
}

// ensemble = NgT;
// ensemble = NPrT;
// ensemble = NPTxyz;
ensemble = NgT;
forceField = "DUFF";
forceFieldFileName = "Field.frc";
cutoffRadius = 35.0;
switchingRadius =32.0;
 
targetTemp = 250.0;
targetPressure = 1.0;
surfaceTension = 0.0;
 
tauThermostat = 1e3;
tauBarostat = 1e5;
 
dt = 20;
runTime = 20;
// tempSet = "true";
sampleTime = 20;
statusTime = 20;
// thermalTime = 1;
//resetTime = 5e2;
// useInitialTime = "false";
// useInitialExtendedSystemState = false;
 
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 30000000
        Hmat: {{ 167.4117949, 0, 0 }, { 0, 192.2101739, 0 }, { 0, 0, 168.7655036 }}
  Thermostat: 6.561577322e-06 , -2.415346931
    Barostat: {{ -1.151721885e-07, 0, 0 }, { 0, 5.240323006e-08, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0    pvqj        0.0       0.0        0.0  7.093100e-05 -1.381681e-03 -1.519481e-03  1.0  0.0  0.0  0.0 -3.765982e-01 -2.299070e-01 -1.026296e-03
         1    pvqj        6.0       0.0        0.0 -4.213290e-05  1.255466e-03  8.277276e-04  0.5  0.0  0.5  0.0 -8.275834e-02  1.091971e+00  2.505469e-05
    </StuntDoubles>
  </Snapshot>
</OpenMD>
