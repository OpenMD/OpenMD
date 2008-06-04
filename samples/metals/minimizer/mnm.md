<OOPSE version=4>
  <MetaData>
#include "water.md"
#include "metals.md"

component{
  type = "SSD_E";
  nMol = 1;
}
component{
  type = "Au";
  nMol = 1;
}

//minimizer = CG;
//minimizerMaxIter = 5000;
//minimizerWriteFrq = 1;
//minimizerStepSize = 0.1;

ensemble = NVT;
targetTemp = 10.0;
tauThermostat = 1000.0;
//tempSet = "true";
//thermalTime = 100.0;
forceField = "MnM";
//electrostaticSummationMethod = "shifted_Force";
//dielectric = 80.0;
cutoffRadius = 25.0;
switchingRadius = 20.0;
dt = 1;
runTime = 1e4;
sampleTime = 10;
statusTime = 10;

  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 10000
        Hmat: {{ 60, 0, 0 }, { 0, 60, 0 }, { 0, 0, 60 }}
  Thermostat: 0.136949 , 1115.21
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0    pvqj       -2.821584235        -3.45911795       -2.048658245 -1.493914e-03 -9.573049e-04  9.328191e-04  7.128266e-02 -4.725419e-01  8.552692e-01  2.003438e-01  1.633433e-03 -7.549648e-03  5.201663e-04
         1      pv       -5.173054631       -4.490862284        2.929414839  1.366366e-04  8.755717e-05 -8.531764e-05
    </StuntDoubles>
  </Snapshot>
</OOPSE>
