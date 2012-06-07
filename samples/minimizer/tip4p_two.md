<OpenMD version=1>
  <MetaData>
#include "water.md"


component{
  type = "TIP4P";
  nMol = 2;
}

minimizer {
  useMinimizer = true;
  method = SD;
  maxIterations = 5000;
}

forceField = "DUFF2";
cutoffMethod = "shifted_force";
dampingAlpha = 0.185;
cutoffRadius = 12.0;

  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 10000
        Hmat: {{ 60, 0, 0 }, { 0, 60, 0 }, { 0, 0, 60 }}
  Thermostat: 0.136949 , 1115.21
    </FrameData>
    <StuntDoubles>
         0    pvqj            -2.2274             -2.572            -2.8204  0.000000e+00  0.000000e+00 -0.000000e+00 -1.195100e-01  8.567900e-01  2.295800e-01  4.460100e-01  2.200000e-05 -2.700000e-05  7.000000e-06
         1    pvqj            -2.2274             -2.572                  0 -0.000000e+00 -0.000000e+00  0.000000e+00  3.998900e-01  7.681700e-01 -4.435100e-01 -2.308700e-01 -4.500000e-05 -7.000000e-06 -7.000000e-06
    </StuntDoubles>
  </Snapshot>
</OpenMD>
