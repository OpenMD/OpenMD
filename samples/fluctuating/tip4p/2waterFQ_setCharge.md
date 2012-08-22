<OpenMD version=1>
  <MetaData>

molecule {
  name = "TIP4P_FQ";

  atom[0]{
    type = "O_TIP4P_FQ";
    position( 0.0, 0.0, -0.6556 );
  }
  atom[1]{
    type = "H_TIP4P_FQ";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP4P_FQ";
    position( 0.0, -0.75695, 0.52032 );
  }
  atom[3]{
    type = "M_TIP4P_FQ";
    position( 0.0, 0.0, 0.08444 );
  }

  rigidBody[0]{
    members(0,1,2,3);
  }
}
component{
  type = "TIP4P_FQ";
  nMol = 2;
}


ensemble = NVE;
forceField = "FlucQ";
forceFieldFileName = "FQ.frc";
cutoffMethod = "shifted_force";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.7;
outputFluctuatingCharges = true;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;

dt = 1.0;
runTime = 1e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 10;
sampleTime = 100;
statusTime = 10;
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 10000
        Hmat: {{ 19.713, 0, 0 }, { 0, 19.713, 0 }, { 0, 0, 19.713 }}
  Thermostat: 0 , 0
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0  pvqjcw        23.95523961       -2.701653649       -15.70063919  2.409292e-03  9.908152e-06 -1.574237e-03  2.074455e-01  5.294415e-01  5.318552e-01  6.275253e-01 -6.757743e-04 -5.046488e-03  1.001059e-02  -0.500000e+00   0.000000e+00 
         1  pvqjcw       -18.95523961        2.701653649        15.70063919 -2.409292e-03 -9.908152e-06  1.574237e-03  1.748160e-03  5.106060e-01 -8.393383e-01 -1.865201e-01  5.628458e-03 -8.785288e-03 -1.926888e-02  0.500000e+00   0.000000e+00 
    </StuntDoubles>
  </Snapshot>
</OpenMD>
