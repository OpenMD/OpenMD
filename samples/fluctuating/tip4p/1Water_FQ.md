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
  nMol = 1;
}

flucQ {
 targetTemp = 1.0e-6;
 tauThermostat = 10.0;
}

ensemble = NVE;
forceField = "FlucQ";
forceFieldFileName = "FQ.frc";
cutoffMethod = "shifted_force";
dielectric = 80.0;
cutoffRadius = 9.0;
switchingRadius = 7.7;
outputFluctuatingCharges = true;
runTime = 1e-1;
sampleTime = 1;
statusTime = 1;
resetTime = 1;

targetTemp = 300;
targetPressure = 1.0;

tauThermostat = 1e3;
tauBarostat = 1e4;


dt = 1.0;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

//tempSet = "true";
//thermalTime = 10;
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 10000
        Hmat: {{ 19.713, 0, 0 }, { 0, 19.713, 0 }, { 0, 0, 19.713 }}
  Thermostat: 0 , 0
    Barostat: {{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }}
    </FrameData>
    <StuntDoubles>
         0  pvqjcw        0.0 0.0 0.0   0.0 0.0 0.0  1.0 0.0 0.0 0.0  0.0 0.0 0.0   0.0 0.0
    </StuntDoubles>
  </Snapshot>
</OpenMD>
