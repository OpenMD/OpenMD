<OpenMD version=1>
  <MetaData>

molecule {
  name = "TIP4P_FQ";

  atom[0]{
    type = "O_TIP4P_FQ";
    position( 0.0, 0.0, -0.06556 );
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
  constrainTotalCharge = true;
}
component{
  type = "TIP4P_FQ";
  nMol = 2;
}

statFileFormat = "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY|ELECTRONIC_TEMPERATURE";
ensemble = NVE;
forceField = "FlucQ";
forceFieldFileName = "FQ.frc";
cutoffMethod = "shifted_force";
cutoffRadius = 9.0;
outputFluctuatingCharges = true;

flucQ {
 targetTemp = 1.0;
 tauThermostat = 100.0;
}


dt = 0.1;
runTime = 1e4;
useInitialTime = "false";
useInitialExtendedSystemState = "false";

sampleTime = 1;
statusTime = 1;
  </MetaData>
  <Snapshot>
    <FrameData>
        Time: 20000
        Hmat: {{ 19.713, 0, 0 }, { 0, 19.713, 0 }, { 0, 0, 19.713 }}
    </FrameData>
    <StuntDoubles>
         0    pvqj          0.0 0.0 0.0 -9.930000e-04 -2.397000e-03 -1.573000e-03  5.072750e-01 -7.075410e-01  3.148250e-01 -3.780790e-01  1.549000e-03  3.620000e-03  1.009100e-02
         1    pvqj          5.0 0.0 0.0 -2.782000e-03 -6.600000e-04  4.056000e-03  9.326000e-03  1.858280e-01 -9.793550e-01 -7.903000e-02 -6.290000e-04 -1.134100e-02 -1.901100e-02
    </StuntDoubles>
  </Snapshot>
</OpenMD>
