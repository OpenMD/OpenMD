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
  nMol = 1;
}

flucQ {
 targetTemp = 10.0;
 tauThermostat = 10.0;
 dragCoefficient = 0.0001;
}

statFileFormat = "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY|ELECTRONIC_TEMPERATURE";

ensemble = NVE;
forceField = "FlucQ";
forceFieldFileName = "FQ.frc";
cutoffMethod = "shifted_force";
cutoffRadius = 9.0;
outputFluctuatingCharges = true;
runTime = 1e6;
sampleTime = 100;
statusTime = 10;

dt = 1;
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
         0  pvqj        0.0 0.0 0.0   0.0 0.0 0.0  1.0 0.0 0.0 0.0  0.0 0.0 0.0 
    </StuntDoubles>
  </Snapshot>
</OpenMD>
