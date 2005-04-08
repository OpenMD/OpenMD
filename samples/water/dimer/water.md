#ifndef _WATER_MD_
#define _WATER_MD_

molecule{
  name = "Cl-";
  nAtoms = 1;
  atom[0]{
    type = "Cl-";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "Na+";
  nAtoms = 1;
  atom[0]{
    type = "Na+";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "SSD_E";
  nAtoms = 1;
  atom[0]{
    type = "SSD_E";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "SSD_RF";
  nAtoms = 1;
  atom[0]{
    type = "SSD_RF";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "SSD";
  nAtoms = 1;
  atom[0]{
    type = "SSD";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "SSD1";
  nAtoms = 1;
  atom[0]{
    type = "SSD1";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "TIP3P";
  nAtoms = 3;
  atom[0]{
    type = "O_TIP3P";
    position( 0.0, 0.0, -0.06556 );
  }
  atom[1]{
    type = "H_TIP3P";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP3P";
    position( 0.0, -0.75695, 0.52032 );
  }

  nRigidBodies = 1;
  rigidBody[0]{ 
    nMembers = 3;
    members(0, 1, 2);
  }

  nCutoffGroups = 1;
  cutoffGroup[0]{ 
    nMembers = 3;
    members(0, 1, 2);
  }
}

molecule{
  name = "TIP4P";
  nAtoms = 4;
  atom[0]{
    type = "O_TIP4P";
    position( 0.0, 0.0, -0.06556 );
  }
  atom[1]{
    type = "H_TIP4P";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP4P";
    position( 0.0, -0.75695, 0.52032 );
  }
  atom[3]{
    type = "EP_TIP4P";
    position( 0.0, 0.0, 0.08444 );
  }
  nRigidBodies = 1;
  rigidBody[0]{
    nMembers = 4;
    members(0, 1, 2, 3);
  }

  nCutoffGroups = 1;
  cutoffGroup[0]{ 
    nMembers = 4;
    members(0, 1, 2, 3);
  }
}

molecule{
  name = "TIP5P";
  nAtoms = 5;
  atom[0]{
    type = "O_TIP5P";
    position( 0.0, 0.0, -0.06556 );
  }
  atom[1]{
    type = "H_TIP5P";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP5P";
    position( 0.0, -0.75695, 0.52032 );
  }
  atom[3]{
    type = "EP_TIP5P";
    position( 0.57154, 0.0, -0.46971 );
  }
  atom[4]{
    type = "EP_TIP5P";
    position( -0.57154, 0.0, -0.46971 );
  }
  nRigidBodies = 1;
  rigidBody[0]{
    nMembers = 5;
    members(0, 1, 2, 3, 4);
  }

  nCutoffGroups = 1;
  cutoffGroup[0]{ 
    nMembers = 5;
    members(0, 1, 2, 3, 4);
  }
}

molecule{
  name = "SPCE";
  nAtoms = 3;
  atom[0]{
    type = "O_SPCE";
    position( 0.0, 0.0, -0.06461 );
  }
  atom[1]{
    type = "H_SPCE";
    position( 0.0, 0.81649, 0.51275 );
  }
  atom[2]{
    type = "H_SPCE";
    position( 0.0, -0.81649, 0.51275 );
  }
  nRigidBodies = 1;
  rigidBody[0]{
    nMembers = 3;
    members(0, 1, 2);
  }

  nCutoffGroups = 1;
  cutoffGroup[0]{ 
    nMembers = 3;
    members(0, 1, 2);
  }

}

molecule{
  name = "DPD";
  nAtoms = 1;
  atom[0]{
    type = "DPD";
    position(0.0, 0.0, 0.0);
  }
}

#endif
