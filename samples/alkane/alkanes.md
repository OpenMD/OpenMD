#ifndef __ALKANES_MD__
#define __ALKANES_MD__

molecule{
  name = "methane";
  nAtoms = 1;
  atom[0]{
    type="CH4";
    position( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "ethane";
  nAtoms = 2;
  atom[0]{
    type = "CH3";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "CH3";
    position( 0.0, 1.54, 0.0 );
  }

  nBonds = 1;
  bond[0]{
    members( 0, 1 );
  }
}

molecule{
  name = "ethane2";
  nAtoms = 2;
  atom[0]{
    type = "CH3";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "CH3";
    position( 0.0, 1.54, 0.0 );
  }

  nBonds = 1;
  bond[0]{
    members( 0, 1 );
  }
}



molecule{

  name = "propane";
  nAtoms = 3;
  
  atom[0]{
    type = "CH3";
    position( -0.626, 1.407, 0.0 );
  }

  atom[1]{
    type = "CH2";
    position( 0.0, 0.0, 0.0 );
  }

  atom[2]{
    type = "CH3";
    position( 1.54, 0.0, 0.0 );
  }

  nBonds = 2;

  bond[0]{
    members( 0, 1 );
  }

  bond[1]{
    members( 1, 2 );
  }

  nBends = 1;
  bend[0]{
    members( 0, 1, 2 );
  }
}

molecule{

  name = "butane";
  nAtoms = 4;
  
  atom[0]{
    type = "CH3";
    position( -0.626, 1.407, 0.0 );
  }

  atom[1]{
    type = "CH2";
    position( 0.0, 0.0, 0.0 );
  }

  atom[2]{
    type = "CH2";
    position( 1.54, 0.0, 0.0 );
  }

  atom[3]{
    type = "CH3";
    position( 2.166, -1.407, 0.0 );
  }
  

  nBonds = 3;

  bond[0]{
    members( 0, 1 );
  }

  bond[1]{
    members( 1, 2 );
  }

  bond[2]{
    members( 2, 3 );
  }

  nBends = 2;
  bend[0]{
    members( 0, 1, 2 );
  }

  bend[1]{
    members( 1, 2, 3 );
  }

  nTorsions = 1;
  torsion[0]{
    members( 0, 1, 2, 3 );
  }
}

molecule{

  name = "pentane";
  nAtoms = 5;
  
  atom[0]{
    type = "CH3";
    position( -0.626, 1.407, 0.0 );
  }

  atom[1]{
    type = "CH2";
    position( 0.0, 0.0, 0.0 );
  }

  atom[2]{
    type = "CH2";
    position( 1.54, 0.0, 0.0 );
  }

  atom[3]{
    type = "CH2";
    position( 2.166, -1.407, 0.0 );
  }

  atom[4]{
    type = "CH3";
    position( 3.706, -1.407, 0.0 );
  }
  

  nBonds = 4;

  bond[0]{
    members( 0, 1 );
  }

  bond[1]{
    members( 1, 2 );
  }

  bond[2]{
    members( 2, 3 );
  }

  bond[3]{
    members( 3, 4 );
  }

  nBends = 3;
  bend[0]{
    members( 0, 1, 2 );
  }

  bend[1]{
    members( 1, 2, 3 );
  }

  bend[2]{
    members( 2, 3, 4 );
  }

  nTorsions = 2;
  torsion[0]{
    members( 0, 1, 2, 3 );
  }

  torsion[1]{
    members( 1, 2, 3, 4 );
  }
}

molecule{

  name = "pseudoButane";
  nAtoms = 4;
  
  atom[0]{
    type = "CH3";
    position( -0.626, 1.407, 0.0 );
  }

  atom[1]{
    type = "CH2";
    position( 0.0, 0.0, 0.0 );
  }

  atom[2]{
    type = "CH2";
    position( 1.54, 0.0, 0.0 );
  }

  atom[3]{
    type = "CH3";
    position( 2.166, -1.407, 0.0 );
  }
  

  nBonds = 3;

  bond[0]{
    members( 0, 1 );
  }

  bond[1]{
    members( 1, 2 );
  }

  bond[2]{
    members( 2, 3 );
  }

  nBends = 2;
  bend[0]{
    members( 0, 1, 2 );
  }

  bend[1]{
    members( 1, 2, 3 );
  }

  nTorsions = 1;
  torsion[0]{
    members( 0, 1, 2, 3 );
  }

  nRigidBodies = 4;
  rigidBody[0]{
    nMembers = 1;
    members(0);
  }
  rigidBody[1]{
    nMembers = 1;
    members(1);
  }

  rigidBody[2]{
    nMembers = 1;
    members(2);
  }

  rigidBody[3]{
    nMembers = 1;
    members(3);
  }


}
#endif
