molecule {

  name = "HEAD";
  nAtoms = 3;

  atom[0]{
    type = "HDP";
    position( 0.0,0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 ); //?
  }

  atom[1]{
    type = "NC4";
    position( 0.0, 0.0, 2.315 );
  }

  atom[2]{
    type = "PO4";
    position( 0.0, 0.0, -2.315 );
  }

  nRigidBodies = 1;
  rigidBody[0]{
    nMembers = 3;
    members(0, 1, 2);
  }

}
