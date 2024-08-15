molecule {
  name = "methanol";
  
  atom[0] { type = "C3"; }
  atom[1] { type = "OH"; }
  atom[2] { type = "HC"; }
  atom[3] { type = "HC"; }
  atom[4] { type = "HC"; }
  atom[5] { type = "HO"; }
  
  bond { members(1, 5); }
  bond { members(0, 4); }
  bond { members(0, 1); }
  bond { members(0, 3); }
  bond { members(0, 2); }
}

molecule{
  name = "SPCE";
  
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
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
}

