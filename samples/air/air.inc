#ifndef _AIR_INC_
#define _AIR_INC_

molecule{
  name = "He";
  
  atom[0]{
    type = "He";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "Ne";
  
  atom[0]{
    type = "Ne";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "Ar";
  
  atom[0]{
    type = "Ar";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "Kr";
  
  atom[0]{
    type = "Kr";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "CH4";
  
  atom[0]{
    type = "CH4";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "O2";
  
  atom[0]{
    type = "X_O2";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "O_O2";
    position( 0.0, 0.0, 0.605 );
  }
  atom[2]{
    type = "O_O2";
    position( 0.0, 0.0, -0.605 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
}

molecule{
  name = "N2";
  
  atom[0]{
    type = "X_N2";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "N_N2";
    position( 0.0, 0.0, 0.550 );
  }
  atom[2]{
    type = "N_N2";
    position( 0.0, 0.0, -0.550 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
}

molecule{
  name = "CO2";
  
  atom[0]{
    type = "C_CO2";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "O_CO2";
    position( 0.0, 0.0, 1.160 );
  }
  atom[2]{
    type = "O_CO2";
    position( 0.0, 0.0, -1.160 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
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


#endif
