#ifndef _MADRID2019_INC_
#define _MADRID2019_INC_

molecule{
  name = "TIP4P-2005";
  
  atom[0]{
    type = "O_TIP4P-2005";
    position( 0.0, 0.0, -0.06556 );
  }
  atom[1]{
    type = "H_TIP4P-2005";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP4P-2005";
    position( 0.0, -0.75695, 0.52032 );
  }
  atom[3]{
    type = "EP_TIP4P-2005";
    position( 0.0, 0.0, 0.08904 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2, 3);
  }  
}

molecule {
  name = "H3O+";
  atom[0] { type = "OH3"; position( 0, 0, 0.046805);}
  atom[1] { type = "H3"; position( 0.93473, 0, -0.247609);}
  atom[2] { type = "H3"; position( -0.467365, -0.8095, -0.247609);}
  atom[3] { type = "H3"; position( -0.467365, 0.8095, -0.247609);}

  rigidBody[0] { members(0, 1, 2, 3); }
}

molecule {
  name = "SO42-";
  atom[0] { type = "SSO4"; position( 0, 0, 0);}
  atom[1] { type = "OSO4"; position( 0.860252, 0.860252, 0.860252);}
  atom[2] { type = "OSO4"; position( 0.860252, -0.860252, -0.860252);}
  atom[3] { type = "OSO4"; position( -0.860252, 0.860252, -0.860252);}
  atom[4] { type = "OSO4"; position( -0.860252, -0.860252, 0.860252);}

  rigidBody[0] { members(0, 1, 2, 3, 4); }
}

molecule {
  name = "OH-";
  atom[0] { type = "OHI"; position( 0, 0, -0.057647);}
  atom[1] { type = "HI"; position( 0, 0, 0.922353);}

  rigidBody[0] { members(0, 1); }
}

molecule {
  name = "NO3-";
  atom[0] { type = "N1"; position( 0.0, 0.0, 0.0); }
  atom[1] { type = "ON"; position( 1.256, 0.0, 0.0); }
  atom[2] { type = "ON"; position( -0.628, 1.087728, 0.0); }
  atom[3] { type = "ON"; position( -0.628, -1.087728, 0.0); }

  rigidBody[0] { members(0, 1, 2, 3); }
}

molecule {
  name = "NH4+";
  atom[0] { type = "NNH4"; position( 0, 0, 0);}
  atom[1] { type = "HNH4"; position( 0.583124, 0.583124, 0.583124);}
  atom[2] { type = "HNH4"; position( 0.583124, -0.583124, -0.583124);}
  atom[3] { type = "HNH4"; position( -0.583124, 0.583124, -0.583124);}
  atom[4] { type = "HNH4"; position( -0.583124, -0.583124, 0.583124);}

  rigidBody[0] { members(0, 1, 2, 3, 4); }
}

molecule{
  name = "Na+";
  atom[0] {
    type = "Na+";
    position(0, 0, 0);
  }
}

molecule{
  name = "Mg2+";
  atom[0] {
    type = "Mg2+";
    position(0, 0, 0);
  }
}

molecule{
  name = "Ca2+";
  atom[0] {
    type = "Ca2+";
    position(0, 0, 0);
  }
}

molecule{
  name = "K+";
  atom[0] {
    type = "K+";
    position(0, 0, 0);
  }
}

molecule{
  name = "Li+";
  atom[0] {
    type = "Li+";
    position(0, 0, 0);
  }
}

molecule{
  name = "Rb+";
  atom[0] {
    type = "Rb+";
    position(0, 0, 0);
  }
}

molecule{
  name = "Cs+";
  atom[0] {
    type = "Cs+";
    position(0, 0, 0);
  }
}

molecule{
  name = "F-";
  atom[0] {
    type = "F-";
    position(0, 0, 0);
  }
}

molecule{
  name = "Cl-";
  atom[0] {
    type = "Cl-";
    position(0, 0, 0);
  }
}

molecule{
  name = "Br-";
  atom[0] {
    type = "Br-";
    position(0, 0, 0);
  }
}

molecule{
  name = "I-";
  atom[0] {
    type = "I-";
    position(0, 0, 0);
  }
}

molecule{
  name = "Ba2+";
  atom[0] {
    type = "Ba2+";
    position(0, 0, 0);
  }
}

molecule{
  name = "Sr2+";
  atom[0] {
    type = "Sr2+";
    position(0, 0, 0);
  }
}

#endif
