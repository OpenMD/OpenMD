#ifndef _WATER_INC_
#define _WATER_INC_

molecule{
  name = "Cl-";
  
  atom[0]{
    type = "Cl-";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "Na+";
  
  atom[0]{
    type = "Na+";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "SSD_E";
  
  atom[0]{
    type = "SSD_E";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "SSD_RF";
  
  atom[0]{
    type = "SSD_RF";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "SSD";
  
  atom[0]{
    type = "SSD";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "SSD1";
  
  atom[0]{
    type = "SSD1";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "TRED";
  
  atom[0]{
    type = "TRED";
    position( 0.0, 0.0, 0.0 );
    orientation( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "EP_TRED";
    position( 0.0, 0.0, 0.5 );
  }

  rigidBody[0]{
    members(0, 1);
  }
}

molecule{
  name = "TIP3P";
  
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
  
  rigidBody[0]{     
    members(0, 1, 2);
  }  
}

molecule{
  name = "TIP3P-FB";
  
  atom[0]{
    type = "O_TIP3P-FB";
    position( 0.0, 0.0, -0.066424 );
  }
  atom[1]{
    type = "H_TIP3P-FB";
    position( 0.0, 0.819341, 0.527225 );
  }
  atom[2]{
    type = "H_TIP3P-FB";
    position( 0.0, -0.819341, 0.527225 );
  }
  
  rigidBody[0]{     
    members(0, 1, 2);
  }  
}

molecule{
  name = "TIP4P";
  
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
  
  rigidBody[0]{    
    members(0, 1, 2, 3);
  }  
}

molecule{
  name = "TIP4P-Ice";
  
  atom[0]{
    type = "O_TIP4P-Ice";
    position( 0.0, 0.0, -0.06556 );
  }
  atom[1]{
    type = "H_TIP4P-Ice";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP4P-Ice";
    position( 0.0, -0.75695, 0.52032 );
  }
  atom[3]{
    type = "EP_TIP4P-Ice";
    position( 0.0, 0.0, 0.09214 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2, 3);
  }  
}

molecule{
  name = "TIP4P-Ew";
  
  atom[0]{
    type = "O_TIP4P-Ew";
    position( 0.0, 0.0, -0.06556 );
  }
  atom[1]{
    type = "H_TIP4P-Ew";
    position( 0.0, 0.75695, 0.52032 );
  }
  atom[2]{
    type = "H_TIP4P-Ew";
    position( 0.0, -0.75695, 0.52032 );
  }
  atom[3]{
    type = "EP_TIP4P-Ew";
    position( 0.0, 0.0, 0.05944 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2, 3);
  }
}

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

molecule{
  name = "TIP4P-FB";
  
  atom[0]{
    type = "O_TIP4P-FB";
    position( 0.0, 0.0, -0.0655549 );
  }
  atom[1]{
    type = "H_TIP4P-FB";
    position( 0.0, 0.75695, 0.520327 );
  }
  atom[2]{
    type = "H_TIP4P-FB";
    position( 0.0, -0.75695, 0.520327 );
  }
  atom[3]{
    type = "EP_TIP4P-FB";
    position( 0.0, 0.0, 0.0397151 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2, 3);
  }  
}

molecule{
  name = "TIP5P";
  
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
  
  rigidBody[0]{
    members(0, 1, 2, 3, 4);
  }
}

molecule{
  name = "TIP5P-E";
  
  atom[0]{
    type = "O_TIP5P-E";
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
  
  rigidBody[0]{    
    members(0, 1, 2, 3, 4);
  }
}

molecule{
  name = "NE6";
  
  atom[0]{
    type = "O_NE6";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "H_NE6";
    position( 0.0, 0.576029, 0.79283665 );
  }
  atom[2]{
    type = "H_NE6";
    position( 0.0, -0.576029, 0.79283665 );
  }
  atom[3]{
    type = "EP_NE6";
    position( 0.0, 0.23, 0.0 );
  }
  atom[4]{
    type = "LP_NE6";
    position( 0.732813007, -0.50364843, 0.0  );
  }
  atom[5]{
    type = "LP_NE6";
    position( -0.732813007, -0.50364843, 0.0  );
  }
 
  rigidBody[0]{    
    members(0, 1, 2, 3, 4, 5);
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

molecule{
  name = "SPC";
  
  atom[0]{
    type = "O_SPC";
    position( 0.0, 0.0, -0.06461 );
  }
  atom[1]{
    type = "H_SPC";
    position( 0.0, 0.81649, 0.51275 );
  }
  atom[2]{
    type = "H_SPC";
    position( 0.0, -0.81649, 0.51275 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
}


molecule{
  name = "SPC-HW";
  
  atom[0]{
    type = "O_SPC-HW";
    position( 0.0, 0.0, -0.06461 );
  }
  atom[1]{
    type = "D_SPC-HW";
    position( 0.0, 0.81649, 0.51275 );
  }
  atom[2]{
    type = "D_SPC-HW";
    position( 0.0, -0.81649, 0.51275 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
}

molecule{
  name = "DPD";
  
  atom[0]{
    type = "DPD";
    position(0.0, 0.0, 0.0);
  }
}

molecule{
  name = "OPC";
  
  atom[0]{
    type = "O_OPC";
    position( 0.0, 0.0, -0.0603651 );
  }
  atom[1]{
    type = "H_OPC";
    position( 0.0, 0.685582, 0.479134 );
  }
  atom[2]{
    type = "H_OPC";
    position( 0.0, -0.685582, 0.479134 );
  }
  atom[3]{
    type = "EP_OPC";
    position( 0.0, 0.0, 0.0990349 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2, 3);
  }  
}

molecule{
  name = "OPC3";
  
  atom[0]{
    type = "O_OPC3";
    position( 0.0, 0.0, -0.0632382 );
  }
  atom[1]{
    type = "H_OPC3";
    position( 0.0, 0.799262, 0.501939 );
  }
  atom[2]{
    type = "H_OPC3";
    position( 0.0, -0.799262, 0.501939 );
  }
  
  rigidBody[0]{    
    members(0, 1, 2);
  }
}


#endif
