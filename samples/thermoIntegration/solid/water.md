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

  
  cutoffGroup{
    
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

  
  cutoffGroup{
    
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

  
  cutoffGroup{
    
    members(0, 1, 2, 3, 4);
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

  
  cutoffGroup{
    
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
