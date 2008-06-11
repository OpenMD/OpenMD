#ifndef __ALKANES_MD__
#define __ALKANES_MD__

molecule{
  name = "S";
  
  atom[0]{
    type = "S";
    position( 0.0, 0.0, 0.0 );
  }
}
 
molecule{

  name = "Butanethiol";
  
  
  atom[0]{
    type = "S";
    position( -0.626, 1.709, 0.0 );
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
  

  

  bond{
    members( 0, 1 );
  }

  bond{
    members( 1, 2 );
  }

  bond{
    members( 2, 3 );
  }

  bond{
    members( 3, 4 );
  }

  
  bend{
    members( 0, 1, 2 );
  }

  bend{
    members( 1, 2, 3 );
  }

  bend{
    members( 2, 3, 4 );
  }

  
  torsion{
    members( 0, 1, 2, 3 );
  }

  torsion{
    members( 1, 2, 3, 4 );
  }
   
}

#endif
