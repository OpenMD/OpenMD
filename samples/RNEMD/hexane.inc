#ifndef __HEXANE_MD__
#define __HEXANE_MD__

molecule{

  name = "Hexane";

  atom[0]{
    type = "CH3";
    position( -2.16600000, 1.40700000, 0.00000000 );
  }

  atom[1]{
    type = "CH2";
    position( -0.62600000, 1.40700000, 0.00000000 );
  }

  atom[2]{
    type = "CH2";
    position( 0.00000000, 0.00000000, 0.00000000 );
  }

  atom[3]{
    type = "CH2";
    position( 1.54000000, 0.00000000, 0.00000000 );
  }

  atom[4]{
    type = "CH2";
    position( 2.16600000, -1.40700000, 0.00000000 );
  }

  atom[5]{
    type = "CH3";
    position( 3.70600000, -1.40700000, 0.00000000 );
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

  bond{
    members( 4, 5 );
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

  bend{
    members( 3, 4, 5 );
  }

  
  torsion{
    members( 0, 1, 2, 3 );
  }

  torsion{
    members( 1, 2, 3, 4 );
  }

  torsion{
    members( 2, 3, 4, 5 );
  }
   
}

#endif