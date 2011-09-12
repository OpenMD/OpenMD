#ifndef __ALKANES_MD__
#define __ALKANES_MD__

molecule{
  name = "methane";
  
  atom[0]{
    type="CH4";
    position( 0.0, 0.0, 0.0 );
  }
}

molecule{
  name = "ethane";
  
  atom[0]{
    type = "CH3";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "CH3";
    position( 0.0, 1.54, 0.0 );
  }

  
  bond{
    members( 0, 1 );
  }
}

molecule{
  name = "ethane2";
  
  atom[0]{
    type = "CH3";
    position( 0.0, 0.0, 0.0 );
  }
  atom[1]{
    type = "CH3";
    position( 0.0, 1.54, 0.0 );
  }

  
  bond{
    members( 0, 1 );
  }
}



molecule{

  name = "propane";
  
  
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

  

  bond{
    members( 0, 1 );
  }

  bond{
    members( 1, 2 );
  }

  
  bend{
    members( 0, 1, 2 );
  }
}

molecule{

  name = "butane";
  
  
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
  

  

  bond{
    members( 0, 1 );
  }

  bond{
    members( 1, 2 );
  }

  bond{
    members( 2, 3 );
  }

  
  bend{
    members( 0, 1, 2 );
  }

  bend{
    members( 1, 2, 3 );
  }

  
  torsion{
    members( 0, 1, 2, 3 );
  }

}

molecule{

  name = "pentane";
  
  
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

molecule{

  name = "pseudoButane";
  
  
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
  

  

  bond{
    members( 0, 1 );
  }

  bond{
    members( 1, 2 );
  }

  bond{
    members( 2, 3 );
  }

  
  bend{
    members( 0, 1, 2 );
  }

  bend{
    members( 1, 2, 3 );
  }

  
  torsion{
    members( 0, 1, 2, 3 );
  }

  
  rigidBody[0]{
    
    members(0);
  }
  rigidBody[1]{
    
    members(1);
  }

  rigidBody[2]{
    
    members(2);
  }

  rigidBody[3]{
    
    members(3);
  }


}
#endif
