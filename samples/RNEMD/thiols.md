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


molecule{

  name = "Hexanethiol";
  
  
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
    type = "CH2";
    position( 3.706, -1.407, 0.0 );
  }

  atom[5]{
    type = "CH2";
    position( 4.332, -2.814, 0.0 );
  }

  atom[6]{
    type = "CH3";
    position( 5.872, -2.814, 0.0 );
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

  bond{
    members( 5, 6 );
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

  bend{
    members( 4, 5, 6 );
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

  torsion{
    members( 3, 4, 5, 6 );
  }
   
}


molecule{

  name = "Octanethiol";
  
  
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
    type = "CH2";
    position( 3.706, -1.407, 0.0 );
  }

  atom[5]{
    type = "CH2";
    position( 4.332, -2.814, 0.0 );
  }

  atom[6]{
    type = "CH2";
    position( 5.872, -2.814, 0.0 );
  }
  
  atom[7]{
    type = "CH2";
    position( 6.498, -4.221, 0.0 );
  }

  atom[8]{
    type = "CH3";
    position( 8.038, -4.221, 0.0 );
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

  bond{
    members( 5, 6 );
  }

  bond{
    members( 6, 7 );
  }

  bond{
    members( 7, 8 );
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

  bend{
    members( 4, 5, 6 );
  }

  bend{
    members( 5, 6, 7 );
  }

  bend{
    members( 6, 7, 8 );
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

  torsion{
    members( 3, 4, 5, 6 );
  }

  torsion{
    members( 4, 5, 6, 7 );
  }

  torsion{
    members( 5, 6, 7, 8 );
  }
   
}


molecule{

  name = "Decanethiol";
  
  
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
    type = "CH2";
    position( 3.706, -1.407, 0.0 );
  }

  atom[5]{
    type = "CH2";
    position( 4.332, -2.814, 0.0 );
  }

  atom[6]{
    type = "CH2";
    position( 5.872, -2.814, 0.0 );
  }
  
  atom[7]{
    type = "CH2";
    position( 6.498, -4.221, 0.0 );
  }

  atom[8]{
    type = "CH2";
    position( 8.038, -4.221, 0.0 );
  }

  atom[9]{
    type = "CH2";
    position( 8.664, -5.628, 0.0 );
  }

  atom[10]{
    type = "CH3";
    position( 10.204, -5.628, 0.0 );
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

  bond{
    members( 5, 6 );
  }

  bond{
    members( 6, 7 );
  }

  bond{
    members( 7, 8 );
  }

  bond{
    members( 8, 9 );
  }

  bond{
    members( 9, 10 );
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

  bend{
    members( 4, 5, 6 );
  }

  bend{
    members( 5, 6, 7 );
  }

  bend{
    members( 6, 7, 8 );
  }

  bend{
    members( 7, 8, 9 );
  }

  bend{
    members( 8, 9, 10 );
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

  torsion{
    members( 3, 4, 5, 6 );
  }

  torsion{
    members( 4, 5, 6, 7 );
  }

  torsion{
    members( 5, 6, 7, 8 );
  }

  torsion{
    members( 6, 7, 8, 9 );
  }

  torsion{
    members( 7, 8, 9, 10 );
  }
   
}


molecule{

  name = "Dodecanethiol";
  
  
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
    type = "CH2";
    position( 3.706, -1.407, 0.0 );
  }

  atom[5]{
    type = "CH2";
    position( 4.332, -2.814, 0.0 );
  }

  atom[6]{
    type = "CH2";
    position( 5.872, -2.814, 0.0 );
  }
  
  atom[7]{
    type = "CH2";
    position( 6.498, -4.221, 0.0 );
  }

  atom[8]{
    type = "CH2";
    position( 8.038, -4.221, 0.0 );
  }

  atom[9]{
    type = "CH2";
    position( 8.664, -5.628, 0.0 );
  }

  atom[10]{
    type = "CH2";
    position( 10.204, -5.628, 0.0 );
  }

  atom[11]{
    type = "CH2";
    position( 10.830, -7.035, 0.0 );
  }

  atom[12]{
    type = "CH3";
    position( 12.370, -7.035, 0.0 );
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

  bond{
    members( 5, 6 );
  }

  bond{
    members( 6, 7 );
  }

  bond{
    members( 7, 8 );
  }

  bond{
    members( 8, 9 );
  }

  bond{
    members( 9, 10 );
  }

  bond{
    members( 10, 11 );
  }

  bond{
    members( 11, 12 );
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

  bend{
    members( 4, 5, 6 );
  }

  bend{
    members( 5, 6, 7 );
  }

  bend{
    members( 6, 7, 8 );
  }

  bend{
    members( 7, 8, 9 );
  }

  bend{
    members( 8, 9, 10 );
  }

  bend{
    members( 9, 10, 11 );
  }

  bend{
    members( 10, 11, 12 );
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

  torsion{
    members( 3, 4, 5, 6 );
  }

  torsion{
    members( 4, 5, 6, 7 );
  }

  torsion{
    members( 5, 6, 7, 8 );
  }

  torsion{
    members( 6, 7, 8, 9 );
  }

  torsion{
    members( 7, 8, 9, 10 );
  }

  torsion{
    members( 8, 9, 10, 11 );
  }

  torsion{
    members( 9, 10, 11, 12 );
  }
   
}

#endif
