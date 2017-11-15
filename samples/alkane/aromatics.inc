molecule{
   name = "toluene";
   
   atom[0]{
      type = "CHar";
      position(1.16383, 0.33380, 2.09425 );
    }
   atom[1]{
      type = "CHar";
      position(1.16404, 0.32490, 0.69741 );
    }
   atom[2]{
      type = "CHar";
      position(0.00000, 0.00000, 0.00000);
    }
   atom[3]{
      type = "CHar";
      position(-1.16627, -0.31714, 0.69742);
    }
   atom[4]{
      type = "CHar";
      position(-1.17048, -0.30932, 2.09426);
    }
   atom[5]{
      type = "RCHar";
      position(0.00000, 0.00000, 2.79623);
    }
   atom[6]{
      type = "CH3";
      position(-0.00080, 0.00357, 4.29930);
    }
   
    rigidBody[0]{
       members(0, 1, 2, 3, 4, 5, 6);
    }

    cutoffGroup{
       members(0, 1, 2, 3, 4, 5, 6);    
    }
}

molecule{
   name = "tolueneAA";
   
   atom[0]{
      type = "CA";
      position(0.000000, 0.000000, -1.857901);
    }
   atom[1]{
      type = "HA";
      position(0.000000, 0.000000, -2.937901);
   }
   atom[2]{
      type = "CA";
      position(0.000000, 1.212440, -1.157901);
    }
   atom[3]{
      type = "HA";
      position(0.000000, 2.147740, -1.697901);
   }
   atom[4]{
      type = "CA";
      position(0.000000, 1.212440, 0.242099);
    }
   atom[5]{
      type = "HA";
      position(0.000000, 2.147740, 0.782099);
   }
   atom[6]{
      type = "CA";
      position(0.000000, -1.212440, -1.157901);
    }
   atom[7]{
      type = "HA";
      position(0.000000, -2.147740, -1.697901);
   }
   atom[8]{
      type = "CA";
      position(0.000000, -1.212440, 0.242099);
    }
   atom[9]{
      type = "HA";
      position(0.000000, -2.147740, 0.782099);
   }
   atom[10]{
      type = "CA";
      position(0.000000, 0.000000, 0.942099);
    }
   atom[11]{
      type = "CTT";
      position(0.000000, 0.000000, 2.442099);
    }
   atom[12]{
      type = "HC";
      position(1.018050, 0.000000, 2.802599);
   }
   atom[13]{
      type = "HC";
      position(-0.509025, 0.881660, 2.802599);
   }
   atom[14]{
      type = "HC";
      position(-0.509025, -0.881660, 2.802599);
   }

   bond{
    members(0, 1);
   }
   bond{
    members(0, 2);
   }
   bond{
    members(0, 6);
   }
   bond{
    members(2, 3);
   }
   bond{
    members(2, 4);
   }
   bond{
    members(4, 5);
   }
   bond{
    members(4, 10);
   }
   bond{
    members(6, 7);
   }
   bond{
    members(6, 8);
   }
   bond{
    members(8, 9);
   }
   bond{
    members(8, 10);
   }
   bond{
    members(10, 11);
   }
   bond{
    members(11, 12);
   }
   bond{
    members(11, 13);
   }
   bond{
    members(11, 14);
   }
}
