  molecule{
    name = "CarbonMonoxide";
    atom[0]{
      type = "CO";
      position(0.0, 0.0, 0.0);
    }
    atom[1]{
      type = "OC";
      position(0.0, 0.0, 1.13);
    }
    atom[2]{
      type = "M_CO";
      position(0.0, 0.0, 0.6457);
    }
    rigidBody[0]{
      members(0,1,2);
    }
  }

