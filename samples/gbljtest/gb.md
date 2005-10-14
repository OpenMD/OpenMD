#ifndef __GB_MD__
#define __GB_MD__

molecule{
  name = "GB";
  nAtoms = 1;
  atom[0]{
     type="GB";
     position( 0.0, 0.0, 0.0 );
     orientation(0.0,0.0,0.0);
  }
}
molecule{
  name = "linear";
  nAtoms = 1;
  atom[0]{
     type="linear";
     position( 0.0, 0.0, 0.0 );
     orientation(0.0,0.0,0.0);
  }
}
#endif
