#ifndef _EULER3_H_
#define _EULER3_H_

class Quaternion;
class Mat3x3d;
class Vector3d;

class Euler3{
  public:
    Euler3(){
      phi = 0;
      theta = 0;
      psi = 0;
    }
    Euler3( double phi, double theta, double psi){
      this->phi = phi;
      this->theta = theta;
      this->psi = psi;
    }

    Euler3(double e[3]){
      phi = e[0];
      theta = e[1];
      psi = e[2];
    }

    Euler3(const Vector3d& v);
    
    Euler3(Mat3x3d& m);
    
    Euler3(Quaternion& q);
    
  public:
    union{
      struct{
        double phi;
        double theta;
        double psi;
      };
      double angle[3];
    };
};

#endif //endif idndef _EULER3_H_
