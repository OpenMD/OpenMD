#include "math/Euler3.hpp"
#include "Mat3x3d.hpp"
#include "Quaternion.hpp"
#include "Vector3d.hpp"

Euler3::Euler3(const Vector3d& v){
  this->phi = v.x;
  this->theta = v.y;
  this->psi = v.z;
}

Euler3::Euler3(Mat3x3d& m){
  *this = m.toEuler();
}

Euler3::Euler3(Quaternion& q){
  *this = q.toEuler();
}
