#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include <cmath>
#include <iostream>

using namespace std;

class Vector3d{

  public:

     Vector3d(){
      this->x = double();
      this->y = double();
      this->z = double();
      
    }

     Vector3d( double x, double y, double z){
      this->x = x;
      this->y = y;
      this->z = z;
    }

     Vector3d(const Vector3d& v1){
      this->x = v1.x;
      this->y = v1.y;  
      this->z = v1.z;
    }

    double& operator[](unsigned int index){
      switch (index){
        case 0 :
          return x;

        case 1 :
          return y;

        case 2 :
          return z;

        default:
          cerr << index <<" is invalid index" << endl;
    	  exit(1);
      }
    }

    const double& operator[](unsigned int index) const{
      switch (index){
        case 0 :
          return x;

        case 1 :
          return y;

        case 2 :
          return z;

        default:
          cerr << index <<" is invalid index" << endl;
    	  exit(1);
      }
    }

    Vector3d& operator=( const Vector3d& v1 ){
      if(this == & v1)
        return *this;
      
      this->x = v1.x;
      this->y = v1.y;  
      this->z = v1.z;

      return *this;
    }

    bool operator ==( const Vector3d& v1 ){
      return this->x == v1.x && this->y == v1.y && this->z == v1.z;
    }

    bool operator !=( const Vector3d& v1 ){
      return this->x != v1.x || this->y != v1.y || this->z != v1.z;
    }

    void neg(){
      this->x = -this->x;
      this->y = -this->y;
      this->z = -this->z;
    }

    void add( const Vector3d& v1 ){
      this->x += v1.x;
      this->y += v1.y;
      this->z += v1.z;
    }

    void add( const Vector3d& v1, const Vector3d  &v2 ){
      this->x = v1.x + v1.x;
      this->y = v1.y + v2.y;
      this->z = v1.z + v2.z;

    }
        
    void sub( const Vector3d& v1 ){
      this->x -= v1.x;
      this->y -= v1.y;
      this->z -= v1.z;
    }

    void sub( const Vector3d& v1, const Vector3d  &v2 ){
      this->x = v1.x - v1.x;
      this->y = v1.y - v2.y;
      this->z = v1.z - v2.z;
    }

    void mul( double r ){
      this->x *= r;
      this->y *= r;
      this->z *= r;
    }

    void mul( double r, const Vector3d& v1 ){
      this->x = r * v1.x;
      this->y = r * v1.y;
      this->z = r * v1.z;

    }

    void mul( const Vector3d& v1, double r ){
      this->x = v1.x * r;
      this->y = v1.y * r;
      this->z = v1.z * r;
    }

    void div( double r){
      this->x /= r;
      this->y /= r;
      this->z /= r;
    }
         
     void div( const Vector3d& v1, double r ){
      this->x = v1.x/r;
      this->y = v1.y/r;
      this->z = v1.z/r;
    }

    void operator +=( const Vector3d& v1 ){
      this->x += v1.x;
      this->y += v1.y;
      this->z += v1.z;
    }
           
    void operator -=( const Vector3d& v1 ){
      this->x -= v1.x;
      this->y -= v1.y;
      this->z -= v1.z;
    }
          
    void operator *=( double r ){
      this->x *= r;
      this->y *= r;
      this->z *= r;
    }
          
    void operator /=( double r ){
      this->x /= r;
      this->y /= r;
      this->z /= r;
    }

    //vector addition
    friend Vector3d operator+ ( const Vector3d& v1, const Vector3d& v2){
      return Vector3d(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z); 
    }

     //vector subtraction
     friend Vector3d operator- ( const Vector3d& v1, const Vector3d& v2) {
      return Vector3d(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z); 
     }
     //unary minus
     friend Vector3d operator- ( const Vector3d& v) {
      return Vector3d(-v.x, -v.y, -v.z); 
     }
     //multiply by a scalar
     friend Vector3d operator* ( const double& r, const Vector3d& v) {
      return Vector3d( r*v.x, r*v.y, r*v.z); 
     }
     
     //multiply by a scalar
     friend Vector3d operator* ( const Vector3d& v, const double& r) {
      return Vector3d( r*v.x, r*v.y, r*v.z); 
     }
     //divide by a scalar
     friend Vector3d operator/ ( const Vector3d& v, const double& r) {
      return Vector3d( v.x/r, v.y/r, v.z/r); 
     }


    //dot product    
    friend double dotProduct (const Vector3d& v1, const Vector3d& v2){
      return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    //cross product
    friend Vector3d crossProduct (const Vector3d& v1, const Vector3d& v2){
      Vector3d result;
      
      result.x =   v1.y * v2.z - v1.z * v2.y ;
      result.y = -v1.x * v2.z + v1.z * v2.x ;
      result.z =   v1.x * v2.y - v1.y * v2.x;

      return result;
    }
    
    void normalize(){
      double len;
      
      len = length();
      x /= len ;
      y /= len;
      z /= len;
    }

    double length(){
      return sqrt(x*x + y*y + z*z);
    }

    double length2(){
      return x*x + y*y + z*z;
    }

  public:

    double x;
    double y;
    double z;
};

    
#endif //end ifndef _VECTOR3D_H_
