#ifndef _ABSTRACT_CLASSES_H_
#define _ABSTRACT_CLASSES_H_

#include <string>
#include "primitives/Atom.hpp"
#ifdef IS_MPI

#include <mpi.h>
#endif

using namespace std;

class Constraint{
  public:
    Constraint(){
    }
    ~Constraint(){
    }

    int get_a(){
      return a;
    }
    void set_a(int index){
      a = index;
    }
    int get_b(){
      return b;
    }
    void set_b(int index){
      b = index;
    }
    double get_dsqr(){
      return dsqr;
    }
    void set_dsqr(double ds){
      dsqr = ds;
    }

  private:
    int a; /* index of constrained atom a */
    int b; /* index of constrained atom b */
    double dsqr; /* the square of the constraint distance */
};

class SRI{
  public:
    SRI(){
      c_potential_E = 0.0;
    }
    virtual ~SRI(){
    }

    virtual void calc_forces() = 0;

    double get_potential(){
      return c_potential_E;
    }
    virtual int is_constrained() = 0;
    virtual Constraint* get_constraint() = 0;
    virtual void constrain(double bond_distance) = 0;

    virtual void printMe(void) = 0;

  protected:
    double c_potential_E;
};


class BaseIntegrator{
  public:
    BaseIntegrator(){
    }
    virtual ~BaseIntegrator(){
    }

    virtual void integrate(void) = 0;
    virtual double getConservedQuantity(void) = 0;
    virtual string getAdditionalParameters(void) = 0;
};
#endif
