#ifndef _SRI_H_
#define _SRI_H_

#include <iostream>

#include "Atom.hpp"
#include "DirectionalAtom.hpp"
#include "AbstractClasses.hpp"

// a little home-made vector structure

struct vect{
  double x;
  double y;
  double z;
  double length;
};

/************************************************************************
 * 
 * This section describes the base bond, bend, and torsion
 * classes. later these classes will be extended to good/evil ends.
 *
 ************************************************************************/

class Bond : public SRI{

public:
  Bond();
  virtual ~Bond();

  void calc_forces();
  int is_constrained() {return c_is_constrained;}
  Constraint *get_constraint() {return c_constraint;}
  void constrain(double bond_distance);

protected:
  virtual double bond_force(double r_ab) = 0;
  void set_atoms( Atom &, Atom & );

  int c_is_constrained;
  Constraint *c_constraint;
  Atom * c_p_a; /* atom a */
  Atom * c_p_b; /* atom b */
};


class Bend : public SRI{
 
public:
  Bend() {}
  virtual ~Bend() {}
  
  virtual void calc_forces();
  int is_constrained() {return 0;}
  Constraint *get_constraint() {return NULL;}
  void constrain(double bond_distance){} /*meaningless for bends */

protected:
  virtual double bend_force(double theta) = 0;
  void set_atoms( Atom &, Atom &, Atom & );

  Atom * c_p_a; /* atom a */
  Atom * c_p_b; /* atom b */
  Atom * c_p_c; /* atom c */
};

class Torsion : public SRI{

public:
  Torsion() {}
  virtual ~Torsion() {}

  void calc_forces();
  int is_constrained() {return 0;}
  Constraint *get_constraint() {return NULL;}
  void constrain(double bond_distance){} /*meaningless for torsions */

  

protected:

  void set_atoms(Atom &, Atom &, Atom &, Atom &);
  virtual double torsion_force(double cos_phi) = 0;

  Atom * c_p_a;
  Atom * c_p_b;
  Atom * c_p_c;
  Atom * c_p_d;
};

/**********************************************************************
 * 
 * These next classes are extensions of the base classes. These are
 * the actual objects which will be used in the simulation.
 *
 **********************************************************************/

class ConstrainedBond : public Bond{

public:
  ConstrainedBond( Atom &a, Atom &b, double constraint );
  ~ConstrainedBond() {}
  
  void printMe( void ){
    std::cerr << c_p_a->getType() << " - " << c_p_b->getType() 
	      << ": " << c_p_a->getIndex() << " - " 
	      << c_p_b->getIndex()
	      << ", d0 = " << d0 << "\n";
  }

private:
  double bond_force( double r_ab ){ return 0.0; }
  double d0;
};

class HarmonicBond : public Bond{

public:
  HarmonicBond(Atom &a, Atom &b, double theR0, double theK0 );
  ~HarmonicBond(){}

  void printMe( void ){
    std::cerr << c_p_a->getType() << " - " << c_p_b->getType() 
	      << ": " << c_p_a->getIndex() << " - " 
	      << c_p_b->getIndex()
	      << ", d0 = " << d0 << ", k0 = " << k0 <<"\n";
  }

  private:
  double bond_force( double r_ab );
  double d0;
  double k0;

};

class QuadraticBend : public Bend{

public:
  QuadraticBend( Atom &a, Atom &b, Atom &c );
  ~QuadraticBend(){}

  void setConstants( double the_c1, double the_c2, double the_c3, 
		     double the_Th0 );
  void printMe( void ){
    std::cerr << c_p_a->getType() << " - " << c_p_b->getType() << " - "
	      << c_p_c->getType() << " : " 
	      << c_p_a->getIndex() << " - " << c_p_b->getIndex() << " - "
	      << c_p_c->getIndex() 
	      <<", k1 = " << c1 << "; k2 = " << c2 
	      << "; k3 = " << c3 << "; theta0 =" << theta0 << "\n";
  }

private:
  double bend_force( double theta );

  double c1, c2, c3;
  double theta0;
};

class GhostBend : public Bend{

public:
  GhostBend( Atom &a, Atom &b );
  ~GhostBend(){}

  void calc_forces( void );

  void setConstants( double the_c1, double the_c2, double the_c3, 
		     double the_Th0 );
  void printMe( void ){
    std::cerr << c_p_a->getType() << " - " << c_p_b->getType()
	      << " : " 
	      << c_p_a->getIndex() << " - " << c_p_b->getIndex() << " - "
	      <<", k1 = " << c1 << "; k2 = " << c2 
	      << "; k3 = " << c3 << "; theta0 =" << theta0 << "\n";
  }

private:
  double bend_force( double theta );

  double c1, c2, c3;
  double theta0;

  DirectionalAtom* atomB;
};

class CubicTorsion : public Torsion{

public:
  CubicTorsion( Atom &a, Atom &b, Atom &c, Atom &d );
  ~CubicTorsion() {}

  void setConstants( double the_k1, double the_k2, double the_k3,
		     double the_k4 );
  void printMe( void ){
    std::cerr << c_p_a->getType() << " - " << c_p_b->getType() << " - "
	      << c_p_c->getType() << " - " << c_p_d->getType() << ": "
	      << c_p_a->getIndex() << " - " << c_p_b->getIndex() << " - "
	      << c_p_c->getIndex() << " - " << c_p_d->getIndex()
	      << ", k1 = " << k1 << "; k2 = " << k2 
	      << "; k3 = " << k3 << "; k4 =" << k4 << "\n";
  }
  
private:
  
  double torsion_force( double cos_phi );
  
  double k1, k2, k3, k4;
};
  
#endif
