/*
 *  GeometryBuilder.cpp
 *  nanorodBuilder
 *
 *  Created by Charles Vardeman II on 4/4/05.
 *  Copyright 2005 University of Notre Dame. All rights reserved.
 *
 */

/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */

#include "GeometryBuilder.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <math.h>

using namespace std;
using namespace oopse;

//typedef CGAL::Homogeneous<int>              Kernel;
typedef CGAL::Simple_cartesian<double>     Kernel;
//typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

typedef Kernel::Point_3                    Point_3;
typedef Kernel::Vector_3                               Vector_3;
typedef CGAL::Polyhedron_traits_with_normals_3<Kernel> Traits;
//typedef CGAL::Polyhedron_3<Kernel>                     Polyhedron;
typedef CGAL::Polyhedron_3<Traits>                     Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Plane_iterator Plane_iterator;
typedef Polyhedron::Vertex_handle Vertex_handle;

Polyhedron nanoRodPolyhedron;
Polyhedron nanoRodTwinnedPolyhedron1;
Polyhedron nanoRodTwinnedPolyhedron2;
Polyhedron nanoRodTwinnedPolyhedron3;
Polyhedron nanoRodTwinnedPolyhedron4;
Polyhedron nanoRodTwinnedPolyhedron5;




//typedef CGAL::Scaling Scaling;
//typedef Aff_transformation_3<Kernel> A;( const Scaling,
//																Kernel::RT s=RT(20),
//																Kernel::RT hw = RT(1));





// A modifier creating a triangle with the incremental builder.
template <class HDS>
class buildSingleCrystal : public CGAL::Modifier_base<HDS> {
 public:
  Vertex_handle end1;
  Vertex_handle neight1;
  Vertex_handle end2;
  Vertex_handle neight2;
  Vertex_handle neight3;
  
  buildSingleCrystal() {}
  void operator()( HDS& hds) {
    // Postcondition: `hds' is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( 12, 15, 6);
    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;
    
    
    
    
    
    B.add_vertex( Point(-0.7887222926324,    0.4874571845315, -0.2562714077342));
    B.add_vertex( Point(-0.4874571845316,    0.4874571845315,  0.6709272557930));
    B.add_vertex( Point(-0.7887222926324,   -0.4874571845315, -0.2562714077342)); //End vertex
    end1 = B.add_vertex( Point( 0.0000000000000,    1.0000000000000,  0.0000000000000));
    neight3 = B.add_vertex( Point(-0.4874571845315,    -0.4874571845316,  0.6709272557930));
    neight1 = B.add_vertex( Point(-0.0000000000000,    0.4874571845316, -0.8293116961175));
    B.add_vertex( Point( 0.0000000000000,    -0.4874571845316, -0.8293116961175));
    B.add_vertex( Point( 0.4874571845315,    0.4874571845316,  0.6709272557930));
    end2 = B.add_vertex( Point(-0.0000000000000,    -1.0000000000000,  0.0000000000000)); //End Vertex
    B.add_vertex( Point( 0.7887222926324,    0.4874571845315, -0.2562714077342));
    neight2 = B.add_vertex( Point( 0.4874571845316,    -0.4874571845315,  0.6709272557930));
    B.add_vertex( Point( 0.7887222926324,    -0.4874571845315, -0.2562714077342));
    
    B.begin_facet();
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 9);
    B.add_vertex_to_facet( 11);
    B.add_vertex_to_facet( 10);
    B.end_facet();
				
    B.begin_facet();
    B.add_vertex_to_facet( 8);
    B.add_vertex_to_facet( 10);
    B.add_vertex_to_facet( 11);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 9);
    B.add_vertex_to_facet( 7);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 9);
    B.add_vertex_to_facet( 5);
    B.add_vertex_to_facet( 6);
    B.add_vertex_to_facet( 11);
    B.end_facet();
				
    B.begin_facet();
    B.add_vertex_to_facet( 8);
    B.add_vertex_to_facet( 11);
    B.add_vertex_to_facet( 6);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 5);
    B.add_vertex_to_facet( 9);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 5);
    B.add_vertex_to_facet( 0);
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 6);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 8);
    B.add_vertex_to_facet( 6);
    B.add_vertex_to_facet( 2);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 0);
    B.add_vertex_to_facet( 5);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 0);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 2);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 8);
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 4);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 0);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 10);
    B.add_vertex_to_facet( 4);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 8);
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 10);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 1);
    B.end_facet();
    
				
    B.end_surface();
  }
};



/* First Crystal*/

 template <class HDS>
 class buildtwinned1 : public CGAL::Modifier_base<HDS> {
public:
   Vertex_handle end1;
   Vertex_handle neight1;
   Vertex_handle end2;
   Vertex_handle neight2;
   Vertex_handle neight3;
   
   buildtwinned1() {}
   void operator()( HDS& hds) {
     // Postcondition: `hds' is a valid polyhedral surface.
     CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
     B.begin_surface( 12, 5, 6);
     typedef typename HDS::Vertex   Vertex;
     typedef typename Vertex::Point Point;
     
     
     
     
     
     B.add_vertex( Point(-0.7887222926324,    0.4874571845315, -0.2562714077342));
     B.add_vertex( Point(-0.4874571845316,    0.4874571845315,  0.6709272557930));
     B.add_vertex( Point(-0.7887222926324,   -0.4874571845315, -0.2562714077342)); //End vertex
     end1 = B.add_vertex( Point( 0.0000000000000,    1.0000000000000,  0.0000000000000));
     neight3 = B.add_vertex( Point(-0.4874571845315,    -0.4874571845316,  0.6709272557930));
     neight1 = B.add_vertex( Point(-0.0000000000000,    0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.0000000000000,    -0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.4874571845315,    0.4874571845316,  0.6709272557930));
     end2 = B.add_vertex( Point(-0.0000000000000,    -1.0000000000000,  0.0000000000000)); //End Vertex
     B.add_vertex( Point( 0.7887222926324,    0.4874571845315, -0.2562714077342));
     neight2 = B.add_vertex( Point( 0.4874571845316,    -0.4874571845315,  0.6709272557930));
     B.add_vertex( Point( 0.7887222926324,    -0.4874571845315, -0.2562714077342));
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 4);
     B.add_vertex_to_facet( 1);
     B.add_vertex_to_facet( 3);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 0);
     B.add_vertex_to_facet( 2);
     B.add_vertex_to_facet( 8);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 0);
     B.add_vertex_to_facet( 1);
     B.add_vertex_to_facet( 4);
     B.add_vertex_to_facet( 2);
     B.end_facet();    
     
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 1);
     B.add_vertex_to_facet( 0);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 4);
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 2);
     B.end_facet();
     
     B.end_surface();
   }
 };
 
 


/*second crystal*/
 template <class HDS>
 class buildtwinned2 : public CGAL::Modifier_base<HDS> {
public:
   Vertex_handle end1;
   Vertex_handle neight1;
   Vertex_handle end2;
   Vertex_handle neight2;
   Vertex_handle neight3;
   
   buildtwinned2() {}
   void operator()( HDS& hds) {
     // Postcondition: `hds' is a valid polyhedral surface.
     CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
     B.begin_surface( 12, 5, 6);
     typedef typename HDS::Vertex   Vertex;
     typedef typename Vertex::Point Point;
     
     
     
     
     
     B.add_vertex( Point(-0.7887222926324,    0.4874571845315, -0.2562714077342));
     B.add_vertex( Point(-0.4874571845316,    0.4874571845315,  0.6709272557930));
     B.add_vertex( Point(-0.7887222926324,   -0.4874571845315, -0.2562714077342)); //End vertex
     end1 = B.add_vertex( Point( 0.0000000000000,    1.0000000000000,  0.0000000000000));
     neight3 = B.add_vertex( Point(-0.4874571845315,    -0.4874571845316,  0.6709272557930));
     neight1 = B.add_vertex( Point(-0.0000000000000,    0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.0000000000000,    -0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.4874571845315,    0.4874571845316,  0.6709272557930));
     end2 = B.add_vertex( Point(-0.0000000000000,    -1.0000000000000,  0.0000000000000)); //End Vertex
     B.add_vertex( Point( 0.7887222926324,    0.4874571845315, -0.2562714077342));
     neight2 = B.add_vertex( Point( 0.4874571845316,    -0.4874571845315,  0.6709272557930));
     B.add_vertex( Point( 0.7887222926324,    -0.4874571845315, -0.2562714077342));
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 2);
     B.add_vertex_to_facet( 0);
     B.add_vertex_to_facet( 3);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 5);
     B.add_vertex_to_facet( 6);
     B.add_vertex_to_facet( 8);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 5);
     B.add_vertex_to_facet( 0);
     B.add_vertex_to_facet( 2);
     B.add_vertex_to_facet( 6);
     B.end_facet();    
     
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 0);
     B.add_vertex_to_facet( 5);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 6);
     B.add_vertex_to_facet( 2);
     B.end_facet();
     
     B.end_surface();
   }
 };

/* Third Crystal*/

 template <class HDS>
 class buildtwinned3 : public CGAL::Modifier_base<HDS> {
public:
   Vertex_handle end1;
   Vertex_handle neight1;
   Vertex_handle end2;
   Vertex_handle neight2;
   Vertex_handle neight3;
   
   buildtwinned3() {}
   void operator()( HDS& hds) {
     // Postcondition: `hds' is a valid polyhedral surface.
     CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
     B.begin_surface( 12, 5, 6);
     typedef typename HDS::Vertex   Vertex;
     typedef typename Vertex::Point Point;
     
     
     
     
     
     B.add_vertex( Point(-0.7887222926324,    0.4874571845315, -0.2562714077342));
     B.add_vertex( Point(-0.4874571845316,    0.4874571845315,  0.6709272557930));
     B.add_vertex( Point(-0.7887222926324,   -0.4874571845315, -0.2562714077342)); //End vertex
     end1 = B.add_vertex( Point( 0.0000000000000,    1.0000000000000,  0.0000000000000));
     neight3 = B.add_vertex( Point(-0.4874571845315,    -0.4874571845316,  0.6709272557930));
     neight1 = B.add_vertex( Point(-0.0000000000000,    0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.0000000000000,    -0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.4874571845315,    0.4874571845316,  0.6709272557930));
     end2 = B.add_vertex( Point(-0.0000000000000,    -1.0000000000000,  0.0000000000000)); //End Vertex
     B.add_vertex( Point( 0.7887222926324,    0.4874571845315, -0.2562714077342));
     neight2 = B.add_vertex( Point( 0.4874571845316,    -0.4874571845315,  0.6709272557930));
     B.add_vertex( Point( 0.7887222926324,    -0.4874571845315, -0.2562714077342));
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 6);
     B.add_vertex_to_facet( 5);
     B.add_vertex_to_facet( 3);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 9);
     B.add_vertex_to_facet( 11);
     B.add_vertex_to_facet( 8);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 9);
     B.add_vertex_to_facet( 5);
     B.add_vertex_to_facet( 6);
     B.add_vertex_to_facet( 11);
     B.end_facet();    
     
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 5);
     B.add_vertex_to_facet( 9);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 11);
     B.add_vertex_to_facet( 6);
     B.end_facet();
     
     B.end_surface();
   }
 };

/*Fourth Segment*/
 template <class HDS>
 class buildtwinned4 : public CGAL::Modifier_base<HDS> {
public:
   Vertex_handle end1;
   Vertex_handle neight1;
   Vertex_handle end2;
   Vertex_handle neight2;
   Vertex_handle neight3;
   
   buildtwinned4() {}
   void operator()( HDS& hds) {
     // Postcondition: `hds' is a valid polyhedral surface.
     CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
     B.begin_surface( 12, 5, 6);
     typedef typename HDS::Vertex   Vertex;
     typedef typename Vertex::Point Point;
     
     
     
     
     
     B.add_vertex( Point(-0.7887222926324,    0.4874571845315, -0.2562714077342));
     B.add_vertex( Point(-0.4874571845316,    0.4874571845315,  0.6709272557930));
     B.add_vertex( Point(-0.7887222926324,   -0.4874571845315, -0.2562714077342)); //End vertex
     end1 = B.add_vertex( Point( 0.0000000000000,    1.0000000000000,  0.0000000000000));
     neight3 = B.add_vertex( Point(-0.4874571845315,    -0.4874571845316,  0.6709272557930));
     neight1 = B.add_vertex( Point(-0.0000000000000,    0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.0000000000000,    -0.4874571845316, -0.8293116961175));
     B.add_vertex( Point( 0.4874571845315,    0.4874571845316,  0.6709272557930));
     end2 = B.add_vertex( Point(-0.0000000000000,    -1.0000000000000,  0.0000000000000)); //End Vertex
     B.add_vertex( Point( 0.7887222926324,    0.4874571845315, -0.2562714077342));
     neight2 = B.add_vertex( Point( 0.4874571845316,    -0.4874571845315,  0.6709272557930));
     B.add_vertex( Point( 0.7887222926324,    -0.4874571845315, -0.2562714077342));
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 11);
     B.add_vertex_to_facet( 9);
     B.add_vertex_to_facet( 3);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 7);
     B.add_vertex_to_facet( 10);
     B.add_vertex_to_facet( 8);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 7);
     B.add_vertex_to_facet( 9);
     B.add_vertex_to_facet( 11);
     B.add_vertex_to_facet( 10);
     B.end_facet();    
     
     
     B.begin_facet();
     B.add_vertex_to_facet( 3);
     B.add_vertex_to_facet( 9);
     B.add_vertex_to_facet( 7);
     B.end_facet();
     
     B.begin_facet();
     B.add_vertex_to_facet( 8);
     B.add_vertex_to_facet( 10);
     B.add_vertex_to_facet( 11);
     B.end_facet();
     
     B.end_surface();
   }
 };
 

/*Fifth Segment*/
template <class HDS>
class buildtwinned5 : public CGAL::Modifier_base<HDS> {
public:
  Vertex_handle end1;
  Vertex_handle neight1;
  Vertex_handle end2;
  Vertex_handle neight2;
  Vertex_handle neight3;
  
  buildtwinned5() {}
  void operator()( HDS& hds) {
    // Postcondition: `hds' is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( 12, 5, 6);
    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;
    
    
    
    
    
    B.add_vertex( Point(-0.7887222926324,    0.4874571845315, -0.2562714077342));
    B.add_vertex( Point(-0.4874571845316,    0.4874571845315,  0.6709272557930));
    B.add_vertex( Point(-0.7887222926324,   -0.4874571845315, -0.2562714077342)); //End vertex
    end1 = B.add_vertex( Point( 0.0000000000000,    1.0000000000000,  0.0000000000000));
    neight3 = B.add_vertex( Point(-0.4874571845315,    -0.4874571845316,  0.6709272557930));
    neight1 = B.add_vertex( Point(-0.0000000000000,    0.4874571845316, -0.8293116961175));
    B.add_vertex( Point( 0.0000000000000,    -0.4874571845316, -0.8293116961175));
    B.add_vertex( Point( 0.4874571845315,    0.4874571845316,  0.6709272557930));
    end2 = B.add_vertex( Point(-0.0000000000000,    -1.0000000000000,  0.0000000000000)); //End Vertex
    B.add_vertex( Point( 0.7887222926324,    0.4874571845315, -0.2562714077342));
    neight2 = B.add_vertex( Point( 0.4874571845316,    -0.4874571845315,  0.6709272557930));
    B.add_vertex( Point( 0.7887222926324,    -0.4874571845315, -0.2562714077342));
    
    B.begin_facet();
    B.add_vertex_to_facet( 8);
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 3);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 10);
    B.add_vertex_to_facet( 8);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 10);
    B.end_facet();    
    
    
    B.begin_facet();
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 3);
    B.end_facet();
    
    B.begin_facet();
    B.add_vertex_to_facet( 10);
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 8);
    B.end_facet();
    
    B.end_surface();
  }
};






struct Normal_vector {
  template <class Facet>
  typename Facet::Plane_3 operator()( Facet& f) {
    typename Facet::Halfedge_handle h = f.halfedge();
    // Facet::Plane_3 is the normal vector type. We assume the
    // CGAL Kernel here and use its global functions.
    return CGAL::cross_product( 
			       h->next()->vertex()->point() - h->vertex()->point(),
			       h->next()->next()->vertex()->point() - h->next()->vertex()->point());
  }
};


bool GeometryBuilder::isInsidePolyhedron(double x, double y, double z) {
	
  Point_3 point(x,y,z);
  Plane_iterator i;
  Facet_iterator j;
  for ( i =nanoRodPolyhedron.planes_begin(), j = nanoRodPolyhedron.facets_begin(); i != nanoRodPolyhedron.planes_end() && j !=nanoRodPolyhedron.facets_end(); ++i, ++j) {
    Halfedge_facet_circulator k = j->facet_begin();
    
    Vector_3 newVector = point - k->vertex()->point();
    Vector_3 normal = *i;		
    double dot_product = newVector.x() * normal.x() + newVector.y() * normal.y() + newVector.z() * normal.z();
    
    if (dot_product < 0) {
      return false;	
    }
  }
  
  return true;
}


GeometryBuilder::GeometryBuilder(double length,double width) {
  // Create the geometry for nanorod
 buildSingleCrystal<HalfedgeDS> singleCrystalNanorod;
  
  nanoRodPolyhedron.delegate( singleCrystalNanorod);
   
  double y1 = singleCrystalNanorod.end1->point().y() - singleCrystalNanorod.neight1->point().y();
  double y2 = singleCrystalNanorod.end2->point().y() - singleCrystalNanorod.neight2->point().y();
  
  double endDist = sqrt(pow(singleCrystalNanorod.neight2->point().x() - singleCrystalNanorod.neight3->point().x(),2)+
                        pow(singleCrystalNanorod.neight2->point().y() - singleCrystalNanorod.neight3->point().y(),2)+
                        pow(singleCrystalNanorod.neight2->point().z() - singleCrystalNanorod.neight3->point().z(),2));
  
  double endRatio1 = y1/endDist;
  double endRatio2 = y2/endDist;
  
  std::cout << "End dist is " << endDist <<" ratio " << endRatio1 << std::endl;
  
  CGAL::Aff_transformation_3<Kernel> aff_tranformation( width,
							0.0,
							0.0,
							0.0,
							0.0,
							length,
							0.0,
							0.0,
							0.0,
							0.0,
							width,
							0.0);	
  std::transform( nanoRodPolyhedron.points_begin(), nanoRodPolyhedron.points_end(), nanoRodPolyhedron.points_begin(), aff_tranformation);
	
  
  double endDist2 = sqrt(pow(singleCrystalNanorod.neight2->point().x() -singleCrystalNanorod.neight3->point().x(),2)+
                        pow(singleCrystalNanorod.neight2->point().y() - singleCrystalNanorod.neight3->point().y(),2)+
                        pow(singleCrystalNanorod.neight2->point().z() - singleCrystalNanorod.neight3->point().z(),2));
  
  Point_3 point1(singleCrystalNanorod.end1->point().x(), endDist2*endRatio1 + singleCrystalNanorod.neight1->point().y(), singleCrystalNanorod.end1->point().z());
  Point_3 point2(singleCrystalNanorod.end2->point().x(), endDist2*endRatio2 + singleCrystalNanorod.neight2->point().y(), singleCrystalNanorod.end2->point().z());
  singleCrystalNanorod.end1->point() = point1;
  singleCrystalNanorod.end2->point() = point2;
	
  // Construct normal vector for each face.
  std::transform( nanoRodPolyhedron.facets_begin(), nanoRodPolyhedron.facets_end(), nanoRodPolyhedron.planes_begin(),
		  Normal_vector());
} 	



GeometryBuilder::GeometryBuilder(double length,double width, bool twinned) {
  // Create the geometry for nanorod
  
  buildtwinned1<HalfedgeDS> crystal1;  
  buildtwinned2<HalfedgeDS> crystal2;  
  buildtwinned3<HalfedgeDS> crystal3;  
  buildtwinned4<HalfedgeDS> crystal4;  
  buildtwinned5<HalfedgeDS> crystal5;  
  
  
  nanoRodTwinnedPolyhedron1.delegate( crystal1);
  nanoRodTwinnedPolyhedron2.delegate( crystal2);
  nanoRodTwinnedPolyhedron3.delegate( crystal3);
  nanoRodTwinnedPolyhedron4.delegate( crystal4);
  nanoRodTwinnedPolyhedron5.delegate( crystal5);
  
  
  
  
  
  
  double y1 = crystal1.end1->point().y() - crystal1.neight1->point().y();
  double y2 = crystal1.end2->point().y() - crystal1.neight2->point().y();
  
  double endDist = sqrt(pow(crystal1.neight2->point().x() - crystal1.neight3->point().x(),2)+
                        pow(crystal1.neight2->point().y() - crystal1.neight3->point().y(),2)+
                        pow(crystal1.neight2->point().z() - crystal1.neight3->point().z(),2));
  
  double endRatio1 = y1/endDist;
  double endRatio2 = y2/endDist;
  
  std::cout << "End dist is " << endDist <<" ratio " << endRatio1 << std::endl;
  
  CGAL::Aff_transformation_3<Kernel> aff_tranformation( width,
                                                        0.0,
                                                        0.0,
                                                        0.0,
                                                        0.0,
                                                        length,
                                                        0.0,
                                                        0.0,
                                                        0.0,
                                                        0.0,
                                                        width,
                                                        0.0);	
  std::transform(nanoRodTwinnedPolyhedron1.points_begin(), nanoRodTwinnedPolyhedron1.points_end(), nanoRodTwinnedPolyhedron1.points_begin(), aff_tranformation);
	
  
  double endDist2 = sqrt(pow(crystal1.neight2->point().x() - crystal1.neight3->point().x(),2)+
                         pow(crystal1.neight2->point().y() - crystal1.neight3->point().y(),2)+
                         pow(crystal1.neight2->point().z() - crystal1.neight3->point().z(),2));
  
  Point_3 point1(crystal1.end1->point().x(), endDist2*endRatio1 + crystal1.neight1->point().y(), crystal1.end1->point().z());
  Point_3 point2(crystal1.end2->point().x(), endDist2*endRatio2 + crystal1.neight2->point().y(), crystal1.end2->point().z());
  crystal1.end1->point() = point1;
  crystal1.end2->point() = point2;
	
  // Construct normal vector for each face.
  std::transform( nanoRodTwinnedPolyhedron1.facets_begin(), nanoRodTwinnedPolyhedron1.facets_end(), nanoRodTwinnedPolyhedron1.planes_begin(),
                  Normal_vector());
} 	









  void GeometryBuilder::dumpGeometry(const std::string& geomFileName){
     
     std::ofstream newGeomFile;
     
     //create new .md file based on old .md file
     newGeomFile.open(geomFileName.c_str());
     
     // Write polyhedron in Object File Format (OFF).
     CGAL::set_ascii_mode( std::cout);
     newGeomFile << "OFF" << std::endl << nanoRodPolyhedron.size_of_vertices() << ' ' 
        << nanoRodPolyhedron.size_of_facets() << " 0" << std::endl;
     std::copy( nanoRodPolyhedron.points_begin(), nanoRodPolyhedron.points_end(),
                std::ostream_iterator<Point_3>( newGeomFile, "\n"));
     for (  Facet_iterator i = nanoRodPolyhedron.facets_begin(); i != nanoRodPolyhedron.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        newGeomFile << CGAL::circulator_size(j) << ' ';
        do {
           newGeomFile << ' ' << std::distance(nanoRodPolyhedron.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        newGeomFile << std::endl;
     }
     
     newGeomFile.close();
     
     
  
	
  }


