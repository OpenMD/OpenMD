/* Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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
 *
 *
 *  AlphaShape.cpp
 *
 *  Purpose: To calculate convexhull, hull volume and radius
 *  using the CGAL library.
 *
 *  Created by Charles F. Vardeman II on 11 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: AlphaShape.cpp,v 1.2 2008-05-14 14:31:48 chuckv Exp $
 *
 */

#include "math/AlphaShape.hpp"

#include <iostream>
#include <list>
#include <fstream>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/alpha_shape_geomview_ostream_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Tetrahedron_3.h>


 struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

 typedef CGAL::Alpha_shape_vertex_base_3<K>               Vb;
 typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
 typedef CGAL::Alpha_shape_cell_base_3<K>                 Fb;
 typedef CGAL::Triangulation_data_structure_3<Vbh,Fb>     Tds;
 typedef CGAL::Delaunay_triangulation_3<K,Tds>            Delaunay;
 typedef CGAL::Triangulation_hierarchy_3<Delaunay>        Delaunay_hierarchy;
 typedef CGAL::Alpha_shape_3<Delaunay_hierarchy>          Alpha_shape_3;

 typedef K::Point_3                                  Point;
 typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;
 typedef Alpha_shape_3::NT                           NT;
 typedef Alpha_shape_3::Cell_handle          Cell_handle;
 typedef Alpha_shape_3::Vertex_handle        Vertex_handle;
 typedef Alpha_shape_3::Facet                Facet;
 typedef Alpha_shape_3::Edge                 Edge;
 typedef CGAL::Tetrahedron_3<K> Tetrahedron;




using namespace oopse;

AlphaShape::AlphaShape(){}


bool AlphaShape::genHull(std::vector<Vector3d> pos)
{
  Delaunay_hierarchy dt;
  //points.reserve(pos.size());
  // Copy the positon vector into a points vector for cgal.
  for (unsigned int i = 0; i < pos.size(); ++i)
    {
      Point pt(pos[i][0],pos[i][1],pos[i][2]);
      dt.insert(pt);
    }

  /* Generate Alpha Shape */
  std::cout << "Generating alpha shape" << std::endl;
  Alpha_shape_3 ashape(dt);
  
   /*
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv.set_line_width(4);
  gv.set_trace(false);
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();
  */
    Alpha_shape_3::NT alpha_solid = ashape.find_alpha_solid();
    Alpha_iterator opt = ashape.find_optimal_alpha(1);
    std::cout << "Smallest alpha value to get a solid through data points is "
    	    << alpha_solid << std::endl;
    std::cout << "Optimal alpha value to get one connected component is "
    	    <<  *opt    << std::endl;
    ashape.set_alpha(*opt);
    assert(ashape.number_of_solid_components() == 1);
    std::list<Cell_handle>        cells;
    std::list<Facet>              facets;
    std::list<Edge>               edges;
    std::list<Vertex_handle>      vertices;
    
    ashape.get_alpha_shape_cells(std::back_inserter(cells),
    			  Alpha_shape_3::INTERIOR);
    ashape.get_alpha_shape_facets(std::back_inserter(facets),
    			  Alpha_shape_3::REGULAR);
    ashape.get_alpha_shape_facets(std::back_inserter(facets),
    			  Alpha_shape_3::SINGULAR);
    ashape.get_alpha_shape_edges(std::back_inserter(edges),
    			  Alpha_shape_3::SINGULAR);
    ashape.get_alpha_shape_vertices(std::back_inserter(vertices),
            Alpha_shape_3::REGULAR);			   
    std::cout << " The 0-shape has : " << std::endl;
    std::cout << cells.size() << " interior tetrahedra" << std::endl;
    std::cout << facets.size() << " boundary facets" << std::endl;
    std::cout << edges.size()  << " singular edges" << std::endl;
    std::cout << vertices.size() << " singular vertices" << std::endl;
  
    RealType volume_;
    std::list<Cell_handle>::iterator thiscell;
    
    for(Alpha_shape_3::Cell_iterator c_it = ashape.cells_begin(); c_it != ashape.cells_end(); ++c_it)
      {
        volume_ += ashape.tetrahedron(c_it).volume();
      }
  
    //gv << (Delaunay) ashape;
    //std::cout << ashape;
    
}

RealType AlphaShape::getVolume()
{
	/*
  std::list<Point> L;
  L.push_front(Point(0,0,0));
  L.push_front(Point(1,0,0));
  L.push_front(Point(0,1,0));

  Triangulation T(L.begin(), L.end());

  int n = T.number_of_vertices();

  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0,0,1);
  V[1] = Point(1,1,1);
  V[2] = Point(2,2,2);

  n = n + T.insert(V.begin(), V.end());

  assert( n == 6 );       // 6 points have been inserted
  assert( T.is_valid() ); // checking validity of T

  double sum_v = 0;
  for(Triangulation::Cell_iterator c_it = T.cells_begin(); c_it != T.cells_end(); ++c_it)
    {
      sum_v += T.tetrahedron(c_it).volume();
    }
  std::cout << "sum_v " << sum_v << std::endl;
*/
	return 0.0;
}

void AlphaShape::geomviewHull(const std::string& geomFileName)
{
/*
  std::ofstream newGeomFile;

  //create new .md file based on old .md file
  newGeomFile.open(geomFileName.c_str());

  // Write polyhedron in Object File Format (OFF).
  CGAL::set_ascii_mode( std::cout);
  newGeomFile << "OFF" << std::endl << ch_polyhedron.size_of_vertices() << ' '
  << ch_polyhedron.size_of_facets() << " 0" << std::endl;
  std::copy( ch_polyhedron.points_begin(), ch_polyhedron.points_end(),
             std::ostream_iterator<Point_3>( newGeomFile, "\n"));
  for (  Facet_iterator i = ch_polyhedron.facets_begin(); i != ch_polyhedron.facets_end(); ++i)
    {
      Halfedge_facet_circulator j = i->facet_begin();
      // Facets in polyhedral surfaces are at least triangles.
      CGAL_assertion( CGAL::circulator_size(j) >= 3);
      newGeomFile << CGAL::circulator_size(j) << ' ';
      do
        {
          newGeomFile << ' ' << std::distance(ch_polyhedron.vertices_begin(), j->vertex());
        }
      while ( ++j != i->facet_begin());
      newGeomFile << std::endl;
    }

  newGeomFile.close();
*/

}
