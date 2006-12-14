/* Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
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
 *  ConvexHull.cpp
 *
 *  Purpose: To calculate convexhull, hull volume and radius
 *  using the CGAL library.
 *
 *  Created by Charles F. Vardeman II on 11 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: ConvexHull.cpp,v 1.1 2006-12-14 19:32:32 chuckv Exp $
 *
 */

#include "math/ConvexHull.hpp"
#include <iostream>
#include <fstream>


using namespace oopse;




bool ConvexHull::genHull(std::vector<Vector3d> pos)
{

  std::vector<Point_3> points;
  points.reserve(pos.size());
  // Copy the positon vector into a points vector for cgal.
  for (int i = 0; i < pos.size(); ++i)
    {
      Point_3 pt(pos[i][0],pos[i][1],pos[i][2]);
      points.push_back(pt);
    }

  // compute convex hull
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);
  // Make sure hull is a polyhedron
  if ( CGAL::assign(ch_polyhedron, ch_object) )
    {
      return true;
    }
  else
    {
      return false;
    }
}

RealType ConvexHull::getVolume()
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

void ConvexHull::geomviewHull(const std::string& geomFileName)
{

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


}
