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
 *  Purpose: To calculate convexhull, hull volume libqhull.
 *
 *  Created by Charles F. Vardeman II on 11 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: ConvexHull.cpp,v 1.6 2008-05-14 14:31:48 chuckv Exp $
 *
 */

/* Standard includes independent of library */
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <iterator>
#include "math/ConvexHull.hpp"
#include "utils/simError.h"
using namespace oopse;

/* CGAL version of convex hull first then QHULL */
#ifdef HAVE_CGAL

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/double.h>

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>

typedef double RT;
typedef CGAL::Simple_cartesian<RT>                K;
typedef CGAL::Convex_hull_traits_3<K>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;
typedef K::Point_3                                Point;


ConvexHull::ConvexHull(){}

bool ConvexHull::genHull(std::vector<Vector3d> pos)
{
  
  std::vector<Point> points;	
  
  
  // Copy the positon vector into a points vector for cgal.
  for (int i = 0; i < pos.size(); ++i)
    {
      Point pt(pos[i][0],pos[i][1],pos[i][2]);
      points.push_back(pt);
    }
  
  // define object to hold convex hull
  Polyhedron_3 ch_object_;
  // compute convex hull
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object_);
  
  for (Polyhedron_3::Vertex_iterator v = ch_object_.vertices_begin(); ch_object_.vertices_end(); ++v){
    std::cout<< v.point()<<std::endl;
  }
  
  
}



#else
#ifdef HAVE_QHULL
/* Old options Qt Qu Qg QG0 FA */
ConvexHull::ConvexHull() : dim_(3), options_("qhull Qt  Qci Tcv Pp")
			   //ConvexHull::ConvexHull() : dim_(3), options_("qhull d Qbb Qt i")
{}

bool ConvexHull::genHull(std::vector<Vector3d> pos)
{
  
  
  int numpoints = pos.size();
  coordT* pt_array;
  coordT* surfpt_array;
  std::list<int> surface_atoms;
  FILE *outdummy = NULL;
  FILE *errdummy = NULL;
  
  pt_array = (coordT*) malloc(sizeof(coordT) * (numpoints * dim_));
  
  
  for (int i = 0; i < numpoints; i++) {
    pt_array[dim_ * i] = pos[i][0];
    pt_array[dim_ * i + 1] = pos[i][1];
    pt_array[dim_ * i + 2] = pos[i][2];
  }
  
  
  
  
  /*
    qh_initflags(const_cast<char *>(options_.c_str()));
    qh_init_B(pospoints, numpoints, dim_, ismalloc);
    qh_qhull();
    qh_check_output();

    qh_produce_output();
  */
  boolT ismalloc = False;
  
  if (!qh_new_qhull(dim_, numpoints, pt_array, ismalloc,
		    const_cast<char *>(options_.c_str()), NULL, stderr)) {
    
    vertexT *vertex, **vertexp;
    facetT *facet;
    setT *vertices;
    unsigned int nf = qh num_facets;
    
    //Matrix idx(nf, dim);
    /*
      int j, i = 0, id = 0;
      
      int id2 = 0;
      coordT *point,**pointp;
      realT dist;
      FORALLfacets {
      j = 0;
      
      if (!facet->simplicial){
      // should never happen with Qt
      sprintf(painCave.errMsg, "ConvexHull: non-simplicaial facet detected");
      painCave.isFatal = 0;
      simError();
      }
			
      vertices = qh_facet3vertex(facet);
      FOREACHvertex_(vertices){
      id = qh_pointid(vertex->point);
      surface_atoms.push_back(id);
      //std::cout << "Ag  " << pos[id][0] << "    " << pos[id][1] << "    " << pos[id][2]<< std::endl;
      }
      qh_settempfree(&vertices);
      
      FOREACHpoint_(facet->coplanarset){
      vertex= qh_nearvertex (facet, point, &dist);
      //id= qh_pointid (vertex->point);
      id2= qh_pointid (point);
      surface_atoms.push_back(id2);
      //std::cout << "Ag  " << pos[id2][0] << "    " << pos[id2][1] << "    " << pos[id2][2]<< std::endl;
      //printf ("%d %d %d " qh_REAL_1 "\n", id, id2, facet->id, dist);
      //std::cout << "Neighbors are: %d $d %d\n" << id << id2 << facet->id;
					
      }
      
      }
		
*/
		
		
  }




  qh_getarea(qh facet_list);
  volume_ = qh totvol;
  area_ = qh totarea;
  //	FILE *newGeomFile;
  
  
  /*
    FORALLfacets {
    for (int k=0; k < qh hull_dim; k++)
    printf ("%6.2g ", facet->normal[k]);
    printf ("\n");
    }
  */
  
  int curlong,totlong;
  //	geomviewHull("junk.off");
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    std::cerr << "qhull internal warning (main): did not free %d bytes of long memory (%d pieces) "
	      << totlong << curlong << std::endl;
  free(pt_array);
  /*
    int j = 0;
    surface_atoms.sort();
    surface_atoms.unique();
    surfpt_array = (coordT*) malloc(sizeof(coordT) * (surface_atoms.size() * dim_));
    for(std::list<int>::iterator list_iter = surface_atoms.begin(); 
    list_iter != surface_atoms.end(); list_iter++)
    {
    int i = *list_iter;
    //surfpt_array[dim_ * j] = pos[i][0];
    //surfpt_array[dim_ * j + 1] = pos[i][1];
    //surfpt_array[dim_ * j + 2] = pos[i][2];
    std::cout << "Ag  " << pos[i][0] << "  " << pos[i][1] << "  "<< pos[i][2] << std::endl;
    j++;
    }
  */	
  
  /*	
	std::string deloptions_ = "qhull d Qt";
	facetT *facet, *neighbor;
	ridgeT *ridge, **ridgep;
	
	if (!qh_new_qhull(dim_, surface_atoms.size(), surfpt_array, ismalloc,
	const_cast<char *>(deloptions_.c_str()), NULL, stderr)) {
	
	qh visit_id++;
	FORALLfacets {
	if (!facet->upperdelaunay) {
	facet->visitid= qh visit_id;
	qh_makeridges(facet);
	FOREACHridge_(facet->ridges) {
	neighbor= otherfacet_(ridge, facet);
	if (neighbor->visitid != qh visit_id) {
	
	FOREACHvertex_(ridge->vertices)
	int id2 = qh_pointid (vertex->point); 
	std::cout << "Ag  " << pos[id2][0] << "    " << pos[id2][1] << "    " << pos[id2][2]<< std::endl;
	}
	}
	}
	


	
	}

	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
	if (curlong || totlong)
	std::cerr << "qhull internal warning (main): did not free %d bytes of long memory (%d pieces) "
	<< totlong << curlong << std::endl;
	free(surfpt_array);
  */		
  return true;
}



RealType ConvexHull::getVolume()
{
  return volume_;
}

void ConvexHull::geomviewHull(const std::string& geomFileName)
{

  FILE *newGeomFile;
  
  //create new .md file based on old .md file
  newGeomFile = fopen(geomFileName.c_str(), "w");
  qh_findgood_all(qh facet_list);
  for (int i = 0; i < qh_PRINTEND; i++)
    qh_printfacets(newGeomFile, qh PRINTout[i], qh facet_list, NULL, !qh_ALL);
  
  fclose(newGeomFile);
}
#endif //QHULL
#endif //CGAL



