/* Copyright (c) 2008, 2009, 2010 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 *
 *
 *  AlphaHull.cpp
 *
 *  Purpose: To calculate Alpha hull, hull volume libqhull.
 *
 *  Created by Charles F. Vardeman II on 11 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id$
 *
 */

/* Standard includes independent of library */

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <iterator>
#include <utility>
#include "math/AlphaHull.hpp"
#include "utils/simError.h"

#ifdef IS_MPI
#include <mpi.h>
#endif

using namespace OpenMD;

#ifdef HAVE_QHULL
extern "C"
{
#include <qhull/libqhull.h>
#include <qhull/mem.h>
#include <qhull/qset.h>
#include <qhull/geom.h>
#include <qhull/merge.h>
#include <qhull/poly.h>
#include <qhull/io.h>
#include <qhull/stat.h>
}
double calculate_circumradius(pointT* p0,pointT* p1,pointT* p2, int dim);

AlphaHull::AlphaHull(double alpha) : Hull(), dim_(3), alpha_(alpha), options_("qhull d QJ Tcv Pp") {
}

void AlphaHull::computeHull(std::vector<StuntDouble*> bodydoubles) { 
 
  int numpoints = bodydoubles.size();
  bool alphashape=true;
  
  Triangles_.clear();
  
  vertexT *vertex, **vertexp;
  facetT *facet, *neighbor;
  setT *vertices, *verticestop, *verticesbottom;
  int curlong, totlong;
  pointT *interiorPoint;
  
  std::vector<double> ptArray(numpoints*dim_);

  // Copy the positon vector into a points vector for qhull.
  std::vector<StuntDouble*>::iterator SD;
  int i = 0;
  for (SD =bodydoubles.begin(); SD != bodydoubles.end(); ++SD){
    Vector3d pos = (*SD)->getPos();      
    ptArray[dim_ * i] = pos.x();
    ptArray[dim_ * i + 1] = pos.y();
    ptArray[dim_ * i + 2] = pos.z();
    i++;
  }

    /* Clean up memory from previous convex hull calculations*/
  boolT ismalloc = False;

  int ridgesCount=0;
  if (qh_new_qhull(dim_, numpoints, &ptArray[0], ismalloc,
                   const_cast<char *>(options_.c_str()), NULL, stderr)) {

    sprintf(painCave.errMsg, "AlphaHull: Qhull failed to compute convex hull");
    painCave.isFatal = 1;
    simError();
    
  } //qh_new_qhull


#ifdef IS_MPI
  //If we are doing the mpi version, set up some vectors for data communication
  
  int nproc = MPI::COMM_WORLD.Get_size();
  int myrank = MPI::COMM_WORLD.Get_rank();
  int localHullSites = 0;

  std::vector<int> hullSitesOnProc(nproc, 0);
  std::vector<int> coordsOnProc(nproc, 0);
  std::vector<int> displacements(nproc, 0);
  std::vector<int> vectorDisplacements(nproc, 0);

  std::vector<double> coords;
  std::vector<double> vels;
  std::vector<int> indexMap;
  std::vector<double> masses;

  FORALLvertices{
    localHullSites++;
    
    int idx = qh_pointid(vertex->point);

    indexMap.push_back(idx);

    coords.push_back(ptArray[dim_  * idx]);
    coords.push_back(ptArray[dim_  * idx + 1]);
    coords.push_back(ptArray[dim_  * idx + 2]);

    StuntDouble* sd = bodydoubles[idx];

    Vector3d vel = sd->getVel();
    vels.push_back(vel.x());
    vels.push_back(vel.y());
    vels.push_back(vel.z());

    masses.push_back(sd->getMass());
  }

  MPI::COMM_WORLD.Allgather(&localHullSites, 1, MPI::INT, &hullSitesOnProc[0],
                            1, MPI::INT);

  int globalHullSites = 0;
  for (int iproc = 0; iproc < nproc; iproc++){
    globalHullSites += hullSitesOnProc[iproc];
    coordsOnProc[iproc] = dim_ * hullSitesOnProc[iproc];
  }

  displacements[0] = 0;
  vectorDisplacements[0] = 0;
  
  for (int iproc = 1; iproc < nproc; iproc++){
    displacements[iproc] = displacements[iproc-1] + hullSitesOnProc[iproc-1];
    vectorDisplacements[iproc] = vectorDisplacements[iproc-1] + coordsOnProc[iproc-1]; 
  }

  std::vector<double> globalCoords(dim_ * globalHullSites);
  std::vector<double> globalVels(dim_ * globalHullSites);
  std::vector<double> globalMasses(globalHullSites);

  int count = coordsOnProc[myrank];
  
  MPI::COMM_WORLD.Allgatherv(&coords[0], count, MPI::DOUBLE, &globalCoords[0],
                             &coordsOnProc[0], &vectorDisplacements[0], 
                             MPI::DOUBLE);

  MPI::COMM_WORLD.Allgatherv(&vels[0], count, MPI::DOUBLE, &globalVels[0], 
                             &coordsOnProc[0], &vectorDisplacements[0],
                             MPI::DOUBLE);

  MPI::COMM_WORLD.Allgatherv(&masses[0], localHullSites, MPI::DOUBLE,
                             &globalMasses[0], &hullSitesOnProc[0], 
                             &displacements[0], MPI::DOUBLE);

  // Free previous hull
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    std::cerr << "qhull internal warning (main): did not free %d bytes of long memory (%d pieces) "
	      << totlong << curlong << std::endl;
  
  if (qh_new_qhull(dim_, globalHullSites, &globalCoords[0], ismalloc,
                   const_cast<char *>(options_.c_str()), NULL, stderr)){
    
    sprintf(painCave.errMsg, "AlphaHull: Qhull failed to compute global convex hull");
    painCave.isFatal = 1;
    simError();
    
  } //qh_new_qhull


#endif

  //Set facet->center as the Voronoi center
  qh_setvoronoi_all();
  
  
  int convexNumVert = qh_setsize(qh_facetvertices (qh facet_list, NULL, false));
  //Insert all the sample points, because, even with alpha=0, the alpha shape/alpha complex will
  //contain them.

  //  tri::Allocator<CMeshO>::AddVertices(pm.cm,convexNumVert);
  
  /*ivp length is 'qh num_vertices' because each vertex is accessed through its ID whose range is 
    0<=qh_pointid(vertex->point)<qh num_vertices*/
  //  vector<tri::Allocator<CMeshO>::VertexPointer> ivp(qh num_vertices);
  /*i=0;
  FORALLvertices{	
    if ((*vertex).point){
      //  pm.cm.vert[i].P()[0] = (*vertex).point[0];
      // pm.cm.vert[i].P()[1] = (*vertex).point[1];
      //pm.cm.vert[i].P()[2] = (*vertex).point[2];
      // ivp[qh_pointid(vertex->point)] = &pm.cm.vert[i];
      i++;
    }
  }
  */
  //Set of alpha complex triangles for alphashape filtering
  setT* set= qh_settemp(4* qh num_facets); 
  
  qh visit_id++;
  int numFacets=0;
  std::vector<std::vector <int> > facetlist;
  interiorPoint = qh interior_point;
  FORALLfacet_(qh facet_list) {
    numFacets++;
    if (!facet->upperdelaunay) {
      //For all facets (that are tetrahedrons)calculate the radius of the empty circumsphere considering 
      //the distance between the circumcenter and a vertex of the facet
      vertexT* vertex = (vertexT *)(facet->vertices->e[0].p);
      double* center = facet->center;
      double radius =  qh_pointdist(vertex->point,center,dim_);
      
      if (radius>alpha_) // if the facet is not good consider the ridges
        {
          //if calculating the alphashape, unmark the facet ('good' is used as 'marked'). 
          facet->good=false;
          
          //Compute each ridge (triangle) once and test the cironference radius with alpha
          facet->visitid= qh visit_id;
          qh_makeridges(facet);
          ridgeT *ridge, **ridgep;
          int goodTriangles=0;
          FOREACHridge_(facet->ridges) {
            neighbor= otherfacet_(ridge, facet);
            if (( neighbor->visitid != qh visit_id)){ 			
              //Calculate the radius of the circumference 
              pointT* p0 = ((vertexT*) (ridge->vertices->e[0].p))->point;
              pointT* p1 = ((vertexT*) (ridge->vertices->e[1].p))->point;
              pointT* p2 = ((vertexT*) (ridge->vertices->e[2].p))->point;
              
              radius = calculate_circumradius(p0,p1,p2, dim_);
              
              if(radius <=alpha_){
                goodTriangles++;
                //save the triangle (ridge) for subsequent filtering
                qh_setappend(&set, ridge); 
              }
            }
          }

          //If calculating the alphashape, mark the facet('good' is used as 'marked'). 
          //This facet will have some triangles hidden by the facet's neighbor.
          if(goodTriangles==4)
            facet->good=true;
          
        }
      else //the facet is good. Put all the triangles of the tetrahedron in the mesh
        {
          //Compute each ridge (triangle) once
          facet->visitid= qh visit_id;
          //If calculating the alphashape, mark the facet('good' is used as 'marked').
          //This facet will have some triangles hidden by the facet's neighbor.
          facet->good=true;
          qh_makeridges(facet);
          ridgeT *ridge, **ridgep;
          FOREACHridge_(facet->ridges) {
            neighbor= otherfacet_(ridge, facet);
            if ((neighbor->visitid != qh visit_id)){
                 qh_setappend(&set, ridge);
            }	
          }
        }
    }
  }
  //assert(numFacets== qh num_facets);
  
  //Filter the triangles (only the ones on the boundary of the alpha complex) and build the mesh


  
  ridgeT *ridge, **ridgep;
  FOREACHridge_(set) {
    if ((!ridge->top->good || !ridge->bottom->good || ridge->top->upperdelaunay || ridge->bottom->upperdelaunay)){
      //        tri::Allocator<CMeshO>::FaceIterator fi=tri::Allocator<CMeshO>::AddFaces(pm.cm,1);
      ridgesCount++;
      int vertex_n, vertex_i;
      Triangle face;
      
      // Vector3d V3dNormal(facet->normal[0], facet->normal[1], facet->normal[2]);
      //face.setNormal(V3dNormal);
      
      
      //coordT *center = qh_getcenter(ridge->vertices);
      //std::cout << "Centers are " << center[0] << "  " <<center[1] << "  " << center[2] << std::endl;
      //Vector3d V3dCentroid(center[0], center[1], center[2]);
      //face.setCentroid(V3dCentroid);

      
      Vector3d faceVel = V3Zero;
      Vector3d p[3];
      RealType faceMass = 0.0;
      
      int ver = 0;
      std::vector<int> virtexlist;
      FOREACHvertex_i_(ridge->vertices){
        int id = qh_pointid(vertex->point);
        p[ver][0] = vertex->point[0];
        p[ver][1] = vertex->point[1];
        p[ver][2] = vertex->point[2];
        Vector3d vel;
        RealType mass;
        ver++;
        virtexlist.push_back(id);
        // std::cout << "Ridge: " << ridgesCount << " Vertex " << id << std::endl; 

        vel = bodydoubles[id]->getVel();
        mass = bodydoubles[id]->getMass();
        face.addVertexSD(bodydoubles[id]);


        faceVel = faceVel + vel;
        faceMass = faceMass + mass;
      } //FOREACH Vertex 
      facetlist.push_back(virtexlist);
      face.addVertices(p[0],p[1],p[2]);
      face.setFacetMass(faceMass);
      face.setFacetVelocity(faceVel/3.0);
      
      RealType area = face.getArea();
      area_ += area;
      Vector3d normal = face.getUnitNormal();
      RealType offset =  ((0.0-p[0][0])*normal[0] + (0.0-p[0][1])*normal[1] + (0.0-p[0][2])*normal[2]);
      RealType dist =  normal[0] * interiorPoint[0] + normal[1]*interiorPoint[1] + normal[2]*interiorPoint[2];
      std::cout << "Dist and normal and area are: " << normal << std::endl;
      volume_ += dist *area/qh hull_dim;
      
      Triangles_.push_back(face);
    }
  }

  std::cout << "Volume is: " << volume_ << std::endl; 

//assert(pm.cm.fn == ridgesCount);
/*
  std::cout <<"OFF"<<std::endl; 
  std::cout << bodydoubles.size() << "  " << facetlist.size() << "  " << 3*facetlist.size() << std::endl; 
  for (SD =bodydoubles.begin(); SD != bodydoubles.end(); ++SD){
    Vector3d pos = (*SD)->getPos();      
    std::cout << pos.x() << "  " << pos.y() << "  " << pos.z() << std::endl; 
  }

  
  std::vector<std::vector<int> >::iterator thisfacet;
  std::vector<int>::iterator thisvertex;

  for (thisfacet = facetlist.begin(); thisfacet != facetlist.end(); thisfacet++){
    std::cout << (*thisfacet).size(); 
    for (thisvertex = (*thisfacet).begin(); thisvertex != (*thisfacet).end(); thisvertex++){
      std::cout << "  " <<  *thisvertex;
    }
    std::cout << std::endl; 
  }
*/



/*  
  FORALLfacets {  
    Triangle face;

    Vector3d V3dNormal(facet->normal[0], facet->normal[1], facet->normal[2]);
    face.setNormal(V3dNormal);
    
    RealType faceArea = qh_facetarea(facet);
    face.setArea(faceArea);
    
    vertices = qh_facet3vertex(facet);
      
    coordT *center = qh_getcenter(vertices);
    Vector3d V3dCentroid(center[0], center[1], center[2]);
    face.setCentroid(V3dCentroid);

    Vector3d faceVel = V3Zero;
    Vector3d p[3];
    RealType faceMass = 0.0;

    int ver = 0;

    FOREACHvertex_(vertices){
      int id = qh_pointid(vertex->point);
      p[ver][0] = vertex->point[0];
      p[ver][1] = vertex->point[1];
      p[ver][2] = vertex->point[2];
      
      Vector3d vel;
      RealType mass;

#ifdef IS_MPI
      vel = Vector3d(globalVels[dim_ * id],
                     globalVels[dim_ * id + 1], 
                     globalVels[dim_ * id + 2]);
      mass = globalMasses[id];

      // localID will be between 0 and hullSitesOnProc[myrank] if we
      // own this guy.

      int localID = id - displacements[myrank];

      if (localID >= 0 && localID < hullSitesOnProc[myrank])
        face.addVertexSD(bodydoubles[indexMap[localID]]);
      
#else
      vel = bodydoubles[id]->getVel();
      mass = bodydoubles[id]->getMass();
      face.addVertexSD(bodydoubles[id]);      
#endif
	
      faceVel = faceVel + vel;
      faceMass = faceMass + mass;
      ver++;      
    } //Foreachvertex

    face.addVertices(p[0], p[1], p[2]);
    face.setFacetMass(faceMass);
    face.setFacetVelocity(faceVel/3.0);
    Triangles_.push_back(face);
    qh_settempfree(&vertices);      

  } //FORALLfacets
*/
  // qh_getarea(qh facet_list);
  //volume_ = qh totvol;
  // area_ = qh totarea;

  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    std::cerr << "qhull internal warning (main): did not free %d bytes of long memory (%d pieces) "
              << totlong << curlong << std::endl;    
}

void AlphaHull::printHull(const std::string& geomFileName) {

#ifdef IS_MPI
  if (worldRank == 0)  {
#endif
  FILE *newGeomFile;
  
  //create new .md file based on old .md file
  newGeomFile = fopen(geomFileName.c_str(), "w");
  qh_findgood_all(qh facet_list);
  for (int i = 0; i < qh_PRINTEND; i++)
    qh_printfacets(newGeomFile, qh PRINTout[i], qh facet_list, NULL, !qh_ALL);
  
  fclose(newGeomFile);
#ifdef IS_MPI
  }
#endif  
}

  double calculate_circumradius(pointT* p0,pointT* p1,pointT* p2, int dim){
    coordT a = qh_pointdist(p0,p1,dim);
    coordT b = qh_pointdist(p1,p2,dim);
    coordT c = qh_pointdist(p2,p0,dim);
    
    coordT sum =(a + b + c)*0.5;
    coordT area = sum*(a+b-sum)*(a+c-sum)*(b+c-sum);
    return (double) (a*b*c)/(4*sqrt(area));
  }

#endif //QHULL
