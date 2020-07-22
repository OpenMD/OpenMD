/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifdef IS_MPI
#include <mpi.h>
#endif

/* Standard includes independent of library */

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <iterator>
#include "math/AlphaHull.hpp"
#include "utils/simError.h"

#include "math/qhull.hpp"

#ifdef HAVE_QHULL
using namespace std;
using namespace OpenMD;

double calculate_circumradius(pointT* p0, pointT* p1, pointT* p2, int dim);

AlphaHull::AlphaHull(double alpha) : Hull(), dim_(4), alpha_(alpha), 
                                     options_("qhull d QJ Tcv Pp") {
}

void AlphaHull::computeHull(vector<StuntDouble*> bodydoubles) { 

#ifdef HAVE_QHULL_REENTRANT
  qhT qh_qh;
  qhT *qh= &qh_qh;
  QHULL_LIB_CHECK
#endif
  
  int numpoints = bodydoubles.size();
  
  Triangles_.clear();
  
  vertexT *vertex;
  facetT *facet, *neighbor;
  pointT *interiorPoint;
  int curlong, totlong;
    
  vector<double> ptArray(numpoints*dim_);
  
  // Copy the positon vector into a points vector for qhull.
  vector<StuntDouble*>::iterator SD;
  int i = 0;

  for (SD =bodydoubles.begin(); SD != bodydoubles.end(); ++SD){
    Vector3d pos = (*SD)->getPos();      
    ptArray[dim_ * i    ] = pos.x();
    ptArray[dim_ * i + 1] = pos.y();
    ptArray[dim_ * i + 2] = pos.z();
    ptArray[dim_ * i + 3] = pos.lengthSquare();
    i++;
  }
  /* Clean up memory from previous convex hull calculations*/
  boolT ismalloc = False;
  
  /* compute the hull for our local points (or all the points for single
     processor versions) */
#ifdef HAVE_QHULL_REENTRANT
  qh_init_A(qh, NULL, NULL, stderr, 0, NULL);
  int exitcode= setjmp(qh->errexit);
  if (!exitcode) {
    qh->NOerrexit = False;
    qh_initflags(qh, const_cast<char *>(options_.c_str()));
    qh_init_B(qh, &ptArray[0], numpoints, dim_, ismalloc);
    qh_qhull(qh);
    qh_check_output(qh);
    exitcode= qh_ERRnone;
    qh->NOerrexit= True;
  } else {
    sprintf(painCave.errMsg, "AlphaHull: Qhull failed to compute convex hull");
    painCave.isFatal = 1;
    simError();
  }
#else
  qh_init_A(NULL, NULL, stderr, 0, NULL);
  int exitcode= setjmp(qh errexit);
  if (!exitcode) {
    qh_initflags(const_cast<char *>(options_.c_str()));
    qh_init_B(&ptArray[0], numpoints, dim_, ismalloc);
    qh_qhull();
    qh_check_output();
    exitcode= qh_ERRnone;
    qh NOerrexit= True;
  } else {
    sprintf(painCave.errMsg, "AlphaHull: Qhull failed to compute convex hull");
    painCave.isFatal = 1;
    simError();
  }
#endif

  
#ifdef IS_MPI
  //If we are doing the mpi version, set up some vectors for data communication
  
  int nproc;
  int myrank;
  MPI_Comm_size( MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);

  int localHullSites = 0;

  vector<int> hullSitesOnProc(nproc, 0);
  vector<int> coordsOnProc(nproc, 0);
  vector<int> displacements(nproc, 0);
  vector<int> vectorDisplacements(nproc, 0);

  vector<double> coords;
  vector<double> vels;
  vector<int> indexMap;
  vector<double> masses;

  FORALLvertices{
    localHullSites++;
    
#ifdef HAVE_QHULL_REENTRANT
    int idx = qh_pointid(qh, vertex->point);
#else
    int idx = qh_pointid(vertex->point);
#endif
    
    indexMap.push_back(idx);

    coords.push_back(ptArray[dim_  * idx]);
    coords.push_back(ptArray[dim_  * idx + 1]);
    coords.push_back(ptArray[dim_  * idx + 2]);
    coords.push_back(ptArray[dim_  * idx + 3]);


    StuntDouble* sd = bodydoubles[idx];

    Vector3d vel = sd->getVel();
    vels.push_back(vel.x());
    vels.push_back(vel.y());
    vels.push_back(vel.z());
    vels.push_back(0.0);

    masses.push_back(sd->getMass());
  }

  MPI_Allgather(&localHullSites, 1, MPI_INT, &hullSitesOnProc[0],
                1, MPI_INT, MPI_COMM_WORLD);

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

  vector<double> globalCoords(dim_ * globalHullSites);
  vector<double> globalVels(dim_ * globalHullSites);
  vector<double> globalMasses(globalHullSites);

  int count = coordsOnProc[myrank];
  
  MPI_Allgatherv(&coords[0], count, MPI_DOUBLE, &globalCoords[0],
                 &coordsOnProc[0], &vectorDisplacements[0], 
                 MPI_DOUBLE, MPI_COMM_WORLD);
  
  MPI_Allgatherv(&vels[0], count, MPI_DOUBLE, &globalVels[0], 
                 &coordsOnProc[0], &vectorDisplacements[0],
                 MPI_DOUBLE, MPI_COMM_WORLD);
  
  MPI_Allgatherv(&masses[0], localHullSites, MPI_DOUBLE,
                 &globalMasses[0], &hullSitesOnProc[0], 
                 &displacements[0], MPI_DOUBLE, MPI_COMM_WORLD);
  
  // Free previous hull
#ifdef HAVE_QHULL_REENTRANT
  qh_freeqhull(qh, !qh_ALL);
  qh_memfreeshort(qh, &curlong, &totlong);
#else
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
#endif
  if (curlong || totlong) {
    sprintf(painCave.errMsg, "AlphaHull: qhull internal warning:\n"
            "\tdid not free %d bytes of long memory (%d pieces)", 
            totlong, curlong);
    painCave.isFatal = 1;
    simError();
  }
  
#ifdef HAVE_QHULL_REENTRANT
  qh_init_A(qh, NULL, NULL, stderr, 0, NULL);
  exitcode= setjmp(qh->errexit);
  if (!exitcode) {
    qh->NOerrexit = False;
    qh_initflags(qh, const_cast<char *>(options_.c_str()));
    qh_init_B(qh, &globalCoords[0], globalHullSites, dim_, ismalloc);
    qh_qhull(qh);
    qh_check_output(qh);
    exitcode= qh_ERRnone;
    qh->NOerrexit= True;
  } else {
    sprintf(painCave.errMsg, "AlphaHull: Qhull failed to compute convex hull");
    painCave.isFatal = 1;
    simError();
  }
#else   
  qh_init_A(NULL, NULL, stderr, 0, NULL);
  exitcode= setjmp(qh errexit);
  if (!exitcode) {
    qh NOerrexit = False;
    qh_initflags(const_cast<char *>(options_.c_str()));
    qh_init_B(&globalCoords[0], globalHullSites, dim_, ismalloc);
    qh_qhull();
    qh_check_output();
    exitcode= qh_ERRnone;
    qh NOerrexit= True;
  } else {
    sprintf(painCave.errMsg, "AlphaHull: Qhull failed to compute convex hull");
    painCave.isFatal = 1;
    simError();
  }
#endif    

#endif

  //Set facet->center as the Voronoi center
#ifdef HAVE_QHULL_REENTRANT
  qh_setvoronoi_all(qh);
#else
  qh_setvoronoi_all();
#endif
  
  //Set of alpha complex triangles for alphashape filtering
  vector<vector <int> > facetlist;
  int numFacets=0;

#ifdef HAVE_QHULL_REENTRANT
  setT* set= qh_settemp(qh, 4* qh->num_facets); 
  qh->visit_id++;
  interiorPoint = qh->interior_point;  
#else
  setT* set= qh_settemp(4* qh num_facets); 
  qh visit_id++;
  interiorPoint = qh interior_point;
#endif

#ifdef HAVE_QHULL_REENTRANT
  FORALLfacet_(qh->facet_list) {
#else  
  FORALLfacet_(qh facet_list) {
#endif
    numFacets++;
    if (!facet->upperdelaunay) {
      //For all facets (that are tetrahedrons)calculate the radius of
      //the empty circumsphere considering the distance between the
      //circumcenter and a vertex of the facet
      vertexT* vertex = (vertexT *)(facet->vertices->e[0].p);

#ifdef HAVE_QHULL_REENTRANT
      double* center = qh_facetcenter(qh, facet->vertices);      
#else
      double* center = qh_facetcenter(facet->vertices);
#endif
      double radius = qh_pointdist(center, vertex->point, dim_-1);      
      //If radius is bigger than alpha, remove the tetrahedron

      if (radius > alpha_) {
          //if calculating the alphashape, unmark the facet ('good' is
          //used as 'marked').
          facet->good=false;
          
          //Compute each ridge (triangle) once and test the
          //cironference radius with alpha
#ifdef HAVE_QHULL_REENTRANT
          facet->visitid= qh->visit_id;
          qh_makeridges(qh, facet);
#else
          facet->visitid= qh visit_id;
          qh_makeridges(facet);
#endif
          ridgeT *ridge, **ridgep;
          int goodTriangles=0;
          FOREACHridge_(facet->ridges) {
            neighbor= otherfacet_(ridge, facet);
#ifdef HAVE_QHULL_REENTRANT
            if (( neighbor->visitid != qh->visit_id)){
#else
            if (( neighbor->visitid != qh visit_id)){
#endif
              //Calculate the radius of the circumference 
              pointT* p0 = ((vertexT*) (ridge->vertices->e[0].p))->point;
              pointT* p1 = ((vertexT*) (ridge->vertices->e[1].p))->point;
              pointT* p2 = ((vertexT*) (ridge->vertices->e[2].p))->point;
              
              radius = calculate_circumradius(p0,p1,p2, dim_-1);
              
              if(radius <=alpha_){
                goodTriangles++;
                //save the triangle (ridge) for subsequent filtering
#ifdef HAVE_QHULL_REENTRANT
                qh_setappend(qh, &set, ridge);
#else
                qh_setappend(&set, ridge);
#endif
              }
            }
          }

          //If calculating the alphashape, mark the facet('good' is
          //used as 'marked').  This facet will have some triangles
          //hidden by the facet's neighbor.
          if(goodTriangles==4)
            facet->good=true;
          
        }
      else //the facet is good. Put all the triangles of the
           //tetrahedron in the mesh
        {
          //Compute each ridge (triangle) once
#ifdef HAVE_QHULL_REENTRANT
          facet->visitid= qh->visit_id;
#else          
          facet->visitid= qh visit_id;
#endif
          //Mark the facet('good' is used as 'marked').  This facet
          //will have some triangles hidden by the facet's neighbor.
          facet->good=true;
#ifdef HAVE_QHULL_REENTRANT
          qh_makeridges(qh, facet);
#else          
          qh_makeridges(facet);
#endif
          ridgeT *ridge, **ridgep;
          FOREACHridge_(facet->ridges) {
            neighbor= otherfacet_(ridge, facet);
#ifdef HAVE_QHULL_REENTRANT
            if ((neighbor->visitid != qh->visit_id)){
              qh_setappend(qh, &set, ridge);
            }
#else            
            if ((neighbor->visitid != qh visit_id)){
              qh_setappend(&set, ridge);
            }
#endif
          }
        }
    }
  }

  int ridgesCount=0;
  
  ridgeT *ridge, **ridgep;
  FOREACHridge_(set) {
    if ((!ridge->top->good || !ridge->bottom->good ||
         ridge->top->upperdelaunay || ridge->bottom->upperdelaunay)){
      ridgesCount++;
      int vertex_n, vertex_i;
      Triangle face;
      
      Vector3d faceVel = V3Zero;
      Vector3d p[3];
      RealType faceMass = 0.0;
      
      int ver = 0;
      vector<int> vertexlist;
      
#ifdef HAVE_QHULL_REENTRANT
      FOREACHvertex_i_(qh, ridge->vertices){
#else
      FOREACHvertex_i_(ridge->vertices){
#endif
#ifdef HAVE_QHULL_REENTRANT
        int id = qh_pointid(qh, vertex->point);
#else
        int id = qh_pointid(vertex->point);
#endif
        p[ver][0] = vertex->point[0];
        p[ver][1] = vertex->point[1];
        p[ver][2] = vertex->point[2];
        Vector3d vel;
        RealType mass;
        ver++;
        vertexlist.push_back(id);

        vel = bodydoubles[id]->getVel();
        mass = bodydoubles[id]->getMass();
        face.addVertexSD(bodydoubles[id]);


        faceVel = faceVel + vel;
        faceMass = faceMass + mass;
      } //FOREACH Vertex 
      facetlist.push_back(vertexlist);
      face.addVertices(p[0],p[1],p[2]);
      face.setFacetMass(faceMass);
      face.setFacetVelocity(faceVel / RealType(3.0));
      
      RealType area = face.getArea();
      area_ += area;
      Vector3d normal = face.getUnitNormal();
      RealType dist =  normal[0] * interiorPoint[0] +
        normal[1]*interiorPoint[1] + normal[2]*interiorPoint[2];
#ifdef HAVE_QHULL_REENTRANT
      volume_ += dist *area/qh->hull_dim;
#else      
      volume_ += dist *area/qh hull_dim;
#endif
      
      Triangles_.push_back(face);
    }
  }


#ifdef HAVE_QHULL_REENTRANT
  qh_freeqhull(qh, !qh_ALL);
  qh_memfreeshort(qh, &curlong, &totlong);  
#else
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
#endif
  if (curlong || totlong) {
    sprintf(painCave.errMsg, "AlphaHull: qhull internal warning:\n"
            "\tdid not free %d bytes of long memory (%d pieces)", 
            totlong, curlong);
    painCave.isFatal = 1;
    simError();
  }
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
