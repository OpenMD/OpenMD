/* Copyright (c) 2010 The University of Notre Dame. All Rights Reserved.
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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). 
 *
 *  ConvexHull.cpp
 *
 *  Purpose: To calculate a convex hull.
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
#include "math/ConvexHull.hpp"
#include "utils/simError.h"

#include "math/qhull.hpp"

#ifdef HAVE_QHULL
using namespace OpenMD;
using namespace std;

ConvexHull::ConvexHull() : Hull(), dim_(3), options_("qhull FA Qt Pp") {
}

void ConvexHull::computeHull(vector<StuntDouble*> bodydoubles) { 
  
  int numpoints = bodydoubles.size();
  
  Triangles_.clear();
  
  vertexT *vertex, **vertexp;
  facetT *facet;
  setT *vertices;
  int curlong, totlong;

  vector<double> ptArray(numpoints*dim_);

  // Copy the positon vector into a points vector for qhull.
  vector<StuntDouble*>::iterator SD;
  int i = 0;
 
  for (SD =bodydoubles.begin(); SD != bodydoubles.end(); ++SD){
    Vector3d pos = (*SD)->getPos();      
    ptArray[dim_ * i] = pos.x();
    ptArray[dim_ * i + 1] = pos.y();
    ptArray[dim_ * i + 2] = pos.z();
    i++;
  }
  
  /* Clean up memory from previous convex hull calculations */
  boolT ismalloc = False;
  
  /* compute the hull for our local points (or all the points for single
     processor versions) */
  if (qh_new_qhull(dim_, numpoints, &ptArray[0], ismalloc,
                   const_cast<char *>(options_.c_str()), NULL, stderr)) {
    
    sprintf(painCave.errMsg, "ConvexHull: Qhull failed to compute convex hull");
    painCave.isFatal = 1;
    simError();
    
  } //qh_new_qhull


#ifdef IS_MPI
  //If we are doing the mpi version, set up some vectors for data communication
  
  int nproc = MPI::COMM_WORLD.Get_size();
  int myrank = MPI::COMM_WORLD.Get_rank();
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

  vector<double> globalCoords(dim_ * globalHullSites);
  vector<double> globalVels(dim_ * globalHullSites);
  vector<double> globalMasses(globalHullSites);

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
  if (curlong || totlong) {
    sprintf(painCave.errMsg, "ConvexHull: qhull internal warning:\n"
            "\tdid not free %d bytes of long memory (%d pieces)", 
            totlong, curlong);
    painCave.isFatal = 1;
    simError();
  }
  
  if (qh_new_qhull(dim_, globalHullSites, &globalCoords[0], ismalloc,
                   const_cast<char *>(options_.c_str()), NULL, stderr)){
    
    sprintf(painCave.errMsg, 
            "ConvexHull: Qhull failed to compute global convex hull");
    painCave.isFatal = 1;
    simError();
    
  } //qh_new_qhull

#endif
  // commented out below, so comment out here also.
  // intPoint = qh interior_point;
  // RealType calcvol = 0.0;
  
  qh_triangulate ();

  FORALLfacets {  
    Triangle face;
    //Qhull sets the unit normal in facet->normal
    Vector3d V3dNormal(facet->normal[0], facet->normal[1], facet->normal[2]);
    face.setUnitNormal(V3dNormal);
    
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


      if (localID >= 0 && localID < hullSitesOnProc[myrank]){
        face.addVertexSD(bodydoubles[indexMap[localID]]);
      }else{
        face.addVertexSD(NULL);
      }
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
    face.setFacetVelocity(faceVel / RealType(3.0));
    /*
    RealType comparea = face.computeArea();
    realT calcarea = qh_facetarea (facet);
    Vector3d V3dCompNorm = -face.computeUnitNormal();
    RealType thisOffset = ((0.0-p[0][0])*V3dCompNorm[0] + (0.0-p[0][1])*V3dCompNorm[1] + (0.0-p[0][2])*V3dCompNorm[2]);
    RealType dist = facet->offset + intPoint[0]*V3dNormal[0] + intPoint[1]*V3dNormal[1] + intPoint[2]*V3dNormal[2];
    cout << "facet offset and computed offset: " << facet->offset << "  " << thisOffset <<  endl;
    calcvol +=  -dist*comparea/qh hull_dim;
    */
    Triangles_.push_back(face);
    qh_settempfree(&vertices);      

  } //FORALLfacets
  
  qh_getarea(qh facet_list);
  volume_ = qh totvol;
  area_ = qh totarea;
     
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong) {
    sprintf(painCave.errMsg, "ConvexHull: qhull internal warning:\n"
            "\tdid not free %d bytes of long memory (%d pieces)", 
            totlong, curlong);
    painCave.isFatal = 1;
    simError();
  }
}

void ConvexHull::printHull(const string& geomFileName) {
  
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
#endif //QHULL
