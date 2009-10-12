/* Copyright (c) 2008, 2009 The University of Notre Dame. All Rights Reserved.
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
 *  @version $Id: ConvexHull.cpp,v 1.14 2009-10-12 20:11:29 chuckv Exp $
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
//#include <CGAL/Homogeneous.h>
#include <CGAL/basic.h>
//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Origin.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/double.h>
#include <CGAL/number_utils.h>


//#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
//#include <CGAL/Lazy_exact_nt.h>



typedef CGAL::MP_Float RT;
//typedef double RT;
//typedef CGAL::Homogeneous<RT>                     K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Vector_3                               Vector_3;
//typedef CGAL::Convex_hull_traits_3<K>             Traits;
typedef CGAL::Polyhedron_traits_with_normals_3<K> Traits;
//typedef Traits::Polyhedron_3                      Polyhedron_3;
typedef CGAL::Polyhedron_3<Traits>                     Polyhedron_3;
typedef K::Point_3                                Point_3;


typedef Polyhedron_3::HalfedgeDS             HalfedgeDS;
typedef Polyhedron_3::Facet_iterator                   Facet_iterator;
typedef Polyhedron_3::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef Polyhedron_3::Halfedge_handle Halfedge_handle;
typedef Polyhedron_3::Facet_iterator Facet_iterator;
typedef Polyhedron_3::Plane_iterator Plane_iterator;
typedef Polyhedron_3::Vertex_iterator Vertex_iterator;
typedef Polyhedron_3::Vertex_handle Vertex_handle;
typedef Polyhedron_3::Point_iterator Point_iterator;
 


class Enriched_Point_3 : public K::Point_3{
public:
  Enriched_Point_3(double x,double y,double z) : K::Point_3(x,y,z), yupMyPoint(false), mySD(NULL) {}

  bool isMyPoint() const{ return yupMyPoint; }
  void myPoint(){ yupMyPoint = true; }
  void setSD(StuntDouble* SD){mySD = SD;}
  StuntDouble* getStuntDouble(){return mySD;}
private:
  bool yupMyPoint;
  StuntDouble* mySD;

};





    // compare Point_3's... used in setting up the STL map from points to indices
template <typename Pt3>
struct Point_3_comp {
  bool operator() (const Pt3 & p, const Pt3 & q) const {
    return CGAL::lexicographically_xyz_smaller(p,q); // this is defined inline & hence we had to create fn object & not ptrfun
  }
};

// coordinate-based hashing inefficient but can we do better if pts are copied?
typedef std::map<Point_3, StuntDouble* ,Point_3_comp<Point_3> > ptMapType;

#ifdef IS_MPI
struct {
  double x,y,z;
} surfacePt;
#endif

ConvexHull::ConvexHull() : Hull(){
  //If we are doing the mpi version, set up some vectors for data communication
#ifdef IS_MPI


 nproc_ = MPI::COMM_WORLD.Get_size();
 myrank_ = MPI::COMM_WORLD.Get_rank();
 NstoProc_ = new int[nproc_];
 vecdispls_   = new int[nproc_];
 displs_ = new int[nproc_];
 // Create a surface point type in MPI to send
 surfacePtType = MPI::DOUBLE.Create_contiguous(3);
 surfacePtType.Commit();
 

#endif
}

void ConvexHull::computeHull(std::vector<StuntDouble*> bodydoubles)
{
 
  std::vector<Enriched_Point_3> points;
  ptMapType myMap;
  Point_iterator   hc;
  
  // Copy the positon vector into a points vector for cgal.
  std::vector<StuntDouble*>::iterator SD;

    for (SD =bodydoubles.begin(); SD != bodydoubles.end(); ++SD)
    {
      Vector3d pos = (*SD)->getPos();
      Enriched_Point_3* pt = new Enriched_Point_3(pos.x(),pos.y(),pos.z());
      pt->setSD(*SD);     
      points.push_back(*pt);
      // myMap[pt]=(*SD);
    }
  
  // define object to hold convex hull
  CGAL::Object ch_object_;
  Polyhedron_3 polyhedron;

  // compute convex hull
  
  std::vector<Enriched_Point_3>::iterator testpt;
  
  
 
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron);
 

 
  Ns_ = polyhedron.size_of_vertices();

#ifdef IS_MPI
  /* Gather an array of the number of verticies on each processor */
  

  surfacePtsGlobal_.clear();
  surfacePtsLocal_.clear();

  MPI::COMM_WORLD.Allgather(&Ns_,1,MPI::INT,&NstoProc_[0],1,MPI::INT);

  for (int i = 0; i < nproc_; i++){
    Nsglobal_ += NstoProc_[i];
  }
  /*Reminder ideally, we would like to reserve size for the vectors here*/
  surfacePtsLocal_.reserve(Ns_);
  surfacePtsGlobal_.resize(Nsglobal_);
  //  std::fill(surfacePtsGlobal_.begin(),surfacePtsGlobal_.end(),0);

  /* Build a displacements array */
  for (int i = 1; i < nproc_; i++){
    vecdispls_[i] = vecdispls_[i-1] + NstoProc_[i-1];
  }
  
  int noffset = vecdispls_[myrank_];
  /* gather the potential hull */
  
  
  for (hc =polyhedron.points_begin();hc != polyhedron.points_end(); ++hc){
    Point_3 mypoint = *hc;
    surfacePt_ mpiSurfacePt;
    mpiSurfacePt.x = CGAL::to_double(mypoint.x());
    mpiSurfacePt.y = CGAL::to_double(mypoint.y());
    mpiSurfacePt.z = CGAL::to_double(mypoint.z());
    surfacePtsLocal_.push_back(mpiSurfacePt);
  }

  MPI::COMM_WORLD.Allgatherv(&surfacePtsLocal_[0],Ns_,surfacePtType,&surfacePtsGlobal_[0],NstoProc_,vecdispls_,surfacePtType);
  std::vector<surfacePt_>::iterator spt;
  std::vector<Enriched_Point_3> gblpoints;

  int mine = 0;
  int pointidx = 0;
  for (spt = surfacePtsGlobal_.begin(); spt != surfacePtsGlobal_.end(); ++spt)
    {     
      surfacePt_ thispos = *spt;
      Enriched_Point_3 ept(thispos.x,thispos.y,thispos.z);
      if (mine >= noffset && mine < noffset + Ns_){
	ept.myPoint();
	ept.setSD(points[pointidx].getStuntDouble());
	pointidx++;
      }
      gblpoints.push_back(ept);

      mine++;
    }

  /* Compute the global hull */
  polyhedron.clear();
  CGAL::convex_hull_3(gblpoints.begin(), gblpoints.end(), polyhedron);


#endif


  
  /* Loop over all of the surface triangles and build data structures for atoms and normals*/
  Facet_iterator j;
  area_ = 0;
  for ( j = polyhedron.facets_begin(); j !=polyhedron.facets_end(); ++j) {
    Halfedge_handle h = j->halfedge();

    Point_3 r0=h->vertex()->point();
    Point_3 r1=h->next()->vertex()->point();
    Point_3 r2=h->next()->next()->vertex()->point();

    Point_3* pr0 = &r0;
    Point_3* pr1 = &r1;
    Point_3* pr2 = &r2;

    Enriched_Point_3* er0 = static_cast<Enriched_Point_3*>(pr0);
    Enriched_Point_3* er1 = static_cast<Enriched_Point_3*>(pr1);
    Enriched_Point_3* er2 = static_cast<Enriched_Point_3*>(pr2);

    // StuntDouble* sd = er0->getStuntDouble();
    std::cerr << "sd globalIndex = " << to_double(er0->x()) << "\n";
   
    Point_3 thisCentroid = CGAL::centroid(r0,r1,r2);

    Vector_3 normal = CGAL::cross_product(r1-r0,r2-r0);

    Triangle* face = new Triangle();
    Vector3d V3dNormal(CGAL::to_double(normal.x()),CGAL::to_double(normal.y()),CGAL::to_double(normal.z()));
    Vector3d V3dCentroid(CGAL::to_double(thisCentroid.x()),CGAL::to_double(thisCentroid.y()),CGAL::to_double(thisCentroid.z()));
    face->setNormal(V3dNormal);
    face->setCentroid(V3dCentroid);
    RealType faceArea = 0.5*V3dNormal.length();
    face->setArea(faceArea);
    area_ += faceArea;
    Triangles_.push_back(face);
    //    ptMapType::const_iterator locn=myMap.find(mypoint);
    //    int myIndex = locn->second;

  }
  
  

 
}
void ConvexHull::printHull(const std::string& geomFileName)
{
  /*
  std::ofstream newGeomFile;
  
  //create new .md file based on old .md file
  newGeomFile.open("testhull.off");
  
  // Write polyhedron in Object File Format (OFF).
  CGAL::set_ascii_mode( std::cout);
  newGeomFile << "OFF" << std::endl << polyhedron.size_of_vertices() << ' '
	      << polyhedron.size_of_facets() << " 0" << std::endl;
  std::copy( polyhedron.points_begin(), polyhedron.points_end(),
	     std::ostream_iterator<Point_3>( newGeomFile, "\n"));
  for (  Facet_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i) {
    Halfedge_facet_circulator j = i->facet_begin();
    // Facets in polyhedral surfaces are at least triangles.
    CGAL_assertion( CGAL::circulator_size(j) >= 3);
    newGeomFile << CGAL::circulator_size(j) << ' ';
    do {
      newGeomFile << ' ' << std::distance(polyhedron.vertices_begin(), j->vertex());
    } while ( ++j != i->facet_begin());
    newGeomFile << std::endl;
  }
  
  newGeomFile.close();
  */
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







#else
#ifdef HAVE_QHULL
/* Old options Qt Qu Qg QG0 FA */
/* More old opts Qc Qi Pp*/
ConvexHull::ConvexHull() : Hull(), dim_(3), options_("qhull Qt Pp"), Ns_(200), nTriangles_(0) {
  //If we are doing the mpi version, set up some vectors for data communication
#ifdef IS_MPI


 nproc_ = MPI::COMM_WORLD.Get_size();
 myrank_ = MPI::COMM_WORLD.Get_rank();
 NstoProc_ = new int[nproc_];
 vecdispls_   = new int[nproc_];
 vecNstoProc_ = new int[nproc_];
 displs_ = new int[nproc_];

 // Create a surface point type in MPI to send
 //surfacePtType = MPI::DOUBLE.Create_contiguous(3);
 // surfacePtType.Commit();
 

#endif
}



void ConvexHull::computeHull(std::vector<StuntDouble*> bodydoubles)
{
  
  std::vector<int> surfaceIDs;
  std::vector<int> surfaceIDsGlobal;
  std::vector<int> localPtsMap;
  int numpoints = bodydoubles.size();

  //coordT* pt_array;
  coordT* surfpt_array;
  vertexT *vertex, **vertexp;
  facetT *facet;
  setT *vertices;
  int curlong,totlong;
  int id;
  
  coordT *point,**pointp;


  FILE *outdummy = NULL;
  FILE *errdummy = NULL;
  
  //pt_array = (coordT*) malloc(sizeof(coordT) * (numpoints * dim_));

//  double* ptArray = new double[numpoints * 3];
  std::vector<double> ptArray(numpoints*3);
  std::vector<bool> isSurfaceID(numpoints); 

  // Copy the positon vector into a points vector for qhull.
  std::vector<StuntDouble*>::iterator SD;
  int i = 0;
  for (SD =bodydoubles.begin(); SD != bodydoubles.end(); ++SD)
    {
      Vector3d pos = (*SD)->getPos();
      
      ptArray[dim_ * i] = pos.x();
      ptArray[dim_ * i + 1] = pos.y();
      ptArray[dim_ * i + 2] = pos.z();
      i++;
    }
  

  
  
  
  
  boolT ismalloc = False;
  /* Clean up memory from previous convex hull calculations*/
  
  Triangles_.clear();
  surfaceSDs_.clear();
  surfaceSDs_.reserve(Ns_);

  if (qh_new_qhull(dim_, numpoints, &ptArray[0], ismalloc,
   		    const_cast<char *>(options_.c_str()), NULL, stderr)) {

      sprintf(painCave.errMsg, "ConvexHull: Qhull failed to compute convex hull");
      painCave.isFatal = 1;
      simError();
      
  } //qh_new_qhull


#ifdef IS_MPI
  std::vector<double> localPts;
  std::vector<double> localVel;
  std::vector<double> localMass;
  int localPtArraySize;
  
 
  std::fill(isSurfaceID.begin(),isSurfaceID.end(),false);
 

  FORALLfacets {
    
    if (!facet->simplicial){
      // should never happen with Qt
      sprintf(painCave.errMsg, "ConvexHull: non-simplicaial facet detected");
      painCave.isFatal = 1;
      simError();
    }
    
    
    vertices = qh_facet3vertex(facet);
    FOREACHvertex_(vertices){
      id = qh_pointid(vertex->point);

      if( !isSurfaceID[id] ){
	isSurfaceID[id] = true;
      }
    }      
    qh_settempfree(&vertices);      
      
  } //FORALLfacets

 


  int idx = 0;
  int nIsIts = 0;
  FORALLvertices {
    idx = qh_pointid(vertex->point);
    localPts.push_back(ptArray[dim_ * idx]);     
    localPts.push_back(ptArray[dim_ * idx + 1]); 
    localPts.push_back(ptArray[dim_ * idx + 2]);

    Vector3d vel = bodydoubles[idx]->getVel();
    localVel.push_back(vel.x());
    localVel.push_back(vel.y());
    localVel.push_back(vel.z());


    RealType bdmass = bodydoubles[idx]->getMass();
    localMass.push_back(bdmass);

    localPtsMap.push_back(idx); 

  }


  localPtArraySize = int(localPts.size()/3.0);
 
  MPI::COMM_WORLD.Allgather(&localPtArraySize,1,MPI::INT,&NstoProc_[0],1,MPI::INT);
  
  Nsglobal_=0;
  for (int i = 0; i < nproc_; i++){
    Nsglobal_ += NstoProc_[i];
    vecNstoProc_[i] = NstoProc_[i]*3;
  }
  
 
  int nglobalPts = Nsglobal_*3;
 

  std::vector<double> globalPts(nglobalPts);
  std::vector<double> globalVel(nglobalPts);
  std::vector<double> globalMass(Nsglobal_);


  
  isSurfaceID.resize(nglobalPts);


  std::fill(globalPts.begin(),globalPts.end(),0.0);
 
  vecdispls_[0] = 0;
  /* Build a displacements array */
  for (int i = 1; i < nproc_; i++){
    vecdispls_[i] = vecdispls_[i-1] + vecNstoProc_[i-1];
  }
  
  displs_[0] = 0;
  for (int i = 1; i < nproc_; i++){
    displs_[i] = displs_[i-1] + NstoProc_[i-1];
  }
   
  int noffset = vecdispls_[myrank_];
  /* gather the potential hull */
  
  MPI::COMM_WORLD.Allgatherv(&localPts[0],localPtArraySize*3,MPI::DOUBLE,&globalPts[0],&vecNstoProc_[0],&vecdispls_[0],MPI::DOUBLE);
  MPI::COMM_WORLD.Allgatherv(&localVel[0],localPtArraySize*3,MPI::DOUBLE,&globalVel[0],&vecNstoProc_[0],&vecdispls_[0],MPI::DOUBLE);
  MPI::COMM_WORLD.Allgatherv(&localMass[0],localPtArraySize,MPI::DOUBLE,&globalMass[0],&NstoProc_[0],&displs_[0],MPI::DOUBLE);

  /*
  int tmpidx = 0;
  
  if (myrank_ == 0){
    for (i = 0; i < nglobalPts-3; i++){      
      std::cout << "Au   " << globalPts[tmpidx] << "  " << globalPts[tmpidx+1] << "  " << globalPts[tmpidx +2] << std::endl;
      tmpidx = tmpidx + 3;
    }
  }
  */
  
  // Free previous hull
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    std::cerr << "qhull internal warning (main): did not free %d bytes of long memory (%d pieces) "
	      << totlong << curlong << std::endl;

  if (qh_new_qhull(dim_, Nsglobal_, &globalPts[0], ismalloc,
   		    const_cast<char *>(options_.c_str()), NULL, stderr)){

      sprintf(painCave.errMsg, "ConvexHull: Qhull failed to compute global convex hull");
      painCave.isFatal = 1;
      simError();
      
  } //qh_new_qhull

#endif






    unsigned int nf = qh num_facets;
     
    /* Build Surface SD list first */

    std::fill(isSurfaceID.begin(),isSurfaceID.end(),false);

    FORALLfacets {
      
      if (!facet->simplicial){
      // should never happen with Qt
	sprintf(painCave.errMsg, "ConvexHull: non-simplicaial facet detected");
	painCave.isFatal = 1;
	simError();
      } //simplicical
      
      Triangle face;
      Vector3d  V3dNormal(facet->normal[0],facet->normal[1],facet->normal[2]);
      face.setNormal(V3dNormal);
 
      

      //RealType faceArea = 0.5*V3dNormal.length();
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
	id = qh_pointid(vertex->point);
	p[ver][0] = vertex->point[0];
	p[ver][1] = vertex->point[1];
	p[ver][2] = vertex->point[2];
	int localindex = id;
#ifdef IS_MPI
	Vector3d velVector(globalVel[dim_ * id],globalVel[dim_ * id + 1], globalVel[dim_ * id + 1]);
	
	faceVel = faceVel + velVector;
	faceMass = faceMass + globalMass[id];
	if (id >= noffset/3 && id < (noffset + localPtArraySize)/3 ){
	  localindex = localPtsMap[id-noffset/3];
#else
	  faceVel = faceVel + bodydoubles[localindex]->getVel();
	  faceMass = faceMass + bodydoubles[localindex]->getMass();
#endif
	  face.addVertexSD(bodydoubles[localindex]);
	  if( !isSurfaceID[id] ){
	    isSurfaceID[id] = true;
#ifdef IS_MPI	    
	    
#endif
	    
	    surfaceSDs_.push_back(bodydoubles[localindex]);
	    
	  } //IF isSurfaceID

#ifdef IS_MPI
	 
	}else{
	  face.addVertexSD(NULL);
	  }
#endif
	ver++;
      } //Foreachvertex
      /*
      if (!SETempty_(facet->coplanarset)){
	FOREACHpoint_(facet->coplanarset){
	  id = qh_pointid(point);
	  surfaceSDs_.push_back(bodydoubles[id]);
	}
      }
      */
      face.addVertices(p[0],p[1],p[2]);
      face.setFacetMass(faceMass);
      face.setFacetVelocity(faceVel/3.0);
      Triangles_.push_back(face);
      qh_settempfree(&vertices);      

    } //FORALLfacets

    /*    
    std::cout << surfaceSDs_.size() << std::endl;
    for (SD = surfaceSDs_.begin(); SD != surfaceSDs_.end(); ++SD){
      Vector3d thisatom = (*SD)->getPos();
      std::cout << "Au " << thisatom.x() << "  " << thisatom.y() << " " << thisatom.z() << std::endl;
    }
    */



    Ns_ = surfaceSDs_.size();
    nTriangles_ = Triangles_.size();
    
    qh_getarea(qh facet_list);
    volume_ = qh totvol;
    area_ = qh totarea;
    
    
    
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
    if (curlong || totlong)
      std::cerr << "qhull internal warning (main): did not free %d bytes of long memory (%d pieces) "
		<< totlong << curlong << std::endl;
    
    
    
}



void ConvexHull::printHull(const std::string& geomFileName)
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



