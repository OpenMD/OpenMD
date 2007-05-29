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
 *  @version $Id: ConvexHull.cpp,v 1.2 2007-05-29 22:50:14 chuckv Exp $
 *
 */

#include "math/ConvexHull.hpp"
#include <iostream>
#include <fstream>

char options[] = "qhull Qt FA"; 
int dim_ = 3;

using namespace oopse;

ConvexHull::ConvexHull(){}


bool ConvexHull::genHull(std::vector<Vector3d> pos)
{
	FILE *outfile = stdout;
	FILE *errfile = stderr;
	facetT *facet;
	int exitcode;
	boolT ismalloc = False;
	int curlong,totlong;
	
	int numpoints = pos.size();
	
	coordT points[numpoints][dim_];
	
	for (int i=0; i<numpoints; i++)
	{
     points[i][0] = pos[i][0];
     points[i][1] = pos[i][1];
     points[i][2] = pos[i][2];		
	}
   


	qh_initflags (options);
    qh_init_B (points[0], numpoints, dim_, ismalloc);
    qh_qhull();
    qh_check_output();
							
							
							
    qh_getarea(qh facet_list);
    volume_ = qh totvol;							
	area_ = qh totarea;						
							
							
							
    qh_freeqhull(!qh_ALL);
	qh_memfreeshort (&curlong, &totlong);
	if (curlong || totlong)
		fprintf (errfile, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n",
				 totlong, curlong);
	

	return true;
}


RealType ConvexHull::getVolume()
{
	return volume_;
}

void ConvexHull::geomviewHull(const std::string& geomFileName)
{

  std::ofstream newGeomFile;

  //create new .md file based on old .md file
  newGeomFile.open(geomFileName.c_str());


  newGeomFile.close();


}
