/*
 * Copyright (c) 2013 The University of Notre Dame. All Rights Reserved.
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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

/*! \file Decahedron.hpp
    \brief Decahedron cluster structure generator
*/

/* Original copyright & license text:

Copyright (c) 2011, Dmitry
Copyright (c) 2009, Richard Brown
Copyright (c) 2011, Evgeny Pr
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the distribution
      
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

*/


#ifndef CLUSTERS_DECAHEDRON_HPP
#define CLUSTERS_DECAHEDRON_HPP

#include <vector>
#include "math/Vector3.hpp"

using namespace std;
namespace OpenMD{

    //! Creates the regular decahedron, Ino decahedron, or truncated
    //! (Marks) decahedron structures (depending on the parameters).
    /*!
        (Heavily modified from Matlab code from: 
        Dmitry, Richard Brown, and Evgeny Pr)

    */
  class Decahedron {
  public:
    //! Default constructor
    Decahedron(int columnAtoms, int shells, int twinAtoms);
    virtual ~Decahedron();
    
    //! Get the generated points in the cluster.
    virtual vector<Vector3d> getPoints();
    
  protected:
    //! Generate the rings of the Decahedron.
    vector<Vector3d> truncatedRing( int n, int k );
    
    int N_; // number of atoms in the central atomic column
    int M_; // number of {002} shells
    int K_; // number of atomic columns along twin boundary

    vector<Vector3d> Points;
    vector<Vector3d> Basis; // Basis vectors of the Decahedron
  };

  class RegularDecahedron : public Decahedron {
  public:
    RegularDecahedron(int shells) : Decahedron(shells-1, shells, shells) {}
  };

  class InoDecahedron : public Decahedron {
  public:
    InoDecahedron(int columnAtoms, int shells) : 
      Decahedron(columnAtoms, shells, columnAtoms) {}
  };

  class CurlingStoneDecahedron : public Decahedron {
  public:
    CurlingStoneDecahedron(int columnAtoms, int shells, int twinAtoms, int truncatedPlanes);
    vector<Vector3d> getPoints();
  private:
    int T_; // number of truncated planes snipped from top and bottom
  };
    
}

#endif
