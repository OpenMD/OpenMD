#include "math/RMSD.hpp"
#include "math/SVD.hpp"

using namespace OpenMD;
using namespace JAMA;

RealType RMSD::calculate_rmsd(std::vector<Vector3d> mov,
                              Vector3d mov_com,
                              Vector3d mov_to_ref) {
  
  assert(mov.size() == ref_.size());
  int n;
  int n_vec = ref_.size();

  /* calculate the centre of mass */
  mov_com = V3Zero;
  
  for (n=0; n < n_vec; n++) {
    mov_com += mov[n];
  }
  
  mov_com /= (RealType)n_vec;
  
  mov_to_ref = ref_com - mov_com;
  
  /* shift mov to center of mass */

  for (n=0; n < n_vec; n++) {
    mov[n] -= mov_com;
  }
  
  /* initialize */
  Mat3x3d R(0.0);
  RealType E0 = 0.0;
  
  for (n=0; n < n_vec; n++) {
    
    /* 
     * E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n) 
     */
    E0 += dot(mov[n], mov[n]) + dot(ref_[n], ref_[n]);
    
    /*
     * correlation matrix R:   
     *   R(i,j) = sum(over n): y(n,i) * x(n,j)  
     *   where x(n) and y(n) are two vector sets   
     */

    R += outProduct(mov[n], ref_[n]);

  }
  E0 *= 0.5;

  // Pack everything into Dynamic arrays for the SVD:
  DynamicRectMatrix<RealType> Rtmp(3,3,0.0);
  DynamicRectMatrix<RealType> vtmp(3, 3);
  DynamicVector<RealType> stmp(3);
  DynamicRectMatrix<RealType> wtmp(3, 3);

  Rtmp.setSubMatrix(0, 0, R);
  SVD<RealType> svd(Rtmp);

  svd.getU(vtmp);
  svd.getSingularValues(stmp);
  svd.getV(wtmp);

  Mat3x3d v;
  Vector3d s;
  Mat3x3d w;

  vtmp.getSubMatrix(0,0,v);
  stmp.getSubVector(0,s);
  wtmp.getSubMatrix(0,0,w);
    
  int is_reflection = (v.determinant() * w.determinant()) < RealType(0.0);
  if (is_reflection)
    s(2) = -s(2);

  RealType rmsd_sq = (E0 - 2.0 * s.sum() )/ (RealType)n_vec;
  rmsd_sq = max(rmsd_sq, RealType(0.0));
  RealType rmsd = sqrt(rmsd_sq);
  return rmsd;
}
                            
RotMat3x3d RMSD::optimal_superposition(std::vector<Vector3d> mov,
				       Vector3d mov_com,
				       Vector3d mov_to_ref) {
    
  assert(mov.size() == ref_.size());
  int n;
  int n_vec = ref_.size();

  /* calculate the centre of mass */
  mov_com = V3Zero;
  
  for (n=0; n < n_vec; n++) {
    mov_com += mov[n];
  }
  
  mov_com /= (RealType)n_vec;
  
  mov_to_ref = ref_com - mov_com;
  
  /* shift mov to center of mass */

  for (n=0; n < n_vec; n++) {
    mov[n] -= mov_com;
  }
  
  /* initialize */
  Mat3x3d R(0.0);
  RealType E0 = 0.0;
  
  for (n=0; n < n_vec; n++) {
    
    /*
     * correlation matrix R:   
     *   R(i,j) = sum(over n): y(n,i) * x(n,j)  
     *   where x(n) and y(n) are two vector sets   
     */

    R += outProduct(mov[n], ref_[n]);
  }

  // Pack everything into Dynamic arrays for the SVD:
  DynamicRectMatrix<RealType> Rtmp(3,3,0.0);
  DynamicRectMatrix<RealType> vtmp(3, 3);
  DynamicVector<RealType> stmp(3);
  DynamicRectMatrix<RealType> wtmp(3, 3);

  Rtmp.setSubMatrix(0, 0, R);
  SVD<RealType> svd(Rtmp);

  svd.getU(vtmp);
  svd.getSingularValues(stmp);
  svd.getV(wtmp);

  Mat3x3d v;
  Vector3d s;
  Mat3x3d w;

  vtmp.getSubMatrix(0,0,v);
  stmp.getSubVector(0,s);
  wtmp.getSubMatrix(0,0,w);
    
  int is_reflection = (v.determinant() * w.determinant()) < 0.0;
  if (is_reflection)
    s(2) = -s(2);

  RotMat3x3d U = v*w;
  return U;
}

