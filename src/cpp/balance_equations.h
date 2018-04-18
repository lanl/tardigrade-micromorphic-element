/*!============================================================
  |                                                           |
  |                    balance_equations.h                    |
  |                                                           |
  -------------------------------------------------------------
  | The header file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#ifndef BALANCE_EQUATIONS_H
#define BALANCE_EQUATIONS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deformation_measures.h>

namespace balance_equations{

    //Forces from the balance of linear momentum
    void compute_internal_force(const int &component, const double (&dNdx)[3], const Vector_9 &cauchy, double &f_i);

    void compute_body_force(const int &component, const double &N, const double &density, const double (&b)[3], double &fb_i);

    void compute_kinematic_force(const int &component, const double &N, const double &density, const double (&a)[3], double &fkin_i);

    //Stresses from the balance of first moment of momentum
    void compute_internal_couple(const int component_i, const int component_j, const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double &cint_ij);

    void compute_body_couple(const int &component_i, const int &component_j, const double &N, const double &density, const double (&l)[9], double &cb_ij);
    
    void compute_kinematic_couple(const int &component_i, const int &component_j, const double &N, const double &density, const double (&omega)[9], double &ckin_ij);
    
    //Compute the jacobians of the balance of linear momentum w.r.t. the indicated displacement dof
    void compute_internal_force_jacobian(const int &component, const int &dof_num, const double &N, const double (&dNdx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, double &dfdU_iK);
    
    
    //Compute the jacobians of the balance of first moment of momentum
}

#endif
