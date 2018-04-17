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
    void compute_internal_force(const double (&dNdx)[3], const Vector_9 &cauchy, double (&fint)[3]);

    void compute_body_force(const double &N, const double &density, const (&b)[3], double (&fb)[3]);

    void compute_kinematic_force(const double &N, const double &density, const (&a)[3], double (&fkin)[3]);

    //Stresses from the balance of first moment of momentum
    void compute_internal_couple(const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double (&cint)[9]);

    void compute_body_couple(const double &N, const double &density, const double (&l)[9], double (&cb)[9]);
    
    void compute_kinematic_couple(const double &N, const double &density, const double (&omega)[9], double (&ckin)[9]);
    
    //Compute the jacobians of the balance of linear momentum w.r.t. the indicated displacement dof
    void compute_internal_force_jacobian(const int &dof_num, const double &N, const double (&dNdx)[3], const Matrix_9x9 &dcauchydF, const Matrix_9x9 &dcauchydchi, const Matrix_9x27 &dcauchydgrad_chi, const Matrix_27x9 &dgrad_chidF, double (&dfintdU)[3,12]);
    
    
    //Compute the jacobians of the balance of first moment of momentum
}

#endif
