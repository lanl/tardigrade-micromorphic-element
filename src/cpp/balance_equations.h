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

namespace balance_equations{

    //Forces from the balance of linear momentum
    void compute_internal_force(const double (&dNdx)[3], const double (&cauchy)[9], double (&fint)[3]);

    void compute_body_force(const double &N, const double &density, const (&b)[3], double (&fb)[3]);

    void compute_kinematic_force(const double &N, const double &density, const (&a)[3], double (&fkin)[3]);

    //Stresses from the balance of first moment of momentum
    void compute_internal_stress(const double &N, const double (&dNdx)[3], const double (&cauchy)[9], const double (&s)[9], const double (&m)[27], double (&sint)[9]);

}

#endif
