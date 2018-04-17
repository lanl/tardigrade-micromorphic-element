/*!============================================================
  |                                                           |
  |                   balance_equations.cpp                   |
  |                                                           |
  -------------------------------------------------------------
  | The source file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deformation_measures.h>

namespace balance_equations{

    void compute_internal_force(const double (&dNdx)[3], const double (&cauchy)[9], double (&fint)[3]){
        /*!================================
        |    compute_internal_force    |
        ================================

        Compute the internal force given the gradient 
        of the shape function and the cauchy stress.

        */

        fint[0] = dNdx[0]*cauchy[1] + dNdx[1]*cauchy[8] + dNdx[2]*cauchy[7];
        find[1] = dNdx[0]*cauchy[5] + dNdx[1]*cauchy[1] + dNdx[2]*cauchy[6];
        find[2] = dNdx[0]*cauchy[4] + dNdx[1]*cauchy[3] + dNdx[2]*cauchy[2];

        return;
    }

    void compute_body_force(const double &N, const double &density, const (&b)[3], double (&fb)[3]){
        /*!============================
        |    compute_body_force    |
        ============================

        Compute the body force given the body force per unit density.

        */

        fb[0] = N*density*b[0];
        fb[1] = N*density*b[1];
        fb[2] = N*density*b[2];

        return;
    }

    void compute_kinematic_force(const double &N, const double &density, const (&a)[3], double (&fkin)[3]){
        /*!=================================
        |    compute_kinematic_force    |
        =================================

        Compute the kinimatic force given the shape 
        function, density, and acceleration.

        */

        fkin[0] = N*density*a[0];
        fkin[1] = N*density*a[1];
        fkin[2] = N*density*a[2];

        return;
    }

    void compute_internal_stress(const double &N, const double (&dNdx)[3], const double (&cauchy)[9], const double (&s)[9], const double (&m)[27], double (&sint)[9]){
        /*!=================================
        |    compute_internal_stress    |
        =================================

        Compute the internal stress

        */

        return;
    }

}
