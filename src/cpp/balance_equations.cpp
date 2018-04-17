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
#include <balance_equations.h>

namespace balance_equations{

    void compute_internal_force(const double (&dNdx)[3], const Vector_9 &cauchy, double (&fint)[3]){
        /*!================================
        |    compute_internal_force    |
        ================================

        Compute the internal force given the gradient 
        of the shape function and the cauchy stress.

        */

        fint[0] = dNdx[0]*cauchy[0] + dNdx[1]*cauchy[8] + dNdx[2]*cauchy[7];
        fint[1] = dNdx[0]*cauchy[5] + dNdx[1]*cauchy[1] + dNdx[2]*cauchy[6];
        fint[2] = dNdx[0]*cauchy[4] + dNdx[1]*cauchy[3] + dNdx[2]*cauchy[2];

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

    void compute_internal_couple(const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double (&cint)[9]){
        /*!=================================
        |    compute_internal_couple    |
        =================================

        Compute the internal couple.

        */

        //Add the contributions of the cauchy stress and symmetric stress
        for (int i=0; i<9; i++){
            cint[i] = N*(cauchy[i]-s[i]);
        }
        
        //Add the contributions of the gradient of the higher order couple stress
        cint[0] -= dNdx[0]*m[0]+dNdx[1]*m[9]+dNdx[2]*m[18];
        cint[5] -= dNdx[0]*m[8]+dNdx[1]*m[17]+dNdx[2]*m[26];
        cint[4] -= dNdx[0]*m[7]+dNdx[1]*m[16]+dNdx[2]*m[25];
        cint[8] -= dNdx[0]*m[5]+dNdx[1]*m[14]+dNdx[2]*m[23];
        cint[1] -= dNdx[0]*m[1]+dNdx[1]*m[10]+dNdx[2]*m[19];
        cint[3] -= dNdx[0]*m[6]+dNdx[1]*m[15]+dNdx[2]*m[24];
        cint[7] -= dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22];
        cint[6] -= dNdx[0]*m[3]+dNdx[1]*m[12]+dNdx[2]*m[21];
        cint[2] -= dNdx[0]*m[2]+dNdx[1]*m[11]+dNdx[2]*m[20];
        
        return;
    }
    
    void compute_body_couple(const double &N, const double &density, const double (&l)[9], double (&cb)[9]){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        cb[0] = N*density*l[0];
        cb[1] = N*density*l[1];
        cb[2] = N*density*l[2];
        cb[3] = N*density*l[6];
        cb[4] = N*density*l[7];
        cb[5] = N*density*l[8];
        cb[6] = N*density*l[3];
        cb[7] = N*density*l[4];
        cb[8] = N*density*l[5];
        
        return;
    }
    
    void compute_kinematic_couple(const double &N, const double &density, const double (&omega)[9], double (&ckin)[9]){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        cb[0] = N*density*omega[0];
        cb[1] = N*density*omega[1];
        cb[2] = N*density*omega[2];
        cb[3] = N*density*omega[6];
        cb[4] = N*density*omega[7];
        cb[5] = N*density*omega[8];
        cb[6] = N*density*omega[3];
        cb[7] = N*density*omega[4];
        cb[8] = N*density*omega[5];
        
        return;
    }
    
    void compute_internal_force_jacobian(const int &dof_num, const double &N, const double (&dNdx)[3], const Matrix_9x9 &DcauchyDF, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, const Matrix_27x9 &dFdgrad_u){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the internal force jacobian for a given DOF.
        
        Note that the DOF are expected to be organized
        
        dof_num:  0,  1,  2,     3,     4,     5,     6,     7,     8,     9,    10,    11
            dof: u1, u2, u3, phi11, phi22, phi33, phi23, phi13, phi12, phi32, phi31, phi21
            
        
        
        */
        
        //!Total derivatives

    }
}
