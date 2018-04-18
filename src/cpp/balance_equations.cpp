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

    void compute_internal_force(const int &component, const double (&dNdx)[3], const Vector_9 &cauchy, double &f_i){
        /*!================================
        |    compute_internal_force    |
        ================================

        Compute the internal force given the gradient 
        of the shape function and the cauchy stress on 
        one of the components.

        */
        
        if (component==0){
            f_i = -(dNdx[0]*cauchy[0] + dNdx[1]*cauchy[8] + dNdx[2]*cauchy[7]);
        }
        else if (component==1){
            f_i = -(dNdx[0]*cauchy[5] + dNdx[1]*cauchy[1] + dNdx[2]*cauchy[6]);
        }
        
        else if (component==2){
            f_i = -(dNdx[0]*cauchy[4] + dNdx[1]*cauchy[3] + dNdx[2]*cauchy[2]);
        }
        
        else{
            std::cout << "Error: Component " << component << "not recognized\n";
        }

        return;
    }

    void compute_body_force(const int &component, const double &N, const double &density, const double (&b)[3], double &fb_i){
        /*!============================
        |    compute_body_force    |
        ============================

        Compute the body force given the body force per unit density
        on one of the components.

        */
        
        if (component==0){
            fb_i = N*density*b[0];
        }
        else if (component==1){
            fb_i = N*density*b[1];
        }
        
        else if (component==2){
            fb_i = N*density*b[2];
        }
        
        else{
            std::cout << "Error: Component " << component << "not recognized\n";
        }

        return;
    }

    void compute_kinematic_force(const int &component, const double &N, const double &density, const double (&a)[3], double &fkin_i){
        /*!=================================
        |    compute_kinematic_force    |
        =================================

        Compute the kinimatic force given the shape 
        function, density, and acceleration.

        */
        
        if (component==0){
            fkin_i = -N*density*a[0];
        }
        else if (component==1){
            fkin_i = -N*density*a[1];
        }
        
        else if (component==2){
            fkin_i = -N*density*a[2];
        }
        
        else{
            std::cout << "Error: Component " << component << "not recognized\n";
        }

        return;
    }

    void compute_internal_couple(const int component_i, const int component_j, const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double &cint_ij){
        /*!=================================
        |    compute_internal_couple    |
        =================================

        Compute the internal couple at a given 
        component pair.

        */

        if ((component_i == 0) && (component_j == 0)){
            cint_ij = N*(cauchy[0] - s[0]) - (dNdx[0]*m[0]+dNdx[1]*m[9]+dNdx[2]*m[18]);
        }
        else if ((component_i == 1) && (component_j == 1)){
            cint_ij = N*(cauchy[1] - s[1]) - (dNdx[0]*m[1]+dNdx[1]*m[10]+dNdx[2]*m[19]);
        }
        else if ((component_i == 2) && (component_j == 2)){
            cint_ij = N*(cauchy[2] - s[2]) - (dNdx[0]*m[2]+dNdx[1]*m[11]+dNdx[2]*m[20]);
        }
        else if ((component_i == 1) && (component_j == 2)){
            cint_ij = N*(cauchy[3] - s[3]) - (dNdx[0]*m[6]+dNdx[1]*m[15]+dNdx[2]*m[24]);
        }
        else if ((component_i == 0) && (component_j == 2)){
            cint_ij = N*(cauchy[4] - s[4]) - (dNdx[0]*m[7]+dNdx[1]*m[16]+dNdx[2]*m[25]);
        }
        else if ((component_i == 0) && (component_j == 1)){
            cint_ij = N*(cauchy[5] - s[5]) - (dNdx[0]*m[8]+dNdx[1]*m[17]+dNdx[2]*m[26]);
        }
        else if ((component_i == 2) && (component_j == 1)){
            cint_ij = N*(cauchy[6] - s[6]) - (dNdx[0]*m[3]+dNdx[1]*m[12]+dNdx[2]*m[21]);
        }
        else if ((component_i == 2) && (component_j == 0)){
            cint_ij = N*(cauchy[7] - s[7]) - (dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22]);
        }
        else if ((component_i == 1) && (component_j == 0)){
            cint_ij = N*(cauchy[8] - s[8]) - (dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22]);
        }
        else {
            std::cout << "Error: Indices " << component_i << ", " << component_j << " not recognized.";
        }
        
        return;
    }
    
    void compute_body_couple(const int &component_i, const int &component_j, const double &N, const double &density, const double (&l)[9], double &cb_ij){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        if ((component_i == 0) && (component_j == 0)){
            cb_ij = N*density*l[0];
        }
        else if ((component_i == 1) && (component_j == 1)){
            cb_ij = N*density*l[1];;
        }
        else if ((component_i == 2) && (component_j == 2)){
            cb_ij = N*density*l[2];;
        }
        else if ((component_i == 1) && (component_j == 2)){
            cb_ij = N*density*l[6];;
        }
        else if ((component_i == 0) && (component_j == 2)){
            cb_ij = N*density*l[7];;
        }
        else if ((component_i == 0) && (component_j == 1)){
            cb_ij = N*density*l[8];;
        }
        else if ((component_i == 2) && (component_j == 1)){
            cb_ij = N*density*l[3];;
        }
        else if ((component_i == 2) && (component_j == 0)){
            cb_ij = N*density*l[4];;
        }
        else if ((component_i == 1) && (component_j == 0)){
            cb_ij = N*density*l[5];;
        }
        else {
            std::cout << "Error: Indices " << component_i << ", " << component_j << " not recognized.";
        }
        
        return;
    }
    
    void compute_kinematic_couple(const int &component_i, const int &component_j, const double &N, const double &density, const double (&omega)[9], double &ckin_ij){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        if ((component_i == 0) && (component_j == 0)){
            ckin_ij = -N*density*omega[0];
        }
        else if ((component_i == 1) && (component_j == 1)){
            ckin_ij = -N*density*omega[1];
        }
        else if ((component_i == 2) && (component_j == 2)){
            ckin_ij = -N*density*omega[2];
        }
        else if ((component_i == 1) && (component_j == 2)){
            ckin_ij = -N*density*omega[6];
        }
        else if ((component_i == 0) && (component_j == 2)){
            ckin_ij = -N*density*omega[7];
        }
        else if ((component_i == 0) && (component_j == 1)){
            ckin_ij = -N*density*omega[8];
        }
        else if ((component_i == 2) && (component_j == 1)){
            ckin_ij = -N*density*omega[3];
        }
        else if ((component_i == 2) && (component_j == 0)){
            ckin_ij = -N*density*omega[4];
        }
        else if ((component_i == 1) && (component_j == 0)){
            ckin_ij = -N*density*omega[5];
        }
        else {
            std::cout << "Error: Indices " << component_i << ", " << component_j << " not recognized.";
        }
        
        return;
    }
    
    void compute_internal_force_jacobian(const int &component, const int &dof_num, const double &N, const double (&dNdx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, double &dfdU_iK){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the internal force jacobian for a given DOF.
        
        Note that the DOF are expected to be organized
        
        dof_num:  0,  1,  2,     3,     4,     5,     6,     7,     8,     9,    10,    11
            dof: u1, u2, u3, phi11, phi22, phi33, phi23, phi13, phi12, phi32, phi31, phi21
        */
        
        double term1; //Terms associated with DcauchyDgrad_u
        double term2; //Terms associated with DcauchyDphi
        double term3; //Terms associated with DcauchyDgrad_phi
        
        //Assign term1
        if (dof_num<3){

            if (component == 0){
                if( dof_num == 0){
                    term1 =  dNdx[0]*DcauchyDgrad_u(0,0)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(0,5)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(0,4)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(8,0)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(8,5)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(8,4)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(7,0)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(7,5)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(7,4)*dNdx[2];
                }
                else if( dof_num == 1){
                    term1 =  dNdx[0]*DcauchyDgrad_u(0,8)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(0,1)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(0,3)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(8,8)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(8,1)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(8,3)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(7,8)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(7,1)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(7,3)*dNdx[2];
                }
                else if( dof_num == 2){
                    term1 =  dNdx[0]*DcauchyDgrad_u(0,7)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(0,6)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(0,2)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(8,7)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(8,6)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(8,2)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(7,7)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(7,6)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(7,2)*dNdx[2];
                }
            }
            else if (component == 1){
                if( dof_num == 0){
                    term1 =  dNdx[0]*DcauchyDgrad_u(5,0)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(5,5)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(5,4)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(1,0)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(1,5)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(1,4)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(6,0)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(6,5)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(6,4)*dNdx[2];
                }
                else if( dof_num == 1){
                    term1 =  dNdx[0]*DcauchyDgrad_u(5,8)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(5,1)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(5,3)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(1,8)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(1,1)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(1,3)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(6,8)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(6,1)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(6,3)*dNdx[2];
                }
                else if( dof_num == 2){
                    term1 =  dNdx[0]*DcauchyDgrad_u(5,7)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(5,6)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(5,2)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(1,7)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(1,6)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(1,2)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(6,7)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(6,6)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(6,2)*dNdx[2];
                }
            }
            else if (component == 2){
                if( dof_num == 0){
                    term1 =  dNdx[0]*DcauchyDgrad_u(4,0)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(4,5)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(4,4)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(3,0)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(3,5)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(3,4)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(2,0)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(2,5)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(2,4)*dNdx[2];
                }
                else if( dof_num == 1){
                    term1 =  dNdx[0]*DcauchyDgrad_u(4,8)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(4,1)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(4,3)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(3,8)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(3,1)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(3,3)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(2,8)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(2,1)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(2,3)*dNdx[2];
                }
                else if( dof_num == 2){
                    term1 =  dNdx[0]*DcauchyDgrad_u(4,7)*dNdx[0] + dNdx[0]*DcauchyDgrad_u(4,6)*dNdx[1] + dNdx[0]*DcauchyDgrad_u(4,2)*dNdx[2]
                            +dNdx[1]*DcauchyDgrad_u(3,7)*dNdx[0] + dNdx[1]*DcauchyDgrad_u(3,6)*dNdx[1] + dNdx[1]*DcauchyDgrad_u(3,2)*dNdx[2]
                            +dNdx[2]*DcauchyDgrad_u(2,7)*dNdx[0] + dNdx[2]*DcauchyDgrad_u(2,6)*dNdx[1] + dNdx[2]*DcauchyDgrad_u(2,2)*dNdx[2];
                }
            }
            else {
                std::cout << "Error: Component value of" << component << " out of range.";
            }
            
        }
        else {
            term1 = term1 = term1 = 0;
        }
        
        //Assign term2
        if (dof_num<3){
            term2 = term2 = term2 = 0;
        }
        else {
            if (component == 0){
                if (dof_num == 3){term2 = (dNdx[0]*DcauchyDphi(0,0) + dNdx[1]*DcauchyDphi(8,0) + dNdx[2]*DcauchyDphi(7,0))*N;}
                else if (dof_num ==  8){term2 = (dNdx[0]*DcauchyDphi(0,5) + dNdx[1]*DcauchyDphi(8,5) + dNdx[2]*DcauchyDphi(7,5))*N;}
                else if (dof_num ==  7){term2 = (dNdx[0]*DcauchyDphi(0,4) + dNdx[1]*DcauchyDphi(8,4) + dNdx[2]*DcauchyDphi(7,4))*N;}
                else if (dof_num == 11){term2 = (dNdx[0]*DcauchyDphi(0,8) + dNdx[1]*DcauchyDphi(8,8) + dNdx[2]*DcauchyDphi(7,8))*N;}
                else if (dof_num ==  4){term2 = (dNdx[0]*DcauchyDphi(0,1) + dNdx[1]*DcauchyDphi(8,1) + dNdx[2]*DcauchyDphi(7,1))*N;}
                else if (dof_num ==  6){term2 = (dNdx[0]*DcauchyDphi(0,3) + dNdx[1]*DcauchyDphi(8,3) + dNdx[2]*DcauchyDphi(7,3))*N;}
                else if (dof_num == 10){term2 = (dNdx[0]*DcauchyDphi(0,7) + dNdx[1]*DcauchyDphi(8,7) + dNdx[2]*DcauchyDphi(7,7))*N;}
                else if (dof_num ==  9){term2 = (dNdx[0]*DcauchyDphi(0,6) + dNdx[1]*DcauchyDphi(8,6) + dNdx[2]*DcauchyDphi(7,6))*N;}
                else if (dof_num ==  5){term2 = (dNdx[0]*DcauchyDphi(0,2) + dNdx[1]*DcauchyDphi(8,2) + dNdx[2]*DcauchyDphi(7,2))*N;}
            }
            else if (component == 1){
                if (dof_num == 3){term2 = (dNdx[0]*DcauchyDphi(5,0) + dNdx[1]*DcauchyDphi(1,0) + dNdx[2]*DcauchyDphi(6,0))*N;}
                else if (dof_num ==  8){term2 = (dNdx[0]*DcauchyDphi(5,5) + dNdx[1]*DcauchyDphi(1,5) + dNdx[2]*DcauchyDphi(6,5))*N;}
                else if (dof_num ==  7){term2 = (dNdx[0]*DcauchyDphi(5,4) + dNdx[1]*DcauchyDphi(1,4) + dNdx[2]*DcauchyDphi(6,4))*N;}
                else if (dof_num == 11){term2 = (dNdx[0]*DcauchyDphi(5,8) + dNdx[1]*DcauchyDphi(1,8) + dNdx[2]*DcauchyDphi(6,8))*N;}
                else if (dof_num ==  4){term2 = (dNdx[0]*DcauchyDphi(5,1) + dNdx[1]*DcauchyDphi(1,1) + dNdx[2]*DcauchyDphi(6,1))*N;}
                else if (dof_num ==  6){term2 = (dNdx[0]*DcauchyDphi(5,3) + dNdx[1]*DcauchyDphi(1,3) + dNdx[2]*DcauchyDphi(6,3))*N;}
                else if (dof_num == 10){term2 = (dNdx[0]*DcauchyDphi(5,7) + dNdx[1]*DcauchyDphi(1,7) + dNdx[2]*DcauchyDphi(6,7))*N;}
                else if (dof_num ==  9){term2 = (dNdx[0]*DcauchyDphi(5,6) + dNdx[1]*DcauchyDphi(1,6) + dNdx[2]*DcauchyDphi(6,6))*N;}
                else if (dof_num ==  5){term2 = (dNdx[0]*DcauchyDphi(5,2) + dNdx[1]*DcauchyDphi(1,2) + dNdx[2]*DcauchyDphi(6,2))*N;}
            }
            else if (component == 2){
                if (dof_num == 3){term2 = (dNdx[0]*DcauchyDphi(4,0) + dNdx[1]*DcauchyDphi(3,0) + dNdx[2]*DcauchyDphi(2,0))*N;}
                else if (dof_num ==  8){term2 = (dNdx[0]*DcauchyDphi(4,5) + dNdx[1]*DcauchyDphi(3,5) + dNdx[2]*DcauchyDphi(2,5))*N;}
                else if (dof_num ==  7){term2 = (dNdx[0]*DcauchyDphi(4,4) + dNdx[1]*DcauchyDphi(3,4) + dNdx[2]*DcauchyDphi(2,4))*N;}
                else if (dof_num == 11){term2 = (dNdx[0]*DcauchyDphi(4,8) + dNdx[1]*DcauchyDphi(3,8) + dNdx[2]*DcauchyDphi(2,8))*N;}
                else if (dof_num ==  4){term2 = (dNdx[0]*DcauchyDphi(4,1) + dNdx[1]*DcauchyDphi(3,1) + dNdx[2]*DcauchyDphi(2,1))*N;}
                else if (dof_num ==  6){term2 = (dNdx[0]*DcauchyDphi(4,3) + dNdx[1]*DcauchyDphi(3,3) + dNdx[2]*DcauchyDphi(2,3))*N;}
                else if (dof_num == 10){term2 = (dNdx[0]*DcauchyDphi(4,7) + dNdx[1]*DcauchyDphi(3,7) + dNdx[2]*DcauchyDphi(2,7))*N;}
                else if (dof_num ==  9){term2 = (dNdx[0]*DcauchyDphi(4,6) + dNdx[1]*DcauchyDphi(3,6) + dNdx[2]*DcauchyDphi(2,6))*N;}
                else if (dof_num ==  5){term2 = (dNdx[0]*DcauchyDphi(4,2) + dNdx[1]*DcauchyDphi(3,2) + dNdx[2]*DcauchyDphi(2,2))*N;}
            }
            else {
                std::cout << "Error: Component value of" << component << " out of range.";
            }
        }
        
        //Assign term3
        if (dof_num<3){
            term3 = term3 = term3 = 0;
        }
        else {
            if (component == 0){
                if (dof_num ==  3){term3 =  dNdx[0]*DcauchyDgrad_phi(0,0)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,5)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,4)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,0)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,5)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,4)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,0)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,5)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,4)*dNdx[2];
                }
                else if (dof_num ==  8){term3 =  dNdx[0]*DcauchyDgrad_phi(0,8)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,1)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,3)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,8)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,1)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,3)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,8)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,1)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,3)*dNdx[2];
                }
                else if (dof_num ==  7){term3 =  dNdx[0]*DcauchyDgrad_phi(0,7)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,6)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,2)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,7)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,6)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,2)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,7)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,6)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,2)*dNdx[2];
                }
                else if (dof_num == 11){term3 =  dNdx[0]*DcauchyDgrad_phi(0,9)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,14)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,13)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,9)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,14)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,13)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,9)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,14)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,13)*dNdx[2];
                }
                else if (dof_num ==  4){term3 =  dNdx[0]*DcauchyDgrad_phi(0,17)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,10)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,12)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,17)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,10)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,12)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,17)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,10)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,12)*dNdx[2];
                }
                else if (dof_num ==  6){term3 =  dNdx[0]*DcauchyDgrad_phi(0,16)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,15)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,11)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,16)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,15)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,11)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,16)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,15)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,11)*dNdx[2];
                }
                else if (dof_num == 10){term3 =  dNdx[0]*DcauchyDgrad_phi(0,18)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,23)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,22)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,18)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,23)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,22)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,18)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,23)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,22)*dNdx[2];
                }
                else if (dof_num ==  9){term3 =  dNdx[0]*DcauchyDgrad_phi(0,26)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,19)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,21)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,26)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,19)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,21)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,26)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,19)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,21)*dNdx[2];
                }
                else if (dof_num ==  5){term3 =  dNdx[0]*DcauchyDgrad_phi(0,25)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(0,24)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(0,20)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(8,25)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(8,24)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(8,20)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(7,25)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(7,24)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(7,20)*dNdx[2];
                }
            }
            else if (component == 1){
                if (dof_num ==  3){term3 =  dNdx[0]*DcauchyDgrad_phi(5,0)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,5)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,4)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,0)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,5)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,4)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,0)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,5)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,4)*dNdx[2];
                }
                else if (dof_num ==  8){term3 =  dNdx[0]*DcauchyDgrad_phi(5,8)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,1)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,3)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,8)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,1)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,3)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,8)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,1)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,3)*dNdx[2];
                }
                else if (dof_num ==  7){term3 =  dNdx[0]*DcauchyDgrad_phi(5,7)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,6)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,2)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,7)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,6)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,2)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,7)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,6)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,2)*dNdx[2];
                }
                else if (dof_num == 11){term3 =  dNdx[0]*DcauchyDgrad_phi(5,9)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,14)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,13)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,9)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,14)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,13)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,9)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,14)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,13)*dNdx[2];
                }
                else if (dof_num ==  4){term3 =  dNdx[0]*DcauchyDgrad_phi(5,17)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,10)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,12)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,17)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,10)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,12)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,17)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,10)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,12)*dNdx[2];
                }
                else if (dof_num ==  6){term3 =  dNdx[0]*DcauchyDgrad_phi(5,16)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,15)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,11)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,16)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,15)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,11)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,16)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,15)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,11)*dNdx[2];
                }
                else if (dof_num == 10){term3 =  dNdx[0]*DcauchyDgrad_phi(5,18)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,23)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,22)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,18)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,23)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,22)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,18)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,23)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,22)*dNdx[2];
                }
                else if (dof_num ==  9){term3 =  dNdx[0]*DcauchyDgrad_phi(5,26)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,19)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,21)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,26)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,19)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,21)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,26)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,19)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,21)*dNdx[2];
                }
                else if (dof_num ==  5){term3 =  dNdx[0]*DcauchyDgrad_phi(5,25)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(5,24)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(5,20)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(1,25)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(1,24)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(1,20)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(6,25)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(6,24)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(6,20)*dNdx[2];
                }
            }
            else if (component == 2){
                if (dof_num ==  3){term3 =  dNdx[0]*DcauchyDgrad_phi(4,0)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,5)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,4)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,0)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,5)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,4)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,0)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,5)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,4)*dNdx[2];
                }
                else if (dof_num ==  8){term3 =  dNdx[0]*DcauchyDgrad_phi(4,8)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,1)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,3)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,8)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,1)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,3)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,8)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,1)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,3)*dNdx[2];
                }
                else if (dof_num ==  7){term3 =  dNdx[0]*DcauchyDgrad_phi(4,7)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,6)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,2)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,7)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,6)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,2)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,7)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,6)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,2)*dNdx[2];
                }
                else if (dof_num == 11){term3 =  dNdx[0]*DcauchyDgrad_phi(4,9)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,14)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,13)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,9)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,14)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,13)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,9)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,14)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,13)*dNdx[2];
                }
                else if (dof_num ==  4){term3 =  dNdx[0]*DcauchyDgrad_phi(4,17)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,10)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,12)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,17)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,10)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,12)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,17)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,10)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,12)*dNdx[2];
                }
                else if (dof_num ==  6){term3 =  dNdx[0]*DcauchyDgrad_phi(4,16)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,15)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,11)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,16)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,15)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,11)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,16)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,15)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,11)*dNdx[2];
                }
                else if (dof_num == 10){term3 =  dNdx[0]*DcauchyDgrad_phi(4,18)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,23)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,22)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,18)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,23)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,22)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,18)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,23)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,22)*dNdx[2];
                }
                else if (dof_num ==  9){term3 =  dNdx[0]*DcauchyDgrad_phi(4,26)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,19)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,21)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,26)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,19)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,21)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,26)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,19)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,21)*dNdx[2];
                }
                else if (dof_num ==  5){term3 =  dNdx[0]*DcauchyDgrad_phi(4,25)*dNdx[0] + dNdx[0]*DcauchyDgrad_phi(4,24)*dNdx[1] + dNdx[0]*DcauchyDgrad_phi(4,20)*dNdx[2]
                                              +dNdx[1]*DcauchyDgrad_phi(3,25)*dNdx[0] + dNdx[1]*DcauchyDgrad_phi(3,24)*dNdx[1] + dNdx[1]*DcauchyDgrad_phi(3,20)*dNdx[2]
                                              +dNdx[2]*DcauchyDgrad_phi(2,25)*dNdx[0] + dNdx[2]*DcauchyDgrad_phi(2,24)*dNdx[1] + dNdx[2]*DcauchyDgrad_phi(2,20)*dNdx[2];
                }
            }
        }
        
        //Assemble the jacobian
        dfdU_iK = -term1;// - term2 - term3;}
    }
}
