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
        of the shape function and the cauchy stress on 
        one of the components.

        */
        
        fint[0] = -(dNdx[0]*cauchy[0] + dNdx[1]*cauchy[8] + dNdx[2]*cauchy[7]);
        fint[1] = -(dNdx[0]*cauchy[5] + dNdx[1]*cauchy[1] + dNdx[2]*cauchy[6]);
        fint[2] = -(dNdx[0]*cauchy[4] + dNdx[1]*cauchy[3] + dNdx[2]*cauchy[2]);

        return;
    }
    
    void compute_internal_force(const int &i, const double (&dNdx)[3], const Vector_9 &cauchy, double &fint_i){
        /*!================================
        |    compute_internal_force    |
        ================================

        Compute the internal force given the gradient 
        of the shape function and the cauchy stress on 
        one of the components.

        */
        
        if( i == 0 ){
            fint_i = -(dNdx[0]*cauchy[0] + dNdx[1]*cauchy[8] + dNdx[2]*cauchy[7]);
        }
        else if ( i == 1 ){
            fint_i = -(dNdx[0]*cauchy[5] + dNdx[1]*cauchy[1] + dNdx[2]*cauchy[6]);
        }
        else if ( i == 2 ){
            fint_i = -(dNdx[0]*cauchy[4] + dNdx[1]*cauchy[3] + dNdx[2]*cauchy[2]);
        }
        else{
            std::cout << "Error: index beyond appropriate range.\n";
        }

        return;
    }

    void compute_body_force(const double &N, const double &density, const double (&b)[3], double (&fb)[3]){
        /*!============================
        |    compute_body_force    |
        ============================

        Compute the body force given the body force per unit density
        on one of the components.

        */
        
        for (int i=0; i<3; i++){
            fb[i] = N*density*b[i];
        }

        return;
    }
    
    void compute_body_force(const int &i, const double &N, const double &density, const double (&b)[3], double &fb_i){
        /*!============================
        |    compute_body_force    |
        ============================

        Compute the body force given the body force per unit density
        on one of the components.

        */
        
        fb_i = N*density*b[i];

        return;
    }

    void compute_kinematic_force(const double &N, const double &density, const double (&a)[3], double (&fkin)[3]){
        /*!=================================
        |    compute_kinematic_force    |
        =================================

        Compute the kinimatic force given the shape 
        function, density, and acceleration.

        */
        
        for (int i=0; i<3; i++){
            fkin[i] = -N*density*a[i];
        }

        return;
    }
    
    void compute_kinematic_force(const int &i, const double &N, const double &density, const double (&a)[3], double &fkin_i){
        /*!=================================
        |    compute_kinematic_force    |
        =================================

        Compute the kinimatic force given the shape 
        function, density, and acceleration.

        */
        
        fkin_i = -N*density*a[i];

        return;
    }

    void compute_internal_couple(const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double (&cint)[9]){
        /*!=================================
        |    compute_internal_couple    |
        =================================

        Compute the internal couple for all 
        indices.

        */
        
        cint[0] = N*(cauchy[0] - s[0]) - (dNdx[0]*m[0]+dNdx[1]*m[9]+dNdx[2]*m[18]);
        cint[1] = N*(cauchy[1] - s[1]) - (dNdx[0]*m[1]+dNdx[1]*m[10]+dNdx[2]*m[19]);
        cint[2] = N*(cauchy[2] - s[2]) - (dNdx[0]*m[2]+dNdx[1]*m[11]+dNdx[2]*m[20]);
        cint[3] = N*(cauchy[3] - s[3]) - (dNdx[0]*m[6]+dNdx[1]*m[15]+dNdx[2]*m[24]);
        cint[4] = N*(cauchy[4] - s[4]) - (dNdx[0]*m[7]+dNdx[1]*m[16]+dNdx[2]*m[25]);
        cint[5] = N*(cauchy[5] - s[5]) - (dNdx[0]*m[8]+dNdx[1]*m[17]+dNdx[2]*m[26]);
        cint[6] = N*(cauchy[6] - s[6]) - (dNdx[0]*m[3]+dNdx[1]*m[12]+dNdx[2]*m[21]);
        cint[7] = N*(cauchy[7] - s[7]) - (dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22]);
        cint[8] = N*(cauchy[8] - s[8]) - (dNdx[0]*m[5]+dNdx[1]*m[14]+dNdx[2]*m[23]);
        
        return;
    }
    
    void compute_internal_couple(const int &i, const int &j, const double &N, const double (&dNdx)[3], const Vector_9 &cauchy, const Vector_9 &s, const Vector_27 &m, double &cint_ij){
        /*!=================================
        |    compute_internal_couple    |
        =================================

        Compute the internal couple at a given 
        component pair.

        */
        
        if ( ( i == 0 ) && ( j == 0 ) ){
            cint_ij = N*(cauchy[0] - s[0]) - (dNdx[0]*m[0]+dNdx[1]*m[ 9]+dNdx[2]*m[18]);
        }
        else if ( ( i == 1 ) && ( j == 1 ) ){
            cint_ij = N*(cauchy[1] - s[1]) - (dNdx[0]*m[1]+dNdx[1]*m[10]+dNdx[2]*m[19]);
        }
        else if ( ( i == 2 ) && ( j == 2 ) ){
            cint_ij = N*(cauchy[2] - s[2]) - (dNdx[0]*m[2]+dNdx[1]*m[11]+dNdx[2]*m[20]);
        }
        else if ( ( i == 1 ) && ( j == 2 ) ){
            cint_ij = N*(cauchy[3] - s[3]) - (dNdx[0]*m[6]+dNdx[1]*m[15]+dNdx[2]*m[24]);
        }
        else if ( ( i == 0 ) && ( j == 2 ) ){
            cint_ij = N*(cauchy[4] - s[4]) - (dNdx[0]*m[7]+dNdx[1]*m[16]+dNdx[2]*m[25]);
        }
        else if ( ( i == 0 ) && ( j == 1 ) ){
            cint_ij = N*(cauchy[5] - s[5]) - (dNdx[0]*m[8]+dNdx[1]*m[17]+dNdx[2]*m[26]);
        }
        else if ( ( i == 2 ) && ( j == 1 ) ){
            cint_ij = N*(cauchy[6] - s[6]) - (dNdx[0]*m[3]+dNdx[1]*m[12]+dNdx[2]*m[21]);
        }
        else if ( ( i == 2 ) && ( j == 0 ) ){
            cint_ij = N*(cauchy[7] - s[7]) - (dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22]);
        }
        else if ( ( i == 1 ) && ( j == 0 ) ){
            cint_ij = N*(cauchy[8] - s[8]) - (dNdx[0]*m[5]+dNdx[1]*m[14]+dNdx[2]*m[23]);
        }
        else{
            std::cout << "Error: Index out of range\n";
        }
        
        return;
    }
    
    void compute_body_couple(const double &N, const double &density, const double (&l)[9], double (&cb)[9]){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        cb[0] = N*density*l[0];
        cb[1] = N*density*l[1];;
        cb[2] = N*density*l[2];;
        cb[3] = N*density*l[6];;
        cb[4] = N*density*l[7];;
        cb[5] = N*density*l[8];;
        cb[6] = N*density*l[3];;
        cb[7] = N*density*l[4];;
        cb[8] = N*density*l[5];;
        
        return;
    }
    
    void compute_body_couple(const int &i, const int &j, const double &N, const double &density, const double (&l)[9], double &cb_ij){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        if ( ( i == 0 ) && ( j == 0 ) ){
            cb_ij = N*density*l[0];
        }
        else if ( ( i == 1 ) && ( j == 1 ) ){
            cb_ij = N*density*l[1];
        }
        else if ( ( i == 2 ) && ( j == 2 ) ){
            cb_ij = N*density*l[2];
        }
        else if ( ( i == 1 ) && ( j == 2 ) ){
            cb_ij = N*density*l[6];
        }
        else if ( ( i == 0 ) && ( j == 2 ) ){
            cb_ij = N*density*l[7];
        }
        else if ( ( i == 0 ) && ( j == 1 ) ){
            cb_ij = N*density*l[8];
        }
        else if ( ( i == 2 ) && ( j == 1 ) ){
            cb_ij = N*density*l[3];
        }
        else if ( ( i == 2 ) && ( j == 0 ) ){
            cb_ij = N*density*l[4];
        }
        else if ( ( i == 1 ) && ( j == 0 ) ){
            cb_ij = N*density*l[5];
        }
        else{
            std::cout << "Error: Index out of range\n";
        }
        
        return;
    }
    
    void compute_kinematic_couple(const double &N, const double &density, const double (&omega)[9], double (&ckin)[9]){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        ckin[0] = -N*density*omega[0];
        ckin[1] = -N*density*omega[1];;
        ckin[2] = -N*density*omega[2];;
        ckin[3] = -N*density*omega[6];;
        ckin[4] = -N*density*omega[7];;
        ckin[5] = -N*density*omega[8];;
        ckin[6] = -N*density*omega[3];;
        ckin[7] = -N*density*omega[4];;
        ckin[8] = -N*density*omega[5];;
        
        return;
    }
    
    void compute_kinematic_couple(const int &i, const int &j, const double &N, const double &density, const double (&omega)[9], double &ckin_ij){
        /*!=============================
        |    compute_body_couple    |
        =============================
        
        Compute the body couple term.
        
        */
        
        if ( ( i == 0 ) && ( j == 0 ) ){
            ckin_ij = -N*density*omega[0];
        }
        else if ( ( i == 1 ) && ( j == 1 ) ){
            ckin_ij = -N*density*omega[1];
        }
        else if ( ( i == 2 ) && ( j == 2 ) ){
            ckin_ij = -N*density*omega[2];
        }
        else if ( ( i == 1 ) && ( j == 2 ) ){
            ckin_ij = -N*density*omega[6];
        }
        else if ( ( i == 0 ) && ( j == 2 ) ){
            ckin_ij = -N*density*omega[7];
        }
        else if ( ( i == 0 ) && ( j == 1 ) ){
            ckin_ij = -N*density*omega[8];
        }
        else if ( ( i == 2 ) && ( j == 1 ) ){
            ckin_ij = -N*density*omega[3];
        }
        else if ( ( i == 2 ) && ( j == 0 ) ){
            ckin_ij = -N*density*omega[4];
        }
        else if ( ( i == 1 ) && ( j == 0 ) ){
            ckin_ij = -N*density*omega[5];
        }
        else{
            std::cout << "Error: Index out of range\n";
        }
        
        return;
    }
    
    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, Matrix_3x12 &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.
        
        */
        
        //Compute required DOF jacobians
        SpMat dgrad_udU(9,12);
        SpMat dphidU(9,12);
        SpMat dgrad_phidU(27,12);
        
        construct_dgrad_udU(detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(detadx, dgrad_phidU);

        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute the divergence of DcauchyDU
        DfintDU.row(0) = -(dNdx[0]*DcauchyDU.row(0) + dNdx[1]*DcauchyDU.row(8) + dNdx[2]*DcauchyDU.row(7));
        DfintDU.row(1) = -(dNdx[0]*DcauchyDU.row(5) + dNdx[1]*DcauchyDU.row(1) + dNdx[2]*DcauchyDU.row(6));
        DfintDU.row(2) = -(dNdx[0]*DcauchyDU.row(4) + dNdx[1]*DcauchyDU.row(3) + dNdx[2]*DcauchyDU.row(2));
        
    }
    
    void compute_internal_force_jacobian(const int &i, const int &dof_num, const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, double &DfintDU_iA){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.
        
        The degrees of freedom are organized:
        
        dif_num:   0   1    2       3       4       5       6       7       8       9      10      11 
        DOF:     u_1 u_2, u_3, phi_11, phi_22, phi_33, phi_23, phi_13, phi_12, phi_32, phi_31, phi_21 
        
        and the jacobian is organized
        
        f_{i}{dof_num}
        
        */
        
        //Compute required DOF jacobians
        SpMat _dgrad_udU(9,12);
        SpMat _dphidU(9,12);
        SpMat _dgrad_phidU(27,12);
        
        construct_dgrad_udU(detadx, _dgrad_udU);
        construct_dphidU(eta, _dphidU);
        construct_dgrad_phidU(detadx, _dgrad_phidU);
        
        //TODO: There is probably a more efficient way to do this.
        Matrix_9x12 dgrad_udU    = _dgrad_udU;
        Matrix_9x12 dphidU       = _dphidU;
        Matrix_27x12 dgrad_phidU = _dgrad_phidU;
        
        if ( i == 0 ){
            DfintDU_iA = -(dNdx[0]*(DcauchyDgrad_u.row(0).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(0).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(0).dot(dgrad_phidU.col(dof_num)))
                         + dNdx[1]*(DcauchyDgrad_u.row(8).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(8).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(8).dot(dgrad_phidU.col(dof_num)))
                         + dNdx[2]*(DcauchyDgrad_u.row(7).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(7).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(7).dot(dgrad_phidU.col(dof_num))));
        }
        else if (i == 1){
            DfintDU_iA = -(dNdx[0]*(DcauchyDgrad_u.row(5).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(5).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(5).dot(dgrad_phidU.col(dof_num)))
                         + dNdx[1]*(DcauchyDgrad_u.row(1).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(1).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(1).dot(dgrad_phidU.col(dof_num)))
                         + dNdx[2]*(DcauchyDgrad_u.row(6).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(6).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(6).dot(dgrad_phidU.col(dof_num))));
        }
        else if (i == 2){
            DfintDU_iA = -(dNdx[0]*(DcauchyDgrad_u.row(4).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(4).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(4).dot(dgrad_phidU.col(dof_num)))
                         + dNdx[1]*(DcauchyDgrad_u.row(3).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(3).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(3).dot(dgrad_phidU.col(dof_num)))
                         + dNdx[2]*(DcauchyDgrad_u.row(2).dot(dgrad_udU.col(dof_num))+DcauchyDphi.row(2).dot(dphidU.col(dof_num))+DcauchyDgrad_phi.row(2).dot(dgrad_phidU.col(dof_num))));
        }
        else{
            std::cout << "Error: Index out of range\n";
        }
        
        return;
    }
    
    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          Matrix_9x12 &DcintDU){
        /*!==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple.
        
        */
        
        //Compute required DOF jacobians
        SpMat dgrad_udU(9,12);
        SpMat dphidU(9,12);
        SpMat dgrad_phidU(27,12);
        
        construct_dgrad_udU(detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(detadx, dgrad_phidU);
        
        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute DsDU
        Matrix_9x12 DsDU;
        DsDU = (DsDgrad_u*dgrad_udU + DsDphi*dphidU + DsDgrad_phi*dgrad_phidU);
        
        //Compute DmDU
        Matrix_27x12 DmDU;
        DmDU = (DmDgrad_u*dgrad_udU + DmDphi*dphidU + DmDgrad_phi*dgrad_phidU);
        
        //Transpose the second and third indices of DmDU
        for (int i=0; i<3; i++){
            DmDU.row(9*i+3).swap(DmDU.row(9*i+6));
            DmDU.row(9*i+4).swap(DmDU.row(9*i+7));
            DmDU.row(9*i+5).swap(DmDU.row(9*i+8));
        }
        
        //Compute N_,k DmDU_kji,alpha
        Matrix_9x12 divm;
        for (int i=0; i<9; i++){
            divm.row(i) = dNdx[0]*DmDU.row(0+i) + dNdx[1]*DmDU.row(9+i) + dNdx[2]*DmDU.row(18+i);
        }
        
        DcintDU = N*(DcauchyDU - DsDU) - divm;
    
    }
    
    void compute_internal_couple_jacobian(const int &i, const int &j, const int &dof_num, const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          double &DcintDU_ijA){
        /*!==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple.
        
        The degrees of freedom are organized:
        
        dif_num:   0   1    2       3       4       5       6       7       8       9      10      11 
        DOF:     u_1 u_2, u_3, phi_11, phi_22, phi_33, phi_23, phi_13, phi_12, phi_32, phi_31, phi_21 
        
        and the jacobian is organized
        
        f_{i}{dof_num}
        
        */
        
        //Compute required DOF jacobians
        SpMat _dgrad_udU(9,12);
        SpMat _dphidU(9,12);
        SpMat _dgrad_phidU(27,12);
        
        construct_dgrad_udU(detadx, _dgrad_udU);
        construct_dphidU(eta, _dphidU);
        construct_dgrad_phidU(detadx, _dgrad_phidU);
        
        //TODO: There is probably a more efficient way to do this.
        Matrix_9x12 dgrad_udU    = _dgrad_udU;
        Matrix_9x12 dphidU       = _dphidU;
        Matrix_27x12 dgrad_phidU = _dgrad_phidU;
        
        int Ihat;
        int Jhat;
        if( i==0 && j==0){Ihat = Jhat = 0;}
        else if( i==1 && j==1){Ihat = Jhat = 1;}
        else if( i==2 && j==2){Ihat = Jhat = 2;}
        else if( i==1 && j==2){Ihat = 3; Jhat = 6;}
        else if( i==0 && j==2){Ihat = 4; Jhat = 7;}
        else if( i==0 && j==1){Ihat = 5; Jhat = 8;}
        else if( i==2 && j==1){Ihat = 6; Jhat = 3;}
        else if( i==2 && j==0){Ihat = 7; Jhat = 4;}
        else if( i==1 && j==0){Ihat = 8; Jhat = 5;}
        
        //Compute required term from DcauchyDU;
        double DcauchyDU_ijA;
        DcauchyDU_ijA =   DcauchyDgrad_u.row(Ihat).dot(dgrad_udU.col(dof_num))
                        + DcauchyDphi.row(Ihat).dot(dphidU.col(dof_num))
                        + DcauchyDgrad_phi.row(Ihat).dot(dgrad_phidU.col(dof_num));
        
        //Compute DsDU
        double DsDU_ijA;
        DsDU_ijA =   DsDgrad_u.row(Ihat).dot(dgrad_udU.col(dof_num))
                   + DsDphi.row(Ihat).dot(dphidU.col(dof_num))
                   + DsDgrad_phi.row(Ihat).dot(dgrad_phidU.col(dof_num));
                   
        double divm_ijA;
        divm_ijA =  ( dNdx[0]*DmDgrad_u.row(0+Jhat)   + dNdx[1]*DmDgrad_u.row(9+Jhat)   + dNdx[2]*DmDgrad_u.row(18+Jhat)  ).dot(dgrad_udU.col(dof_num))
                   +( dNdx[0]*DmDphi.row(0+Jhat)      + dNdx[1]*DmDphi.row(9+Jhat)      + dNdx[2]*DmDphi.row(18+Jhat)     ).dot(dphidU.col(dof_num))
                   +( dNdx[0]*DmDgrad_phi.row(0+Jhat) + dNdx[1]*DmDgrad_phi.row(9+Jhat) + dNdx[2]*DmDgrad_phi.row(18+Jhat)).dot(dgrad_phidU.col(dof_num));
        
        DcintDU_ijA = N*(DcauchyDU_ijA - DsDU_ijA) - divm_ijA;
    
    }
    
    void construct_dgrad_udU(const double (&detadx)[3], SpMat &dgrad_udU){
        /*!==========================
        |    construct_dgrad_udU    |
        =============================
        
        Construct the derivative of the gradient of u w.r.t. the DOF vector.
        
        dgrad_udU is a 9x12 matrix
        
        */
        
        std::vector<T> tripletList;
        tripletList.reserve(9);

        //Assemble dgrad_udU
        tripletList.push_back(T(0,0,detadx[0]));
        tripletList.push_back(T(5,0,detadx[1]));
        tripletList.push_back(T(4,0,detadx[2]));
        tripletList.push_back(T(8,1,detadx[0]));
        tripletList.push_back(T(1,1,detadx[1]));
        tripletList.push_back(T(3,1,detadx[2]));
        tripletList.push_back(T(7,2,detadx[0]));
        tripletList.push_back(T(6,2,detadx[1]));
        tripletList.push_back(T(2,2,detadx[2]));
        
        dgrad_udU.setFromTriplets(tripletList.begin(), tripletList.end()); 

        return;
    }
    
    void construct_dphidU(const double &eta, SpMat &dphidU){
        /*!==========================
        |    construct_dphidU    |
        ==========================
        
        Construct the jacobian of phi_iI w.r.t. the degree of freedom vector.
        
        dphidU is a 9x12 matrix
        
        */
        
        std::vector<T> tripletList;
        tripletList.reserve(9);

        //Assemble dphidU
        tripletList.push_back(T(0, 3,eta));
        tripletList.push_back(T(5, 8,eta));
        tripletList.push_back(T(4, 7,eta));
        tripletList.push_back(T(8,11,eta));
        tripletList.push_back(T(1, 4,eta));
        tripletList.push_back(T(3, 6,eta));
        tripletList.push_back(T(7,10,eta));
        tripletList.push_back(T(6, 9,eta));
        tripletList.push_back(T(2, 5,eta));
        
        dphidU.setFromTriplets(tripletList.begin(), tripletList.end());

        return;
    }
    
    void construct_dgrad_phidU(const double (&detadx)[3], SpMat &dgrad_phidU){
        /*!===============================
        |    construct_dgrad_phidU    |
        ===============================
        
        Construct the jacobian of phi_iI,l w.r.t. the degree of freedom vector.
        
        dphidU is a 27x12 matrix
        
        */
        
        std::vector<T> tripletList;
        tripletList.reserve(27);
        
        tripletList.push_back(T( 0, 3,detadx[0]));
        tripletList.push_back(T( 5, 3,detadx[1]));
        tripletList.push_back(T( 4, 3,detadx[2]));
        tripletList.push_back(T( 9,11,detadx[0]));
        tripletList.push_back(T(14,11,detadx[1]));
        tripletList.push_back(T(13,11,detadx[2]));
        tripletList.push_back(T(18,10,detadx[0]));
        tripletList.push_back(T(23,10,detadx[1]));
        tripletList.push_back(T(22,10,detadx[2]));
        tripletList.push_back(T( 8, 8,detadx[0]));
        tripletList.push_back(T( 1, 8,detadx[1]));
        tripletList.push_back(T( 3, 8,detadx[2]));
        tripletList.push_back(T(17, 4,detadx[0]));
        tripletList.push_back(T(10, 4,detadx[1]));
        tripletList.push_back(T(12, 4,detadx[2]));
        tripletList.push_back(T(26, 9,detadx[0]));
        tripletList.push_back(T(19, 9,detadx[1]));
        tripletList.push_back(T(21, 9,detadx[2]));
        tripletList.push_back(T( 7, 7,detadx[0]));
        tripletList.push_back(T( 6, 7,detadx[1]));
        tripletList.push_back(T( 2, 7,detadx[2]));
        tripletList.push_back(T(16, 6,detadx[0]));
        tripletList.push_back(T(15, 6,detadx[1]));
        tripletList.push_back(T(11, 6,detadx[2]));
        tripletList.push_back(T(25, 5,detadx[0]));
        tripletList.push_back(T(24, 5,detadx[1]));
        tripletList.push_back(T(20, 5,detadx[2]));
        
        dgrad_phidU.setFromTriplets(tripletList.begin(), tripletList.end());

        return;
    }
}
