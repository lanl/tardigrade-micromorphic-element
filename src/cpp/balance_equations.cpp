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

    void compute_internal_force(const double (&dNdx)[3], const std::vector<double> &cauchy, double (&fint)[3]){
        /*!================================
        |    compute_internal_force    |
        ================================

        Compute the internal force given the gradient 
        of the shape function and the cauchy stress on 
        one of the components.

        */
        Vector_9 _cauchy;
        map_vector_to_eigen(cauchy,_cauchy);
        compute_internal_force(dNdx, _cauchy, fint);        

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
            assert(1==0); //TODO: Replace with better error handling
        }

        return;
    }

    void compute_internal_force(const int &i, const double (&dNdx)[3], const std::vector<double> &cauchy, double &fint_i){
        /*!================================
        |    compute_internal_force    |
        ================================

        Compute the internal force given the gradient 
        of the shape function and the cauchy stress on 
        one of the components.

        */

        Vector_9 _cauchy;
        map_vector_to_eigen(cauchy,_cauchy);
        compute_internal_force(i, dNdx, _cauchy, fint_i);        

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
        
        cint[0] = -(N*(cauchy[0] - s[0]) - (dNdx[0]*m[0]+dNdx[1]*m[9]+dNdx[2]*m[18]));
        cint[1] = -(N*(cauchy[1] - s[1]) - (dNdx[0]*m[1]+dNdx[1]*m[10]+dNdx[2]*m[19]));
        cint[2] = -(N*(cauchy[2] - s[2]) - (dNdx[0]*m[2]+dNdx[1]*m[11]+dNdx[2]*m[20]));
        cint[3] = -(N*(cauchy[3] - s[3]) - (dNdx[0]*m[6]+dNdx[1]*m[15]+dNdx[2]*m[24]));
        cint[4] = -(N*(cauchy[4] - s[4]) - (dNdx[0]*m[7]+dNdx[1]*m[16]+dNdx[2]*m[25]));
        cint[5] = -(N*(cauchy[5] - s[5]) - (dNdx[0]*m[8]+dNdx[1]*m[17]+dNdx[2]*m[26]));
        cint[6] = -(N*(cauchy[6] - s[6]) - (dNdx[0]*m[3]+dNdx[1]*m[12]+dNdx[2]*m[21]));
        cint[7] = -(N*(cauchy[7] - s[7]) - (dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22]));
        cint[8] = -(N*(cauchy[8] - s[8]) - (dNdx[0]*m[5]+dNdx[1]*m[14]+dNdx[2]*m[23]));
        
        return;
    }
    
    void compute_internal_couple(const double &N, const double (&dNdx)[3], const std::vector<double> &cauchy, const std::vector<double> &s, const std::vector<double> &m, double (&cint)[9]){
        /*!=================================
        |    compute_internal_couple    |
        =================================

        Compute the internal couple for all 
        indices.

        */

        Vector_9  _cauchy;
        Vector_9  _s;
        Vector_27 _m;
        
        map_vector_to_eigen(cauchy,_cauchy);
        map_vector_to_eigen(s,_s);
        map_vector_to_eigen(m,_m);
        compute_internal_couple(N, dNdx, _cauchy, _s, _m, cint);
        
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
            cint_ij = -(N*(cauchy[0] - s[0]) - (dNdx[0]*m[0]+dNdx[1]*m[ 9]+dNdx[2]*m[18]));
        }
        else if ( ( i == 1 ) && ( j == 1 ) ){
            cint_ij = -(N*(cauchy[1] - s[1]) - (dNdx[0]*m[1]+dNdx[1]*m[10]+dNdx[2]*m[19]));
        }
        else if ( ( i == 2 ) && ( j == 2 ) ){
            cint_ij = -(N*(cauchy[2] - s[2]) - (dNdx[0]*m[2]+dNdx[1]*m[11]+dNdx[2]*m[20]));
        }
        else if ( ( i == 1 ) && ( j == 2 ) ){
            cint_ij = -(N*(cauchy[3] - s[3]) - (dNdx[0]*m[6]+dNdx[1]*m[15]+dNdx[2]*m[24]));
        }
        else if ( ( i == 0 ) && ( j == 2 ) ){
            cint_ij = -(N*(cauchy[4] - s[4]) - (dNdx[0]*m[7]+dNdx[1]*m[16]+dNdx[2]*m[25]));
        }
        else if ( ( i == 0 ) && ( j == 1 ) ){
            cint_ij = -(N*(cauchy[5] - s[5]) - (dNdx[0]*m[8]+dNdx[1]*m[17]+dNdx[2]*m[26]));
        }
        else if ( ( i == 2 ) && ( j == 1 ) ){
            cint_ij = -(N*(cauchy[6] - s[6]) - (dNdx[0]*m[3]+dNdx[1]*m[12]+dNdx[2]*m[21]));
        }
        else if ( ( i == 2 ) && ( j == 0 ) ){
            cint_ij = -(N*(cauchy[7] - s[7]) - (dNdx[0]*m[4]+dNdx[1]*m[13]+dNdx[2]*m[22]));
        }
        else if ( ( i == 1 ) && ( j == 0 ) ){
            cint_ij = -(N*(cauchy[8] - s[8]) - (dNdx[0]*m[5]+dNdx[1]*m[14]+dNdx[2]*m[23]));
        }
        else{
            std::cout << "Error: Index out of range\n";
            assert(1==2); //TODO: Replace with better error handling
        }
        
        return;
    }

    void compute_internal_couple(const int &i, const int &j, const double &N, const double (&dNdx)[3], const std::vector<double> &cauchy, const std::vector<double> &s, const std::vector<double> &m, double &cint_ij){
        /*!=================================
        |    compute_internal_couple    |
        =================================

        Compute the internal couple at a given 
        component pair.

        */

        Vector_9  _cauchy;
        Vector_9  _s;
        Vector_27 _m;
        
        map_vector_to_eigen(cauchy,_cauchy);
        map_vector_to_eigen(s,_s);
        map_vector_to_eigen(m,_m);
        compute_internal_couple(i, j, N, dNdx, _cauchy, _s, _m, cint_ij);

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
            assert(2==3); //TODO: Replace with better error handling
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
            assert(3==4); //TODO: Replace with better error handling
        }
        
        return;
    }
    
    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3], const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi, Matrix_3x12 &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.

        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */
        
        //Compute required DOF jacobians
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat dgrad_udU(9,12);
        //SpMat dphidU(9,12);
        //SpMat dgrad_phidU(27,12);
        
        Matrix_9x12  dgrad_udU;
        Matrix_9x12  dphidU;
        Matrix_27x12 dgrad_phidU;

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


    void compute_internal_force_jacobian(const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3],
                                         const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi,
                                         const std::vector<std::vector<double>> &DcauchyDgrad_phi, std::vector<std::vector<double>> &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */

        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;
        Matrix_3x12 _DfintDU;

        map_vector_to_eigen(DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);
        
        compute_internal_force_jacobian(N, dNdx, eta, detadx, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi, _DfintDU);

        map_eigen_to_vector(_DfintDU, DfintDU);

        return;
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
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */
        
        //Compute required DOF jacobians
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat _dgrad_udU(9,12);
        //SpMat _dphidU(9,12);
        //SpMat _dgrad_phidU(27,12);

        Matrix_9x12  dgrad_udU;
        Matrix_9x12  dphidU;
        Matrix_27x12 dgrad_phidU;
        
        construct_dgrad_udU(detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(detadx, dgrad_phidU);
        
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
            assert(4==5); //TODO: Replace with better error handling
        }
        
        return;
    }

    void compute_internal_force_jacobian(const int &i, const int &dof_num, const double &N, const double(&dNdx)[3], const double &eta, const double (&detadx)[3],
                                         const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi,
                                         const std::vector<std::vector<double>> &DcauchyDgrad_phi, double &DfintDU_iA){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */

        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;

        map_vector_to_eigen(DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);
        
        compute_internal_force_jacobian(i, dof_num, N, dNdx, eta, detadx, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi, DfintDU_iA);

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
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */
        
        //Compute required DOF jacobians
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat dgrad_udU(9,12);
        //SpMat dphidU(9,12);
        //SpMat dgrad_phidU(27,12);
        
        Matrix_9x12  dgrad_udU;
        Matrix_9x12  dphidU;
        Matrix_27x12 dgrad_phidU;
        
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
        
        DcintDU = -(N*(DcauchyDU - DsDU) - divm);
    
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
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */
        
        //Compute required DOF jacobians
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat _dgrad_udU(9,12);
        //SpMat _dphidU(9,12);
        //SpMat _dgrad_phidU(27,12);

        Matrix_9x12  _dgrad_udU;
        Matrix_9x12  _dphidU;
        Matrix_27x12 _dgrad_phidU;
        
        construct_dgrad_udU(detadx, _dgrad_udU);
        construct_dphidU(eta, _dphidU);
        construct_dgrad_phidU(detadx, _dgrad_phidU);
        
        //TODO: There is probably a more efficient way to do this.
        Matrix_9x12 dgrad_udU    = _dgrad_udU;
        Matrix_9x12 dphidU       = _dphidU;
        Matrix_27x12 dgrad_phidU = _dgrad_phidU;
        
        int Ihat = 0;
        int Jhat = 0;
        if( i==0 && j==0){Ihat = Jhat = 0;}
        else if( i==1 && j==1){Ihat = Jhat = 1;}
        else if( i==2 && j==2){Ihat = Jhat = 2;}
        else if( i==1 && j==2){Ihat = 3; Jhat = 6;}
        else if( i==0 && j==2){Ihat = 4; Jhat = 7;}
        else if( i==0 && j==1){Ihat = 5; Jhat = 8;}
        else if( i==2 && j==1){Ihat = 6; Jhat = 3;}
        else if( i==2 && j==0){Ihat = 7; Jhat = 4;}
        else if( i==1 && j==0){Ihat = 8; Jhat = 5;}
        else {
            std::cout << "Error: (i,j) = ("<<i<<","<<j<<") not recognized.";
            assert(-1==-2);
        }
        
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
        
        DcintDU_ijA = -(N*(DcauchyDU_ijA - DsDU_ijA) - divm_ijA);
    
    }
    
    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          std::vector<std::vector<double>> &DcintDU){
        /*!==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple.
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */

        Matrix_9x9   _DcauchyDgrad_u;
        Matrix_9x9   _DcauchyDphi;
        Matrix_9x27  _DcauchyDgrad_phi;

        Matrix_9x9   _DsDgrad_u;
        Matrix_9x9   _DsDphi;
        Matrix_9x27  _DsDgrad_phi;
        
        Matrix_27x9  _DmDgrad_u;
        Matrix_27x9  _DmDphi;
        Matrix_27x27 _DmDgrad_phi;

        Matrix_9x12  _DcintDU;

        map_vector_to_eigen(DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);

        map_vector_to_eigen(DsDgrad_u,   _DsDgrad_u);
        map_vector_to_eigen(DsDphi,      _DsDphi);
        map_vector_to_eigen(DsDgrad_phi, _DsDgrad_phi);

        map_vector_to_eigen(DmDgrad_u,   _DmDgrad_u);
        map_vector_to_eigen(DmDphi,      _DmDphi);
        map_vector_to_eigen(DmDgrad_phi, _DmDgrad_phi);

        compute_internal_couple_jacobian(N, dNdx, eta, detadx,
                                         _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                         _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                                         _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                                         _DcintDU);

        map_eigen_to_vector(_DcintDU, DcintDU);
        return;
    }

    void compute_internal_couple_jacobian(const int &i, const int &j, const int &dof_num, const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          double &DcintDU_ijA){

        /*!
        ==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================

        Compute the jacobian of the internal body couple.
        
        The degrees of freedom are organized:
        
        dif_num:   0   1    2       3       4       5       6       7       8       9      10      11 
        DOF:     u_1 u_2, u_3, phi_11, phi_22, phi_33, phi_23, phi_13, phi_12, phi_32, phi_31, phi_21 
        
        and the jacobian is organized
        
        f_{i}{dof_num}
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */
        Matrix_9x9   _DcauchyDgrad_u;
        Matrix_9x9   _DcauchyDphi;
        Matrix_9x27  _DcauchyDgrad_phi;

        Matrix_9x9   _DsDgrad_u;
        Matrix_9x9   _DsDphi;
        Matrix_9x27  _DsDgrad_phi;
        
        Matrix_27x9  _DmDgrad_u;
        Matrix_27x9  _DmDphi;
        Matrix_27x27 _DmDgrad_phi;

        map_vector_to_eigen(DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);

        map_vector_to_eigen(DsDgrad_u,   _DsDgrad_u);
        map_vector_to_eigen(DsDphi,      _DsDphi);
        map_vector_to_eigen(DsDgrad_phi, _DsDgrad_phi);

        map_vector_to_eigen(DmDgrad_u,   _DmDgrad_u);
        map_vector_to_eigen(DmDphi,      _DmDphi);
        map_vector_to_eigen(DmDgrad_phi, _DmDgrad_phi);

        compute_internal_couple_jacobian(i, j, dof_num, N, dNdx, eta, detadx,
                                         _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                         _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                                         _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                                         DcintDU_ijA);
        return;
    }

    /*
    =============================================
    |    Jacobians for current configuration    |
    =============================================
    */
    void compute_internal_force_jacobian(const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&phi)[9],        const Matrix_3x3 &F,
                                         const Vector_9 &cauchy, const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                         Matrix_3x12 &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force in the 
        current configuration.

        */
        
        //Compute required DOF jacobians
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat dgrad_udU(9,12);
        //SpMat dphidU(9,12);
        //SpMat dgrad_phidU(27,12);
        
        Matrix_9x12  dgrad_udU;
        Matrix_9x12  dphidU;
        Matrix_27x12 dgrad_phidU;

        //Compute the test and interpolation function gradients in the 
        //reference configuration.
        double dNdX[3]   = {0,0,0};
        double detadX[3] = {0,0,0};
        for (int I=0; I<3; I++){
            for (int i=0; i<3; i++){
                dNdX[I] += dNdx[i]*F(i,I);
                detadX[I] += detadx[i]*F(i,I);
            }
        }
        //Compute the inverse deformation gradient
        Matrix_3x3 Finv;
        Finv = F.inverse();

        construct_dgrad_udU(Finv, detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(phi, Finv, detadX, dgrad_udU, dgrad_phidU);

        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute the divergence of DcauchyDU
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        //Note: dfintdU = -D/DU(N,j sigma_ji) = -N,J (DFinvdU_JjIhat sigma_ji + Finv_Jj DcauchyDU jiIhat) = -N,J (-Dgrad_udU_JjIhat sigma_ji + Finv_Jj DsigmaDU_jiIhat) 

        double tmp;
        int Jhat;
        int Khat;

        for (int i=0; i<3; i++){
            for (int Ihat=0; Ihat<12; Ihat++){
                tmp = 0;
                for (int J=0; J<3; J++){
                    for (int j=0; j<3; j++){
                        Jhat = sot_to_voigt_map[J][j];
                        Khat = sot_to_voigt_map[j][i];
                        tmp += -dNdX[J]*(Finv(J,j) * DcauchyDU(Khat,Ihat) - dgrad_udU(Jhat,Ihat)*cauchy(Khat));

                    }
                }
                DfintDU(i,Ihat) = tmp;
            }
        }
        return;
    }

    void compute_internal_force_jacobian(const int &component,   const int &dof_num,
                                         const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&phi)[9],        const Matrix_3x3 &F,
                                         const Vector_9 &cauchy, const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                         double &DfintDU_iA){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force in the 
        current configuration.

        */
        
        //Compute required DOF jacobians
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat dgrad_udU(9,12);
        //SpMat dphidU(9,12);
        //SpMat dgrad_phidU(27,12);
        
        Matrix_9x12  dgrad_udU;
        Matrix_9x12  dphidU;
        Matrix_27x12 dgrad_phidU;

        //Compute the test and interpolation function gradients in the 
        //reference configuration.
        double dNdX[3]   = {0,0,0};
        double detadX[3] = {0,0,0};
        for (int I=0; I<3; I++){
            for (int i=0; i<3; i++){
                dNdX[I] += dNdx[i]*F(i,I);
                detadX[I] += detadx[i]*F(i,I);
            }
        }
        //Compute the inverse deformation gradient
        Matrix_3x3 Finv;
        Finv = F.inverse();

        construct_dgrad_udU(Finv,detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(phi, Finv, detadX, dgrad_udU, dgrad_phidU);

        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute the divergence of DcauchyDU
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        //Note: dfintdU = -D/DU(N,j sigma_ji) = -N,J (DFinvdU_JjIhat sigma_ji + Finv_Jj DcauchyDU jiIhat) = -N,J (-Dgrad_udU_JjIhat sigma_ji + Finv_Jj DsigmaDU_jiIhat) 

        int Jhat;
        int Khat;

        DfintDU_iA = 0;

        for (int J=0; J<3; J++){
            for (int j=0; j<3; j++){
                Jhat = sot_to_voigt_map[J][j];
                Khat = sot_to_voigt_map[j][component];
                DfintDU_iA += -dNdX[J]*(Finv(J,j) * DcauchyDU(Khat,dof_num) - dgrad_udU(Jhat,dof_num)*cauchy(Khat));
            }
        }

        return;
    }

    void compute_internal_force_jacobian(const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&phi)[9],        const std::vector<std::vector<double>> &F,
                                         const std::vector<double> &cauchy, const std::vector<std::vector<double>> &DcauchyDgrad_u, 
                                         const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                         std::vector<std::vector<double>> &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force in the 
        current configuration.

        */

        //Map the vectors to eigen matrices and vectors
        Matrix_3x3  _F;
        Vector_9    _cauchy;
        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;
        Matrix_3x12 _DfintDU;

        map_vector_to_eigen(               F,                _F);
        map_vector_to_eigen(          cauchy,           _cauchy);
        map_vector_to_eigen(  DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(     DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);

        compute_internal_force_jacobian(       N,            dNdx,          eta,            detadx,
                                             phi,              _F,
                                         _cauchy, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                        _DfintDU);

        map_eigen_to_vector(_DfintDU,DfintDU);        
        return;
    }

    void compute_internal_force_jacobian(const int &component,   const int &dof_num,
                                         const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&phi)[9],        const std::vector<std::vector<double>> &F,
                                         const std::vector<double> &cauchy, const std::vector<std::vector<double>> &DcauchyDgrad_u, 
                                         const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                         double &DfintDU_iA){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force in the 
        current configuration.

        */

        //Map the vectors to eigen matrices and vectors
        Matrix_3x3  _F;
        Vector_9    _cauchy;
        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;

        map_vector_to_eigen(               F,                _F);
        map_vector_to_eigen(          cauchy,           _cauchy);
        map_vector_to_eigen(  DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(     DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);

        compute_internal_force_jacobian(  component,         dof_num,
                                                  N,            dNdx,          eta,            detadx,
                                                phi,              _F,
                                            _cauchy, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                         DfintDU_iA);

        return;
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

    void construct_dgrad_udU(const double (&detadx)[3], Matrix_9x12 &dgrad_udU){
        /*!==========================
        |    construct_dgrad_udU    |
        =============================
        
        Construct the derivative of the gradient of u w.r.t. the DOF vector.
        
        dgrad_udU is a 9x12 matrix
        
        */

        dgrad_udU = Matrix_9x12::Zero();

        dgrad_udU(0,0) = detadx[0];
        dgrad_udU(5,0) = detadx[1];
        dgrad_udU(4,0) = detadx[2];
        dgrad_udU(8,1) = detadx[0];
        dgrad_udU(1,1) = detadx[1];
        dgrad_udU(3,1) = detadx[2];
        dgrad_udU(7,2) = detadx[0];
        dgrad_udU(6,2) = detadx[1];
        dgrad_udU(2,2) = detadx[2];

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
    
    void construct_dphidU(const double &eta, Matrix_9x12 &dphidU){
        /*!==========================
        |    construct_dphidU    |
        ==========================
        
        Construct the jacobian of phi_iI w.r.t. the degree of freedom vector.
        
        dphidU is a 9x12 matrix
        
        */

        dphidU = Matrix_9x12::Zero();

        dphidU( 0, 3) = eta;
        dphidU( 5, 8) = eta;
        dphidU( 4, 7) = eta;
        dphidU( 8,11) = eta;
        dphidU( 1, 4) = eta;
        dphidU( 3, 6) = eta;
        dphidU( 7,10) = eta;
        dphidU( 6, 9) = eta;
        dphidU( 2, 5) = eta;

        return;
    }

    void construct_dgrad_phidU(const double (&detadx)[3], SpMat &dgrad_phidU){
        /*!===============================
        |    construct_dgrad_phidU    |
        ===============================
        
        Construct the jacobian of phi_iI,l w.r.t. the degree of freedom vector.
        
        dgrad_phidU is a 27x12 matrix
        
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

    void construct_dgrad_phidU(const double (&detadx)[3], Matrix_27x12 &dgrad_phidU){
        /*!===============================
        |    construct_dgrad_phidU    |
        ===============================
        
        Construct the jacobian of phi_iI,l w.r.t. the degree of freedom vector.
        
        dgrad_phidU is a 27x12 matrix
        
        */

        dgrad_phidU = Matrix_27x12::Zero();

        dgrad_phidU( 0, 3) = detadx[0];
        dgrad_phidU( 5, 3) = detadx[1];
        dgrad_phidU( 4, 3) = detadx[2];
        dgrad_phidU( 9,11) = detadx[0];
        dgrad_phidU(14,11) = detadx[1];
        dgrad_phidU(13,11) = detadx[2];
        dgrad_phidU(18,10) = detadx[0];
        dgrad_phidU(23,10) = detadx[1];
        dgrad_phidU(22,10) = detadx[2];
        dgrad_phidU( 8, 8) = detadx[0];
        dgrad_phidU( 1, 8) = detadx[1];
        dgrad_phidU( 3, 8) = detadx[2];
        dgrad_phidU(17, 4) = detadx[0];
        dgrad_phidU(10, 4) = detadx[1];
        dgrad_phidU(12, 4) = detadx[2];
        dgrad_phidU(26, 9) = detadx[0];
        dgrad_phidU(19, 9) = detadx[1];
        dgrad_phidU(21, 9) = detadx[2];
        dgrad_phidU( 7, 7) = detadx[0];
        dgrad_phidU( 6, 7) = detadx[1];
        dgrad_phidU( 2, 7) = detadx[2];
        dgrad_phidU(16, 6) = detadx[0];
        dgrad_phidU(15, 6) = detadx[1];
        dgrad_phidU(11, 6) = detadx[2];
        dgrad_phidU(25, 5) = detadx[0];
        dgrad_phidU(24, 5) = detadx[1];
        dgrad_phidU(20, 5) = detadx[2];

        return;
    }
    
    void construct_dgrad_udU(const Matrix_3x3 &Finv, const double (&detadx)[3], Matrix_9x12 &dgrad_udU){
        /*!==========================
        |    construct_dgrad_udU    |
        =============================
        
        Construct the derivative of the gradient of u w.r.t. the DOF vector.
        
        dgrad_udU is a 9x12 matrix
        
        This is for the balance equation in the current configuration.
        
        */
        
        dgrad_udU = Matrix_9x12::Zero();

        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};


        int Ihat;
        for (int k=0; k<3; k++){
            for (int i=0; i<3; i++){
                Ihat = sot_to_voigt_map[k][i];
                for (int Jhat=0; Jhat<3; Jhat++){
                    dgrad_udU(Ihat,Jhat) = Finv(k,Jhat)*detadx[i];
                }
            }
        }

        return;
    }
    
    void construct_dgrad_phidU(const double (&phi)[9], const Matrix_3x3 &Finv, const double (&detadX)[3], const Matrix_9x12 &dgrad_udU, Matrix_27x12 &dgrad_phidU){
        /*!
        ===============================
        |    construct_dgrad_phidU    |
        ===============================
        
        Construct the derivative of the gradient of grad_phi in the current configuration.
        */
        
        dgrad_phidU = Matrix_27x12::Zero();
        
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
                                  
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
                                  
        double tmp;
        double tmp1;
        int Jhat;
        int Khat;
        int Lhat;
        
        for (int i=0; i<3; i++){
            for (int I=0; I<3; I++){
                Jhat = sot_to_voigt_map[i][I];
                tmp1 = phi[Jhat];
                for (int j=0; j<3; j++){
                    Khat = tot_to_voigt_map[i][I][j];
                    for (int Ihat=0; Ihat<12; Ihat++){
                        tmp = 0;
                        
                        if(Jhat == (Ihat-3)){
                        
                            for (int J=0; J<3; J++){
                                Lhat = sot_to_voigt_map[J][j];
                                tmp += detadX[J]*(-dgrad_udU(Lhat,Ihat)*tmp1 + Finv(J,j));
                            }
                        
                        }
                        else{
                            
                            for (int J=0; J<3; J++){
                                Lhat = sot_to_voigt_map[J][j];
                                tmp += detadX[J]*(-dgrad_udU(Lhat,Ihat)*tmp1);
                            }
                            
                        }
                        
                        
                        dgrad_phidU(Khat,Ihat) = tmp;
                    }
                }
            }
        }
        return;
    }

    void map_eigen_to_vector(const Vector_9 &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        int A = V.size();
        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void map_eigen_to_vector(const Vector_27 &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        int A = V.size();
        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void map_eigen_to_vector(const Eigen::VectorXd &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        int A = V.size();
        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void map_eigen_to_vector(const Matrix_3x3 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.
        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();
        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

   void map_eigen_to_vector(const Matrix_9x9 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_eigen_to_vector(const Matrix_9x27 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_eigen_to_vector(const Matrix_27x9 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_eigen_to_vector(const Matrix_27x27 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */
        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_eigen_to_vector(const Matrix_3x12 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_eigen_to_vector(const Matrix_9x12 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_eigen_to_vector(const Eigen::MatrixXd &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void map_vector_to_eigen(const std::vector<double> &v, Vector_9 &V){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = V.size();
        unsigned int a = v.size();

        if (A != a){std::cout << "Error: Vectors are of different sizes!\n";assert(5==6);}
        for (unsigned int i=0; i<a; i++){
            V[i] = v[i];
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<double> &v, Vector_27 &V){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = V.size();
        unsigned int a = v.size();

        if (A != a){std::cout << "Error: Vectors are of different sizes!\n";assert(6==7);}
        for (unsigned int i=0; i<a; i++){
            V[i] = v[i];
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<double> &v, Eigen::VectorXd &V){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = V.size();
        unsigned int a = v.size();

        if (A != a){std::cout << "Error: Vectors are of different sizes!\n";assert(7==8);}
        for (unsigned int i=0; i<a; i++){
            V[i] = v[i];
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_3x3 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(8==9);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(9==10);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_3x12 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(8==9);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(9==10);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_9x12 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(10==11);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(11==12);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_9x9 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(12==13);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(13==14);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_9x27 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(14==15);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(15==16);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_27x9 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(16==17);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(17==18);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Matrix_27x27 &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(18==19);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(19==20);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }


    void map_vector_to_eigen(const std::vector<std::vector<double>> &v, Eigen::MatrixXd &M){
        /*!
        =============================
        |    map_vector_to_eigen    |
        =============================

        Map a standard vector to an eigen matrix.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a != A){
            std::cout << "Error: Vector and matrix do not have the same number of rows!\n";
            assert(20==21);
        }

        for (unsigned int i=0; i<a; i++){
            if(v[i].size() != B){
                std::cout << "Error Vector and matrix do not have the same number of columns!\n";
                assert(21==22);
            }
            for (unsigned int j=0; j<v[i].size(); j++){
                M(i,j) = v[i][j];
            }
        }
        return;
    }

}
