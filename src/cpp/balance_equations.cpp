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

#include <balance_equations.h>

namespace balance_equations{

    void compute_internal_force( const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &PK2, double ( &fint )[ 3 ]){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress.
         *
         * Note: This has already accounted for the volume change term.
         *
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double ( &fint )[ 3 ]: The internal force.
         */

        for ( unsigned int i = 0; i < 3; i++ ){
            compute_internal_force( i, dNdX, F, PK2, fint[ i ] );
        }

        return;
    }
    
    void compute_internal_force( const int &i,  const double (&dNdX)[3], const variableVector &F, const variableVector &PK2,
                                 double &fint ){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress on a single component.
         *
         * :param const int &i: The component to compute the stress on.
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double ( &fint )[ 3 ]: The internal force.
         */

        fint_i = 0;
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                fint_i -= dNdX[ I ] * PK2[ dim * I + J ] * F[ dim * i + J ];
            }
        }

        return;
    }


    void compute_body_force( const double &N, const double &density, const double ( &b )[ 3 ], double ( &fb )[ 3 ] ){
        /*!
         * Compute the body force
         *
         * :param const double &N: The value of the shape function
         * :param const double &density: The density in the current configuration
         * :param const double ( &b )[ 3 ]: The body force in the current configuration
         * :param const double ( &fb )[ 3 ]: The body force in the current configuration
         */
        
        for ( int i = 0; i < 3; i++ ){
            fb[ i ] = compute_body_force( i, N, density, b, fb[ i ] );
        }

        return;
    }
    
    void compute_body_force( const int &i, const double &N, const double &density, const double ( &b )[ 3 ], double &fb_i ){
        /*!
         * Compute a single component of the body force
         *
         * :param const int &i: The component to compute the body force on
         * :param const double &N: The value of the shape function
         * :param const double &density: The density in the current configuration
         * :param const double ( &b )[ 3 ]: The body force in the current configuration
         * :param const double fb_i: The body force in the current configuration in direction i
         */
        
        fb_i = N * density * b[ i ];

        return;
    }

    void compute_inertial_force( const double &N, const double &density, const double ( &a )[ 3 ], double ( &finertial )[ 3 ] ){
        /*!
         * Compute the inertial force
         *
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const double ( &finertial )[ 3 ]: The inertial force in the current configuration
         */
        
        for (int i=0; i<3; i++){
            compute_inertial_force( N, density, a, finertial[ i ];
        }

        return;
    }

    void compute_inertial_force( const int &i, const double &N, const double &density, const double ( &a )[ 3 ], double finertial_i ){
        /*!
         * Compute the inertial force for the ith component.
         *
         * :param const int &i: The component to compute the inertial force on
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const double ( &finertial )[ 3 ]: The inertial force in the current configuration in direction i
         */
       
        finertial_i = N * density * a[ i ]; 

        return;
    }
    
    void compute_internal_couple( const double &N, const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &chi,
                                  const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                  double ( &cint )[ 9 ] ){
        /*!
         * Compute the internal couple defined as
         * cint_{ ij } = N F_{ iI } ( PK2_{ IJ } - SIGMA_{ IJ } ) F_{ jJ } - N_{ ,K } F_{ jJ } \chi_{ iI } M_{ KJI }
         *
         * :param const double &N: The shape-function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference configuration
         * :param const variableVector &F: The deformation gradient
         * :param const variableVector &chi: The micro-deformation
         * :param const variableVector &PK2: The second Piola-Kirchoff stress
         * :param const variableVector &SIGMA: The symmetric micro stress in the reference configuration.
         * :param const variableVector &M: The higher order stress in the reference configuration.
         * :param double ( &cint )[ 9 ] ): The internal couple stored as
         *     [ cint_{ 11 }, cint_{ 12 }, cint_{ 13 }, cint_{ 21 }, cint_{ 22 }, cint_{ 23 }, cint_{ 31 }, cint_{ 32 }, cint_{ 33 } ]
         */

        //Assume 3d
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                compute_internal_couple( i, j, N, dNdX, F, chi, PK2, SIGMA, M, cint[ dim * i + j ] );
            }
        }

        return;
    }

    void compute_internal_couple( const int &i, const int &j, const double &N, const double ( &dNdX )[ 3 ],
                                  const variableVector &F, const variableVector &chi,
                                  const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                  double cint_ij ){
        /*!
         * Compute the internal couple at index ij defined as
         * cint_{ ij } = N F_{ iI } ( PK2_{ IJ } - SIGMA_{ IJ } ) F_{ jJ } - N_{ ,K } F_{ jJ } \chi_{ iI } M_{ KJI }
         *
         * :param const int &i: The first index.
         * :param const int &j: The second index.
         * :param const double &N: The shape-function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference configuration
         * :param const variableVector &F: The deformation gradient
         * :param const variableVector &chi: The micro-deformation
         * :param const variableVector &PK2: The second Piola-Kirchoff stress
         * :param const variableVector &SIGMA: The symmetric micro stress in the reference configuration.
         * :param const variableVector &M: The higher order stress in the reference configuration.
         * :param double ( &cint )[ 9 ] ): The internal couple stored as
         *     [ cint_{ 11 }, cint_{ 12 }, cint_{ 13 }, cint_{ 21 }, cint_{ 22 }, cint_{ 23 }, cint_{ 31 }, cint_{ 32 }, cint_{ 33 } ]
         */

        //Assume 3D
        const unsigned int dim = 3;

        cint_ij = 0;
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                cint_ij += N * F[ dim * i + I ] * ( PK2[ dim * I + J ] - SIGMA[ dim * I + J ] ) * F[ dim * j + J ];

                for ( unsigned int K = 0; K < dim; K++ ){
                    cint_ij -= dNdX[ K ] * F[ dim * j + J ] * chi[ dim * i + I ] * M[ dim * dim * K + dim * J + I ];
                }
            }
        }

        return;
    }
    
    void compute_body_couple( const double &N, const double &density, const double ( &l )[ 9 ], double ( &cb )[ 9 ] ){
        /*!
         * Compute the body couple term
         *
         * couple_body_{ ij } = N * density * l_{ ji }
         *
         * :param const double &N: The shape-function value.
         * :param const double &density: The density in the current configuration.
         * :param const double ( %l )[ 9 ]: The body-force couple in the current configuration.
         * :param const double ( &cb )[ 9 ]: The couple term.
         */
        
        //Assume 3D
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                compute_body_couple( N, density, l, cb[ dim * i + j ] );
            }
        }
        
        return;
    }

    void compute_body_couple( const int &i, const int &j, const double &N, const double &density, const double ( &l )[ 9 ], 
                              double cb_ij ){
        /*!
         * Compute the body couple term for the indices i and j
         *
         * couple_body_{ ij } = N * density * l_{ ji }
         *
         * :param const int &i: The first index.
         * :param const int &j: The second index.
         * :param const double &N: The shape-function value.
         * :param const double &density: The density in the current configuration.
         * :param const double ( %l )[ 9 ]: The body-force couple in the current configuration.
         * :param const double ( &cb )[ 9 ]: The couple term.
         */

        //Assume 3D
        const unsigned int dim = 3;
        
        cb_ij = N * density * l[ dim * j + i ];
        
        return;
    }
    
    void compute_inertial_couple( const double &N, const double &density, const double ( &omega )[ 9 ], double ( &cinertial )[ 9 ] ){
        /*!
         * Compute the inertial couple in the current configuration
         *
         * cinertial_{ ij } = -N * density * omega_{ ji }
         *
         * :param const double &N: The shape-function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &omega )[ 9 ]: The micro-inertia tensor in the current configuration.
         * :param const double ( &cinertial )[ 9 ]: The inertial couple.
         */

        //Assume 3D
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                compute_internal_couple( i, j, N, density, omega, cinertial[ dim * i + j ] );
            }
        }
        
        return;
    }
    
    void compute_inertial_couple( const int &i, const int &j, const double &N, const double &density, const double ( &omega )[ 9 ], 
                                  double cinertial_ij ){
        /*!
         * Compute the inertial couple in the current configuration for the indices i and j
         *
         * cinertial_{ ij } = -N * density * omega_{ ji }
         *
         * :param const int &i: The first index.
         * :param const int &j: The second index.
         * :param const double &N: The shape-function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &omega )[ 9 ]: The micro-inertia tensor in the current configuration.
         * :param const double ( &cinertial )[ 9 ]: The inertial couple.
         */

        //Assume 3D
        const unsigned int dim = 3;

        cinertial_ij = - N * density * omega[ dim * j + i ];
        
        return;
    }
    
    void compute_internal_force_jacobian(const double &N, const double(&dNdX)[3], const double &eta, const double (&detadX)[3], const Matrix_3x3 &F,
                                         const Vector_9 &PK2, const Matrix_9x9 &DPK2Dgrad_u, const Matrix_9x9 &DPK2Dphi, const Matrix_9x27 &DPK2Dgrad_phi,
                                         Matrix_3x12 &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.

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

        construct_dgrad_udU(detadX, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(detadX, dgrad_phidU);

        //Compute DcauchyDU;
        Matrix_9x12 DPK2DU;
        DPK2DU = (DPK2Dgrad_u*dgrad_udU + DPK2Dphi*dphidU + DPK2Dgrad_phi*dgrad_phidU);

        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        //Must compute DfintDU_jIhat = -dNdX[I]*(DPK2DU_IJIhat F_jJ + PK2_IJ DFDU_jJIhat)
        //NOTE: DFDU_jJIhat = Dgrad_udU_jJIhat

        double tmp;
        int Jhat;
        int Khat;

        for (int j=0; j<3; j++){
            for (int Ihat=0; Ihat<12; Ihat++){
                tmp = 0;
            
                for (int I=0; I<3; I++){
                    for (int J=0; J<3; J++){
                        Jhat = sot_to_voigt_map[I][J];
                        Khat = sot_to_voigt_map[j][J];
                        tmp += -dNdX[I]*(DPK2DU(Jhat,Ihat) * F(j,J) + PK2(Jhat)*dgrad_udU(Khat,Ihat));
                    }
                }
                DfintDU(j,Ihat) = tmp;

            }
        }
        
    }


    void compute_internal_force_jacobian(const double &N, const double(&dNdX)[3], const double &eta, const double (&detadX)[3], const std::vector<std::vector<double>> &F,
                                         const std::vector<double> &PK2, const std::vector<std::vector<double>> &DPK2Dgrad_u, const std::vector<std::vector<double>> &DPK2Dphi,
                                         const std::vector<std::vector<double>> &DPK2Dgrad_phi, std::vector<std::vector<double>> &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.
        
        Note: Currently an incorrect implementation which is retained to be the 
              foundation of a total-lagrangian implementation.        
        */

        Matrix_3x3  _F;
        Vector_9    _PK2;
        Matrix_9x9  _DPK2Dgrad_u;
        Matrix_9x9  _DPK2Dphi;
        Matrix_9x27 _DPK2Dgrad_phi;
        Matrix_3x12 _DfintDU;

        map_vector_to_eigen(F,             _F);
        map_vector_to_eigen(PK2,           _PK2);
        map_vector_to_eigen(DPK2Dgrad_u,   _DPK2Dgrad_u);
        map_vector_to_eigen(DPK2Dphi,      _DPK2Dphi);
        map_vector_to_eigen(DPK2Dgrad_phi, _DPK2Dgrad_phi);
        
        compute_internal_force_jacobian(N, dNdX, eta, detadX, _F, _PK2, _DPK2Dgrad_u, _DPK2Dphi, _DPK2Dgrad_phi, _DfintDU);

        map_eigen_to_vector(_DfintDU, DfintDU);

        return;
    }

    void compute_internal_force_jacobian(const int &component_i, const int &dof_num, const double &N, const double(&dNdX)[3], const double &eta, const double (&detadX)[3], const Matrix_3x3 &F, const Vector_9 &PK2, const Matrix_9x9 &DPK2Dgrad_u, const Matrix_9x9 &DPK2Dphi, const Matrix_9x27 &DPK2Dgrad_phi, double &DfintDU_iA){
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
        //Note: Removing sparse matrices because there 
        //      is a segmentation fault when being created 
        //      while in a moose run.
        //SpMat _dgrad_udU(9,12);
        //SpMat _dphidU(9,12);
        //SpMat _dgrad_phidU(27,12);

        Vector_9  dgrad_udU;
        Vector_9  dphidU;
        Vector_27 dgrad_phidU;

        construct_dgrad_udU(detadX, dof_num, dgrad_udU);
        construct_dphidU(eta, dof_num, dphidU);
        construct_dgrad_phidU(detadX, dof_num, dgrad_phidU);

        //Compute DcauchyDU;
        Vector_9 DPK2DU;
        DPK2DU = (  DPK2Dgrad_u*dgrad_udU
                  + DPK2Dphi*dphidU
                  + DPK2Dgrad_phi*dgrad_phidU);

        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};

        //Must compute DfintDU_jIhat = -dNdX[I]*(DPK2DU_IJIhat F_jJ + PK2_IJ DFDU_jJIhat)
        //NOTE: DFDU_jJIhat = Dgrad_udU_jJIhat

        //double tmp;
        int Jhat;
        int Khat;

        DfintDU_iA = 0;

        for (int I=0; I<3; I++){
            for (int J=0; J<3; J++){
                Jhat = sot_to_voigt_map[I][J];
                Khat = sot_to_voigt_map[component_i][J];
                DfintDU_iA += -dNdX[I]*(DPK2DU(Jhat) * F(component_i,J) + PK2(Jhat)*dgrad_udU(Khat));
            }
        }

        return;
    }

    void compute_internal_force_jacobian(const int &component_i, const int &dof_num, const double &N, const double(&dNdX)[3], const double &eta, const double (&detadX)[3], const std::vector<std::vector<double>> &F,
                                         const std::vector<double> &PK2, const std::vector<std::vector<double>> &DPK2Dgrad_u, const std::vector<std::vector<double>> &DPK2Dphi,
                                         const std::vector<std::vector<double>> &DPK2Dgrad_phi, double &DfintDU_iA){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force.
        
        Total Lagrangian formulation.  
        */

        std::vector<double> dgrad_udU(9,0);
        std::vector<double> dphidU(9,0);
        std::vector<double> dgrad_phidU(27,0);

        std::vector<int> non_zero_terms_dgrad_udU(3,0);
        std::vector<int> non_zero_terms_dgrad_phidU(3,0);

        construct_dgrad_udU(detadX, dof_num, non_zero_terms_dgrad_udU, dgrad_udU);
        construct_dphidU(eta, dof_num, dphidU);
        construct_dgrad_phidU(detadX, dof_num, non_zero_terms_dgrad_phidU, dgrad_phidU);

        //Compute DcauchyDU;
        std::vector<double> DPK2DU(9,0);

//        double tmp;
//        int    nzt;
//        for (int i=0; i<9; i++){
//            tmp = DPK2Dphi[i][dof_num-3]*eta;
//            for (int j=0; j<3; j++){
//                nzt  = non_zero_terms_dgrad_udU[j];
//                tmp +=   DPK2Dgrad_u[i][nzt]*dgrad_udU[nzt];
////                       + DPK2Dphi[i][j]*dphidU[j]
////                       + DPK2Dgrad_phi[i][j]*dgrad_phidU[j];
//            }
//
//            for (int j=0; j<3; j++){
//                nzt  = non_zero_terms_dgrad_phidU[j];
//                tmp += DPK2Dgrad_phi[i][nzt]*dgrad_phidU[nzt];
//            }
//
//            DPK2DU[i] = tmp;
//        }

        double tmp;
        int    nzt;
        if(dof_num<3){

            for (int i=0; i<9; i++){
                tmp = 0;
                for (int j=0; j<3; j++){
                    nzt   =   non_zero_terms_dgrad_udU[j];
                    tmp +=   DPK2Dgrad_u[i][nzt]*dgrad_udU[nzt];
                }

                DPK2DU[i]   = tmp;
            }
        }
        else{

            for (int i=0; i<9; i++){
                tmp = DPK2Dphi[i][dof_num-3]*eta;//0;

                for (int j=0; j<3; j++){
                    nzt   = non_zero_terms_dgrad_phidU[j];
                    tmp += DPK2Dgrad_phi[i][nzt]*dgrad_phidU[nzt];
                }

                DPK2DU[i]   = tmp;
            }
        }

        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        //Must compute DfintDU_jIhat = -dNdX[I]*(DPK2DU_IJIhat F_jJ + PK2_IJ DFDU_jJIhat)
        //NOTE: DFDU_jJIhat = Dgrad_udU_jJIhat

        //double tmp;
        int Jhat;
        int Khat;

        DfintDU_iA = 0;

        for (int I=0; I<3; I++){
            for (int J=0; J<3; J++){
                Jhat = sot_to_voigt_map[I][J];
                Khat = sot_to_voigt_map[component_i][J];
                DfintDU_iA += -dNdX[I]*(DPK2DU[Jhat] * F[component_i][J] + PK2[Jhat]*dgrad_udU[Khat]);
            }
        }

        return;
    }
    
    void compute_internal_couple_jacobian(const double &N, const double (&dNdX)[3], const double &eta, const double (&detadx)[3],
                                          const Matrix_3x3  &F,             const Matrix_3x3  &chi,
                                          const Vector_9    &PK2,           const Vector_9    &SIGMA,      const Vector_27    &M,
                                          const Matrix_9x9  &DPK2Dgrad_u,   const Matrix_9x9  &DPK2Dphi,   const Matrix_9x27  &DPK2Dgrad_phi,
                                          const Matrix_9x9  &DSIGMADgrad_u, const Matrix_9x9  &DSIGMADphi, const Matrix_9x27  &DSIGMADgrad_phi,
                                          const Matrix_27x9 &DMDgrad_u,     const Matrix_27x9 &DMDphi,     const Matrix_27x27 &DMDgrad_phi,
                                          Matrix_9x12 &DcintDU){
        /*!==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple.
        
        Total-Lagrangian formulation
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
        
        //Compute DPK2DU;
        Matrix_9x12 DPK2DU;
        DPK2DU = (DPK2Dgrad_u*dgrad_udU + DPK2Dphi*dphidU + DPK2Dgrad_phi*dgrad_phidU);
        
        //Compute DSIGMADU
        Matrix_9x12 DSIGMADU;
        DSIGMADU = (DSIGMADgrad_u*dgrad_udU + DSIGMADphi*dphidU + DSIGMADgrad_phi*dgrad_phidU);
        
        //Compute DmDU
        Matrix_27x12 DMDU;
        DMDU = (DMDgrad_u*dgrad_udU + DMDphi*dphidU + DMDgrad_phi*dgrad_phidU);
        
        double tmp;
        double F_iI;
        double F_jJ;
        double chi_iI;
        double M_KJI;
        double tmp2;
        double tmp3;
        //int    A;
        int    Ihat;
        int    Jhat;
        int    Khat;
        int    Lhat;
        int    Mhat;
        
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
                                      
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
        
        //Note: dFdU = dgrad_udU and dchidU = dphidU
        
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                Ihat = sot_to_voigt_map[i][j];
                for (int A=0; A<12; A++){
                    tmp = 0;
                    
                    for (int I=0; I<3; I++){
                        Jhat = sot_to_voigt_map[i][I];
                        F_iI = F(i,I);
                        for (int J=0; J<3; J++){
                            Khat = sot_to_voigt_map[I][J];
                            Lhat = sot_to_voigt_map[j][J];
                            tmp2 = PK2(Khat) - SIGMA(Khat);
                            tmp += N*( dgrad_udU(Jhat,A)*tmp2*F(j,J)
                                      +F_iI*(DPK2DU(Khat,A) - DSIGMADU(Khat,A))*F(j,J)
                                      +F_iI*tmp2*dgrad_udU(Lhat,A));
                        }
                    }
                    
                    for (int I=0; I<3; I++){
                        Jhat = sot_to_voigt_map[i][I];
                        chi_iI = chi(i,I);
                        tmp3 = dphidU(Jhat,A);
                        for (int J=0; J<3; J++){
                            Lhat = sot_to_voigt_map[j][J];
                            F_jJ = F(j,J);
                            for (int K=0; K<3; K++){
                                Mhat = tot_to_voigt_map[K][J][I];
                                M_KJI = M(Mhat);
                                tmp -= dNdX[K]*( dgrad_udU(Lhat,A) * chi_iI * M_KJI
                                                +F_jJ * chi_iI * DMDU(Mhat,A)
                                                +F_jJ * M_KJI * tmp3);
                            }
                        }
                    }
                                
                    DcintDU(Ihat,A) = tmp;
                }
            }
        }
        return;
    }
    
    void compute_internal_couple_jacobian(const int &component_i, const int &component_j, const int &dof_num,
                                          const double &N, const double (&dNdX)[3], const double &eta, const double (&detadx)[3],
                                          const Matrix_3x3  &F,             const Matrix_3x3  &chi,
                                          const Vector_9    &PK2,           const Vector_9    &SIGMA,      const Vector_27    &M,
                                          const Matrix_9x9  &DPK2Dgrad_u,   const Matrix_9x9  &DPK2Dphi,   const Matrix_9x27  &DPK2Dgrad_phi,
                                          const Matrix_9x9  &DSIGMADgrad_u, const Matrix_9x9  &DSIGMADphi, const Matrix_9x27  &DSIGMADgrad_phi,
                                          const Matrix_27x9 &DMDgrad_u,     const Matrix_27x9 &DMDphi,     const Matrix_27x27 &DMDgrad_phi,
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
        
        Total-Lagrangian formulation.
        */
        
        //Identify the voigt notation component
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};

        Vector_9  dgrad_udU;
        Vector_9  dphidU;
        Vector_27 dgrad_phidU;
        
        construct_dgrad_udU(detadx, dof_num, dgrad_udU);
        construct_dphidU(eta, dof_num, dphidU);
        construct_dgrad_phidU(detadx, dof_num, dgrad_phidU);
        
        //Compute DPK2DU;
        Vector_9 DPK2DU;
        DPK2DU = (  DPK2Dgrad_u*dgrad_udU
                  + DPK2Dphi*dphidU
                  + DPK2Dgrad_phi*dgrad_phidU);
        
        //Compute DSIGMADU
        Vector_9 DSIGMADU;
        DSIGMADU = (  DSIGMADgrad_u*dgrad_udU
                    + DSIGMADphi*dphidU
                    + DSIGMADgrad_phi*dgrad_phidU);
        
        //Compute DmDU
        Vector_27 DMDU;
        DMDU = (  DMDgrad_u*dgrad_udU
                + DMDphi*dphidU
                + DMDgrad_phi*dgrad_phidU);
        
        //Temporary variables

        Vector_9 PK2mSIGMA;
        PK2mSIGMA = PK2 - SIGMA;

        Vector_9 DPK2DUmDSIGMADU;
        DPK2DUmDSIGMADU = DPK2DU - DSIGMADU;

        double F_iI;
        double F_jJ;
        double chi_iI;
        double M_KJI;

        double dFdU_iIA;
        double dFdU_jJA;
        double dchidU_iIA;
        double PK2mSIGMA_IJ;

        //double tmp2;
        //double tmp3;
        int    Jhat;
        int    Khat;
        int    Lhat;
        int    Mhat;
                                      
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
        
        //Note: dFdU = dgrad_udU and dchidU = dphidU
        
        
        
        
//        for (int i=0; i<3; i++){
//            for (int j=0; j<3; j++){
//                Ihat = sot_to_voigt_map[i][j];
//                for (int A=0; A<12; A++){
        DcintDU_ijA = 0;
        
        for (int I=0; I<3; I++){
            Jhat = sot_to_voigt_map[component_i][I];
            F_iI = F(component_i,I);
            dFdU_iIA = dgrad_udU(Jhat);
            for (int J=0; J<3; J++){
                Khat = sot_to_voigt_map[I][J];
                Lhat = sot_to_voigt_map[component_j][J];
                PK2mSIGMA_IJ = PK2mSIGMA(Khat);
                DcintDU_ijA += N*( dFdU_iIA*PK2mSIGMA_IJ*F(component_j,J)
                                  +F_iI*DPK2DUmDSIGMADU(Khat)*F(component_j,J)
                                  +F_iI*PK2mSIGMA_IJ*dgrad_udU(Lhat));
            }
        }
        
        for (int I=0; I<3; I++){
            Jhat = sot_to_voigt_map[component_i][I];
            chi_iI = chi(component_i,I);
            dchidU_iIA = dphidU(Jhat);
            for (int J=0; J<3; J++){
                Lhat = sot_to_voigt_map[component_j][J];
                dFdU_jJA = dgrad_udU(Lhat);
                F_jJ = F(component_j,J);
                for (int K=0; K<3; K++){
                    Mhat = tot_to_voigt_map[K][J][I];
                    M_KJI = M(Mhat);
                    DcintDU_ijA -= dNdX[K]*( dFdU_jJA * chi_iI * M_KJI
                                            +    F_jJ * chi_iI * DMDU(Mhat)
                                            +    F_jJ * M_KJI  * dchidU_iIA);
                }
            }
        }

        return;
    }
    
    void compute_internal_couple_jacobian(const int &component_i, const int &component_j, const int &dof_num,
                                          const double &N, const double (&dNdX)[3], const double &eta, const double (&detadX)[3],
                                          const std::vector<std::vector<double>> &F,           const std::vector<std::vector<double>> &chi,
                                          const std::vector<double>              &PK2,           const std::vector<double>              &SIGMA,      const std::vector<double>         &M,
                                          const std::vector<std::vector<double>> &DPK2Dgrad_u, const std::vector<std::vector<double>> &DPK2Dphi, const std::vector<std::vector<double>> &DPK2Dgrad_phi,
                                          const std::vector<std::vector<double>> &DSIGMADgrad_u,      const std::vector<std::vector<double>> &DSIGMADphi,      const std::vector<std::vector<double>> &DSIGMADgrad_phi,
                                          const std::vector<std::vector<double>> &DMDgrad_u,      const std::vector<std::vector<double>> &DMDphi,      const std::vector<std::vector<double>> &DMDgrad_phi,
                                          double &DcintDU_ijA){
        /*!==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple.
        
        Total-Lagrangian formulation.
        */

        //Identify the voigt notation component
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};

        std::vector<double>  dgrad_udU(9,0);
        std::vector<double>  dphidU(9,0);
        std::vector<double> dgrad_phidU(27,0);

        std::vector<int> non_zero_terms_dgrad_udU(3,0);
        std::vector<int> non_zero_terms_dgrad_phidU(3,0);
        
        construct_dgrad_udU(detadX, dof_num, non_zero_terms_dgrad_udU, dgrad_udU);
        construct_dphidU(eta, dof_num, dphidU);
        construct_dgrad_phidU(detadX, dof_num, non_zero_terms_dgrad_phidU, dgrad_phidU);
        
        //Compute DPK2DU;
        std::vector<double> DPK2DU(9,0);
        std::vector<double> DSIGMADU(9,0);
        std::vector<double> DMDU(27,0);

        double tmp1;
        double tmp2;

        int nzt;

        if(dof_num<3){

            for (int i=0; i<9; i++){
                tmp1 = 0;
                tmp2 = 0;
                for (int j=0; j<3; j++){
                    nzt   =   non_zero_terms_dgrad_udU[j];
                    tmp1 +=   DPK2Dgrad_u[i][nzt]*dgrad_udU[nzt];

                    tmp2 +=   DSIGMADgrad_u[i][nzt]*dgrad_udU[nzt];
                }

                DPK2DU[i]   = tmp1;
                DSIGMADU[i] = tmp2;
            }

            for (int i=0; i<27; i++){
                tmp1 = 0;
                for (int j=0; j<3; j++){
                    nzt   =   non_zero_terms_dgrad_udU[j];
                    tmp1 +=   DMDgrad_u[i][nzt]*dgrad_udU[nzt];
                }
                DMDU[i] = tmp1;
           }

        }
        else{

            for (int i=0; i<9; i++){
                tmp1 = DPK2Dphi[i][dof_num-3]*eta;//0;
                tmp2 = DSIGMADphi[i][dof_num-3]*eta;//0;
//                for (int j=0; j<9; j++){
//                    tmp1 += DPK2Dphi[i][j]*dphidU[j];
//
//                    tmp2 += DSIGMADphi[i][j]*dphidU[j];
//                }

                for (int j=0; j<3; j++){
                    nzt   = non_zero_terms_dgrad_phidU[j];
                    tmp1 += DPK2Dgrad_phi[i][nzt]*dgrad_phidU[nzt];
                    tmp2 += DSIGMADgrad_phi[i][nzt]*dgrad_phidU[nzt];
                }

                DPK2DU[i]   = tmp1;
                DSIGMADU[i] = tmp2;
            }
        
            for (int i=0; i<27; i++){
                tmp1 = DMDphi[i][dof_num-3]*eta;//0;
//                for (int j=0; j<9; j++){
//                    tmp1 += DMDphi[i][j]*dphidU[j];
//                }
                DMDU[i] = tmp1;
           }

           for (int i=0; i<27; i++){
                tmp1 = 0;
                for (int j=0; j<3; j++){
                    nzt   = non_zero_terms_dgrad_phidU[j];
                    tmp1 += DMDgrad_phi[i][nzt]*dgrad_phidU[nzt];
                }

                DMDU[i]  += tmp1;
            }
        }
        
        //Temporary variables

//        std::vector<double> PK2mSIGMA(9,0);
//        for (int i=0; i<9; i++){PK2mSIGMA[i] = PK2[i] - SIGMA[i];}

        std::vector<double> DPK2DUmDSIGMADU(9,0);
        for (int i=0; i<9; i++){DPK2DUmDSIGMADU[i] = DPK2DU[i] - DSIGMADU[i];}

        double F_iI;
        double F_jJ;
        double chi_iI;
        double M_KJI;

        double dFdU_iIA;
        double dFdU_jJA;
        double dchidU_iIA;
        double PK2mSIGMA_IJ;

        //double tmp3;
        int    Jhat;
        int    Khat;
        int    Lhat;
        int    Mhat;
                                      
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
        
        //Note: dFdU = dgrad_udU and dchidU = dphidU
        //
        DcintDU_ijA = 0;
        
        for (int I=0; I<3; I++){
            Jhat = sot_to_voigt_map[component_i][I];
            F_iI = F[component_i][I];
            dFdU_iIA = dgrad_udU[Jhat];
            for (int J=0; J<3; J++){
                Khat = sot_to_voigt_map[I][J];
                Lhat = sot_to_voigt_map[component_j][J];
                PK2mSIGMA_IJ = PK2[Khat] - SIGMA[Khat];//PK2mSIGMA[Khat];
                DcintDU_ijA += N*( dFdU_iIA*PK2mSIGMA_IJ*F[component_j][J]
                                  +F_iI*DPK2DUmDSIGMADU[Khat]*F[component_j][J]
                                  +F_iI*PK2mSIGMA_IJ*dgrad_udU[Lhat]);
            }
        }
        
        for (int I=0; I<3; I++){
            Jhat = sot_to_voigt_map[component_i][I];
            chi_iI = chi[component_i][I];
            dchidU_iIA = dphidU[Jhat];
            for (int J=0; J<3; J++){
                Lhat = sot_to_voigt_map[component_j][J];
                dFdU_jJA = dgrad_udU[Lhat];
                F_jJ = F[component_j][J];
                for (int K=0; K<3; K++){
                    Mhat = tot_to_voigt_map[K][J][I];
                    M_KJI = M[Mhat];
                    DcintDU_ijA -= dNdX[K]*( dFdU_jJA * chi_iI * M_KJI
                                            +    F_jJ * chi_iI * DMDU[Mhat]
                                            +    F_jJ * M_KJI  * dchidU_iIA);
                }
            }
        }

        return;
    }

    void compute_internal_couple_jacobian(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const std::vector<std::vector<double>> &F,           const std::vector<std::vector<double>> &chi,
                                          const std::vector<double>              &PK2,           const std::vector<double>              &SIGMA,      const std::vector<double>         &M,
                                          const std::vector<std::vector<double>> &DPK2Dgrad_u, const std::vector<std::vector<double>> &DPK2Dphi, const std::vector<std::vector<double>> &DPK2Dgrad_phi,
                                          const std::vector<std::vector<double>> &DSIGMADgrad_u,      const std::vector<std::vector<double>> &DSIGMADphi,      const std::vector<std::vector<double>> &DSIGMADgrad_phi,
                                          const std::vector<std::vector<double>> &DMDgrad_u,      const std::vector<std::vector<double>> &DMDphi,      const std::vector<std::vector<double>> &DMDgrad_phi,
                                          std::vector<std::vector<double>> &DcintDU){

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
        
        Total-Lagrangian formulation.
        */
        Matrix_3x3   _F;
        Matrix_3x3   _chi;
        
        Vector_9     _PK2;
        Vector_9     _SIGMA;
        Vector_27    _M;
        
        Matrix_9x9   _DPK2Dgrad_u;
        Matrix_9x9   _DPK2Dphi;
        Matrix_9x27  _DPK2Dgrad_phi;

        Matrix_9x9   _DSIGMADgrad_u;
        Matrix_9x9   _DSIGMADphi;
        Matrix_9x27  _DSIGMADgrad_phi;
        
        Matrix_27x9  _DMDgrad_u;
        Matrix_27x9  _DMDphi;
        Matrix_27x27 _DMDgrad_phi;

        Matrix_9x12  _DcintDU;
        
        map_vector_to_eigen(F,             _F);
        map_vector_to_eigen(chi,           _chi);
        
        map_vector_to_eigen(PK2,           _PK2);
        map_vector_to_eigen(SIGMA,         _SIGMA);
        map_vector_to_eigen(M,             _M);
        
        map_vector_to_eigen(DPK2Dgrad_u,   _DPK2Dgrad_u);
        map_vector_to_eigen(DPK2Dphi,      _DPK2Dphi);
        map_vector_to_eigen(DPK2Dgrad_phi, _DPK2Dgrad_phi);

        map_vector_to_eigen(DSIGMADgrad_u,   _DSIGMADgrad_u);
        map_vector_to_eigen(DSIGMADphi,      _DSIGMADphi);
        map_vector_to_eigen(DSIGMADgrad_phi, _DSIGMADgrad_phi);

        map_vector_to_eigen(DMDgrad_u,   _DMDgrad_u);
        map_vector_to_eigen(DMDphi,      _DMDphi);
        map_vector_to_eigen(DMDgrad_phi, _DMDgrad_phi);
        
        compute_internal_couple_jacobian( N, dNdx, eta, detadx,
                                         _F,             _chi,
                                         _PK2,           _SIGMA,      _M,
                                         _DPK2Dgrad_u,   _DPK2Dphi,   _DPK2Dgrad_phi,
                                         _DSIGMADgrad_u, _DSIGMADphi, _DSIGMADgrad_phi,
                                         _DMDgrad_u,     _DMDphi,     _DMDgrad_phi,
                                         _DcintDU);

        map_eigen_to_vector(_DcintDU,DcintDU);
                                         
        return;
    }

    /*
    =============================================
    |    Jacobians for current configuration    |
    =============================================
    */
    void compute_internal_force_jacobian_current(const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                         const Vector_9 &cauchy, const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                         Matrix_3x12 &DfintDU){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force in the 
        current configuration.

        */

        //Compute the deformation gradient (grad_u is assumed to be w.r.t. the current configuration)
        Matrix_3x3 grad_u_mat;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                grad_u_mat(i,j) = grad_u[i][j];
            }
        }
        Matrix_3x3 F = (Matrix_3x3::Identity() - grad_u_mat).inverse();
        
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
//        double detadX[3] = {0,0,0};
        for (int I=0; I<3; I++){
            for (int i=0; i<3; i++){
                dNdX[I] += dNdx[i]*F(i,I);
//                detadX[I] += detadx[i]*F(i,I);
            }
        }
        //Compute the inverse deformation gradient
        Matrix_3x3 Finv;
        Finv = F.inverse();

        construct_dgrad_udU(Finv, detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(grad_phi, Finv, detadx, dgrad_phidU);

        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute the divergence of DcauchyDU
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        //Note: dfintdU = -D/DU(N,j sigma_ji) = -N,J (DFinvdU_JjIhat sigma_ji + Finv_Jj DcauchyDU jiIhat) = -N,J (-Dgrad_udU_JjIhat sigma_ji + Finv_Jj DsigmaDU_jiIhat) 

        double tmp;
        //int Jhat;
        int Khat;

        for (int i=0; i<3; i++){
            for (int Ihat=0; Ihat<12; Ihat++){
                tmp = 0;
                for (int J=0; J<3; J++){
                    for (int j=0; j<3; j++){
                        //Jhat = sot_to_voigt_map[J][j];
                        Khat = sot_to_voigt_map[j][i];
                        tmp += -dNdX[J]*(Finv(J,j) * DcauchyDU(Khat,Ihat));// - dgrad_udU(Jhat,Ihat)*cauchy(Khat));

                    }
                }

//                for (int j=0; j<3; j++){
//                    Khat = sot_to_voigt_map[j][i];
//                    tmp += -dNdx[j]*DcauchyDU(Khat,Ihat);
//                }
//
//                if( Ihat < 3){
//                    for (int j=0; j<3; j++){
//                        tmp += dNdx[Ihat]*detadx[j]*cauchy(Khat);
//                    }
//                }
                
                DfintDU(i,Ihat) = tmp;
            }
        }
        return;
    }

    void compute_internal_force_jacobian_current(const int &component,   const int &dof_num,
                                         const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                         const Vector_9 &cauchy, const Matrix_9x9 &DcauchyDgrad_u, const Matrix_9x9 &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                         double &DfintDU_iA){
        /*!=========================================
        |    compute_internal_force_jacobian    |
        =========================================
        
        Compute the jacobian of the internal force in the 
        current configuration.

        */

        //Compute the deformation gradient (grad_u is assumed to be w.r.t. the current configuration)
        Matrix_3x3 grad_u_mat;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                grad_u_mat(i,j) = grad_u[i][j];
            }
        }
        Matrix_3x3 F = (Matrix_3x3::Identity() - grad_u_mat).inverse();
        
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
//        double detadX[3] = {0,0,0};
        for (int I=0; I<3; I++){
            for (int i=0; i<3; i++){
                dNdX[I] += dNdx[i]*F(i,I);
//                detadX[I] += detadx[i]*F(i,I);
            }
        }
        //Compute the inverse deformation gradient
        Matrix_3x3 Finv;
        Finv = F.inverse();

        construct_dgrad_udU(Finv,detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(grad_phi, Finv, detadx, dgrad_phidU);

        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute the divergence of DcauchyDU
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        
        //Note: dfintdU = -D/DU(N,j sigma_ji) = -N,J (DFinvdU_JjIhat sigma_ji + Finv_Jj DcauchyDU jiIhat) = -N,J (-Dgrad_udU_JjIhat sigma_ji + Finv_Jj DsigmaDU_jiIhat) 

        //int Jhat;
        int Khat;

        DfintDU_iA = 0;

        for (int J=0; J<3; J++){
            for (int j=0; j<3; j++){
                //Jhat = sot_to_voigt_map[J][j];
                Khat = sot_to_voigt_map[j][component];
                DfintDU_iA += -dNdX[J]*(Finv(J,j) * DcauchyDU(Khat,dof_num));// - dgrad_udU(Jhat,dof_num)*cauchy(Khat));
            }
        }

        return;
    }

    void compute_internal_force_jacobian_current(const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
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
        Vector_9    _cauchy;
        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;
        Matrix_3x12 _DfintDU;

        map_vector_to_eigen(          cauchy,           _cauchy);
        map_vector_to_eigen(  DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(     DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);

        compute_internal_force_jacobian_current(       N,            dNdx,          eta,            detadx,
                                                  grad_u,        grad_phi,
                                                 _cauchy, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                                _DfintDU);

        map_eigen_to_vector(_DfintDU,DfintDU);        
        return;
    }

    void compute_internal_force_jacobian_current(const int &component,   const int &dof_num,
                                         const double &N,        const double(&dNdx)[3],           const double &eta,             const double(&detadx)[3],
                                         const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
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
        Vector_9    _cauchy;
        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;

        map_vector_to_eigen(          cauchy,           _cauchy);
        map_vector_to_eigen(  DcauchyDgrad_u,   _DcauchyDgrad_u);
        map_vector_to_eigen(     DcauchyDphi,      _DcauchyDphi);
        map_vector_to_eigen(DcauchyDgrad_phi, _DcauchyDgrad_phi);

        compute_internal_force_jacobian_current(  component,         dof_num,
                                                          N,            dNdx,          eta,            detadx,
                                                     grad_u,        grad_phi,
                                                    _cauchy, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                                 DfintDU_iA);

        return;
    }
    
    void compute_internal_couple_jacobian_current(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          Matrix_9x12 &DcintDU){
        /*!
        ==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple for the 
        internal couple in the current configuration.        
        */
        
        DcintDU = Matrix_9x12::Zero();
        
        //Compute the deformation gradient (grad_u is assumed to be w.r.t. the current configuration)
        Matrix_3x3 grad_u_mat;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                grad_u_mat(i,j) = grad_u[i][j];
            }
        }
        Matrix_3x3 F = (Matrix_3x3::Identity() - grad_u_mat).inverse();
        
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
        
        //Compute the inverse deformation gradient
        Matrix_3x3 Finv;
        Finv = F.inverse();

        construct_dgrad_udU(Finv, detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(grad_phi, Finv, detadx, dgrad_phidU);
        
        //Compute DcauchyDU;
        Matrix_9x12 DcauchyDU;
        DcauchyDU = (DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU);
        
        //Compute DsDU
        Matrix_9x12 DsDU;
        DsDU = (DsDgrad_u*dgrad_udU + DsDphi*dphidU + DsDgrad_phi*dgrad_phidU);
        
        //Compute DmDU
        Matrix_27x12 DmDU;
        DmDU = (DmDgrad_u*dgrad_udU + DmDphi*dphidU + DmDgrad_phi*dgrad_phidU);
        
        //Compute N_,k DmDU_kji,alpha
        
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
        
        Matrix_9x12 DdivmDU;
        DdivmDU = Matrix_9x12::Zero();
        
        int Ihat;
        int Khat;
        double tmp;
        
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                Ihat = sot_to_voigt_map[i][j];
                
                for (int Jhat=0; Jhat<12; Jhat++){
                    
                    tmp = 0;
                    
                    for (int k=0; k<3; k++){
                        Khat = tot_to_voigt_map[k][j][i];
                        tmp += dNdx[k]*DmDU(Khat,Jhat);
                        
                    }
                    
                    DdivmDU(Ihat,Jhat) = tmp;
                }
            }
        }
        
        DcintDU = -(N*(DcauchyDU - DsDU) - DdivmDU);
        return;
    }
    
    void compute_internal_couple_jacobian_current(const int &component_i, const int &component_j, const int &dof_num,
                                          const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const Matrix_9x9  &DcauchyDgrad_u, const Matrix_9x9  &DcauchyDphi, const Matrix_9x27 &DcauchyDgrad_phi,
                                          const Matrix_9x9  &DsDgrad_u,      const Matrix_9x9  &DsDphi,      const Matrix_9x27 &DsDgrad_phi,
                                          const Matrix_27x9 &DmDgrad_u,      const Matrix_27x9 &DmDphi,      const Matrix_27x27 &DmDgrad_phi,
                                          double &DcintDU_ijA){
        /*!
        ==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple for the 
        internal couple in the current configuration.        
        */
        
        //Compute the deformation gradient (grad_u is assumed to be w.r.t. the current configuration)
        Matrix_3x3 grad_u_mat;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                grad_u_mat(i,j) = grad_u[i][j];
            }
        }
        Matrix_3x3 F = (Matrix_3x3::Identity() - grad_u_mat).inverse();
        
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
        
        //Compute the inverse deformation gradient
        Matrix_3x3 Finv;
        Finv = F.inverse();

        construct_dgrad_udU(Finv, detadx, dgrad_udU);
        construct_dphidU(eta, dphidU);
        construct_dgrad_phidU(grad_phi, Finv, detadx, dgrad_phidU);
        
        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};
                                      
        int Ihat;
        Ihat = sot_to_voigt_map[component_i][component_j];
        
        DcintDU_ijA  = -N*(  DcauchyDgrad_u.row(Ihat).dot(dgrad_udU.col(dof_num))
                           + DcauchyDphi.row(Ihat).dot(dphidU.col(dof_num))
                           + DcauchyDgrad_phi.row(Ihat).dot(dgrad_phidU.col(dof_num)));
        DcintDU_ijA +=  N*(  DsDgrad_u.row(Ihat).dot(dgrad_udU.col(dof_num))
                           + DsDphi.row(Ihat).dot(dphidU.col(dof_num))
                           + DsDgrad_phi.row(Ihat).dot(dgrad_phidU.col(dof_num)));
        
        //Compute DmDU
        Matrix_27x12 DmDU;
        DmDU = (DmDgrad_u*dgrad_udU + DmDphi*dphidU + DmDgrad_phi*dgrad_phidU);
        
        //Compute N_,k DmDU_kji,alpha
        
        
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
        
        for (int k=0; k<3; k++){
            
            Ihat = tot_to_voigt_map[k][component_j][component_i];
            DcintDU_ijA += dNdx[k]*(  DmDgrad_u.row(Ihat).dot(dgrad_udU.col(dof_num))
                                    + DmDphi.row(Ihat).dot(dphidU.col(dof_num))
                                    + DmDgrad_phi.row(Ihat).dot(dgrad_phidU.col(dof_num)));
            
        }
        
        return;
    }
    
    void compute_internal_couple_jacobian_current(const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          std::vector<std::vector<double>> &DcintDU){
        /*!
        ==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple for the 
        internal couple in the current configuration.        
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

        compute_internal_couple_jacobian_current(N, dNdx, eta, detadx,
                                                          grad_u, grad_phi,
                                                 _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                                      _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                                                      _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                                                        _DcintDU);
        
        map_eigen_to_vector(_DcintDU, DcintDU);
        return;
    }
    
    void compute_internal_couple_jacobian_current(const int &component_i, const int &component_j, const int &dof_num,
                                          const double &N, const double (&dNdx)[3], const double &eta, const double (&detadx)[3],
                                          const double (&grad_u)[3][3], const double (&grad_phi)[9][3],
                                          const std::vector<std::vector<double>> &DcauchyDgrad_u, const std::vector<std::vector<double>> &DcauchyDphi, const std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                          const std::vector<std::vector<double>> &DsDgrad_u,      const std::vector<std::vector<double>> &DsDphi,      const std::vector<std::vector<double>> &DsDgrad_phi,
                                          const std::vector<std::vector<double>> &DmDgrad_u,      const std::vector<std::vector<double>> &DmDphi,      const std::vector<std::vector<double>> &DmDgrad_phi,
                                          double &DcintDU_ijA){
        /*!
        ==========================================
        |    compute_internal_couple_jacobian    |
        ==========================================
        
        Compute the jacobian of the internal body couple for the 
        internal couple in the current configuration.        
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

        compute_internal_couple_jacobian_current(component_i, component_j, dof_num,
                                                 N, dNdx, eta, detadx,
                                                          grad_u, grad_phi,
                                                 _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                                      _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                                                      _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                                                     DcintDU_ijA);
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

    void construct_dgrad_udU(const double (&detadx)[3], const int dof_num, Vector_9 &dgrad_udU){
        /*!
        =============================
        |    construct_dgrad_udU    |
        =============================

        Construct dgrad_udU for a single degree of freedom.

        */

        dgrad_udU = Vector_9::Zero();
        if (dof_num==0){
            dgrad_udU(0) = detadx[0];
            dgrad_udU(5) = detadx[1];
            dgrad_udU(4) = detadx[2];
        }
        else if (dof_num==1){
            dgrad_udU(8) = detadx[0];
            dgrad_udU(1) = detadx[1];
            dgrad_udU(3) = detadx[2];
        }
        else if (dof_num==2){
            dgrad_udU(7) = detadx[0];
            dgrad_udU(6) = detadx[1];
            dgrad_udU(2) = detadx[2];
        }
    }
    
    void construct_dgrad_udU(const double (&detadx)[3], const int dof_num, std::vector<double> &dgrad_udU){
        /*!
        =============================
        |    construct_dgrad_udU    |
        =============================

        Construct dgrad_udU for a single degree of freedom.

        NOTE: IT IS ASSUMED THAT dgrad_uDU IS A ZERO VECTOR OF LENGTH 9!

        */

        if (dof_num==0){
            dgrad_udU[0] = detadx[0];
            dgrad_udU[5] = detadx[1];
            dgrad_udU[4] = detadx[2];
        }
        else if (dof_num==1){
            dgrad_udU[8] = detadx[0];
            dgrad_udU[1] = detadx[1];
            dgrad_udU[3] = detadx[2];
        }
        else if (dof_num==2){
            dgrad_udU[7] = detadx[0];
            dgrad_udU[6] = detadx[1];
            dgrad_udU[2] = detadx[2];
        }
    }

    void construct_dgrad_udU(const double (&detadx)[3], const int dof_num, std::vector<int> &non_zero_terms, std::vector<double> &dgrad_udU){
        /*!
        =============================
        |    construct_dgrad_udU    |
        =============================

        Construct dgrad_udU for a single degree of freedom.

        NOTE: IT IS ASSUMED THAT dgrad_uDU IS A ZERO VECTOR OF LENGTH 9!

        */

        if (dof_num==0){
            dgrad_udU[0] = detadx[0];
            dgrad_udU[5] = detadx[1];
            dgrad_udU[4] = detadx[2];
            non_zero_terms[0] = 0;
            non_zero_terms[1] = 5;
            non_zero_terms[2] = 4;
        }
        else if (dof_num==1){
            dgrad_udU[8] = detadx[0];
            dgrad_udU[1] = detadx[1];
            dgrad_udU[3] = detadx[2];
            non_zero_terms[0] = 8;
            non_zero_terms[1] = 1;
            non_zero_terms[2] = 3;
        }
        else if (dof_num==2){
            dgrad_udU[7] = detadx[0];
            dgrad_udU[6] = detadx[1];
            dgrad_udU[2] = detadx[2];
            non_zero_terms[0] = 7;
            non_zero_terms[1] = 6;
            non_zero_terms[2] = 2;
        }
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

    void construct_dphidU(const double &eta, const int dof_num, Vector_9 &dphidU){
        /*!
        ==========================
        |    construct_dphidU    |
        ==========================

        Construct dphidU for a single degree of freedom.
        */

        dphidU = Vector_9::Zero();

        if(dof_num-3>=0){dphidU(dof_num-3) = eta;}
    }

    void construct_dphidU(const double &eta, const int dof_num, std::vector<double> &dphidU){
        /*!
        ==========================
        |    construct_dphidU    |
        ==========================

        Construct dphidU for a single degree of freedom.

        NOTE: dphidU IS ASSUMED TO BE A ZERO VECTOR OF LENGTH 9!
        */

        if(dof_num-3>=0){dphidU[dof_num-3] = eta;}
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

    void construct_dgrad_phidU(const double (&detadx)[3], const int dof_num, Vector_27 &dgrad_phidU){
        /*!
        ===============================
        |    construct_dgrad_phidU    |
        ===============================

        Construct the jacobian of phi_iI,J w.r.t. a single degree of freedom.

        */

        dgrad_phidU = Vector_27::Zero();

        if (dof_num == 3){
            dgrad_phidU( 0) = detadx[0];
            dgrad_phidU( 5) = detadx[1];
            dgrad_phidU( 4) = detadx[2];
        }
        else if (dof_num == 4){
            dgrad_phidU(17) = detadx[0];
            dgrad_phidU(10) = detadx[1];
            dgrad_phidU(12) = detadx[2];
        }

        else if (dof_num == 5){
            dgrad_phidU(25) = detadx[0];
            dgrad_phidU(24) = detadx[1];
            dgrad_phidU(20) = detadx[2];
        }

        else if (dof_num == 6){
            dgrad_phidU(16) = detadx[0];
            dgrad_phidU(15) = detadx[1];
            dgrad_phidU(11) = detadx[2];
        }

        else if (dof_num == 7){
            dgrad_phidU( 7) = detadx[0];
            dgrad_phidU( 6) = detadx[1];
            dgrad_phidU( 2) = detadx[2];
        }

        else if (dof_num == 8){
            dgrad_phidU( 8) = detadx[0];
            dgrad_phidU( 1) = detadx[1];
            dgrad_phidU( 3) = detadx[2];
        }

        else if (dof_num == 9){
            dgrad_phidU(26) = detadx[0];
            dgrad_phidU(19) = detadx[1];
            dgrad_phidU(21) = detadx[2];
        }

        else if (dof_num == 10){
            dgrad_phidU(18) = detadx[0];
            dgrad_phidU(23) = detadx[1];
            dgrad_phidU(22) = detadx[2];
        }

        else if (dof_num == 11){
            dgrad_phidU( 9) = detadx[0];
            dgrad_phidU(14) = detadx[1];
            dgrad_phidU(13) = detadx[2];

        }

    }
    
    void construct_dgrad_phidU(const double (&detadx)[3], const int dof_num, std::vector<double> &dgrad_phidU){
        /*!
        ===============================
        |    construct_dgrad_phidU    |
        ===============================

        Construct the jacobian of phi_iI,J w.r.t. a single degree of freedom.

        NOTE: dgrad_phidU IS ASSUMED TO BE A ZERO VECTOR OF LENGTH 27!

        */

        if (dof_num == 3){
            dgrad_phidU[ 0] = detadx[0];
            dgrad_phidU[ 5] = detadx[1];
            dgrad_phidU[ 4] = detadx[2];
        }
        else if (dof_num == 4){
            dgrad_phidU[17] = detadx[0];
            dgrad_phidU[10] = detadx[1];
            dgrad_phidU[12] = detadx[2];
        }

        else if (dof_num == 5){
            dgrad_phidU[25] = detadx[0];
            dgrad_phidU[24] = detadx[1];
            dgrad_phidU[20] = detadx[2];
        }

        else if (dof_num == 6){
            dgrad_phidU[16] = detadx[0];
            dgrad_phidU[15] = detadx[1];
            dgrad_phidU[11] = detadx[2];
        }

        else if (dof_num == 7){
            dgrad_phidU[ 7] = detadx[0];
            dgrad_phidU[ 6] = detadx[1];
            dgrad_phidU[ 2] = detadx[2];
        }

        else if (dof_num == 8){
            dgrad_phidU[ 8] = detadx[0];
            dgrad_phidU[ 1] = detadx[1];
            dgrad_phidU[ 3] = detadx[2];
        }

        else if (dof_num == 9){
            dgrad_phidU[26] = detadx[0];
            dgrad_phidU[19] = detadx[1];
            dgrad_phidU[21] = detadx[2];
        }

        else if (dof_num == 10){
            dgrad_phidU[18] = detadx[0];
            dgrad_phidU[23] = detadx[1];
            dgrad_phidU[22] = detadx[2];
        }

        else if (dof_num == 11){
            dgrad_phidU[ 9] = detadx[0];
            dgrad_phidU[14] = detadx[1];
            dgrad_phidU[13] = detadx[2];

        }

    }

    void construct_dgrad_phidU(const double (&detadx)[3], const int dof_num, std::vector<int> &non_zero_terms, std::vector<double> &dgrad_phidU){
        /*!
        ===============================
        |    construct_dgrad_phidU    |
        ===============================

        Construct the jacobian of phi_iI,J w.r.t. a single degree of freedom.

        NOTE: dgrad_phidU IS ASSUMED TO BE A ZERO VECTOR OF LENGTH 27!

        */

        if (dof_num == 3){
            dgrad_phidU[ 0] = detadx[0];
            dgrad_phidU[ 5] = detadx[1];
            dgrad_phidU[ 4] = detadx[2];
            non_zero_terms[0] = 0;
            non_zero_terms[1] = 5;
            non_zero_terms[2] = 4;
        }
        else if (dof_num == 4){
            dgrad_phidU[17] = detadx[0];
            dgrad_phidU[10] = detadx[1];
            dgrad_phidU[12] = detadx[2];
            non_zero_terms[0] = 17;
            non_zero_terms[1] = 10;
            non_zero_terms[2] = 12;
        }

        else if (dof_num == 5){
            dgrad_phidU[25] = detadx[0];
            dgrad_phidU[24] = detadx[1];
            dgrad_phidU[20] = detadx[2];
            non_zero_terms[0] = 25;
            non_zero_terms[1] = 24;
            non_zero_terms[2] = 20;
        }

        else if (dof_num == 6){
            dgrad_phidU[16] = detadx[0];
            dgrad_phidU[15] = detadx[1];
            dgrad_phidU[11] = detadx[2];
            non_zero_terms[0] = 16;
            non_zero_terms[1] = 15;
            non_zero_terms[2] = 11;
        }

        else if (dof_num == 7){
            dgrad_phidU[ 7] = detadx[0];
            dgrad_phidU[ 6] = detadx[1];
            dgrad_phidU[ 2] = detadx[2];
            non_zero_terms[0] = 7;
            non_zero_terms[1] = 6;
            non_zero_terms[2] = 2;
        }

        else if (dof_num == 8){
            dgrad_phidU[ 8] = detadx[0];
            dgrad_phidU[ 1] = detadx[1];
            dgrad_phidU[ 3] = detadx[2];
            non_zero_terms[0] = 8;
            non_zero_terms[1] = 1;
            non_zero_terms[2] = 3;
        }

        else if (dof_num == 9){
            dgrad_phidU[26] = detadx[0];
            dgrad_phidU[19] = detadx[1];
            dgrad_phidU[21] = detadx[2];
            non_zero_terms[0] = 26;
            non_zero_terms[1] = 19;
            non_zero_terms[2] = 21;
        }

        else if (dof_num == 10){
            dgrad_phidU[18] = detadx[0];
            dgrad_phidU[23] = detadx[1];
            dgrad_phidU[22] = detadx[2];
            non_zero_terms[0] = 18;
            non_zero_terms[1] = 23;
            non_zero_terms[2] = 22;
        }

        else if (dof_num == 11){
            dgrad_phidU[ 9] = detadx[0];
            dgrad_phidU[14] = detadx[1];
            dgrad_phidU[13] = detadx[2];
            non_zero_terms[0] = 9;
            non_zero_terms[1] = 14;
            non_zero_terms[2] = 13;

        }

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
    
    void construct_dgrad_phidU(const double (&grad_phi)[9][3], const Matrix_3x3 &Finv, const double (&detadx)[3], Matrix_27x12 &dgrad_phidU){
        /*!
        ===============================
        |    construct_dgrad_phidU    |
        ===============================
        
        Construct the derivative of the gradient of grad_phi in the current configuration.
        */
        
        dgrad_phidU = Matrix_27x12::Zero();

        int Ihat;
        int Khat;

        double tmp;
        double tmp1;

        int sot_to_voigt_map[3][3] = {{0,5,4},
                                      {8,1,3},
                                      {7,6,2}};

        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);

        for (int i=0; i<3; i++){
            for (int I=0; I<3; I++){
                Ihat = sot_to_voigt_map[i][I];
                for (int j=0; j<3; j++){
                    Khat = tot_to_voigt_map[i][I][j];
                    tmp1 = grad_phi[Ihat][j];
                    for (int Jhat=0; Jhat<12; Jhat++){
                        tmp = 0;
                        Khat = tot_to_voigt_map[i][I][j];

                        if (Jhat<3){
                            tmp = -detadx[Jhat]*tmp1;
                        }

                        if (Ihat==(Jhat-3)){tmp += detadx[j];}

                        dgrad_phidU(Khat,Jhat) = tmp;
                    }
                }
            }
        }
        return;
    }
}
