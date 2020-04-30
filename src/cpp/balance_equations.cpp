/*!============================================================
  |                                                           |
  |                   balance_equations.cpp                   |
  |                                                           |
  -------------------------------------------------------------
  | The source file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
  |                                                           |
  | Note: Compile with SKIP_ERROR_HANDLING to ignore the      |
  |       error handling which adds some overhead. This is    |
  |       only recommended if the code is known to be         |
  |       functioning correctly and performance is the        |
  |       primary concern.                                    |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#include <balance_equations.h>

namespace balance_equations{

    int compute_internal_force( const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &PK2, double ( &fint )[ 3 ]){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress.
         *
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
         *
         * Note: This has already accounted for the volume change term.
         *
         * Function returns 0 for no errors, 1 if F is the wrong size, and 2 if PK2 is the wrong size.
         *
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double ( &fint )[ 3 ]: The internal force.
         */

        int errorCode;
        for ( unsigned int i = 0; i < 3; i++ ){
            errorCode = compute_internal_force( i, dNdX, F, PK2, fint[ i ] );

#ifndef SKIP_ERROR_HANDLING
            if ( errorCode != 0 ){
                return errorCode;
            }
#endif
        }

        return 0;
    }
    
    int compute_internal_force( const unsigned int &i,
                                const double (&dNdX)[3], const variableVector &F, const variableVector &PK2, double &fint_i ){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress on a single component.
         *
         * fint_i = - N_{ , I } * PK2_{ I J } * F_{ i J }
         *
         * Function returns 0 for no errors, 1 if F is the wrong size, and 2 if PK2 is the wrong size.
         *
         * :param const unsigned int &i: The component to compute the stress on.
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double fint_i: The i'th component of the internal force.
         */

        //Assume 3D
        const unsigned int dim = 3;

#ifndef SKIP_ERROR_HANDLING
        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( PK2.size() != dim * dim ){
            return 2;
        }
#endif

        fint_i = 0;
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                fint_i -= dNdX[ I ] * PK2[ dim * I + J ] * F[ dim * i + J ];
            }
        }

        return 0;
    }

    void compute_body_force( const double &N, const double &density, const double ( &b )[ 3 ], double ( &fb )[ 3 ] ){
        /*!
         * Compute the body force
         *
         * fb_i = N * \rho * b_i
         *
         * :param const double &N: The value of the shape function
         * :param const double &density: The density in the current configuration
         * :param const double ( &b )[ 3 ]: The body force in the current configuration
         * :param const double ( &fb )[ 3 ]: The body force in the current configuration
         */

        //Assume 3D
        const unsigned int dim = 3;
        
        for ( unsigned int i = 0; i < dim; i++ ){
            compute_body_force( i, N, density, b, fb[ i ] );
        }

        return;
    }
    
    void compute_body_force( const unsigned int &i, const double &N, const double &density, const double ( &b )[ 3 ], double &fb_i ){
        /*!
         * Compute a single component of the body force
         *
         * fb_i = N * \rho * b_i
         *
         * :param const unsigned int &i: The component to compute the body force on
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
         * finertial_i = -N * \rho * a_i
         *
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const double ( &finertial )[ 3 ]: The inertial force in the current configuration
         */

        //Assume 3D
        const unsigned int dim = 3;
        
        for ( unsigned int i = 0; i < dim; i++ ){
            compute_inertial_force( i, N, density, a, finertial[ i ] );
        }

        return;
    }

    void compute_inertial_force( const unsigned int &i,
                                 const double &N, const double &density, const double ( &a )[ 3 ], double &finertial_i ){
        /*!
         * Compute the inertial force for the ith component.
         *
         * finertial_i = -N * \rho * a_i
         *
         * :param const unsigned int &i: The component to compute the inertial force on
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const double finertial_i: The inertial force in the current configuration in direction i
         */
       
        finertial_i = -N * density * a[ i ]; 

        return;
    }
    
    int compute_internal_couple( const double &N, const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &chi,
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

        int errorCode;
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                errorCode = compute_internal_couple( i, j, N, dNdX, F, chi, PK2, SIGMA, M, cint[ dim * i + j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){
                    return errorCode;
                }
#endif
            }
        }

        return 0;
    }

    int compute_internal_couple( const unsigned int &i, const unsigned int &j,
                                 const double &N, const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &chi,
                                 const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                 double &cint_ij ){
        /*!
         * Compute the internal couple at index ij defined as
         * cint_{ ij } = N F_{ iI } ( PK2_{ IJ } - SIGMA_{ IJ } ) F_{ jJ } - N_{ ,K } F_{ jJ } \chi_{ iI } M_{ KJI }
         *
         * Function returns 0 for no errors, 1 if F is the wrong size, 2 if chi is the wrong size,
         * 3 if PK2 is the wrong size, 4 if SIGMA is the wrong size, and 5 if M is the wrong size
         *
         * :param const unsigned int &i: The first index.
         * :param const unsigned int &j: The second index.
         * :param const double &N: The shape-function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference configuration
         * :param const variableVector &F: The deformation gradient
         * :param const variableVector &chi: The micro-deformation
         * :param const variableVector &PK2: The second Piola-Kirchoff stress
         * :param const variableVector &SIGMA: The symmetric micro stress in the reference configuration.
         * :param const variableVector &M: The higher order stress in the reference configuration.
         * :param double cint_ij: The ij'th internal couple
         */

        //Assume 3D
        const unsigned int dim = 3;

#ifndef SKIP_ERROR_HANDLING
        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( chi.size() != dim * dim ){
            return 2;
        }

        if ( PK2.size() != dim * dim ){
            return 3;
        }

        if ( SIGMA.size() != dim * dim ){
            return 4;
        }

        if ( M.size() != dim * dim * dim ){
            return 5;
        }
#endif

        cint_ij = 0;
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                cint_ij += N * F[ dim * i + I ] * ( PK2[ dim * I + J ] - SIGMA[ dim * I + J ] ) * F[ dim * j + J ];

                for ( unsigned int K = 0; K < dim; K++ ){
                    cint_ij -= dNdX[ K ] * F[ dim * j + J ] * chi[ dim * i + I ] * M[ dim * dim * K + dim * J + I ];
                }
            }
        }

        return 0;
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
         * :param const double ( &cb )[ 9 ]: The body force couple term.
         */
        
        //Assume 3D
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                compute_body_couple( i, j, N, density, l, cb[ dim * i + j ] );
            }
        }
        
        return;
    }

    void compute_body_couple( const unsigned int &i, const unsigned int &j,
                              const double &N, const double &density, const double ( &l )[ 9 ], double &cb_ij ){
        /*!
         * Compute the body couple term for the indices i and j
         *
         * couple_body_{ ij } = N * density * l_{ ji }
         *
         * :param const unsigned int &i: The first index.
         * :param const unsigned int &j: The second index.
         * :param const double &N: The shape-function value.
         * :param const double &density: The density in the current configuration.
         * :param const double ( %l )[ 9 ]: The body-force couple in the current configuration.
         * :param const double &cb_ij: The ij'th body force couple term.
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
                compute_inertial_couple( i, j, N, density, omega, cinertial[ dim * i + j ] );
            }
        }
        
        return;
    }
    
    void compute_inertial_couple( const unsigned int &i, const unsigned int &j,
                                  const double &N, const double &density, const double ( &omega )[ 9 ], double &cinertial_ij ){
        /*!
         * Compute the inertial couple in the current configuration for the indices i and j
         *
         * cinertial_{ ij } = -N * density * omega_{ ji }
         *
         * :param const unsigned int &i: The first index.
         * :param const unsigned int &j: The second index.
         * :param const double &N: The shape-function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &omega )[ 9 ]: The micro-inertia tensor in the current configuration.
         * :param const double &cinternal_ij: The ij'th inertial couple.
         */

        //Assume 3D
        const unsigned int dim = 3;

        cinertial_ij = - N * density * omega[ dim * j + i ];
        
        return;
    }

    int compute_internal_force_jacobian( const double &N, const double ( &dNdX )[ 3 ], const double &eta, const double ( &detadX )[ 3 ],
                                         const variableVector &F, const variableVector &PK2,
                                         const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                         const variableMatrix &DPK2Dgrad_phi,
                                         variableMatrix &DfintDU ){
        /*!
         * Compute the jacobian of the internal force
         *
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
         *
         * returns 0 if there are no errors, 1 if F is not the right size, 2 if PK2 is not the right size, 3 if DPK2Dgrad_u is
         * not the right size, 4 if DPK2Dphi is not the right size, and 5 if DPK2Dgrad_phi is not the right size.
         *
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. X in the reference configuration.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. X in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The second Piola Kirchoff stress.
         * :param const variableMatrix &DPK2Dgrad_u: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient of the 
         *     macro displacement w.r.t. X in the reference configuration.
         * :param const variableVector &DPK2Dphi: The Jacobian of the second Piola Kirchoff stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient 
         *     of the micro displacement.
         * :param variableMatrix DfintDU: The Jacobian of the internal force w.r.t. the degree of freedom vector which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         */

        //Assume 3D
        const unsigned int dim = 3;

        //Assume 12 degrees of freedom
        const unsigned int NDOF = 12;

        int errorCode;
        DfintDU = variableMatrix( dim, variableVector( NDOF, 0 ) );
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < NDOF; j++ ){
                errorCode = compute_internal_force_jacobian( i, j, N, dNdX, eta, detadX, F, PK2,
                                                             DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                                             DfintDU[ i ][ j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){
                    return errorCode;
                }
#endif
            }
        }

        return 0;
    }

    int compute_internal_force_jacobian( const unsigned int &i, const unsigned int &j,
                                         const double &N, const double ( &dNdX )[ 3 ], const double &eta, const double ( &detadX )[ 3 ],
                                         const variableVector &F, const variableVector &PK2,
                                         const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                         const variableMatrix &DPK2Dgrad_phi,
                                         variableType &DfintDU_ij ){
        /*!
         * Compute the jacobian of the internal force
         *
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
         *
         * returns 0 if there are no errors, 1 if F is not the right size, 2 if PK2 is not the right size, 3 if DPK2Dgrad_u is
         * not the right size, 4 if DPK2Dphi is not the right size, 5 if DPK2Dgrad_phi is not the right size, 6 if 
         * i is out of range, 7 if j is out of range.
         *
         * :param const unsigned int &i: The row index of the Jacobian ( the internal force vector index )
         * :param const unsigned int &j: The column index of the Jacobian ( the degree of freedom vector index )
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. X in the reference configuration.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. X in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The second Piola Kirchoff stress.
         * :param const variableMatrix &DPK2Dgrad_u: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient of the 
         *     macro displacement w.r.t. X in the reference configuration.
         * :param const variableVector &DPK2Dphi: The Jacobian of the second Piola Kirchoff stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient 
         *     of the micro displacement.
         * :param variableType DfintDU_ij: The ij'th component of the internal force jacobian w.r.t. the degree of freedom vector
         *     which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         *     where it's position in this vector is j
         */

        //Assume 3D
        const unsigned int dim = 3;

        //Assume 12 DOF
        const unsigned int NDOF = 12;

#ifndef SKIP_ERROR_HANDLING
        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( PK2.size() != dim * dim ){
            return 2;
        }

        if ( DPK2Dgrad_u.size() != dim * dim ){
            return 3;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_u.size(); i++ ){
            if ( DPK2Dgrad_u[ i ].size() != dim * dim ){
                return 3;
            }
        }

        if ( DPK2Dphi.size() != dim * dim ){
            return 4;
        }

        for ( unsigned int i = 0; i < DPK2Dphi.size(); i++ ){
            if ( DPK2Dphi[ i ].size() != dim * dim ){
                return 4;
            }
        }

        if ( DPK2Dgrad_phi.size() != dim * dim ){
            return 5;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_phi.size(); i++ ){
            if ( DPK2Dgrad_phi[ i ].size() != dim * dim * dim ){
                return 5;
            }
        }

        if ( i >= dim ){
            return 6;
        }

        if ( j >= NDOF ){
            return 7;
        }
#endif

        DfintDU_ij = 0;

        //Assemble the Jacobians w.r.t. the macro displacement
        if ( j < 3 ){
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        DfintDU_ij -= dNdX[ I ] * DPK2Dgrad_u[ dim * I + J ][ dim * j + K ] * F[ dim * i + J ] * detadX[ K ];
                    }

                    if ( i == j ){
                        DfintDU_ij -= dNdX[ I ] * PK2[ dim * I + J ] * detadX[ J ];
                    }
                }
            }
        }
        else if ( ( j < 12 ) && ( j >= 3 ) ){

            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    DfintDU_ij -= dNdX[ I ] * DPK2Dphi[ dim * I + J ][ j - 3 ] * eta * F[ dim * i + J ];

                    for ( unsigned int K = 0; K < dim; K++ ){
                        DfintDU_ij -= dNdX[ I ] * DPK2Dgrad_phi[ dim * I + J ][ dim * ( j - 3 ) + K ] * detadX[ K ] * F[ dim * i + J ];
                    }
                }
            }
        }

        return 0;
    }

    int compute_internal_couple_jacobian( const double &N, const double ( &dNdX )[ 3 ],
                                          const double &eta, const double ( &detadX )[ 3 ],
                                          const variableVector &F, const variableVector &chi,
                                          const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                          const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                          const variableMatrix &DPK2Dgrad_phi,
                                          const variableMatrix &DSIGMADgrad_u, const variableMatrix &DSIGMADphi,
                                          const variableMatrix &DSIGMADgrad_phi,
                                          const variableMatrix &DMDgrad_u, const variableMatrix &DMDphi,
                                          const variableMatrix &DMDgrad_phi,
                                          variableMatrix &DcintDU ){
        /*!
         * Compute the jacobian of the internal couple
         * cint_{ ij } = N F_{ iI } ( PK2_{ IJ } - SIGMA_{ IJ } ) F_{ jJ } - N_{ ,K } F_{ jJ } \chi_{ iI } M_{ KJI }
         *
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference coordinates.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. the reference coordinates.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &chi: The micro deformation tensor.
         * :param const variableVector &PK2: The second Piola Kirchoff stress tensor.
         * :param const variableVector &SIGMA: The symmetric micro stress tensor in the reference configuration.
         * :param const variableVector &M: The higher order stress tensor in the reference configuration.
         * :param const variableMatrix &DPK2Dgrad_u: The derivative of the PK2 stress w.r.t. the gradient of the
         *     macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DPK2Dphi: The derivative of the PK2 stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The derivative of the PK2 stress w.r.t. the gradient of the 
         *     micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADgrad_u: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADphi: The derivative of the symmetric micro stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DSIGMADgrad_phi: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDgrad_u: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDphi: The derivative of the higher order stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DMDgrad_phi: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param variableMatrix &DcintDU: The Jacobian of the internal couple stress w.r.t. the degree of freedom
         *     vector which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         */

        //Assume 3D
        const unsigned int dim = 3;

        //Assume 12 degrees of freedom
        const unsigned int NDOF = 12;

        int errorCode;

        //Compute the Jacobian terms
        DcintDU = variableMatrix( dim * dim, variableVector( NDOF, 0 ) );
        for ( unsigned int i = 0; i < dim * dim; i++ ){
            for ( unsigned int j = 0; j < NDOF; j++ ){
                errorCode = compute_internal_couple_jacobian( i, j, N, dNdX, eta, detadX, F, chi, PK2, SIGMA, M,
                                                              DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                                              DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                                              DMDgrad_u, DMDphi, DMDgrad_phi,
                                                              DcintDU[ i ][ j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){
                    return errorCode;
                }
#endif
            }
        }

        return 0;
    }

    int compute_internal_couple_jacobian( const unsigned int &i, const unsigned int &j,
                                          const double &N, const double ( &dNdX )[ 3 ],
                                          const double &eta, const double ( &detadX )[ 3 ],
                                          const variableVector &F, const variableVector &chi,
                                          const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                          const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                          const variableMatrix &DPK2Dgrad_phi,
                                          const variableMatrix &DSIGMADgrad_u, const variableMatrix &DSIGMADphi,
                                          const variableMatrix &DSIGMADgrad_phi,
                                          const variableMatrix &DMDgrad_u, const variableMatrix &DMDphi,
                                          const variableMatrix &DMDgrad_phi,
                                          variableType &DcintDU_ij ){
        /*!
         * Compute the jacobian of the internal couple
         * cint_{ ij } = N F_{ iI } ( PK2_{ IJ } - SIGMA_{ IJ } ) F_{ jJ } - N_{ ,K } F_{ jJ } \chi_{ iI } M_{ KJI }
         *
         * :param const unsigned int &i: The row index of the Jacobian ( the internal couple vector index )
         * :param const unsigned int &j: The column index of the Jacobian ( the degree of freedom vector index )
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference coordinates.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. the reference coordinates.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &chi: The micro deformation tensor.
         * :param const variableVector &PK2: The second Piola Kirchoff stress tensor.
         * :param const variableVector &SIGMA: The symmetric micro stress tensor in the reference configuration.
         * :param const variableVector &M: The higher order stress tensor in the reference configuration.
         * :param const variableMatrix &DPK2Dgrad_u: The derivative of the PK2 stress w.r.t. the gradient of the
         *     macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DPK2Dphi: The derivative of the PK2 stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The derivative of the PK2 stress w.r.t. the gradient of the 
         *     micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADgrad_u: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADphi: The derivative of the symmetric micro stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DSIGMADgrad_phi: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDgrad_u: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDphi: The derivative of the higher order stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DMDgrad_phi: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param variableMatrix &DcintDU: The Jacobian of the internal couple stress w.r.t. the degree of freedom
         *     vector which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         */

        //Assume 3D
        unsigned int dim = 3;

        //Assume 12 DOF
        unsigned int NDOF = 12;

#ifndef SKIP_ERROR_HANDING

        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( chi.size() != dim * dim ){
            return 2;
        }

        if ( PK2.size() != dim * dim ){
            return 3;
        }

        if ( SIGMA.size() != dim * dim ){
            return 4;
        }

        if ( M.size() != dim * dim * dim ){
            return 5;
        }

        if ( DPK2Dgrad_u.size() != dim * dim ){
            return 6;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_u.size(); i++ ){
            if ( DPK2Dgrad_u[ i ].size() != dim * dim ){
                return 6;
            }
        }

        if ( DPK2Dphi.size() != dim * dim ){
            return 7;
        }

        for ( unsigned int i = 0; i < DPK2Dphi.size(); i++ ){
            if ( DPK2Dphi[ i ].size() != dim * dim ){
                return 7;
            }
        }

        if ( DPK2Dgrad_phi.size() != dim * dim ){
            return 8;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_phi.size(); i++ ){
            if ( DPK2Dgrad_phi[ i ].size() != dim * dim * dim ){
                return 8;
            }
        }

        if ( DSIGMADgrad_u.size() != dim * dim ){
            return 9;
        }

        for ( unsigned int i = 0; i < DSIGMADgrad_u.size(); i++ ){
            if ( DSIGMADgrad_u[ i ].size() != dim * dim ){
                return 9;
            }
        }

        if ( DSIGMADphi.size() != dim * dim ){
            return 10;
        }

        for ( unsigned int i = 0; i < DSIGMADphi.size(); i++ ){
            if ( DSIGMADphi[ i ].size() != dim * dim ){
                return 10;
            }
        }

        if ( DSIGMADgrad_phi.size() != dim * dim ){
            return 11;
        }

        for ( unsigned int i = 0; i < DSIGMADgrad_phi.size(); i++ ){
            if ( DSIGMADgrad_phi[ i ].size() != dim * dim * dim ){
                return 11;
            }
        }

        if ( DMDgrad_u.size() != dim * dim * dim ){
            return 12;
        }

        for ( unsigned int i = 0; i < DMDgrad_u.size(); i++ ){
            if ( DMDgrad_u[ i ].size() != dim * dim ){
                return 12;
            }
        }

        if ( DMDphi.size() != dim * dim * dim ){
            return 13;
        }

        for ( unsigned int i = 0; i < DMDphi.size(); i++ ){
            if ( DMDphi[ i ].size() != dim * dim ){
                return 13;
            }
        }

        if ( DMDgrad_phi.size() != dim * dim * dim ){
            return 14;
        }

        for ( unsigned int i = 0; i < DMDgrad_phi.size(); i++ ){
            if ( DMDgrad_phi[ i ].size() != dim * dim * dim ){
                return 14;
            }
        }

#endif

        if ( j < 3 ){

            DcintDU_ij = 0;
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        DcintDU_ij += N * F[ dim * ( i / 3 ) + I ] * ( DPK2Dgrad_u[ dim * I + J ][ dim * j + K ]
                                                                     - DSIGMADgrad_u[ dim * I + J ][ dim * j + K ]
                                                                     ) * detadX[ K ] * F[ dim * ( i % 3 ) + J ];

                        for ( unsigned int L = 0; L < dim; L++ ){
                            DcintDU_ij -= dNdX[ K ] * F[ dim * ( i % 3 ) + J ] * chi[ dim * ( i / 3 ) + I ]
                                        * DMDgrad_u[ dim * dim * K + dim * J + I ][ dim * j + L ] * detadX[ L ];
                        }
                        
                        if ( ( i % 3) == j ){
                            DcintDU_ij -= dNdX[ K ] * detadX[ J ] * chi[ dim * ( i / 3 ) + I ] * M[ dim * dim * K + dim * J + I ];
                        }
                    }

                    if ( ( i / 3 ) == j ){
                        DcintDU_ij += N * detadX[ I ] * ( PK2[ dim * I + J ] - SIGMA[ dim * I + J ] ) * F[ dim * ( i % 3 ) + J ];
                    }

                    if ( ( i % 3 ) == j ){
                        DcintDU_ij += N * F[ dim * ( i / 3 ) + I ] * ( PK2[ dim * I + J ] - SIGMA[ dim * I + J ] ) * detadX[ J ];
                    }
                }
            }
        }

        return 0;
    }
}
