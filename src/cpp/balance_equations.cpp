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
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
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
                                 double &fint_i ){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress on a single component.
         *
         * fint_i = - N_{ , I } * PK2_{ I J } * F_{ i J }
         *
         * :param const int &i: The component to compute the stress on.
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double fint_i: The i'th component of the internal force.
         */

        //Assume 3D
        const unsigned int dim = 3;

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
    
    void compute_body_force( const int &i, const double &N, const double &density, const double ( &b )[ 3 ], double &fb_i ){
        /*!
         * Compute a single component of the body force
         *
         * fb_i = N * \rho * b_i
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

    void compute_inertial_force( const int &i, const double &N, const double &density, const double ( &a )[ 3 ], double &finertial_i ){
        /*!
         * Compute the inertial force for the ith component.
         *
         * finertial_i = -N * \rho * a_i
         *
         * :param const int &i: The component to compute the inertial force on
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const double finertial_i: The inertial force in the current configuration in direction i
         */
       
        finertial_i = -N * density * a[ i ]; 

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
                                  double &cint_ij ){
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
         * :param double cint_ij: The ij'th internal couple
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

    void compute_body_couple( const int &i, const int &j, const double &N, const double &density, const double ( &l )[ 9 ], 
                              double &cb_ij ){
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
    
    void compute_inertial_couple( const int &i, const int &j, const double &N, const double &density, const double ( &omega )[ 9 ], 
                                  double &cinertial_ij ){
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
         * :param const double &cinternal_ij: The ij'th inertial couple.
         */

        //Assume 3D
        const unsigned int dim = 3;

        cinertial_ij = - N * density * omega[ dim * j + i ];
        
        return;
    }
}
