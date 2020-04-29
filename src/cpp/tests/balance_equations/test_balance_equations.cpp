/*!============================================================================
   |                                                                          |
   |                       test_balance_equations.cpp                         |
   |                                                                          |
   ----------------------------------------------------------------------------
   | The unit test file for balance_equations.h/cpp. This file tests the      |
   | classes and functions defined in balance_equations.h/cpp.                |
   |                                                                          |
   | Generated files:                                                         |
   |    results.tex:  A LaTeX file which contains the results as they will be |
   |                  included in the generated report.                       |
   ============================================================================
   | Dependencies:                                                            |
   | Eigen:  An implementation of various matrix commands. Available at       |
   |         eigen.tuxfamily.org                                              |
   ============================================================================*/

#include<iostream>
#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<balance_equations.h>

typedef balance_equations::variableType variableType;
typedef balance_equations::variableVector variableVector;
typedef balance_equations::variableMatrix variableMatrix;

int test_compute_internal_force( std::ofstream &results ){
    /*!
     * Test the computation of the internal force.
     *
     * :param std::ofstream &results: The output file
     */

    variableVector F = { 1, 2, 3,
                         4, 5, 6,
                         7, 8, 9 };

    variableVector PK2 = { 10, 11, 12,
                           13, 14, 15,
                           16, 17, 18 };

    double dNdX[ 3 ] = { 0.1, 0.2, 0.3 };

    variableVector answer = { -55.2, -136.2, -217.2 };

    double result[ 3 ];

    balance_equations::compute_internal_force( dNdX, F, PK2, result );

    for ( unsigned int i = 0; i < 3; i++ ){
        if ( !vectorTools::fuzzyEquals( result[ i ], answer[ i ] ) ){
            results << "test_compute_internal_force (test 1) & False\n";
            return 1;
        }

        double temp;
        balance_equations::compute_internal_force( i, dNdX, F, PK2, temp );

        if ( !vectorTools::fuzzyEquals( answer[ i ], temp ) ){
            results << "test_compute_internal_force (test 2) & False\n";
            return 1;
        }
    }

    results << "test_compute_internal_force & True\n";
    return 0;
}

int test_compute_body_force( std::ofstream &results ){
    /*!
     * Test the computation of the body force term.
     *
     * :param std::ofstream &results: The output file.
     */

    double b[ 3 ] = { 1, 2, 3 };
    variableType density = 4.15;
    variableType N = 0.271;

    double answer[ 3 ] = { 1.12465, 2.2493 , 3.37395 };
    double result[ 3 ];

    balance_equations::compute_body_force( N, density, b, result );

    for ( unsigned int i = 0; i < 3; i++ ){
        if ( !vectorTools::fuzzyEquals( result[ i ], answer[ i ] ) ){
            results << "test_compute_body_force (test 1) & False\n";
            return 1;
        }

        double temp;
        balance_equations::compute_body_force( i, N, density, b, temp );

        if ( !vectorTools::fuzzyEquals( answer[ i ], temp ) ){
            results << "test_compute_body_force (test 2) & False\n";
            return 1;
        }
    }

    results << "test_compute_body_force & True\n";
    return 0;
}

int test_compute_inertial_force( std::ofstream &results ){
    /*!
     * Test the computation of the intertial force term
     *
     * :param std::ofstream &results: The output file.
     */

    variableType N = 0.1262;
    variableType density = 2.137;
    double a[ 3 ] = { .261, .781, .512 };

    double answer[ 3 ] = { -0.07038893, -0.21062742, -0.13808097 };
    double result[ 3 ];

    balance_equations::compute_inertial_force( N, density, a, result );

    for ( unsigned int i = 0; i < 3; i++ ){
        if ( !vectorTools::fuzzyEquals( result[ i ], answer[ i ] ) ){
            results << "test_compute_inertial_force (test 1) & False\n";
            return 1;
        }

        double temp;
        balance_equations::compute_inertial_force( i, N, density, a, temp );

        if ( !vectorTools::fuzzyEquals( answer[ i ], temp ) ){
            results << "test_compute_inertial_force (test 2) & False\n";
            return 1;
        }
    }

    results << "test_compute_inertial_force & True\n";
    return 0;
}

int test_compute_internal_couple( std::ofstream &results ){
    /*!
     * Test the computation of the internal couple
     *
     * :param std::ofstream &results: The output file.
     */

    variableType N = 0.14940184445609295;
    variableVector F = { 0.12460812, 0.31971128, 0.67550862,
                         0.61095383, 0.87972732, 0.30872424,
                         0.32158187, 0.25480281, 0.45570006 };

    variableVector PK2 = { 0.69090695, 0.72388584, 0.14880964,
                           0.67520596, 0.15106516, 0.77810545,
                           0.07641724, 0.09367471, 0.15905979 };

    variableVector SIGMA = { 0.0651695 , 0.52150417, 0.91873444,
                             0.5622355 , 0.50199447, 0.26729942,
                             0.89858519, 0.09043229, 0.47643382 };

    double dNdX[ 3 ] = { 0.47156478, 0.70952355, 0.91143106 };

    variableVector chi = { 0.39978213, 0.75363408, 0.36656435,
                           0.5132958 , 0.70216153, 0.28540018,
                           0.78854133, 0.1480301 , 0.67998555 };

    variableVector M = { 0.56108707, 0.96515216, 0.75422717, 0.05993671, 0.31106832,
                         0.56857958, 0.137549  , 0.42936996, 0.69791426, 0.71632842,
                         0.93923384, 0.96486191, 0.02510584, 0.08724847, 0.35031414,
                         0.34463859, 0.54019298, 0.13329819, 0.31331963, 0.56039227,
                         0.86780486, 0.86123368, 0.44166717, 0.63443824, 0.74293527,
                         0.19282243, 0.05489031 };

    double answer[ 9 ] = { -1.46993722, -2.88148157, -1.60662834,
                           -1.4602623 , -2.75116191, -1.54919786,
                           -1.61520049, -3.13204731, -1.70457928 };

    double result[ 9 ];

    balance_equations::compute_internal_couple( N, dNdX, F, chi, PK2, SIGMA, M, result );

    for ( unsigned int i = 0; i < 9; i++ ){
        if ( !vectorTools::fuzzyEquals( result[ i ], answer[ i ] ) ){
            results << "test_compute_internal_couple (test 1) & False\n";
            return 1;
        }

        double temp;
        balance_equations::compute_internal_couple( i / 3, i % 3, N, dNdX, F, chi, PK2, SIGMA, M, temp );
        if ( !vectorTools::fuzzyEquals( answer[ i ], temp ) ){
            results << "test_compute_internal_couple (test 2) & False\n";
            return 1;
        }
    }

    results << "test_compute_internal_couple & True\n";
    return 0;
}

int test_compute_body_couple( std::ofstream &results ){
    /*!
     * Test the computation of the body force couple
     *
     * :param std::ofstream &results: The output file
     */

    variableType N = 0.24498934519109938;
    variableType density = 0.5124212990408172;
    double l[ 9 ] = { 0.37607047, 0.28205641, 0.84962159,
                      0.37917927, 0.89360894, 0.17644279,
                      0.72736826, 0.53971981, 0.22072536 };
    double answer[ 9 ] = { 0.04721104, 0.04760132, 0.09131218,
                           0.03540873, 0.11218166, 0.06775522,
                           0.10665959, 0.02215023, 0.02770937 };

    double result[ 9 ];

    balance_equations::compute_body_couple( N, density, l, result );

    for ( unsigned int i = 0; i < 9; i++ ){
        if ( !vectorTools::fuzzyEquals( result[ i ], answer[ i ] ) ){
            results << "test_compute_body_couple (test 1) & False\n";
            return 1;
        }

        double temp;
        balance_equations::compute_body_couple( i / 3, i % 3, N, density, l, temp );
        if ( !vectorTools::fuzzyEquals( temp, answer[ i ] ) ){
            results << "test_compute_body_couple (test 2) & False\n";
            return 1;
        }
    }

    results << "test_compute_body_couple & True\n";
    return 0;
}

int test_compute_inertial_couple( std::ofstream &results ){
    /*!
     * Test the computation of the inertial couple term
     *
     * :param std::ofstream &results: The output file.
     */

    variableType N = 0.5619258626089889;
    variableType density = 0.9913690562264066;
    double omega[ 9 ] = { 0.25013355, 0.27169262, 0.23286954,
                          0.90426712, 0.10991883, 0.58144111,
                          0.35675061, 0.75618734, 0.02762653 };
    double answer[ 9 ] = { -0.13934337, -0.50374543, -0.19873717,
                           -0.15135341, -0.06123313, -0.42125375,
                           -0.12972601, -0.32390684, -0.01539008 };
    double result[ 9 ];

    balance_equations::compute_inertial_couple( N, density, omega, result );

    for ( unsigned int i = 0; i < 9; i++ ){
        if ( !vectorTools::fuzzyEquals( result[ i ], answer[ i ] ) ){
            results << "test_compute_inertial_couple (test 1) & False\n";
            return 1;
        }

        double temp;
        balance_equations::compute_inertial_couple( i / 3, i % 3, N, density, omega, temp );
        if ( !vectorTools::fuzzyEquals( temp, answer[ i ] ) ){
            results << "test_compute_inertial_couple (test 2) & False\n";
            return 1;
        }
    }

    results << "test_compute_inertial_couple & True\n";
    return 0;
}

int main(){
    /*!==========================
    |         main            |
    ===========================
    
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or 
    False if the test passes or fails respectively.*/

    std::ofstream results;
    //Open the results file
    results.open ("results.tex");

    //Tests of the terms in the balance of linear momentum
    test_compute_internal_force( results );
    test_compute_body_force( results );
    test_compute_inertial_force( results );

    //Tests of the terms in the balance of the first moment of momentum
    test_compute_internal_couple( results );
    test_compute_body_couple( results );
    test_compute_inertial_couple( results );

    //Close the results file
    results.close();
}
