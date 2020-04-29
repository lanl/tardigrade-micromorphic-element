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

    test_compute_internal_force( results );
    test_compute_body_force( results );
    test_compute_inertial_force( results );

    //Close the results file
    results.close();
}
