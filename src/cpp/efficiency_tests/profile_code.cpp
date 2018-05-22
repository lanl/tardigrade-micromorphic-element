/*!============================================================================
   |                                                                          |
   |              test_micromorphic_linear_elasticity_voigt.cpp               |
   |                                                                          |
   ----------------------------------------------------------------------------
   | The unit test file for micromorphic_linear_elasticity_voigt.h/cpp. This  |
   | file tests the classes and functions defined in                          |
   | micromorphic_linear_elasticity_voigt.h/cpp.                              |
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
#include<fstream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
#include<finite_difference.h>
#include<deformation_measures.h>
#include<balance_equations.h>
#include<micromorphic_linear_elasticity_voigt.h>
#include<micromorphic_material_library.h>

bool compare_std_vectors(const std::vector<double> v1, const std::vector<double> v2, double tol = 1e-6){
    /*!
    =============================
    |    compare_std_vectors    |
    =============================

    Compare two standard vectors to make sure their values are close

    */

    if(v1.size() != v2.size()){return false;}
    else{
        for (int i=0; i<v1.size(); i++){
            if(1e-6<fabs(v1[i]-v2[i])){return false;}
        }
    }

    return true;
}

bool compare_matrix(const std::vector<std::vector<double>> m1, const std::vector<std::vector<double>> m2, double tol = 1e-6){
    /*!
    ========================
    |    compare_matrix    |
    ========================

    Compare matrices formed of std::vectors

    */

    if(m1.size() != m2.size()){return false;}
    else{
        for (int i=0; i<m1.size(); i++){
            if(!compare_std_vectors(m1[i],m2[i])){return false;}
        }
    }
    return true;
}

void print_vector(const std::vector<double> V){
    /*!
    ======================
    |    print_vector    |
    ======================

    Print a std::vector
    */

    for (int i=0; i<V.size(); i++){
        std::cout << V[i] << " ";
    }
    std::cout << "\n";
}

void print_matrix(const std::vector<std::vector<double>> &M){
    /*!
    ======================
    |    print_matrix    |
    ======================

    Print a matrix stored in std::vectors

    */

    for (int i=0; i<M.size(); i++){
        print_vector(M[i]);
    }
}

void define_grad_u_MOOSE(double (&grad_u)[3][3]){
    /*!
    =============================
    |    define_grad_u_MOOSE    |
    =============================

    Define the gradient of u for comparison to outputs from MOOSE.
    */

    for (int i=0; i<3; i++){

        for (int j=0; j<3; j++){

            grad_u[i][j] = 0.;

        }

    }
}

void define_phi_MOOSE(double (&phi)[9]){
    /*!
    ==========================
    |    define_phi_MOOSE    |
    ==========================

    Define phi for comparision to outputs from MOOSE.
    */

    phi[0] = -2.3;
    phi[1] =  0.124;
    phi[2] =  7.2;
    phi[3] = -8.2;
    phi[4] =  0.28;
    phi[5] =  7.21;
    phi[6] =  2.1;
    phi[7] = -9.2;
    phi[8] =  3.1;
}

void define_grad_phi_data_MOOSE(double (&grad_phi_data)[9][3]){
    /*!
    =============================
    |    define_grad_u_MOOSE    |
    =============================

    Define the gradient of u for comparison to outputs from MOOSE.
    */

    for (int i=0; i<9; i++){

        for (int j=0; j<3; j++){

            grad_phi_data[i][j] = 0.;

        }

    }
}

void voigt_3x3(const Matrix_3x3 &A, Vector_9 &v){
    /*!===================
       |    voigt_3x3    |
       ===================

       Put the 3x3 matrix in voigt notation.

    */

    v(0) = A(0,0);
    v(1) = A(1,1);
    v(2) = A(2,2);
    v(3) = A(1,2);
    v(4) = A(0,2);
    v(5) = A(0,1);
    v(6) = A(2,1);
    v(7) = A(2,0);
    v(8) = A(1,0);

}

void define_grad_u(Matrix_3x3 &grad_u){
    /*!=======================
       |    define_grad_u    |
       =======================

       Define the gradient of u to be used

    */

    grad_u << 0.69646919,  0.28613933,  0.22685145,  0.55131477,  0.71946897,
              0.42310646,  0.9807642 ,  0.68482974,  0.4809319;

    return;
}

void define_phi(Matrix_3x3 &phi){
    /*!====================
       |    define_phi    |
       ====================

       Define the values of phi to be used.

    */

    phi << 0.39211752,  0.34317802,  0.72904971,  0.43857224,  0.0596779,
           0.39804426,  0.73799541,  0.18249173,  0.17545176;
}

void define_grad_phi(Matrix_3x9 &grad_phi){
    /*!=========================
       |    define_grad_phi    |
       =========================

       Define the gradient of phi to be used.

    */

    grad_phi << -1.81005245, -1.29847083, -0.48170751, -0.75470999, -0.4763441 ,
                -1.11329654, -0.95632783, -0.9133694 , -1.99856773, -1.49033497,
                -0.53589066, -0.44965118, -0.24378446, -0.18042363, -0.91358478,
                -0.70051865, -1.09881086, -1.45762653, -3.00984646, -0.42927004,
                -0.25701678, -0.61784346, -0.60307115, -1.35759442, -1.25541793,
                -2.06541739, -1.12273603;

}

void define_PK2(Vector_9 &PK2){
    /*!==================
    |    define_PK2    |
    ====================
    
    Define the expected value of the PK2 stress.
    
    */
    
    PK2 << 59.03427247, 160.14321551, 105.82494785, 214.53178439,
           237.37311424, 249.55639324, 110.81529   ,  25.70187797,
           17.88329254;
}

void define_SIGMA(Vector_9 &SIGMA){
    /*!======================
    |    define_SIGMA    |
    ======================
    
    Define the expected value of the symmetric stress.
    
    */
    
    SIGMA << 119.79705742647155, 323.7978637411792, 216.13272323537976,
             326.0970440545271, 263.5407883419666, 267.61699074207695,
             326.0970440545271, 263.5407883419666, 267.61699074207695;
}

void define_M(Vector_27 &M){
    /*!==================
    |    define_M    |
    ==================
    
    Define the expected value of the higher-order stress.
    
    */
    
    M << 81.31981706487286, 31.798586060251207, 27.355416705438905,
         4.4605220584734795, 15.239752824275838, 22.917719613671604,
         4.761661534444574, 13.617364734132286, 20.631211107480663,
         25.446753288061686, 41.98935229144166, 17.27436660090204,
         12.970633348345526, 4.29549454562416, 21.185457632363434,
         10.684683278867416, 4.658645978608793, 26.552528036554328,
         18.65235152889546, 16.02598269360608, 29.787283731069458,
         11.17944236411268, 18.36235422268097, 4.559452481399969,
         14.133232412473218, 23.03915084858808, 4.849232622884019;
}

void define_cauchy(Vector_9 &cauchy){
    /*!=======================
    |    define_cauchy    |
    =======================

    Define the expected value of the cauchy stress.

    */

    cauchy << -48.84329314813607, -167.47094085138679, -318.99997557253636,
              -284.3961266524248, -34.71321545072979, -30.136059986911263,
              -205.35807451733308, -282.42756314420484, -216.4214058205378;

}

void define_s(Vector_9 &s){
    /*!==================
    |    define_s    |
    ==================

    Define the expected value of the symmetric stress in the 
    current configuration.

    */

    s << -99.10582206835683, -338.1967416137846, -642.1118225175197,
         -491.7929923630433, -318.2450527392582, -246.8057390623807,
         -491.7929923630433, -318.24505273925826, -246.80573906238072;

}

void define_m(Vector_27 &m){
    /*!==================
    |    define_m    |
    ==================

    Define the expected value of the higher-order stress in the 
    current configuration.

    */

    m << -17.153265392410482, -12.480771016951493, -19.147225245448556,
         -11.26628746512982,  -14.179072065222964, -15.619789973641112,
         -27.06888669123007,  -21.847082786416394, -9.37069422181913,
         -7.149691649496841,  -97.42799562868105,  -125.8651897371535,
         -139.98361743799694, -12.071479260463587, -18.679837246649882,
         -118.73866356078591, -197.65579976538393, -200.62313337544705,
         -20.39625460432565,  -109.58264684993866, -145.89097968348062,
         -132.89656976443473, -21.83692364109874,  -32.86501911273255, 
         -165.24571159030918, -239.6065692494737,  -201.04091236148963;
}

void define_parameters(double (&params)[18],bool MOOSE=false){
    /*!========================
    |    define_parameters    |
    ===========================

    Define the parameters to be used in the 
    test functions.

    */

    params[ 0] = 0.191519450379;
    params[ 1] = 0.62210877104;
    params[ 2] = 0.437727739007;
    params[ 3] = 0.785358583714;
    params[ 4] = 0.779975808119;
    params[ 5] = 0.272592605283;
    params[ 6] = 0.276464255143;
    params[ 7] = 0.801872177535;
    params[ 8] = 0.958139353684;
    params[ 9] = 0.875932634742;
    params[10] = 0.357817269958;
    params[11] = 0.500995125523;
    params[12] = 0.683462935172;
    params[13] = 0.712702026983;
    params[14] = 0.37025075479;
    params[15] = 0.561196186066;
    params[16] = 0.503083165308;
    params[17] = 0.0137684495907;


    if(MOOSE){
    params[ 0] =  8e9;
    params[ 1] =  11e9;
    params[ 2] =  2e9;
    params[ 3] =  1.538e9;
    params[ 4] = -1e9;
    params[ 5] = -1.39e9;
    params[ 6] = -2.11e9;
    params[ 7] =  0.;
    params[ 8] =  0.;
    params[ 9] =  0.;
    params[10] =  0.;
    params[11] =  0.;
    params[12] =  0.;
    params[13] =  0.769e6;
    params[14] =  0.;
    params[15] =  0.;
    params[16] =  0.;
    params[17] =  0.;
    }

    return;
}

void define_deformation_gradient(Matrix_3x3 &F){
    /*!=====================================
       |    define_deformation_gradient    |
       =====================================

       Define the deformation gradient to be used.

    */

    Matrix_3x3 grad_u;

    define_grad_u(grad_u);

    F = (Matrix_3x3::Identity() - grad_u).inverse();

    return;
}

void define_chi(Matrix_3x3 &chi){
    /*!====================
      |    define_chi    |
      ====================

      Define the micro-deformation tensor to be used.

    */

    Matrix_3x3 phi;
    define_phi(phi);

    chi = phi + Matrix_3x3::Identity();
}

void define_A(SpMat &A){
    /*!===============
    |    get_A    |
    ===============

    Get the fourth order A stiffness 
    tensor in voigt notation.

    */

    std::vector<T> tripletList;
    tripletList.reserve(21);

    tripletList.push_back(T(0,0,1.43573699245856));
    tripletList.push_back(T(0,1,0.191519450378892));
    tripletList.push_back(T(0,2,0.191519450378892));
    tripletList.push_back(T(1,0,0.191519450378892));
    tripletList.push_back(T(1,1,1.43573699245856));
    tripletList.push_back(T(1,2,0.191519450378892));
    tripletList.push_back(T(2,0,0.191519450378892));
    tripletList.push_back(T(2,1,0.191519450378892));
    tripletList.push_back(T(2,2,1.43573699245856));
    tripletList.push_back(T(3,3,0.622108771039832));
    tripletList.push_back(T(3,6,0.622108771039832));
    tripletList.push_back(T(4,4,0.622108771039832));
    tripletList.push_back(T(4,7,0.622108771039832));
    tripletList.push_back(T(5,5,0.622108771039832));
    tripletList.push_back(T(5,8,0.622108771039832));
    tripletList.push_back(T(6,3,0.622108771039832));
    tripletList.push_back(T(6,6,0.622108771039832));
    tripletList.push_back(T(7,4,0.622108771039832));
    tripletList.push_back(T(7,7,0.622108771039832));
    tripletList.push_back(T(8,5,0.622108771039832));
    tripletList.push_back(T(8,8,0.622108771039832));

    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

void define_B(SpMat &B){
    /*!===============
    |    get_B    |
    ===============

    Get the forth order B stiffness
    tensor in voigt notation.

    */

    std::vector<T> tripletList;
    tripletList.reserve(21);

    tripletList.push_back(T(0,0,0.152009058408597));
    tripletList.push_back(T(0,1,-0.347630844706655));
    tripletList.push_back(T(0,2,-0.347630844706655));
    tripletList.push_back(T(1,0,-0.347630844706655));
    tripletList.push_back(T(1,1,0.152009058408597));
    tripletList.push_back(T(1,2,-0.347630844706655));
    tripletList.push_back(T(2,0,-0.347630844706655));
    tripletList.push_back(T(2,1,-0.347630844706655));
    tripletList.push_back(T(2,2,0.152009058408597));
    tripletList.push_back(T(3,3,0.503511552975707));
    tripletList.push_back(T(3,6,-0.00387164986045507));
    tripletList.push_back(T(4,4,0.503511552975707));
    tripletList.push_back(T(4,7,-0.00387164986045507));
    tripletList.push_back(T(5,5,0.503511552975707));
    tripletList.push_back(T(5,8,-0.00387164986045507));
    tripletList.push_back(T(6,3,-0.00387164986045507));
    tripletList.push_back(T(6,6,0.503511552975707));
    tripletList.push_back(T(7,4,-0.00387164986045507));
    tripletList.push_back(T(7,7,0.503511552975707));
    tripletList.push_back(T(8,5,-0.00387164986045507));
    tripletList.push_back(T(8,8,0.503511552975707));

    B.setFromTriplets(tripletList.begin(), tripletList.end());
}

void define_C(SpMat &C){
    /*!==================
    |    define_C    |
    ==================

    Get the sixth order C stiffness 
    tensor in voigt notation.

    */

    std::vector<T> tripletList;
    tripletList.reserve(183);

    tripletList.push_back(T(0,0,8.97047749088427));
    tripletList.push_back(T(0,1,1.66068457301634));
    tripletList.push_back(T(0,2,1.66068457301634));
    tripletList.push_back(T(0,14,2.14259741437930));
    tripletList.push_back(T(0,17,2.63594416596082));
    tripletList.push_back(T(0,22,2.14259741437930));
    tripletList.push_back(T(0,25,2.63594416596082));
    tripletList.push_back(T(1,0,1.66068457301634));
    tripletList.push_back(T(1,1,1.63171548300639));
    tripletList.push_back(T(1,2,0.357817269957867));
    tripletList.push_back(T(1,14,1.37432904562166));
    tripletList.push_back(T(1,17,1.18589138191610));
    tripletList.push_back(T(1,22,0.500995125523459));
    tripletList.push_back(T(1,25,0.801872177535019));
    tripletList.push_back(T(2,0,1.66068457301634));
    tripletList.push_back(T(2,1,0.357817269957867));
    tripletList.push_back(T(2,2,1.63171548300639));
    tripletList.push_back(T(2,14,0.500995125523459));
    tripletList.push_back(T(2,17,0.801872177535019));
    tripletList.push_back(T(2,22,1.37432904562166));
    tripletList.push_back(T(2,25,1.18589138191610));
    tripletList.push_back(T(3,3,0.712702026982900));
    tripletList.push_back(T(3,6,0.561196186065625));
    tripletList.push_back(T(3,13,0.503083165307810));
    tripletList.push_back(T(3,16,0.370250754790395));
    tripletList.push_back(T(3,23,0.370250754790395));
    tripletList.push_back(T(3,26,0.0137684495906822));
    tripletList.push_back(T(4,4,2.09171782703280));
    tripletList.push_back(T(4,7,1.88958629453973));
    tripletList.push_back(T(4,12,0.875932634742095));
    tripletList.push_back(T(4,15,0.958139353683705));
    tripletList.push_back(T(4,18,1.18589138191610));
    tripletList.push_back(T(4,19,0.801872177535019));
    tripletList.push_back(T(4,20,2.63594416596082));
    tripletList.push_back(T(5,5,2.09171782703280));
    tripletList.push_back(T(5,8,1.88958629453973));
    tripletList.push_back(T(5,9,1.18589138191610));
    tripletList.push_back(T(5,10,2.63594416596082));
    tripletList.push_back(T(5,11,0.801872177535019));
    tripletList.push_back(T(5,21,0.958139353683705));
    tripletList.push_back(T(5,24,0.875932634742095));
    tripletList.push_back(T(6,3,0.561196186065625));
    tripletList.push_back(T(6,6,0.712702026982900));
    tripletList.push_back(T(6,13,0.370250754790395));
    tripletList.push_back(T(6,16,0.0137684495906822));
    tripletList.push_back(T(6,23,0.503083165307810));
    tripletList.push_back(T(6,26,0.370250754790395));
    tripletList.push_back(T(7,4,1.88958629453973));
    tripletList.push_back(T(7,7,1.40993341174572));
    tripletList.push_back(T(7,12,0.958139353683705));
    tripletList.push_back(T(7,15,0.683462935172136));
    tripletList.push_back(T(7,18,1.37432904562166));
    tripletList.push_back(T(7,19,0.500995125523459));
    tripletList.push_back(T(7,20,2.14259741437930));
    tripletList.push_back(T(8,5,1.88958629453973));
    tripletList.push_back(T(8,8,1.40993341174572));
    tripletList.push_back(T(8,9,1.37432904562166));
    tripletList.push_back(T(8,10,2.14259741437930));
    tripletList.push_back(T(8,11,0.500995125523459));
    tripletList.push_back(T(8,21,0.683462935172136));
    tripletList.push_back(T(8,24,0.958139353683705));
    tripletList.push_back(T(9,5,1.18589138191610));
    tripletList.push_back(T(9,8,1.37432904562166));
    tripletList.push_back(T(9,9,1.63171548300639));
    tripletList.push_back(T(9,10,1.66068457301634));
    tripletList.push_back(T(9,11,0.357817269957867));
    tripletList.push_back(T(9,21,0.500995125523459));
    tripletList.push_back(T(9,24,0.801872177535019));
    tripletList.push_back(T(10,5,2.63594416596082));
    tripletList.push_back(T(10,8,2.14259741437930));
    tripletList.push_back(T(10,9,1.66068457301634));
    tripletList.push_back(T(10,10,8.97047749088427));
    tripletList.push_back(T(10,11,1.66068457301634));
    tripletList.push_back(T(10,21,2.14259741437930));
    tripletList.push_back(T(10,24,2.63594416596082));
    tripletList.push_back(T(11,5,0.801872177535019));
    tripletList.push_back(T(11,8,0.500995125523459));
    tripletList.push_back(T(11,9,0.357817269957867));
    tripletList.push_back(T(11,10,1.66068457301634));
    tripletList.push_back(T(11,11,1.63171548300639));
    tripletList.push_back(T(11,21,1.37432904562166));
    tripletList.push_back(T(11,24,1.18589138191610));
    tripletList.push_back(T(12,4,0.875932634742095));
    tripletList.push_back(T(12,7,0.958139353683705));
    tripletList.push_back(T(12,12,2.09171782703280));
    tripletList.push_back(T(12,15,1.88958629453973));
    tripletList.push_back(T(12,18,0.801872177535019));
    tripletList.push_back(T(12,19,1.18589138191610));
    tripletList.push_back(T(12,20,2.63594416596082));
    tripletList.push_back(T(13,3,0.503083165307810));
    tripletList.push_back(T(13,6,0.370250754790395));
    tripletList.push_back(T(13,13,0.712702026982900));
    tripletList.push_back(T(13,16,0.561196186065625));
    tripletList.push_back(T(13,23,0.0137684495906822));
    tripletList.push_back(T(13,26,0.370250754790395));
    tripletList.push_back(T(14,0,2.14259741437930));
    tripletList.push_back(T(14,1,1.37432904562166));
    tripletList.push_back(T(14,2,0.500995125523459));
    tripletList.push_back(T(14,14,1.40993341174572));
    tripletList.push_back(T(14,17,1.88958629453973));
    tripletList.push_back(T(14,22,0.683462935172136));
    tripletList.push_back(T(14,25,0.958139353683705));
    tripletList.push_back(T(15,4,0.958139353683705));
    tripletList.push_back(T(15,7,0.683462935172136));
    tripletList.push_back(T(15,12,1.88958629453973));
    tripletList.push_back(T(15,15,1.40993341174572));
    tripletList.push_back(T(15,18,0.500995125523459));
    tripletList.push_back(T(15,19,1.37432904562166));
    tripletList.push_back(T(15,20,2.14259741437930));
    tripletList.push_back(T(16,3,0.370250754790395));
    tripletList.push_back(T(16,6,0.0137684495906822));
    tripletList.push_back(T(16,13,0.561196186065625));
    tripletList.push_back(T(16,16,0.712702026982900));
    tripletList.push_back(T(16,23,0.370250754790395));
    tripletList.push_back(T(16,26,0.503083165307810));
    tripletList.push_back(T(17,0,2.63594416596082));
    tripletList.push_back(T(17,1,1.18589138191610));
    tripletList.push_back(T(17,2,0.801872177535019));
    tripletList.push_back(T(17,14,1.88958629453973));
    tripletList.push_back(T(17,17,2.09171782703280));
    tripletList.push_back(T(17,22,0.958139353683705));
    tripletList.push_back(T(17,25,0.875932634742095));
    tripletList.push_back(T(18,4,1.18589138191610));
    tripletList.push_back(T(18,7,1.37432904562166));
    tripletList.push_back(T(18,12,0.801872177535019));
    tripletList.push_back(T(18,15,0.500995125523459));
    tripletList.push_back(T(18,18,1.63171548300639));
    tripletList.push_back(T(18,19,0.357817269957867));
    tripletList.push_back(T(18,20,1.66068457301634));
    tripletList.push_back(T(19,4,0.801872177535019));
    tripletList.push_back(T(19,7,0.500995125523459));
    tripletList.push_back(T(19,12,1.18589138191610));
    tripletList.push_back(T(19,15,1.37432904562166));
    tripletList.push_back(T(19,18,0.357817269957867));
    tripletList.push_back(T(19,19,1.63171548300639));
    tripletList.push_back(T(19,20,1.66068457301634));
    tripletList.push_back(T(20,4,2.63594416596082));
    tripletList.push_back(T(20,7,2.14259741437930));
    tripletList.push_back(T(20,12,2.63594416596082));
    tripletList.push_back(T(20,15,2.14259741437930));
    tripletList.push_back(T(20,18,1.66068457301634));
    tripletList.push_back(T(20,19,1.66068457301634));
    tripletList.push_back(T(20,20,8.97047749088427));
    tripletList.push_back(T(21,5,0.958139353683705));
    tripletList.push_back(T(21,8,0.683462935172136));
    tripletList.push_back(T(21,9,0.500995125523459));
    tripletList.push_back(T(21,10,2.14259741437930));
    tripletList.push_back(T(21,11,1.37432904562166));
    tripletList.push_back(T(21,21,1.40993341174572));
    tripletList.push_back(T(21,24,1.88958629453973));
    tripletList.push_back(T(22,0,2.14259741437930));
    tripletList.push_back(T(22,1,0.500995125523459));
    tripletList.push_back(T(22,2,1.37432904562166));
    tripletList.push_back(T(22,14,0.683462935172136));
    tripletList.push_back(T(22,17,0.958139353683705));
    tripletList.push_back(T(22,22,1.40993341174572));
    tripletList.push_back(T(22,25,1.88958629453973));
    tripletList.push_back(T(23,3,0.370250754790395));
    tripletList.push_back(T(23,6,0.503083165307810));
    tripletList.push_back(T(23,13,0.0137684495906822));
    tripletList.push_back(T(23,16,0.370250754790395));
    tripletList.push_back(T(23,23,0.712702026982900));
    tripletList.push_back(T(23,26,0.561196186065625));
    tripletList.push_back(T(24,5,0.875932634742095));
    tripletList.push_back(T(24,8,0.958139353683705));
    tripletList.push_back(T(24,9,0.801872177535019));
    tripletList.push_back(T(24,10,2.63594416596082));
    tripletList.push_back(T(24,11,1.18589138191610));
    tripletList.push_back(T(24,21,1.88958629453973));
    tripletList.push_back(T(24,24,2.09171782703280));
    tripletList.push_back(T(25,0,2.63594416596082));
    tripletList.push_back(T(25,1,0.801872177535019));
    tripletList.push_back(T(25,2,1.18589138191610));
    tripletList.push_back(T(25,14,0.958139353683705));
    tripletList.push_back(T(25,17,0.875932634742095));
    tripletList.push_back(T(25,22,1.88958629453973));
    tripletList.push_back(T(25,25,2.09171782703280));
    tripletList.push_back(T(26,3,0.0137684495906822));
    tripletList.push_back(T(26,6,0.370250754790395));
    tripletList.push_back(T(26,13,0.370250754790395));
    tripletList.push_back(T(26,16,0.503083165307810));
    tripletList.push_back(T(26,23,0.561196186065625));
    tripletList.push_back(T(26,26,0.712702026982900));

    C.setFromTriplets(tripletList.begin(), tripletList.end());
    return;
}

void define_D(SpMat &D){
    /*!===============
    |    get_D    |
    ===============

    Get the forth order D stiffness
    tensor in voigt notation.

    */

    std::vector<T> tripletList;
    tripletList.reserve(21);

    tripletList.push_back(T(0,0,1.33828709399996));
    tripletList.push_back(T(0,1,0.785358583713769));
    tripletList.push_back(T(0,2,0.785358583713769));
    tripletList.push_back(T(1,0,0.785358583713769));
    tripletList.push_back(T(1,1,1.33828709399996));
    tripletList.push_back(T(1,2,0.785358583713769));
    tripletList.push_back(T(2,0,0.785358583713769));
    tripletList.push_back(T(2,1,0.785358583713769));
    tripletList.push_back(T(2,2,1.33828709399996));
    tripletList.push_back(T(3,3,0.276464255143097));
    tripletList.push_back(T(3,6,0.276464255143097));
    tripletList.push_back(T(4,4,0.276464255143097));
    tripletList.push_back(T(4,7,0.276464255143097));
    tripletList.push_back(T(5,5,0.276464255143097));
    tripletList.push_back(T(5,8,0.276464255143097));
    tripletList.push_back(T(6,3,0.276464255143097));
    tripletList.push_back(T(6,6,0.276464255143097));
    tripletList.push_back(T(7,4,0.276464255143097));
    tripletList.push_back(T(7,7,0.276464255143097));
    tripletList.push_back(T(8,5,0.276464255143097));
    tripletList.push_back(T(8,8,0.276464255143097));

    D.setFromTriplets(tripletList.begin(), tripletList.end());
    return;
}

void define_N(double &N){
    /*!==================
    |    define_N    |
    ==================
    
    Define the value of the shape function to be used.
    
    */
    
   // N = 0.261;
    N = 0.490563;
    return;
}

void define_dNdx(double (&dNdx)[3], bool MOOSE = false){
    /*!=====================
    |    define_dNdx    |
    =====================
    
    Define the gradient of the shape function 
    to be used.
    */
    
//    dNdx[0] =  1.42;
//    dNdx[1] =  0.271;
//    dNdx[2] = -2.31;

    dNdx[0] = -0.622008;
    dNdx[1] = -0.622008;
    dNdx[2] = -0.622008;
    return;
}

void define_eta(double &eta){
    /*!====================
    |    define_eta    |
    ====================

    Define the value of the interpolation function to be used.

    */

//    eta = 0.826;
    eta = 0.490563;
    return;
}

void define_detadx(double (&detadx)[3]){
    /*!=======================
    |    define_detadx    |
    =======================

    Define the value of the gradient of the interpolation function.

    */

//    detadx[0] = 0.172;
//    detadx[1] = 3.121;
//    detadx[2] = 0.761;

    detadx[0] = -0.622008;
    detadx[1] = -0.622008;
    detadx[2] = -0.622008;
    return;
}

void define_density(double &density){
    /*!========================
    |    define_density    |
    ========================
    
    Define the density.
    
    */
    
    density = 1412.32;
}

void map_eigen_vector_to_std_vector(const Eigen::MatrixXd &EV, std::vector<double> &V){
    /*!========================================
    |    map_eigen_vector_to_std_vector    |
    ========================================
    
    Map an eigen vector to a standard vector.
    
    */
    
    int nrows = EV.rows();
    int ncols = EV.cols();
    
    V.resize(nrows*ncols);
    
    int n = 0;
    for (int i=0; i<nrows; i++){
        for (int j=0; j<ncols; j++){
            V[n] = EV(i,j);
            n++;
        }
    }
    return;
}

std::vector<double> parse_pk2_stress_RCG(std::vector<double> RCGvec){
    /*!==============================
    |    parse_pk2_stress_RCG    |
    ==============================
    
    Parse the PK2 stress as a function of the right cauchy-green 
    deformation tensor in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 RCG_voigt;
    Matrix_3x3 RCG;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    for (int i=0; i<9; i++){RCG_voigt[i] = RCGvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(RCG_voigt,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    std::vector<double> PK2_vec;
    PK2_vec.resize(9);
    for (int i=0; i<9; i++){PK2_vec[i] = PK2(i);}
    return PK2_vec;
}

std::vector<double> parse_symmetric_stress_RCG(std::vector<double> RCGvec){
    /*!====================================
    |    parse_symmetric_stress_RCG    |
    ====================================
    
    Parse the symmetric stress as a function of the right cauchy-green 
    deformation tensor in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 RCG_voigt;
    Matrix_3x3 RCG;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    for (int i=0; i<9; i++){RCG_voigt[i] = RCGvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(RCG_voigt,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    std::vector<double> SIGMA_vec;
    SIGMA_vec.resize(9);
    for (int i=0; i<9; i++){SIGMA_vec[i] = SIGMA(i);}
    return SIGMA_vec;
}

std::vector<double> parse_pk2_stress_Psi(std::vector<double> Psivec){
    /*!==============================
    |    parse_pk2_stress_Psi    |
    ==============================
    
    Parse the PK2 stress as a function of the deformation measure Psi in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 Psi_voigt;
    Matrix_3x3 Psi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    for (int i=0; i<9; i++){Psi_voigt[i] = Psivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(Psi_voigt,Psi);

    Matrix_3x3 RCG   = F.transpose()*F;
    Matrix_3x3 RCGinv = RCG.inverse();
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    std::vector<double> PK2_vec;
    PK2_vec.resize(9);
    for (int i=0; i<9; i++){PK2_vec[i] = PK2(i);}
    return PK2_vec;
}

std::vector<double> parse_symmetric_stress_Psi(std::vector<double> Psivec){
    /*!=================================
    |    parse_symmetric_stress_Psi    |
    ====================================
    
    Parse the symmetric stress in the reference configuration 
    as a function of the deformation measure Psi in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 Psi_voigt;
    Matrix_3x3 Psi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    for (int i=0; i<9; i++){Psi_voigt[i] = Psivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(Psi_voigt,Psi);

    Matrix_3x3 RCG   = F.transpose()*F;
    Matrix_3x3 RCGinv = RCG.inverse();
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    std::vector<double> SIGMA_vec;
    SIGMA_vec.resize(9);
    for (int i=0; i<9; i++){SIGMA_vec[i] = SIGMA(i);}
    return SIGMA_vec;
}

std::vector<double> parse_pk2_stress_Gamma(std::vector<double> Gammavec){
    /*!================================
    |    parse_pk2_stress_Gamma    |
    ================================
    
    Parse the PK2 stress as a function of the deformation measure Gamma in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27  Gamma_voigt;
    Matrix_3x9 Gamma;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    define_deformation_gradient(F);
    define_chi(chi);
    
    for (int i=0; i<27; i++){Gamma_voigt[i] = Gammavec[i];}
    deformation_measures::undo_voigt_3x9_tensor(Gamma_voigt,Gamma);

    Matrix_3x3 RCG    = F.transpose()*F;
    Matrix_3x3 RCGinv = RCG.inverse();
    Matrix_3x3 Psi    = F.transpose()*chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    std::vector<double> PK2_vec;
    PK2_vec.resize(9);
    for (int i=0; i<9; i++){PK2_vec[i] = PK2(i);}
    return PK2_vec;
}

std::vector<double> parse_symmetric_stress_Gamma(std::vector<double> Gammavec){
    /*!======================================
    |    parse_symmetric_stress_Gamma    |
    ======================================
    
    Parse the symmetric stress as a function of the 
    deformation measure Gamma in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27  Gamma_voigt;
    Matrix_3x9 Gamma;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    define_deformation_gradient(F);
    define_chi(chi);
    
    for (int i=0; i<27; i++){Gamma_voigt[i] = Gammavec[i];}
    deformation_measures::undo_voigt_3x9_tensor(Gamma_voigt,Gamma);

    Matrix_3x3 RCG    = F.transpose()*F;
    Matrix_3x3 RCGinv = RCG.inverse();
    Matrix_3x3 Psi    = F.transpose()*chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    std::vector<double> SIGMA_vec;
    SIGMA_vec.resize(9);
    for (int i=0; i<9; i++){SIGMA_vec[i] = SIGMA(i);}
    return SIGMA_vec;
}

std::vector<double> parse_pk2_stress_F(std::vector<double> Fvec){
    /*!============================
    |    parse_pk2_stress_F    |
    ============================
    
    Parse the PK2 stress as a function of the deformation 
    gradient in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 F_voigt;
    Matrix_3x3 F;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){F_voigt[i] = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    std::vector<double> PK2_vec;
    PK2_vec.resize(9);
    for (int i=0; i<9; i++){PK2_vec[i] = PK2(i);}
    return PK2_vec;
}

std::vector<double> parse_pk2_stress_chi(std::vector<double> chivec){
    /*!==============================
    |    parse_pk2_stress_chi    |
    ==============================
    
    Parse the PK2 stress as a function of the micro-displacement 
    tensor chi in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 chi_voigt;
    Matrix_3x3 chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){chi_voigt[i] = chivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    std::vector<double> PK2_vec;
    PK2_vec.resize(9);
    for (int i=0; i<9; i++){PK2_vec[i] = PK2(i);}
    return PK2_vec;
}

std::vector<double> parse_pk2_stress_grad_chi(std::vector<double> grad_chivec){
    /*!===================================
    |    parse_pk2_stress_grad_chi    |
    ===================================
    
    Parse the PK2 stress as a function of the gradient of 
    the micro-displacement tensor chi in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 grad_chi_voigt;
    Matrix_3x9 grad_chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<27; i++){grad_chi_voigt[i] = grad_chivec[i];}
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    define_deformation_gradient(F);
    define_chi(chi);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    std::vector<double> PK2_vec;
    PK2_vec.resize(9);
    for (int i=0; i<9; i++){PK2_vec[i] = PK2(i);}
    return PK2_vec;
}

std::vector<double> parse_symmetric_stress_F(std::vector<double> Fvec){
    /*!==================================
    |    parse_symmetric_stress_F    |
    ==================================
    
    Parse the symmetric stress as a function of the deformation 
    gradient.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 F_voigt;
    Matrix_3x3 F;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){F_voigt[i] = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    std::vector<double> SIGMA_vec;
    SIGMA_vec.resize(9);
    for (int i=0; i<9; i++){SIGMA_vec[i] = SIGMA(i);}
    return SIGMA_vec;
}

std::vector<double> parse_symmetric_stress_chi(std::vector<double> chivec){
    /*!====================================
    |    parse_symmetric_stress_chi    |
    ====================================
    
    Parse the symmetric stress as a function of the micro-displacement 
    tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 chi_voigt;
    Matrix_3x3 chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){chi_voigt[i] = chivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    std::vector<double> SIGMA_vec;
    SIGMA_vec.resize(9);
    for (int i=0; i<9; i++){SIGMA_vec[i] = SIGMA(i);}
    return SIGMA_vec;
}

std::vector<double> parse_symmetric_stress_grad_chi(std::vector<double> grad_chivec){
    /*!====================================
    |    parse_symmetric_stress_chi    |
    ====================================
    
    Parse the symmetric stress as a function of the micro-displacement 
    tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 grad_chi_voigt;
    Matrix_3x9 grad_chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<27; i++){grad_chi_voigt[i] = grad_chivec[i];}
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    define_deformation_gradient(F);
    define_chi(chi);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    std::vector<double> SIGMA_vec;
    SIGMA_vec.resize(9);
    for (int i=0; i<9; i++){SIGMA_vec[i] = SIGMA(i);}
    return SIGMA_vec;
}

std::vector<double> parse_higher_order_stress_F(std::vector<double> Fvec){
    /*!=====================================
    |    parse_higher_order_stress_F    |
    =====================================
    
    Parse the symmetric stress as a function of the deformation 
    gradient.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 F_voigt;
    Matrix_3x3 F;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x9 grad_chi;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<9; i++){F_voigt[i] = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    //Deformation measures
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures in voigt notation
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    std::vector<double> M_vec;
    M_vec.resize(27);
    for (int i=0; i<27; i++){M_vec[i] = M(i);}
    return M_vec;
}

std::vector<double> parse_higher_order_stress_chi(std::vector<double> chivec){
    /*!=======================================
    |    parse_higher_order_stress_chi    |
    =======================================
    
    Parse the symmetric stress as a function of the micro-deformation
    tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 chi_voigt;
    Matrix_3x3 chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<9; i++){chi_voigt[i] = chivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    //Deformation measures
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures in voigt notation
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    std::vector<double> M_vec;
    M_vec.resize(27);
    for (int i=0; i<27; i++){M_vec[i] = M(i);}
    return M_vec;
}

std::vector<double> parse_higher_order_stress_grad_chi(std::vector<double> grad_chivec){
    /*!============================================
    |    parse_higher_order_stress_grad_chi    |
    ============================================
    
    Parse the symmetric stress as a function of the gradient of 
    the micro-deformation tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 grad_chi_voigt;
    Matrix_3x9 grad_chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<27; i++){grad_chi_voigt[i] = grad_chivec[i];}
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    define_deformation_gradient(F);
    define_chi(chi);
    
    //Deformation measures
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures in voigt notation
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    std::vector<double> M_vec;
    M_vec.resize(27);
    for (int i=0; i<27; i++){M_vec[i] = M(i);}
    return M_vec;
}

std::vector<double> parse_higher_order_stress_gamma(std::vector<double> Gammavec){
    /*!=========================================
    |    parse_higher_order_stress_gamma    |
    =========================================
    
    Parse the higher order stress as a function of the micro 
    deformation measure gamma.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 Gamma_voigt;
    Matrix_3x9 Gamma;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<27; i++){Gamma_voigt[i] = Gammavec[i];}
    deformation_measures::undo_voigt_3x9_tensor(Gamma_voigt,Gamma);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    std::vector<double> M_vec;
    M_vec.resize(27);
    for (int i=0; i<27; i++){M_vec[i] = M(i);}
    return M_vec;
}

std::vector<double> parse_cauchy_stress_F(std::vector<double> Fvec){
    /*!===============================
    |    parse_cauchy_stress_F    |
    ===============================
    
    Parse the cauchy stress as a function of the deformation 
    gradient in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 F_voigt;
    Matrix_3x3 F;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){F_voigt[i] = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);

    //Map the PK2 stress to the cauchy stress
    double Jac = F.determinant();
    Vector_9 cauchy;
    Matrix_3x3 PK2_mat;
    deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);
    deformation_measures::voigt_3x3_tensor(F*PK2_mat*F.transpose()/Jac,cauchy);
    
    std::vector<double> cauchy_vec;
    cauchy_vec.resize(9);
    for (int i=0; i<9; i++){cauchy_vec[i] = cauchy(i);}
    return cauchy_vec;
}

std::vector<double> parse_cauchy_stress_chi(std::vector<double> chivec){
    /*!=================================
    |    parse_cauchy_stress_chi    |
    =================================
    
    Parse the cauchy stress as a function of the micro-displacement 
    tensor chi in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 chi_voigt;
    Matrix_3x3 chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){chi_voigt[i] = chivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    //Map the PK2 stress to the cauchy stress
    double Jac = F.determinant();
    Vector_9 cauchy;
    Matrix_3x3 PK2_mat;
    deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);
    deformation_measures::voigt_3x3_tensor(F*PK2_mat*F.transpose()/Jac,cauchy);
    
    std::vector<double> cauchy_vec;
    cauchy_vec.resize(9);
    for (int i=0; i<9; i++){cauchy_vec[i] = cauchy(i);}
    return cauchy_vec;
}

std::vector<double> parse_cauchy_stress_grad_chi(std::vector<double> grad_chivec){
    /*!======================================
    |    parse_cauchy_stress_grad_chi    |
    ======================================
    
    Parse the cauchy stress as a function of the gradient of 
    the micro-displacement tensor chi in vector form.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 grad_chi_voigt;
    Matrix_3x9 grad_chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<27; i++){grad_chi_voigt[i] = grad_chivec[i];}
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    define_deformation_gradient(F);
    define_chi(chi);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    //Map the PK2 stress to the cauchy stress
    double Jac = F.determinant();
    Vector_9 cauchy;
    Matrix_3x3 PK2_mat;
    deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);
    deformation_measures::voigt_3x3_tensor(F*PK2_mat*F.transpose()/Jac,cauchy);
    
    std::vector<double> cauchy_vec;
    cauchy_vec.resize(9);
    for (int i=0; i<9; i++){cauchy_vec[i] = cauchy(i);}
    return cauchy_vec;
}

std::vector<double> parse_cauchy_stress_grad_u(std::vector<double> grad_uvec){
    /*!===============================
    |    parse_cauchy_stress_F    |
    ===============================
    
    Parse the cauchy stress as a function of the gradient of u in the 
    current configuration.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 _grad_u;
    for (int i=0; i<9; i++){_grad_u[i] = grad_uvec[i];}
    Matrix_3x3 grad_u;
    
    deformation_measures::undo_voigt_3x3_tensor(_grad_u, grad_u);
    
    double _grad_u_data[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u_data[i][j] = grad_u(i,j);
        }
    }
    
    Matrix_3x3 F;
    deformation_measures::get_deformation_gradient(_grad_u_data, F);
 
    Matrix_3x3 chi;
    define_chi(chi);
    
    Matrix_3x9 grad_chi;
    Eigen::Matrix<double,9,3> grad_phi_data; //The incoming data
    grad_phi_data << 0.53155137,  0.53182759,  0.63440096,  0.09210494,  0.43370117,
                     0.43086276,  0.31728548,  0.41482621,  0.86630916,  0.4936851 ,
                     0.42583029,  0.31226122,  0.72244338,  0.32295891,  0.36178866,
                     0.84943179,  0.72445532,  0.61102351,  0.50183668,  0.62395295,
                     0.1156184 ,  0.42635131,  0.89338916,  0.94416002,  0.22826323,
                     0.29371405,  0.63097612;

    double grad_phi_data_array[9][3]; //The format expected by the function

    for (int i=0; i<9; i++){

        for (int j=0; j<3; j++){

            grad_phi_data_array[i][j] = grad_phi_data(i,j);

        }

    }
    
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);
    
    //Map the PK2 stress to the cauchy stress
    double Jac = F.determinant();
    Vector_9 cauchy;
    Matrix_3x3 PK2_mat;
    deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);
    deformation_measures::voigt_3x3_tensor(F*PK2_mat*F.transpose()/Jac,cauchy);
    
    std::vector<double> cauchy_vec;
    cauchy_vec.resize(9);
    for (int i=0; i<9; i++){cauchy_vec[i] = cauchy(i);}
    return cauchy_vec;
}

std::vector<double> parse_cauchy_stress_grad_phi(std::vector<double> grad_phivec){
    /*!======================================
    |    parse_cauchy_stress_grad_phi    |
    ======================================
    
    Parse the cauchy stress as a function of the gradient of u in the 
    current configuration.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 _grad_phi;
    for (int i=0; i<27; i++){_grad_phi[i] = grad_phivec[i];}
    
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    double _grad_u_data[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u_data[i][j] = grad_u(i,j);
        }
    }
    
    Matrix_3x3 F;
    deformation_measures::get_deformation_gradient(_grad_u_data, F);
    
    Matrix_3x3 chi;
    define_chi(chi);
    
    deformation_measures::perform_right_positive_cyclic_permutation(_grad_phi);
    Matrix_3x9 grad_phi;
    Matrix_3x9 grad_chi;
    deformation_measures::undo_voigt_3x9_tensor(_grad_phi, grad_phi);
    grad_phi = F.transpose()*grad_phi.eval();
    deformation_measures::voigt_3x9_tensor(grad_phi,_grad_phi);
    deformation_measures::perform_left_positive_cyclic_permutation(_grad_phi);
    deformation_measures::undo_voigt_3x9_tensor(_grad_phi,grad_chi);
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 PK2;
    micro_material::compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, PK2);

    //Map the PK2 stress to the cauchy stress
    double Jac = F.determinant();
    Vector_9 cauchy;
    Matrix_3x3 PK2_mat;
    deformation_measures::undo_voigt_3x3_tensor(PK2,PK2_mat);
    deformation_measures::voigt_3x3_tensor(F*PK2_mat*F.transpose()/Jac,cauchy);
    
    std::vector<double> cauchy_vec;
    cauchy_vec.resize(9);
    for (int i=0; i<9; i++){cauchy_vec[i] = cauchy(i);}
    return cauchy_vec;
}

std::vector<double> parse_s_stress_F(std::vector<double> Fvec){
    /*!==========================
    |    parse_s_stress_F    |
    ==========================
    
    Parse the symmetric stress in the current configuration 
    as a function of the deformation gradient.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 F_voigt;
    Matrix_3x3 F;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){F_voigt[i] = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);

    //Map the symmetric stress to the reference configuration
    double Jac = F.determinant();
    Vector_9 s;
    Matrix_3x3 SIGMA_mat;
    deformation_measures::undo_voigt_3x3_tensor(SIGMA,SIGMA_mat);
    deformation_measures::voigt_3x3_tensor(F*SIGMA_mat*F.transpose()/Jac,s);
    
    std::vector<double> s_vec;
    s_vec.resize(9);
    for (int i=0; i<9; i++){s_vec[i] = s(i);}
    return s_vec;
}

std::vector<double> parse_s_stress_chi(std::vector<double> chivec){
    /*!============================
    |    parse_s_stress_chi    |
    ============================
    
    Parse the symmetric stress in the current configuration 
    as a function of the micro-displacement tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 chi_voigt;
    Matrix_3x3 chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<9; i++){chi_voigt[i] = chivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    //Map the symmetric stress to the reference configuration
    double Jac = F.determinant();
    Vector_9 s;
    Matrix_3x3 SIGMA_mat;
    deformation_measures::undo_voigt_3x3_tensor(SIGMA,SIGMA_mat);
    deformation_measures::voigt_3x3_tensor(F*SIGMA_mat*F.transpose()/Jac,s);
    
    std::vector<double> s_vec;
    s_vec.resize(9);
    for (int i=0; i<9; i++){s_vec[i] = s(i);}
    return s_vec;
}

std::vector<double> parse_s_stress_grad_chi(std::vector<double> grad_chivec){
    /*!=================================
    |    parse_s_stress_grad_chi    |
    =================================
    
    Parse the symmetric stress in the current configuration 
    as a function of the gradient of the micro-displacement 
    tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 grad_chi_voigt;
    Matrix_3x9 grad_chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    for (int i=0; i<27; i++){grad_chi_voigt[i] = grad_chivec[i];}
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    define_deformation_gradient(F);
    define_chi(chi);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    //Map the symmetric stress to the reference configuration
    double Jac = F.determinant();
    Vector_9 s;
    Matrix_3x3 SIGMA_mat;
    deformation_measures::undo_voigt_3x3_tensor(SIGMA,SIGMA_mat);
    deformation_measures::voigt_3x3_tensor(F*SIGMA_mat*F.transpose()/Jac,s);
    
    std::vector<double> s_vec;
    s_vec.resize(9);
    for (int i=0; i<9; i++){s_vec[i] = s(i);}
    return s_vec;
}

std::vector<double> parse_s_stress_grad_u(std::vector<double> grad_uvec){
    /*!===============================
    |    parse_s_stress_grad_u    |
    ===============================
    
    Parse the symmetric stress in the current configuration as a function
    of the gradient of u in the current configuration.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 _grad_u;
    for (int i=0; i<9; i++){_grad_u[i] = grad_uvec[i];}
    Matrix_3x3 grad_u;
    
    deformation_measures::undo_voigt_3x3_tensor(_grad_u, grad_u);
    
    double _grad_u_data[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u_data[i][j] = grad_u(i,j);
        }
    }
    
    Matrix_3x3 F;
    deformation_measures::get_deformation_gradient(_grad_u_data, F);
    
    Matrix_3x3 chi;
    define_chi(chi);
    
    Matrix_3x9 grad_chi;
    Eigen::Matrix<double,9,3> grad_phi_data; //The incoming data
    grad_phi_data << 0.53155137,  0.53182759,  0.63440096,  0.09210494,  0.43370117,
                     0.43086276,  0.31728548,  0.41482621,  0.86630916,  0.4936851 ,
                     0.42583029,  0.31226122,  0.72244338,  0.32295891,  0.36178866,
                     0.84943179,  0.72445532,  0.61102351,  0.50183668,  0.62395295,
                     0.1156184 ,  0.42635131,  0.89338916,  0.94416002,  0.22826323,
                     0.29371405,  0.63097612;

    double grad_phi_data_array[9][3]; //The format expected by the function

    for (int i=0; i<9; i++){

        for (int j=0; j<3; j++){

            grad_phi_data_array[i][j] = grad_phi_data(i,j);

        }

    }
    
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    
    
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    //Map the symmetric stress to the reference configuration
    double Jac = F.determinant();
    Vector_9 s;
    Matrix_3x3 SIGMA_mat;
    deformation_measures::undo_voigt_3x3_tensor(SIGMA,SIGMA_mat);
    deformation_measures::voigt_3x3_tensor(F*SIGMA_mat*F.transpose()/Jac,s);
    
    std::vector<double> s_vec;
    s_vec.resize(9);
    for (int i=0; i<9; i++){s_vec[i] = s(i);}
    return s_vec;
}

std::vector<double> parse_s_stress_grad_phi(std::vector<double> grad_phivec){
    /*!======================================
    |    parse_s_stress_grad_phi    |
    ======================================
    
    Parse the symmetric stress in the current configuration as a function of the gradient of u in the 
    current configuration.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 _grad_phi;
    for (int i=0; i<27; i++){_grad_phi[i] = grad_phivec[i];}
    
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    double _grad_u_data[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u_data[i][j] = grad_u(i,j);
        }
    }
    
    Matrix_3x3 F;
    deformation_measures::get_deformation_gradient(_grad_u_data, F);
    
    Matrix_3x3 chi;
    define_chi(chi);
    
    deformation_measures::perform_right_positive_cyclic_permutation(_grad_phi);
    Matrix_3x9 grad_phi;
    Matrix_3x9 grad_chi;
    deformation_measures::undo_voigt_3x9_tensor(_grad_phi, grad_phi);
    grad_phi = F.transpose()*grad_phi.eval();
    deformation_measures::voigt_3x9_tensor(grad_phi,_grad_phi);
    deformation_measures::perform_left_positive_cyclic_permutation(_grad_phi);
    deformation_measures::undo_voigt_3x9_tensor(_grad_phi,grad_chi);
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_9 SIGMA;
    micro_material::compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma, A, B, C, D, SIGMA);
    
    //Map the symmetric stress to the reference configuration
    double Jac = F.determinant();
    Vector_9 s;
    Matrix_3x3 SIGMA_mat;
    deformation_measures::undo_voigt_3x3_tensor(SIGMA,SIGMA_mat);
    deformation_measures::voigt_3x3_tensor(F*SIGMA_mat*F.transpose()/Jac,s);
    
    std::vector<double> s_vec;
    s_vec.resize(9);
    for (int i=0; i<9; i++){s_vec[i] = s(i);}
    return s_vec;
}

std::vector<double> parse_m_stress_F(std::vector<double> Fvec){
    /*!==========================
    |    parse_m_stress_F    |
    ==========================
    
    Parse the higher order stress in the current configuration 
    as a function of the deformation gradient.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 F_voigt;
    Matrix_3x3 F;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<9; i++){F_voigt[i] = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    //Deformation measures
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures in voigt notation
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);

    //Map the higher order stress to the current configuration
    double Jac = F.determinant();
    Matrix_3x9 M_mat;
    Vector_27 m;
    deformation_measures::undo_voigt_3x9_tensor(M,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the first index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the second index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(chi*M_mat,m);               //Map the third index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    m = m/Jac;
    
    std::vector<double> m_vec;
    m_vec.resize(27);
    for (int i=0; i<27; i++){m_vec[i] = m(i);}
    return m_vec;
}

std::vector<double> parse_m_stress_chi(std::vector<double> chivec){
    /*!============================
    |    parse_m_stress_chi    |
    ============================
    
    Parse the higher order stress in the current configuration as a function of the micro-deformation
    tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 chi_voigt;
    Matrix_3x3 chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<9; i++){chi_voigt[i] = chivec[i];}
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    
    //Deformation measures
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures in voigt notation
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    //Map the higher order stress to the current configuration
    double Jac = F.determinant();
    Matrix_3x9 M_mat;
    Vector_27 m;
    deformation_measures::undo_voigt_3x9_tensor(M,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the first index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the second index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(chi*M_mat,m);               //Map the third index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    m = m/Jac;
    
    std::vector<double> m_vec;
    m_vec.resize(27);
    for (int i=0; i<27; i++){m_vec[i] = m(i);}
    return m_vec;
}

std::vector<double> parse_m_stress_grad_chi(std::vector<double> grad_chivec){
    /*!=================================
    |    parse_m_stress_grad_chi    |
    =================================
    
    Parse the higher order stress in the current configuration 
    as a function of the gradient of the micro-deformation tensor chi.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 grad_chi_voigt;
    Matrix_3x9 grad_chi;
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    
    SpMat C(27,27);
    define_C(C);
    
    for (int i=0; i<27; i++){grad_chi_voigt[i] = grad_chivec[i];}
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    define_deformation_gradient(F);
    define_chi(chi);
    
    //Deformation measures
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures in voigt notation
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    //Map the higher order stress to the current configuration
    double Jac = F.determinant();
    Matrix_3x9 M_mat;
    Vector_27 m;
    deformation_measures::undo_voigt_3x9_tensor(M,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the first index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the second index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(chi*M_mat,m);               //Map the third index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    m = m/Jac;
    
    std::vector<double> m_vec;
    m_vec.resize(27);
    for (int i=0; i<27; i++){m_vec[i] = m(i);}
    return m_vec;
}

std::vector<double> parse_m_stress_grad_u(std::vector<double> grad_uvec){
    /*!===============================
    |    parse_s_stress_grad_u    |
    ===============================
    
    Parse the higher order stress in the current configuration as a function
    of the gradient of u in the current configuration.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_9 _grad_u;
    for (int i=0; i<9; i++){_grad_u[i] = grad_uvec[i];}
    Matrix_3x3 grad_u;
    
    deformation_measures::undo_voigt_3x3_tensor(_grad_u, grad_u);
    
    double _grad_u_data[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u_data[i][j] = grad_u(i,j);
        }
    }
    
    Matrix_3x3 F;
    deformation_measures::get_deformation_gradient(_grad_u_data, F);
    
    Matrix_3x3 chi;
    define_chi(chi);
    
    Matrix_3x9 grad_chi;
    Eigen::Matrix<double,9,3> grad_phi_data; //The incoming data
    grad_phi_data << 0.53155137,  0.53182759,  0.63440096,  0.09210494,  0.43370117,
                     0.43086276,  0.31728548,  0.41482621,  0.86630916,  0.4936851 ,
                     0.42583029,  0.31226122,  0.72244338,  0.32295891,  0.36178866,
                     0.84943179,  0.72445532,  0.61102351,  0.50183668,  0.62395295,
                     0.1156184 ,  0.42635131,  0.89338916,  0.94416002,  0.22826323,
                     0.29371405,  0.63097612;

    double grad_phi_data_array[9][3]; //The format expected by the function

    for (int i=0; i<9; i++){

        for (int j=0; j<3; j++){

            grad_phi_data_array[i][j] = grad_phi_data(i,j);

        }

    }
    
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    
    
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    //Map the higher order stress to the current configuration
    double Jac = F.determinant();
    Matrix_3x9 M_mat;
    Vector_27 m;
    deformation_measures::undo_voigt_3x9_tensor(M,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the first index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the second index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(chi*M_mat,m);               //Map the third index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    m = m/Jac;
    
    std::vector<double> m_vec;
    m_vec.resize(27);
    for (int i=0; i<27; i++){m_vec[i] = m(i);}
    return m_vec;
}

std::vector<double> parse_m_stress_grad_phi(std::vector<double> grad_phivec){
    /*!======================================
    |    parse_m_stress_grad_phi    |
    ======================================
    
    Parse the higher order stress in the current configuration as a function of the gradient of phi in the 
    current configuration.
    
    */
    
    //Map the deformation gradient from the incoming vector to an eigen matrix
    Vector_27 _grad_phi;
    for (int i=0; i<27; i++){_grad_phi[i] = grad_phivec[i];}
    
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    double _grad_u_data[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u_data[i][j] = grad_u(i,j);
        }
    }
    
    Matrix_3x3 F;
    deformation_measures::get_deformation_gradient(_grad_u_data, F);
    
    Matrix_3x3 chi;
    define_chi(chi);
    
    deformation_measures::perform_right_positive_cyclic_permutation(_grad_phi);
    Matrix_3x9 grad_phi;
    Matrix_3x9 grad_chi;
    deformation_measures::undo_voigt_3x9_tensor(_grad_phi, grad_phi);
    grad_phi = F.transpose()*grad_phi.eval();
    deformation_measures::voigt_3x9_tensor(grad_phi,_grad_phi);
    deformation_measures::perform_left_positive_cyclic_permutation(_grad_phi);
    deformation_measures::undo_voigt_3x9_tensor(_grad_phi,grad_chi);
    
    //Define additional required values
    double t  = 0;
    double dt = 0;
    
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    Matrix_3x3 RCG;
    deformation_measures::get_right_cauchy_green(F,RCG);
    Matrix_3x3 RCGinv = RCG.inverse();

    Matrix_3x3 Psi   = F.transpose()*chi;
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    //Compute the required measures
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    deformation_measures::voigt_3x3_tensor(0.5*(RCG - Matrix_3x3::Identity()),E_voigt);
    deformation_measures::voigt_3x3_tensor(Psi - Matrix_3x3::Identity(),E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    Vector_27 M;
    micro_material::compute_higher_order_stress(Gamma_voigt, C, M);
    
    //Map the higher order stress to the current configuration
    double Jac = F.determinant();
    Matrix_3x9 M_mat;
    Vector_27 m;
    deformation_measures::undo_voigt_3x9_tensor(M,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the first index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(F*M_mat,m);                 //Map the second index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    deformation_measures::undo_voigt_3x9_tensor(m,M_mat);
    deformation_measures::voigt_3x9_tensor(chi*M_mat,m);               //Map the third index
    deformation_measures::perform_left_positive_cyclic_permutation(m); //Cycle the indices
    m = m/Jac;
    
    std::vector<double> m_vec;
    m_vec.resize(27);
    for (int i=0; i<27; i++){m_vec[i] = m(i);}
    return m_vec;
}

std::vector<double> parse_balance_of_linear_momentum_U(std::vector<double> U){
    /*!============================================
    |    parse_balance_of_linear_momentum_U    |
    ============================================
    
    Parse the balance of linear momentum as a function of the 
    degree of freedom vector.
    
    */
    
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);

    double eta;
    define_eta(eta);

    double detadx[3];
    define_detadx(detadx);
    
    //Expand U
    double grad_u[3][3];
    
    grad_u[0][0] = U[0]*detadx[0];
    grad_u[1][1] = U[1]*detadx[1];
    grad_u[2][2] = U[2]*detadx[2];
    grad_u[1][2] = U[1]*detadx[2];
    grad_u[0][2] = U[0]*detadx[2];
    grad_u[0][1] = U[0]*detadx[1];
    grad_u[2][1] = U[2]*detadx[1];
    grad_u[2][0] = U[2]*detadx[0];
    grad_u[1][0] = U[1]*detadx[0];

    double phi[9];
    for (int i=0; i<9; i++){phi[i] = U[i+3]*eta;}
    
    double grad_phi_data[9][3];
    
    for (int I=3; I<12; I++){
        for (int j=0; j<3; j++){            
            grad_phi_data[I-3][j] = U[I]*detadx[j];
        }
    }

    //The required values for the material model
    //Assign required values
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    double t  = 0;
    double dt = 9;
    double params[18];
    std::vector<double> SDVS;
    
    //Output
    Vector_9  PK2;
    Vector_9  SIGMA;
    Vector_27 M;
    
    Vector_9  cauchy;
    Vector_9  s;
    Vector_27 m;
    
    
    define_parameters(params,false); //TODO: remove true!
    
    deformation_measures::get_deformation_gradient(grad_u, F);
    
    deformation_measures::assemble_chi(phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data, F, grad_chi);
    
    Matrix_3x9 grad_phi_mat;
    Vector_27  grad_phi;
    deformation_measures::assemble_grad_chi(grad_phi_data, Matrix_3x3::Identity(), grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    
    double fint[3];
    
    balance_equations::compute_internal_force(dNdx, cauchy, fint);
    
    std::vector<double> result;
    result.resize(3);
    for (int i=0; i<3; i++){result[i] = fint[i];}
    return result;
    
}

std::vector<double> parse_balance_of_first_moment_of_momentum_U(std::vector<double> U){
    /*!=====================================================
    |    parse_balance_of_first_moment_of_momentum_U    |
    =====================================================
    
    Parse the balance of the first moment of momentum as a function of the 
    degree of freedom vector.
    
    */
    
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);

    double eta;
    define_eta(eta);

    double detadx[3];
    define_detadx(detadx);
    
    //Expand U
    double grad_u[3][3];
    
    grad_u[0][0] = U[0]*detadx[0];
    grad_u[1][1] = U[1]*detadx[1];
    grad_u[2][2] = U[2]*detadx[2];
    grad_u[1][2] = U[1]*detadx[2];
    grad_u[0][2] = U[0]*detadx[2];
    grad_u[0][1] = U[0]*detadx[1];
    grad_u[2][1] = U[2]*detadx[1];
    grad_u[2][0] = U[2]*detadx[0];
    grad_u[1][0] = U[1]*detadx[0];

    double phi[9];
    for (int i=0; i<9; i++){phi[i] = U[i+3]*eta;}
    
    double grad_phi_data[9][3];
    
    for (int I=3; I<12; I++){
        for (int j=0; j<3; j++){            
            grad_phi_data[I-3][j] = U[I]*detadx[j];
        }
    }
    
    //The required values for the material model
    //Assign required values
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    double t  = 0;
    double dt = 9;
    double params[18];
    std::vector<double> SDVS;
    
    //Output
    Vector_9  PK2;
    Vector_9  SIGMA;
    Vector_27 M;
    
    Vector_9  cauchy;
    Vector_9  s;
    Vector_27 m;
    
    
    define_parameters(params,false); //TODO: set to false!
    
    deformation_measures::get_deformation_gradient(grad_u, F);
    
    deformation_measures::assemble_chi(phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data, F, grad_chi);
    
    Matrix_3x9 grad_phi_mat;
    Vector_27  grad_phi;
    deformation_measures::assemble_grad_chi(grad_phi_data, Matrix_3x3::Identity(), grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    
    double cint[9];
    
    balance_equations::compute_internal_couple(N, dNdx, cauchy, s, m, cint);
    
    std::vector<double> result;
    result.resize(9);
    for (int i=0; i<9; i++){result[i] = cint[i];}
    return result;
    
}

int test_micromorphic_material_library(){
    /*!
    ============================================
    |    test_micromorphic_material_library    |
    ============================================
    
    Test the version of the model that resides in the micromorphic material library 
    to make sure nothing was lost in translation.
    
    */

    //The DOF vector
    std::vector<double> U = {1.22,2.1,4.1,-2.3,.124,7.2,-8.2,.28,7.21,2.1,-9.2,3.1};

    //The shape function values
    double N;
    define_N(N);

    double dNdX[3];
    define_dNdx(dNdX);

    double eta;
    define_eta(eta);

    double detadX[3];
    define_detadx(detadX);

    //Compute the values explicitly
    double grad_u[3][3];

    grad_u[0][0] = U[0]*detadX[0];
    grad_u[1][1] = U[1]*detadX[1];
    grad_u[2][2] = U[2]*detadX[2];
    grad_u[1][2] = U[1]*detadX[2];
    grad_u[0][2] = U[0]*detadX[2];
    grad_u[0][1] = U[0]*detadX[1];
    grad_u[2][1] = U[2]*detadX[1];
    grad_u[2][0] = U[2]*detadX[0];
    grad_u[1][0] = U[1]*detadX[0];

    double phi[9];
    for (int i=0; i<9; i++){phi[i] = U[i+3]*eta;}

    double grad_phi_data[9][3];

    for (int I=3; I<12; I++){
        for (int j=0; j<3; j++){
            grad_phi_data[I-3][j] = U[I]*detadX[j];
        }
    }

    //The required values for the material model
    //Assign required values
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;

    double t  = 0;
    double dt = 9;
    double params[18];
    std::vector<double> SDVS;
    
    //Output
    double fint[3];
    double cint[9];

    Vector_9  PK2;
    Vector_9  SIGMA;
    Vector_27 M;

    Matrix_9x9   DPK2Dgrad_u;
    Matrix_9x9   DPK2Dphi;
    Matrix_9x27  DPK2Dgrad_phi;
    Matrix_9x9   DSIGMADgrad_u;
    Matrix_9x9   DSIGMADphi;
    Matrix_9x27  DSIGMADgrad_phi;
    Matrix_27x9  DMDgrad_u;
    Matrix_27x9  DMDphi;
    Matrix_27x27 DMDgrad_phi;

    define_parameters(params);

    deformation_measures::get_deformation_gradient(grad_u, F, false);

    deformation_measures::assemble_chi(phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data, grad_chi);

    Matrix_3x9 grad_phi_mat;
    Vector_27  grad_phi;
    deformation_measures::assemble_grad_chi(grad_phi_data, Matrix_3x3::Identity(), grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);

    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               DPK2Dgrad_u,   DPK2Dphi,   DPK2Dgrad_phi,
                               DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                               DMDgrad_u,     DMDphi,     DMDgrad_phi);

    //Compute the values using the version stored in the material library.
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();

    auto material = factory.GetMaterial("LinearElasticity");

    //Output

    Vector_9  _PK2;
    Vector_9  _SIGMA;
    Vector_27 _M;

    Matrix_9x9   _DPK2Dgrad_u;
    Matrix_9x9   _DPK2Dphi;
    Matrix_9x27  _DPK2Dgrad_phi;
    Matrix_9x9   _DSIGMADgrad_u;
    Matrix_9x9   _DSIGMADphi;
    Matrix_9x27  _DSIGMADgrad_phi;
    Matrix_27x9  _DMDgrad_u;
    Matrix_27x9  _DMDphi;
    Matrix_27x27 _DMDgrad_phi;

    //Vector output
    std::vector<std::vector<double>> _F_v;
    std::vector<std::vector<double>> _chi_v;

    balance_equations::map_eigen_to_vector(F,_F_v);
    balance_equations::map_eigen_to_vector(chi,_chi_v);

    std::vector<double> _PK2_v;
    std::vector<double> _SIGMA_v;
    std::vector<double> _M_v;

    std::vector<std::vector<double>> _DPK2Dgrad_u_v;
    std::vector<std::vector<double>> _DPK2Dphi_v;
    std::vector<std::vector<double>> _DPK2Dgrad_phi_v;
    std::vector<std::vector<double>> _DSIGMADgrad_u_v;
    std::vector<std::vector<double>> _DSIGMADphi_v;
    std::vector<std::vector<double>> _DSIGMADgrad_phi_v;
    std::vector<std::vector<double>> _DMDgrad_u_v;
    std::vector<std::vector<double>> _DMDphi_v;
    std::vector<std::vector<double>> _DMDgrad_phi_v;

    std::vector<std::vector<double>> _DfintDU_v;
    std::vector<std::vector<double>> _DcintDU_v;
    
    std::vector<double> time;
    time.resize(2);
    time[0] = t;
    time[1] = dt;

    std::vector<double> fparams;
    fparams.resize(18);
    for (int i=0; i<18; i++){fparams[i] = params[i];}

    std::vector<double> ADDDOF;
    std::vector<std::vector<double>> ADD_grad_DOF;
    std::vector<Eigen::VectorXd> ADD_TERMS;
    std::vector<Eigen::MatrixXd> ADD_JACOBIANS;

    std::vector<std::vector<double>>              ADD_TERMS_v;
    std::vector<std::vector<std::vector<double>>> ADD_JACOBIANS_v;

    //Compare the result from the stress calculation alone
    material->evaluate_model(time, fparams, grad_u, phi, grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _PK2,  _SIGMA,     _M,  ADD_TERMS);

    bool tot_result = true;

    tot_result *= PK2.isApprox(_PK2);
    tot_result *= SIGMA.isApprox(_SIGMA);
    tot_result *= M.isApprox(_M);

//    std::cout << "tot_result: " << tot_result << "\n";

    //Compare the result from the vector form of the stress calculation alone
    material->evaluate_model(time,    fparams,  grad_u, phi,       grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _PK2_v, _SIGMA_v,          _M_v, ADD_TERMS_v);

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(PK2[i]-_PK2_v[i]);
    }

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(SIGMA[i]-_SIGMA_v[i]);
    }

    for (int i=0; i<27; i++){
        tot_result *= 1e-6>fabs(M[i]-_M_v[i]);
    }

//    std::cout << "tot_result: " << tot_result << "\n";

    //Compare the results from the jacobian calculation.
    material->evaluate_model(time, fparams, grad_u, phi, grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _PK2, _SIGMA, _M,
                             _DPK2Dgrad_u,   _DPK2Dphi,   _DPK2Dgrad_phi,
                             _DSIGMADgrad_u, _DSIGMADphi, _DSIGMADgrad_phi,
                             _DMDgrad_u,     _DMDphi,     _DMDgrad_phi,
                             ADD_TERMS,      ADD_JACOBIANS);

    tot_result *= PK2.isApprox(_PK2);
    tot_result *= SIGMA.isApprox(_SIGMA);
    tot_result *= M.isApprox(_M);

    tot_result *= DPK2Dgrad_u.isApprox(_DPK2Dgrad_u);
    tot_result *= DPK2Dphi.isApprox(_DPK2Dphi); //d(x)dchi = D(x)Dphi
    tot_result *= DPK2Dgrad_phi.isApprox(_DPK2Dgrad_phi);

    tot_result *= DSIGMADgrad_u.isApprox(_DSIGMADgrad_u);
    tot_result *= DSIGMADphi.isApprox(_DSIGMADphi); //d(x)dchi = D(x)Dphi
    tot_result *= DSIGMADgrad_phi.isApprox(_DSIGMADgrad_phi);
    
    tot_result *= DMDgrad_u.isApprox(_DMDgrad_u);
    tot_result *= DMDphi.isApprox(_DMDphi); //d(x)dchi = D(x)Dphi
    tot_result *= DMDgrad_phi.isApprox(_DMDgrad_phi);

    //Compare the result from the vector form of the jacobian calculation
    material->evaluate_model(time, fparams, grad_u, phi, grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _PK2_v, _SIGMA_v, _M_v,
                             _DPK2Dgrad_u_v,   _DPK2Dphi_v,   _DPK2Dgrad_phi_v,
                             _DSIGMADgrad_u_v, _DSIGMADphi_v, _DSIGMADgrad_phi_v,
                             _DMDgrad_u_v,     _DMDphi_v,     _DMDgrad_phi_v,
                             ADD_TERMS_v,      ADD_JACOBIANS_v);

    //Now compute the values for the balance equations
    balance_equations::compute_internal_force(dNdX, _F_v, _PK2_v, fint);
    balance_equations::compute_internal_couple(N, dNdX, _F_v, _chi_v, _PK2_v, _SIGMA_v, _M_v, cint);
    
    balance_equations::compute_internal_force_jacobian(N, dNdX, eta, detadX,
                                                       _F_v,
                                                       _PK2_v,
                                                       _DPK2Dgrad_u_v, _DPK2Dphi_v, _DPK2Dgrad_phi_v,
                                                       _DfintDU_v);

    balance_equations::compute_internal_couple_jacobian(N, dNdX, eta, detadX,
                                                        _F_v, _chi_v,
                                                        _PK2_v, _SIGMA_v, _M_v,
                                                        _DPK2Dgrad_u_v, _DPK2Dphi_v, _DPK2Dgrad_phi_v,
                                                        _DSIGMADgrad_u_v, _DSIGMADphi_v, _DSIGMADgrad_phi_v,
                                                        _DMDgrad_u_v, _DMDphi_v, _DMDgrad_phi_v,
                                                        _DcintDU_v);

    double tmp;

    for (int i=0; i<3; i++){
        for (int A=0; A<12; A++){
            balance_equations::compute_internal_force_jacobian(i, A
                                                               N, dNdX, eta, detadX,
                                                               _F_v,
                                                               _PK2_v,
                                                               _DPK2Dgrad_u_v, _DPK2Dphi_v, _DPK2Dgrad_phi_v,
                                                               tmp);

        }
    }

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            for (int A=0; A<12; A++){
                balance_equations::compute_internal_couple_jacobian(i, j, A,
                                                                    N, dNdX, eta, detadX,
                                                                    _F_v, _chi_v,
                                                                    _PK2_v, _SIGMA_v, _M_v,
                                                                    _DPK2Dgrad_u_v, _DPK2Dphi_v, _DPK2Dgrad_phi_v,
                                                                    _DSIGMADgrad_u_v, _DSIGMADphi_v, _DSIGMADgrad_phi_v,
                                                                    _DMDgrad_u_v, _DMDphi_v, _DMDgrad_phi_v,
                                                                    tmp);
            }
        }
    }

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(PK2[i]-_PK2_v[i]);
    }

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(SIGMA[i]-_SIGMA_v[i]);
    }

    for (int i=0; i<27; i++){
        tot_result *= 1e-6>fabs(M[i]-_M_v[i]);
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DPK2Dgrad_u(i,j) - _DPK2Dgrad_u_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DPK2Dphi(i,j) - _DPK2Dphi_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<27; j++){
            tot_result *= 1e-6>fabs(DPK2Dgrad_phi(i,j) - _DPK2Dgrad_phi_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DSIGMADgrad_u(i,j) - _DSIGMADgrad_u_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DSIGMADphi(i,j) - _DSIGMADphi_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<27; j++){
            tot_result *= 1e-6>fabs(DSIGMADgrad_phi(i,j) - _DSIGMADgrad_phi_v[i][j]);
        }
    }

    for (int i=0; i<27; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DMDgrad_u(i,j) - _DMDgrad_u_v[i][j]);
        }
    }

    for (int i=0; i<27; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DMDphi(i,j) - _DMDphi_v[i][j]);
        }
    }

    for (int i=0; i<27; i++){
        for (int j=0; j<27; j++){
            tot_result *= 1e-6>fabs(DMDgrad_phi(i,j) - _DMDgrad_phi_v[i][j]);
        }
    }

}

int main(){
    /*!==========================
    |         main            |
    ===========================
    
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or 
    False if the test passes or fails respectively.*/
    
    //!Test of the instance of the model in the material library
    test_micromorphic_material_library();
}

