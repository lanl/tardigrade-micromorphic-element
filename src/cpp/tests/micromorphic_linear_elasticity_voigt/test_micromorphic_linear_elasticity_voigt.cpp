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

void define_parameters(double (&params)[18]){
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
    
    N = 0.261;
    return;
}

void define_dNdx(double (&dNdx)[3]){
    /*!=====================
    |    define_dNdx    |
    =====================
    
    Define the gradient of the shape function 
    to be used.
    */
    
    dNdx[0] =  1.42;
    dNdx[1] =  0.271;
    dNdx[2] = -2.31;
    return;
}

void define_eta(double &eta){
    /*!====================
    |    define_eta    |
    ====================

    Define the value of the interpolation function to be used.

    */

    eta = 0.826;
    return;
}

void define_detadx(double (&detadx)[3]){
    /*!=======================
    |    define_detadx    |
    =======================

    Define the value of the gradient of the interpolation function.

    */

    detadx[0] = 0.172;
    detadx[1] = 3.121;
    detadx[2] = 0.761;
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
    
    
    define_parameters(params);
    
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
    
    
    define_parameters(params);
    
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

int test_compute_A_voigt(std::ofstream &results){
    /*!==============================
    |    test_compute_A_voigt    |
    ==============================

    Test the computation of the forth order A 
    stiffness tensor in voigt form.

    */

    SpMat  A( 9, 9); //The expected result
    SpMat _A( 9, 9); //The result from the function
    Matrix_9x9 _A_alt; //The alternative calculation.

    define_A(A);     //Set the expected result
    micro_material::compute_A_voigt(0.191519450379, 0.62210877104, _A);
    micro_material::compute_A_voigt(0.191519450379, 0.62210877104, _A_alt);

    bool tot_result = A.isApprox(_A);
    Matrix_9x9 A_dense;
    A_dense = A;
    tot_result *= A_dense.isApprox(_A_alt);

    if (tot_result){
        results << "test_compute_A_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_A_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_B_voigt(std::ofstream &results){
    /*!==============================
    |    test_compute_B_voigt    |
    ==============================

    Test the computation of the forth order B
    stiffness tensor in voigt form.

    */

    SpMat  B( 9, 9); //The expected result
    SpMat _B( 9, 9); //The result from the function
    Matrix_9x9 _B_alt; //The alternative calculation

    define_B(B);     //Set the expected result
    micro_material::compute_B_voigt(0.4377277390071145, 0.7799758081188035,
                                    0.2725926052826416, 0.2764642551430967,
                                    0.7853585837137692, _B);
    micro_material::compute_B_voigt(0.4377277390071145, 0.7799758081188035,
                                    0.2725926052826416, 0.2764642551430967,
                                    0.7853585837137692, _B_alt);

    bool tot_result = B.isApprox(_B);
    Matrix_9x9 B_dense;
    B_dense = B;
    tot_result *= B_dense.isApprox(_B_alt);

    if (tot_result){
        results << "test_compute_B_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_B_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_C_voigt(std::ofstream &results){
    /*!==============================
    |    test_compute_C_voigt    |
    ==============================

    Test the computation of the sixth order C
    stiffness tensor in voigt form.

    */

    SpMat  C(27,27); //The expected result
    SpMat _C(27,27); //The result from the function
    Matrix_27x27 _C_alt; //The alternative calculation

    define_C(C);     //Set the expected result
    micro_material::compute_C_voigt(0.8018721775350193, 0.9581393536837052,
                                    0.8759326347420947, 0.35781726995786667,
                                    0.5009951255234587, 0.6834629351721363,
                                    0.7127020269829002, 0.37025075479039493,
                                    0.5611961860656249, 0.5030831653078097,
                                    0.013768449590682241, _C);
    micro_material::compute_C_voigt(0.8018721775350193, 0.9581393536837052,
                                    0.8759326347420947, 0.35781726995786667,
                                    0.5009951255234587, 0.6834629351721363,
                                    0.7127020269829002, 0.37025075479039493,
                                    0.5611961860656249, 0.5030831653078097,
                                    0.013768449590682241, _C_alt);

    bool tot_result = C.isApprox(_C);
    Matrix_27x27 C_dense;
    C_dense = C;
    tot_result *= C_dense.isApprox(_C_alt);

    if (tot_result){
        results << "test_compute_C_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_C_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_D_voigt(std::ofstream &results){
    /*!==============================
    |    test_compute_D_voigt    |
    ==============================

    Test the computation of the forth order D
    stiffness tensor in voigt form.

    */

    SpMat  D( 9, 9); //The expected result
    SpMat _D( 9, 9); //The result from the function
    Matrix_9x9 _D_alt; //The alternative calculation

    define_D(D);     //Set the expected result
    micro_material::compute_D_voigt(0.2764642551430967, 0.7853585837137692, _D);
    micro_material::compute_D_voigt(0.2764642551430967, 0.7853585837137692, _D_alt);

    bool tot_result = D.isApprox(_D);
    Matrix_9x9 D_dense;
    D_dense = D;
    tot_result *= D_dense.isApprox(_D_alt);

    if (tot_result){
        results << "test_compute_D_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_D_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_PK2_stress(std::ofstream &results){
    /*!=================================
    |    test_compute_PK2_stress    |
    =================================

    Test the computation of the PK2 stress 
    tensor in voigt form.

    */

    Vector_9  PK2; //The expected result
    Vector_9 _PK2; //The result of the function

    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9  SIGMA;          //The symmetric stress
    Vector_27 M;              //The higher order stress  

    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);  //Note: grad_chi = grad_phi

    define_PK2(PK2);
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, _PK2, SIGMA, M);
    
    bool tot_result = PK2.isApprox(_PK2,1e-6);

    if (tot_result){
        results << "test_compute_PK2_stress_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_PK2_stress_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_symmetric_stress(std::ofstream &results){
    /*!=======================================
    |    test_compute_symmetric_stress    |
    =======================================

    Test the computation of the symmetric stress 
    tensor in voigt form.

    */

    Vector_9  SIGMA; //The expected result
    Vector_9 _SIGMA; //The result of the function

    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_27 M;              //The higher order stress  

    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);  //Note: grad_chi = grad_phi

    define_SIGMA(SIGMA);
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, _SIGMA, M);
    
    bool tot_result = SIGMA.isApprox(_SIGMA,1e-6);

    if (tot_result){
        results << "test_compute_symmetric_stress_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_symmetric_stress_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_higher_order_stress(std::ofstream &results){
    /*!==========================================
    |    test_compute_higher_order_stress    |
    ==========================================

    Test the computation of the higher-order stress 
    tensor in voigt form.

    */

    Vector_27  M; //The expected result
    Vector_27 _M; //The result of the function

    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress

    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);  //Note: grad_chi = grad_phi

    define_M(M);
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, _M);
    
    bool tot_result = M.isApprox(_M,1e-6);

    if (tot_result){
        results << "test_compute_higher_order_stress_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_higher_order_stress_voigt & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPK2dRCG(std::ofstream &results){
    /*!===============================
    |    test_compute_dPK2dRCG    |
    ===============================
    
    Test the computation of the derivative of the second 
    Piola-Kirchoff stress w.r.t. the right cauchy-green 
    deformation tensor.
    
    */
    
    Matrix_9x9 dPK2dRCG; //The expected result
    std::vector< std::vector<double> > dPK2dRCG_vec; //The vector form.
    Matrix_9x9 _dPK2dRCG; //The result of the function.
    
    Matrix_3x3 F0; //The base point about which to compute the derivative
    define_deformation_gradient(F0);
    
    Matrix_3x3 RCG0 = F0.transpose()*F0;
    
    Vector_9 RCG0_vec; //The vector form of RCG
    deformation_measures::voigt_3x3_tensor(RCG0,RCG0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = RCG0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "dPK2dRCG\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_RCG, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dRCG_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPK2dRCG_vec.size(); i++){
        for (int j=0; j<dPK2dRCG_vec[i].size(); j++){
            dPK2dRCG(j,i) = dPK2dRCG_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Poputate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    Matrix_3x3 E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());
    Matrix_3x3 E_micro = F.transpose()*chi - Matrix_3x3::Identity();
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    
    deformation_measures::voigt_3x3_tensor(E,E_voigt);
    deformation_measures::voigt_3x3_tensor(E_micro,E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    //Evaluate the function
    t0 = Clock::now();
    micro_material::compute_dPK2dRCG(RCG0, RCG0.inverse(), Gamma, Gamma_voigt,
                                     E, E_micro, E_voigt, E_micro_voigt,
                                     A, B, C, D, _dPK2dRCG);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dPK2dRCG.isApprox(_dPK2dRCG,1e-6);

    Matrix_9x9   A_dense;
    Matrix_9x9   B_dense;
    Matrix_27x27 C_dense;
    Matrix_9x9   D_dense;

    A_dense = A;
    B_dense = B;
    C_dense = C;
    D_dense = D;

    t0 = Clock::now();
    micro_material::compute_dPK2dRCG(RCG0, RCG0.inverse(), Gamma, Gamma_voigt,
                                     E, E_micro, E_voigt, E_micro_voigt,
                                     A_dense, B_dense, C_dense, D_dense, _dPK2dRCG);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (dense): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    tot_result *= dPK2dRCG.isApprox(_dPK2dRCG,1e-6);
    
    if (tot_result){
        results << "test_compute_dPK2dRCG & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPK2dRCG & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dSIGMAdRCG(std::ofstream &results){
    /*!=================================
    |    test_compute_dSIGMAdRCG    |
    =================================
    
    Test the computation of the derivative of the symmetric 
    stress in the reference configuration w.r.t. the right cauchy-green 
    deformation tensor.
    
    */
    
    Matrix_9x9 dSIGMAdRCG; //The expected result
    std::vector< std::vector<double> > dSIGMAdRCG_vec; //The vector form.
    Matrix_9x9 _dSIGMAdRCG; //The result of the function.
    
    Matrix_3x3 F0; //The base point about which to compute the derivative
    define_deformation_gradient(F0);
    
    Matrix_3x3 RCG0 = F0.transpose()*F0;
    
    Vector_9 RCG0_vec; //The vector form of RCG
    deformation_measures::voigt_3x3_tensor(RCG0,RCG0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = RCG0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndSIGMAdRCG\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_RCG, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdRCG_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dSIGMAdRCG_vec.size(); i++){
        for (int j=0; j<dSIGMAdRCG_vec[i].size(); j++){
            dSIGMAdRCG(j,i) = dSIGMAdRCG_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Poputate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    Matrix_3x3 E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());
    Matrix_3x3 E_micro = F.transpose()*chi - Matrix_3x3::Identity();
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    
    deformation_measures::voigt_3x3_tensor(E,E_voigt);
    deformation_measures::voigt_3x3_tensor(E_micro,E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    //Evaluate the function
    Matrix_9x9 terms[4];
    Matrix_9x9 dPK2dRCG;
    t0 = Clock::now();
    micro_material::compute_dPK2dRCG(RCG0, RCG0.inverse(), Gamma, Gamma_voigt,
                                     E, E_micro, E_voigt, E_micro_voigt,
                                     A, B, C, D, terms, dPK2dRCG);
    micro_material::compute_dSIGMAdRCG(terms, _dSIGMAdRCG);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dSIGMAdRCG.isApprox(_dSIGMAdRCG,1e-6);

    Matrix_9x9   A_dense;
    Matrix_9x9   B_dense;
    Matrix_27x27 C_dense;
    Matrix_9x9   D_dense;

    A_dense = A;
    B_dense = B;
    C_dense = C;
    D_dense = D;

    t0 = Clock::now();
    micro_material::compute_dPK2dRCG(RCG0, RCG0.inverse(), Gamma, Gamma_voigt,
                                     E, E_micro, E_voigt, E_micro_voigt,
                                     A_dense, B_dense, C_dense, D_dense, terms, dPK2dRCG);
    micro_material::compute_dSIGMAdRCG(terms, _dSIGMAdRCG);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (dense): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    tot_result *= dSIGMAdRCG.isApprox(_dSIGMAdRCG,1e-6);
    
    if (tot_result){
        results << "test_compute_dSIGMAdRCG & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdRCG & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPK2dPsi(std::ofstream &results){
    /*!===============================
    |    test_compute_dPK2dPsi    |
    ===============================
    
    Test the computation of the derivative of the second 
    Piola-Kirchoff stress w.r.t. the deformation measures 
    Psi.
    
    */
    
    Matrix_9x9 dPK2dPsi; //The expected result
    std::vector< std::vector<double> > dPK2dPsi_vec; //The vector form.
    Matrix_9x9 _dPK2dPsi; //The result of the function.
    
    //Define the fundamental measures
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);
    
    //Define the derived deformation measures
    Matrix_3x3 RCG   = F.transpose()*F;
    Matrix_3x3 Psi0  = F.transpose()*chi; //The base point about which to compute the derivative
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    Vector_9 Phi0_vec; //The vector form of RCG
    deformation_measures::voigt_3x3_tensor(Psi0,Phi0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = Phi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndPK2dPsi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_Psi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dPsi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPK2dPsi_vec.size(); i++){
        for (int j=0; j<dPK2dPsi_vec[i].size(); j++){
            dPK2dPsi(j,i) = dPK2dPsi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Poputate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    
    Matrix_3x3 E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());
    Matrix_3x3 E_micro = F.transpose()*chi - Matrix_3x3::Identity();
    
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    
    deformation_measures::voigt_3x3_tensor(E,E_voigt);
    deformation_measures::voigt_3x3_tensor(E_micro,E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    //Evaluate the function
    t0 = Clock::now();
    micro_material::compute_dPK2dPsi(RCG.inverse(), E_micro, E_voigt, E_micro_voigt,
                                     B, D, _dPK2dPsi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    bool tot_result = dPK2dPsi.isApprox(_dPK2dPsi,1e-6);

    Matrix_9x9   A_dense;
    Matrix_9x9   B_dense;
    Matrix_27x27 C_dense;
    Matrix_9x9   D_dense;

    A_dense = A;
    B_dense = B;
    C_dense = C;
    D_dense = D;

    t0 = Clock::now();
    micro_material::compute_dPK2dPsi(RCG.inverse(), E_micro, E_voigt, E_micro_voigt,
                                     B_dense, D_dense, _dPK2dPsi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (dense): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    tot_result *= dPK2dPsi.isApprox(_dPK2dPsi,1e-6);
    
    if (tot_result){
        results << "test_compute_dPK2dPsi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPK2dPsi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dSIGMAdPsi(std::ofstream &results){
    /*!=================================
    |    test_compute_dSIGMAdPsi    |
    =================================
    
    Test the computation of the derivative of the symmetric 
    stress in the reference configuration w.r.t. the deformation 
    measure Psi.
    
    */
    
    Matrix_9x9 dSIGMAdPsi; //The expected result
    std::vector< std::vector<double> > dSIGMAdPsi_vec; //The vector form.
    Matrix_9x9 _dSIGMAdPsi; //The result of the function.
    
    //Define the fundamental measures
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);
    
    //Define the derived deformation measures
    Matrix_3x3 RCG   = F.transpose()*F;
    Matrix_3x3 Psi0  = F.transpose()*chi; //The base point about which to compute the derivative
    Matrix_3x9 Gamma = F.transpose()*grad_chi;
    
    Vector_9 Phi0_vec; //The vector form of RCG
    deformation_measures::voigt_3x3_tensor(Psi0,Phi0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = Phi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndSIGMAdPsi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_Psi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdPsi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dSIGMAdPsi_vec.size(); i++){
        for (int j=0; j<dSIGMAdPsi_vec[i].size(); j++){
            dSIGMAdPsi(j,i) = dSIGMAdPsi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Poputate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    
    Matrix_3x3 E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());
    Matrix_3x3 E_micro = F.transpose()*chi - Matrix_3x3::Identity();
    
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    
    deformation_measures::voigt_3x3_tensor(E,E_voigt);
    deformation_measures::voigt_3x3_tensor(E_micro,E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    //Evaluate the function
    Matrix_9x9 dPK2dPsi;
    Matrix_9x9 terms[3];
    t0 = Clock::now();
    micro_material::compute_dPK2dPsi(RCG.inverse(), E_micro, E_voigt, E_micro_voigt,
                                     B, D, terms, dPK2dPsi);
    micro_material::compute_dSIGMAdPsi(terms, _dSIGMAdPsi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dSIGMAdPsi.isApprox(_dSIGMAdPsi,1e-6);

    Matrix_9x9   A_dense;
    Matrix_9x9   B_dense;
    Matrix_27x27 C_dense;
    Matrix_9x9   D_dense;

    A_dense = A;
    B_dense = B;
    C_dense = C;
    D_dense = D;

    t0 = Clock::now();
    micro_material::compute_dPK2dPsi(RCG.inverse(), E_micro, E_voigt, E_micro_voigt,
                                     B_dense, D_dense, terms, dPK2dPsi);
    micro_material::compute_dSIGMAdPsi(terms, _dSIGMAdPsi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (dense): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    tot_result *= dSIGMAdPsi.isApprox(_dSIGMAdPsi,1e-6);

    if (tot_result){
        results << "test_compute_dSIGMAdPsi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdPsi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPK2dGamma(std::ofstream &results){
    /*!=================================
    |    test_compute_dPK2dGamma    |
    =================================
    
    Test the computation of the derivative of the second 
    Piola-Kirchoff stress w.r.t. the deformation measure
    Gamma.
    
    */
    
    Matrix_9x27 dPK2dGamma;                            //The expected result
    std::vector< std::vector<double> > dPK2dGamma_vec; //The vector form.
    Matrix_9x27 _dPK2dGamma;                           //The result of the function.
    
    //Define the fundamental measures
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);
    
    //Define the derived deformation measures
    Matrix_3x3 RCG    = F.transpose()*F;
    Matrix_3x3 Psi    = F.transpose()*chi; //The base point about which to compute the derivative
    Matrix_3x9 Gamma0 = F.transpose()*grad_chi;
    
    Vector_27 Gamma0_vec; //The vector form of RCG
    deformation_measures::voigt_3x9_tensor(Gamma0,Gamma0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = Gamma0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndPK2dGamma\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_Gamma, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dGamma_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPK2dGamma_vec.size(); i++){
        for (int j=0; j<dPK2dGamma_vec[i].size(); j++){
            dPK2dGamma(j,i) = dPK2dGamma_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Poputate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    
    Matrix_3x3 E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());
    Matrix_3x3 E_micro = F.transpose()*chi - Matrix_3x3::Identity();
    
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    
    deformation_measures::voigt_3x3_tensor(E,E_voigt);
    deformation_measures::voigt_3x3_tensor(E_micro,E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma0,Gamma_voigt);
    
    //Evaluate the function
    t0 = Clock::now();
    micro_material::compute_dPK2dGamma(RCG.inverse(), Gamma0, Gamma_voigt, C, _dPK2dGamma);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dPK2dGamma.isApprox(_dPK2dGamma,1e-6);
    
    if (tot_result){
        results << "test_compute_dPK2dGamma & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPK2dGamma & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dSIGMAdGamma(std::ofstream &results){
    /*!===================================
    |    test_compute_dSIGMAdGamma    |
    ===================================
    
    Test the computation of the derivative of the symmetric 
    stress in the reference configuration w.r.t. the 
    deformation measure Gamma.
    
    */
    
    Matrix_9x27 dSIGMAdGamma;                            //The expected result
    std::vector< std::vector<double> > dSIGMAdGamma_vec; //The vector form.
    Matrix_9x27 _dSIGMAdGamma;                           //The result of the function.
    
    //Define the fundamental measures
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);
    
    //Define the derived deformation measures
    Matrix_3x3 RCG    = F.transpose()*F;
    Matrix_3x3 Psi    = F.transpose()*chi; //The base point about which to compute the derivative
    Matrix_3x9 Gamma0 = F.transpose()*grad_chi;
    
    Vector_27 Gamma0_vec; //The vector form of RCG
    deformation_measures::voigt_3x9_tensor(Gamma0,Gamma0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = Gamma0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndSIGMAdGamma\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_Gamma, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdGamma_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dSIGMAdGamma_vec.size(); i++){
        for (int j=0; j<dSIGMAdGamma_vec[i].size(); j++){
            dSIGMAdGamma(j,i) = dSIGMAdGamma_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    
    Matrix_3x3 E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());
    Matrix_3x3 E_micro = F.transpose()*chi - Matrix_3x3::Identity();
    
    Vector_9 E_voigt;
    Vector_9 E_micro_voigt;
    Vector_27 Gamma_voigt;
    
    deformation_measures::voigt_3x3_tensor(E,E_voigt);
    deformation_measures::voigt_3x3_tensor(E_micro,E_micro_voigt);
    deformation_measures::voigt_3x9_tensor(Gamma0,Gamma_voigt);
    
    //Evaluate the function
    Matrix_9x27 dPK2dGamma;
    Matrix_9x27 terms[2];
    t0 = Clock::now();
    micro_material::compute_dPK2dGamma(RCG.inverse(), Gamma0, Gamma_voigt, C, terms, dPK2dGamma);
    micro_material::compute_dSIGMAdGamma(terms, _dSIGMAdGamma);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dSIGMAdGamma.isApprox(_dSIGMAdGamma,1e-6);
    
    if (tot_result){
        results << "test_compute_dSIGMAdGamma & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdGamma & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dMdGamma(std::ofstream &results){
    /*!===============================
    |    test_compute_dMdGamma    |
    ===============================
    
    Test the computation of the derivative of the higher order  
    stress in the reference configuration w.r.t. the 
    deformation measure Gamma.
    
    */
    
    Matrix_27x27 dMdGamma;                           //The expected result
    std::vector< std::vector<double> > dMdGamma_vec; //The vector form.
    Matrix_27x27 _dMdGamma;                          //The result of the function.
    
    //Define the fundamental measures
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_3x9 grad_chi;
    
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);
    
    //Define the derived deformation measures
    Matrix_3x3 RCG    = F.transpose()*F;
    Matrix_3x3 Psi    = F.transpose()*chi; //The base point about which to compute the derivative
    Matrix_3x9 Gamma0 = F.transpose()*grad_chi;
    
    Vector_27 Gamma0_vec; //The vector form of RCG
    deformation_measures::voigt_3x9_tensor(Gamma0,Gamma0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = Gamma0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndMdGamma\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_higher_order_stress_gamma, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dMdGamma_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dMdGamma_vec.size(); i++){
        for (int j=0; j<dMdGamma_vec[i].size(); j++){
            dMdGamma(j,i) = dMdGamma_vec[i][j];
        }
    }
    
    //Populate the stiffness matrices
    SpMat C(27,27);
    
    define_C(C);
    
    //Evaluate the function
    t0 = Clock::now();
    micro_material::compute_dMdGamma(Matrix_27x27(C), _dMdGamma);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dMdGamma.isApprox(_dMdGamma,1e-6);
    
    if (tot_result){
        results << "test_compute_dMdGamma & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dMdGamma & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPK2dF(std::ofstream &results){
    /*!=============================
    |    test_compute_dPK2dF    |
    =============================
    
    Test the computation of the derivative of the second 
    Piola-Kirchoff stress w.r.t. the deformation 
    gradient.
    
    */
    
    Matrix_9x9 dPK2dF; //The expected result
    std::vector< std::vector<double> > dPK2dF_vec; //The vector form.
    Matrix_9x9 _dPK2dF; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x3 F0; //The base point about which to compute the derivative
    define_deformation_gradient(F0);
    Vector_9 F0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F0,F0_vec);
    
    Matrix_3x3 RCG = F0.transpose()*F0;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndPK2dF\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dF_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPK2dF_vec.size(); i++){
        for (int j=0; j<dPK2dF_vec[i].size(); j++){
            dPK2dF(j,i) = dPK2dF_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F0,       chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               _dPK2dF,  dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dPK2dF.isApprox(_dPK2dF,1e-6);
    
    if (tot_result){
        results << "test_compute_dPK2dF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPK2dF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPK2dchi(std::ofstream &results){
    /*!===============================
    |    test_compute_dPK2dchi    |
    ===============================
    
    Test the computation of the derivative of the second 
    Piola-Kirchoff stress w.r.t. the micro-deformation tensor 
    chi.
    
    */
    
    Matrix_9x9 dPK2dchi; //The expected result
    std::vector< std::vector<double> > dPK2dchi_vec; //The vector form.
    Matrix_9x9 _dPK2dchi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    //Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x3 chi0; //The base point about which to compute the derivative
    define_chi(chi0);
    Vector_9 chi0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(chi0,chi0_vec);
    
    Matrix_3x3 F;
    define_deformation_gradient(F);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = chi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndPK2dchi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dchi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPK2dchi_vec.size(); i++){
        for (int j=0; j<dPK2dchi_vec[i].size(); j++){
            dPK2dchi(j,i) = dPK2dchi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi0,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   _dPK2dchi,  dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //std::cout << " dPK2dchi:\n" <<  dPK2dchi << "\n";
    //std::cout << "_dPK2dchi:\n" << _dPK2dchi << "\n";
    
    bool tot_result = dPK2dchi.isApprox(_dPK2dchi,1e-6);
    
    if (tot_result){
        results << "test_compute_dPK2dchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPK2dchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPK2dgrad_chi(std::ofstream &results){
    /*!====================================
    |    test_compute_dPK2dgrad_chi    |
    ====================================
    
    Test the computation of the derivative of the second 
    Piola-Kirchoff stress w.r.t. the gradient of the 
    micro-deformation tensor chi.
    
    */
    
    Matrix_9x27 dPK2dgrad_chi; //The expected result
    std::vector< std::vector<double> > dPK2dgrad_chi_vec; //The vector form.
    Matrix_9x27 _dPK2dgrad_chi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    //Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x9 grad_chi0; //The base point about which to compute the derivative
    define_grad_phi(grad_chi0);
    Vector_27 grad_chi0_vec; //The vector form of F
    deformation_measures::voigt_3x9_tensor(grad_chi0,grad_chi0_vec);
    
    Matrix_3x3 F;
    define_deformation_gradient(F);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = grad_chi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndPK2dgrad_chi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dgrad_chi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPK2dgrad_chi_vec.size(); i++){
        for (int j=0; j<dPK2dgrad_chi_vec[i].size(); j++){
            dPK2dgrad_chi(j,i) = dPK2dgrad_chi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof tensor
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi0,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   _dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dPK2dgrad_chi.isApprox(_dPK2dgrad_chi,1e-6);
    
    if (tot_result){
        results << "test_compute_dPK2dgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPK2dgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dSIGMAdF(std::ofstream &results){
    /*!=============================
    |    test_compute_dSIGMAdF    |
    =============================
    
    Test the computation of the derivative of the symmetric 
    stress in the reference configuration w.r.t. the deformation 
    gradient.
    
    */
    
    Matrix_9x9 dSIGMAdF;                             //The expected result
    std::vector< std::vector<double> > dSIGMAdF_vec; //The vector form.
    Matrix_9x9 _dSIGMAdF;                            //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    //Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x3 F0; //The base point about which to compute the derivative
    define_deformation_gradient(F0);
    Vector_9 F0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F0,F0_vec);
    
    Matrix_3x3 RCG = F0.transpose()*F0;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndSIGMAdF\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdF_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dSIGMAdF_vec.size(); i++){
        for (int j=0; j<dSIGMAdF_vec[i].size(); j++){
            dSIGMAdF(j,i) = dSIGMAdF_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,         dt,         params,  
                               F0,        chi,        grad_chi,
                               SDVS,      PK2,        SIGMA,    M,
                               dPK2dF,    dPK2dchi,   dPK2dgrad_chi,
                               _dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,      dMdchi,     dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dSIGMAdF.isApprox(_dSIGMAdF,1e-6);
    
    if (tot_result){
        results << "test_compute_dSIGMAdF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dSIGMAdchi(std::ofstream &results){
    /*!=================================
    |    test_compute_dSIGMAdchi    |
    =================================
    
    Test the computation of the derivative of the symmetric 
    stress in the reference configuration w.r.t. the micro-displacement
    tensor chi.
    
    */
    
    Matrix_9x9 dSIGMAdchi;                             //The expected result
    std::vector< std::vector<double> > dSIGMAdchi_vec; //The vector form.
    Matrix_9x9 _dSIGMAdchi;                            //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    //Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x3 chi0; //The base point about which to compute the derivative
    define_chi(chi0);
    Vector_9 chi0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(chi0,chi0_vec);
    
    Matrix_3x3 F;
    define_deformation_gradient(F);
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = chi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndSIGMAdchih\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdchi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dSIGMAdchi_vec.size(); i++){
        for (int j=0; j<dSIGMAdchi_vec[i].size(); j++){
            dSIGMAdchi(j,i) = dSIGMAdchi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,         dt,          params,  
                               F,         chi0,        grad_chi,
                               SDVS,      PK2,         SIGMA,    M,
                               dPK2dF,    dPK2dchi,    dPK2dgrad_chi,
                               dSIGMAdF,  _dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,      dMdchi,      dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dSIGMAdchi.isApprox(_dSIGMAdchi,1e-6);
    
    if (tot_result){
        results << "test_compute_dSIGMAdchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dSIGMAdgrad_chi(std::ofstream &results){
    /*!======================================
    |    test_compute_dSIGMAdgrad_chi    |
    ======================================
    
    Test the computation of the derivative of the symmetric 
    stress in the reference configuration w.r.t. the gradient of 
    the micro-displacement tensor chi.
    
    */
    
    Matrix_9x27 dSIGMAdgrad_chi;                            //The expected result
    std::vector< std::vector<double> > dSIGMAdgrad_chi_vec; //The vector form.
    Matrix_9x27 _dSIGMAdgrad_chi;                           //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    //Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x9 grad_chi0; //The base point about which to compute the derivative
    define_grad_phi(grad_chi0);
    Vector_27 grad_chi0_vec; //The vector form of F
    deformation_measures::voigt_3x9_tensor(grad_chi0,grad_chi0_vec);
    
    Matrix_3x3 F;
    define_deformation_gradient(F);
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = grad_chi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndSIGMAdgrad_chih\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdgrad_chi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dSIGMAdgrad_chi_vec.size(); i++){
        for (int j=0; j<dSIGMAdgrad_chi_vec[i].size(); j++){
            dSIGMAdgrad_chi(j,i) = dSIGMAdgrad_chi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    
    t0 = Clock::now();
    micro_material::get_stress(t,         dt,          params,  
                               F,         chi,         grad_chi0,
                               SDVS,      PK2,         SIGMA,    M,
                               dPK2dF,    dPK2dchi,    dPK2dgrad_chi,
                               dSIGMAdF,  dSIGMAdchi,  _dSIGMAdgrad_chi,
                               dMdF,      dMdchi,      dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    bool tot_result = dSIGMAdgrad_chi.isApprox(_dSIGMAdgrad_chi,1e-6);
    
    if (tot_result){
        results << "test_compute_dSIGMAdgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dMdF(std::ofstream &results){
    /*!===========================
    |    test_compute_dMdF    |
    ===========================
    
    Test the computation of the derivative of the higher order 
    stress in the reference configuration w.r.t. the deformation 
    gradient.
    
    */
    
    Matrix_27x9 dMdF;                             //The expected result
    std::vector< std::vector<double> > dMdF_vec;  //The vector form.
    Matrix_27x9 _dMdF;                            //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    //Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x3 F0; //The base point about which to compute the derivative
    define_deformation_gradient(F0);
    Vector_9 F0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F0,F0_vec);
    
    Matrix_3x3 RCG = F0.transpose()*F0;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndMdF\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_higher_order_stress_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dMdF_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dMdF_vec.size(); i++){
        for (int j=0; j<dMdF_vec[i].size(); j++){
            dMdF(j,i) = dMdF_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F0,       chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               _dMdF,    dMdchi,     dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dMdF.isApprox(_dMdF,1e-6);
    
    if (tot_result){
        results << "test_compute_dMdF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dMdF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dMdchi(std::ofstream &results){
    /*!=============================
    |    test_compute_dMdchi    |
    =============================
    
    Test the computation of the derivative of the higher order 
    stress in the reference configuration w.r.t. the micro-deformation
    tensor chi.
    
    */
    
    Matrix_27x9 dMdchi;                             //The expected result
    std::vector< std::vector<double> > dMdchi_vec;  //The vector form.
    Matrix_27x9 _dMdchi;                            //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    //Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x3 chi0;
    Matrix_3x3 F; //The base point about which to compute the derivative
    define_deformation_gradient(F);
    define_chi(chi0);
    Vector_9 chi0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(chi0,chi0_vec);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = chi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndMdchi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_higher_order_stress_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dMdchi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dMdchi_vec.size(); i++){
        for (int j=0; j<dMdchi_vec[i].size(); j++){
            dMdchi(j,i) = dMdchi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi0,       grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     _dMdchi,    dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dMdchi.isApprox(_dMdchi,1e-6);
    
    if (tot_result){
        results << "test_compute_dMdchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dMdchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dMdgrad_chi(std::ofstream &results){
    /*!==================================
    |    test_compute_dMdgrad_chi    |
    ==================================
    
    Test the computation of the derivative of the higher order 
    stress in the reference configuration w.r.t. the gradient 
    of the micro-deformation tensor chi.
    
    */
    
    Matrix_27x27 dMdgrad_chi;                            //The expected result
    std::vector< std::vector<double> > dMdgrad_chi_vec;  //The vector form.
    Matrix_27x27 _dMdgrad_chi;                           //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    //Matrix_27x27 dMdgrad_chi;
    
    Matrix_3x9 grad_chi0;
    Matrix_3x3 F; //The base point about which to compute the derivative
    define_deformation_gradient(F);
    define_grad_phi(grad_chi0);
    Vector_27 grad_chi0_vec; //The vector form of F
    deformation_measures::voigt_3x9_tensor(grad_chi0,grad_chi0_vec);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = grad_chi0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndMdgrad_chi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_higher_order_stress_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dMdgrad_chi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dMdgrad_chi_vec.size(); i++){
        for (int j=0; j<dMdgrad_chi_vec[i].size(); j++){
            dMdgrad_chi(j,i) = dMdgrad_chi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi0,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     _dMdgrad_chi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    bool tot_result = dMdgrad_chi.isApprox(_dMdgrad_chi,1e-6);
    
    if (tot_result){
        results << "test_compute_dMdgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dMdgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dcauchydF(std::ofstream &results){
    /*!================================
    |    test_compute_dcauchydF    |
    ================================
    
    Test the computation of the derivative of the cauchy 
    stress w.r.t. the deformation gradient.
    
    */
    
    Matrix_9x9 dcauchydF; //The expected result
    std::vector< std::vector<double> > dcauchydF_vec; //The vector form.
    Matrix_9x9 _dcauchydF; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x3 F; //The base point about which to compute the derivative
    define_deformation_gradient(F);
    Vector_9 x0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F,x0_vec);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndcauchydF\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_cauchy_stress_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dcauchydF_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dcauchydF_vec.size(); i++){
        for (int j=0; j<dcauchydF_vec[i].size(); j++){
            dcauchydF(j,i) = dcauchydF_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                _dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                      dsdF,      dsdchi,      dsdgrad_chi,
                                                      dmdF,      dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    //std::cout << " dcauchydF:\n" <<  dcauchydF << "\n";
    //std::cout << "_dcauchydF:\n" << _dcauchydF << "\n";
    
    bool tot_result = dcauchydF.isApprox(_dcauchydF,1e-6);
    
    if (tot_result){
        results << "test_compute_dcauchydF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dcauchydF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dsdF(std::ofstream &results){
    /*!===========================
    |    test_compute_dsdF    |
    ===========================
    
    Test the computation of the derivative of the symmetric  
    stress in the current configuration w.r.t. the deformation gradient.
    
    */
    
    Matrix_9x9 dsdF; //The expected result
    std::vector< std::vector<double> > dsdF_vec; //The vector form.
    Matrix_9x9 _dsdF; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
//    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x3 F; //The base point about which to compute the derivative
    define_deformation_gradient(F);
    Vector_9 x0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F,x0_vec);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndsdF\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_s_stress_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dsdF_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dsdF_vec.size(); i++){
        for (int j=0; j<dsdF_vec[i].size(); j++){
            dsdF(j,i) = dsdF_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                    _dsdF,      dsdchi,      dsdgrad_chi,
                                                      dmdF,      dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    //std::cout << " dsdF:\n" <<  dsdF << "\n";
    //std::cout << "_dsdF:\n" << _dsdF << "\n";
    
    bool tot_result = dsdF.isApprox(_dsdF,1e-6);
    
    if (tot_result){
        results << "test_compute_dsdF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dsdF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dmdF(std::ofstream &results){
    /*!===========================
    |    test_compute_dmdF    |
    ===========================
    
    Test the computation of the derivative of the higher order 
    stress in the current configuration w.r.t. the deformation gradient.
    
    */
    
    Matrix_27x9 dmdF; //The expected result
    std::vector< std::vector<double> > dmdF_vec; //The vector form.
    Matrix_27x9 _dmdF; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
//    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x3 F; //The base point about which to compute the derivative
    define_deformation_gradient(F);
    Vector_9 x0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F,x0_vec);
    
    Matrix_3x3 RCG = F.transpose()*F;
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndmdF\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_m_stress_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dmdF_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dmdF_vec.size(); i++){
        for (int j=0; j<dmdF_vec[i].size(); j++){
            dmdF(j,i) = dmdF_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 chi;           //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_chi(chi);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                     dsdF,      dsdchi,      dsdgrad_chi,
                                                    _dmdF,      dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    //std::cout << " dmdF:\n" <<  dmdF << "\n";
    //std::cout << "_dmdF:\n" << _dmdF << "\n";
    
    bool tot_result = dmdF.isApprox(_dmdF,1e-6);
    
    if (tot_result){
        results << "test_compute_dmdF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dmdF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dcauchydchi(std::ofstream &results){
    /*!==================================
    |    test_compute_dcauchydchi    |
    ==================================
    
    Test the computation of the derivative of the cauchy 
    stress w.r.t. the micro-displacement tensor.
    
    */
    
    Matrix_9x9 dcauchydchi; //The expected result
    std::vector< std::vector<double> > dcauchydchi_vec; //The vector form.
    Matrix_9x9 _dcauchydchi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
//    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x3 chi;
    define_chi(chi);
    Vector_9 x0_vec;
    deformation_measures::voigt_3x3_tensor(chi,x0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndcauchydchi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_cauchy_stress_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dcauchydchi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dcauchydchi_vec.size(); i++){
        for (int j=0; j<dcauchydchi_vec[i].size(); j++){
            dcauchydchi(j,i) = dcauchydchi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The micro-dof
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                 dcauchydF, _dcauchydchi, dcauchydgrad_chi,
                                                      dsdF,       dsdchi,      dsdgrad_chi,
                                                      dmdF,       dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //std::cout << " dcauchydchi:\n" <<  dcauchydchi << "\n";
    //std::cout << "_dcauchydchi:\n" << _dcauchydchi << "\n";
    
    bool tot_result = dcauchydchi.isApprox(_dcauchydchi,1e-6);
    
    if (tot_result){
        results << "test_compute_dcauchydchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dcauchydchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dsdchi(std::ofstream &results){
    /*!=============================
    |    test_compute_dsdchi    |
    =============================
    
    Test the computation of the derivative of the symmetric  
    stress in the current configuration w.r.t. the deformation gradient.
    
    */
    
    Matrix_9x9 dsdchi; //The expected result
    std::vector< std::vector<double> > dsdchi_vec; //The vector form.
    Matrix_9x9 _dsdchi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
//    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x3 chi; //The base point about which to compute the derivative
    define_chi(chi);
    Vector_9 x0_vec; //The vector form of chi
    deformation_measures::voigt_3x3_tensor(chi,x0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndsdchi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_s_stress_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dsdchi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dsdchi_vec.size(); i++){
        for (int j=0; j<dsdchi_vec[i].size(); j++){
            dsdchi(j,i) = dsdchi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                     dsdF,     _dsdchi,      dsdgrad_chi,
                                                     dmdF,      dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    //std::cout << " dsdchi:\n" <<  dsdchi << "\n";
    //std::cout << "_dsdchi:\n" << _dsdchi << "\n";
    
    bool tot_result = dsdchi.isApprox(_dsdchi,1e-6);
    
    if (tot_result){
        results << "test_compute_dsdchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dsdchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dmdchi(std::ofstream &results){
    /*!=============================
    |    test_compute_dmdchi    |
    =============================
    
    Test the computation of the derivative of the higher order 
    stress in the current configuration w.r.t. the micro-displacement
    tensor.
    
    */
    
    Matrix_27x9 dmdchi; //The expected result
    std::vector< std::vector<double> > dmdchi_vec; //The vector form.
    Matrix_27x9 _dmdchi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
//    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x3 chi; //The base point about which to compute the derivative
    define_chi(chi);
    Vector_9 x0_vec; //The vector form of chi
    deformation_measures::voigt_3x3_tensor(chi,x0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndmdchi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_m_stress_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dmdchi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dmdchi_vec.size(); i++){
        for (int j=0; j<dmdchi_vec[i].size(); j++){
            dmdchi(j,i) = dmdchi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x9 grad_chi;      //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_grad_phi(grad_chi); //Note: grad_phi == grad_chi
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                     dsdF,      dsdchi,      dsdgrad_chi,
                                                     dmdF,     _dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //std::cout << " dmdchi:\n" <<  dmdchi << "\n";
    //std::cout << "_dmdchi:\n" << _dmdchi << "\n";
    
    bool tot_result = dmdchi.isApprox(_dmdchi,1e-6);
    
    if (tot_result){
        results << "test_compute_dmdchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dmdchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dcauchydgrad_chi(std::ofstream &results){
    /*!=======================================
    |    test_compute_dcauchydgrad_chi    |
    =======================================
    
    Test the computation of the derivative of the cauchy 
    stress w.r.t. the gradient of the micro-displacement tensor.
    
    */
    
    Matrix_9x27 dcauchydgrad_chi; //The expected result
    std::vector< std::vector<double> > dcauchydgrad_chi_vec; //The vector form.
    Matrix_9x27 _dcauchydgrad_chi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
//    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x9 grad_chi;
    define_grad_phi(grad_chi);
    Vector_27 x0_vec;
    deformation_measures::voigt_3x9_tensor(grad_chi,x0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndcauchydgrad_chi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_cauchy_stress_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dcauchydgrad_chi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dcauchydgrad_chi_vec.size(); i++){
        for (int j=0; j<dcauchydgrad_chi_vec[i].size(); j++){
            dcauchydgrad_chi(j,i) = dcauchydgrad_chi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The micro-dof
    Matrix_3x3 chi;           //The micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                 dcauchydF,  dcauchydchi, _dcauchydgrad_chi,
                                                      dsdF,       dsdchi,       dsdgrad_chi,
                                                      dmdF,       dmdchi,       dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //std::cout << " dcauchydgrad_chi:\n" <<  dcauchydgrad_chi << "\n";
    //std::cout << "_dcauchydgrad_chi:\n" << _dcauchydgrad_chi << "\n";
    
    bool tot_result = dcauchydgrad_chi.isApprox(_dcauchydgrad_chi,1e-6);
    
    if (tot_result){
        results << "test_compute_dcauchydgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dcauchydgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dsdgrad_chi(std::ofstream &results){
    /*!==================================
    |    test_compute_dsdgrad_chi    |
    ==================================
    
    Test the computation of the derivative of the symmetric  
    stress in the current configuration w.r.t. the gradient 
    of the micro-displacement tensor.
    
    */
    
    Matrix_9x27 dsdgrad_chi; //The expected result
    std::vector< std::vector<double> > dsdgrad_chi_vec; //The vector form.
    Matrix_9x27 _dsdgrad_chi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
//    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x9 grad_chi; //The base point about which to compute the derivative
    define_grad_phi(grad_chi);
    Vector_27 x0_vec; //The vector form of chi
    deformation_measures::voigt_3x9_tensor(grad_chi,x0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndsdgrad_chi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_s_stress_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dsdgrad_chi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dsdgrad_chi_vec.size(); i++){
        for (int j=0; j<dsdgrad_chi_vec[i].size(); j++){
            dsdgrad_chi(j,i) = dsdgrad_chi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                     dsdF,      dsdchi,     _dsdgrad_chi,
                                                     dmdF,      dmdchi,      dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    //std::cout << " dsdgrad_chi:\n" <<  dsdgrad_chi << "\n";
    //std::cout << "_dsdgrad_chi:\n" << _dsdgrad_chi << "\n";
    
    bool tot_result = dsdgrad_chi.isApprox(_dsdgrad_chi,1e-6);
    
    if (tot_result){
        results << "test_compute_dsdgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dsdgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dmdgrad_chi(std::ofstream &results){
    /*!==================================
    |    test_compute_dmdgrad_chi    |
    ==================================
    
    Test the computation of the derivative of the higher order 
    stress in the current configuration w.r.t. the micro-displacement
    tensor.
    
    */
    
    Matrix_27x27 dmdgrad_chi; //The expected result
    std::vector< std::vector<double> > dmdgrad_chi_vec; //The vector form.
    Matrix_27x27 _dmdgrad_chi; //The result of the function.
    
    //Define the other gradients
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
//    Matrix_27x27 dmdgrad_chi;
    
    Matrix_3x9 grad_chi; //The base point about which to compute the derivative
    define_grad_phi(grad_chi);
    Vector_27 x0_vec; //The vector form of chi
    deformation_measures::voigt_3x9_tensor(grad_chi,x0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = x0_vec(i);}
    
    //Initialize the finite difference operator
    std::cout << "\ndmdgrad_chi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_m_stress_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dmdgrad_chi_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dmdgrad_chi_vec.size(); i++){
        for (int j=0; j<dmdgrad_chi_vec[i].size(); j++){
            dmdgrad_chi(j,i) = dmdgrad_chi_vec[i][j];
        }
    }
    
    //Obtain the required values
    double     t = 0.;        //The current time value (unneeded)
    double    dt = 0.;        //The change in time (unneeded)
    double params[18];        //The material parameters
    Matrix_3x3 F;             //The deformation gradient
    Matrix_3x3 chi;           //The gradient of the micro-dof
    std::vector<double> SDVS; //The state variables (unneeded)
    Vector_9 PK2;             //The second piola kirchhoff stress
    Vector_9 SIGMA;           //The symmetric stress
    Vector_27 M;              //The higher order stress

    //The jacobians in the reference configuration
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //The stresses in the current configuation
    Vector_9  cauchy; //The cauchy stress
    Vector_9  s;      //The symmetric stress
    Vector_27 m;      //The higher-order stress
    
    //Populate the stiffness matrices
    SpMat A(9,9);
    SpMat B(9,9);
    SpMat C(27,27);
    SpMat D(9,9);
    
    define_A(A);
    define_B(B);
    define_C(C);
    define_D(D);
    
    //Populate the required values
    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    
    t0 = Clock::now();
    micro_material::get_stress(t,        dt,         params,  
                               F,        chi,        grad_chi,
                               SDVS,     PK2,        SIGMA,    M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
    deformation_measures::map_jacobians_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m,
                                                dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi, dMdF, dMdchi, dMdgrad_chi,
                                                dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                     dsdF,      dsdchi,      dsdgrad_chi,
                                                     dmdF,      dmdchi,     _dmdgrad_chi);

    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //std::cout << " dmdgrad_chi:\n" <<  dmdgrad_chi << "\n";
    //std::cout << "_dmdgrad_chi:\n" << _dmdgrad_chi << "\n";
    
    bool tot_result = dmdgrad_chi.isApprox(_dmdgrad_chi,1e-6);
    
    if (tot_result){
        results << "test_compute_dmdgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dmdgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_total_derivatives_DcauchyDgrad_u(std::ofstream &results){
    /*!=======================================================
    |    test_compute_total_derivatives_DcauchyDgrad_u    |
    =======================================================
 
    Test the computation of the total derivatives of the stress 
    measures w.r.t. the gradients and magnitudes of the degrees of freedom.
 
    testing DcauchyDgrad_u
 
    */
    
    const int mhat = 9;
    const int nhat = 9;
    
    Matrix_9x9  r; //The expected result
    Matrix_9x9 _r; //The function output
    
    //The gradient of the macro-deformation
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    
    double _grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }
    
    //The value of the micro-displacement
    Matrix_3x3 phi;
    define_phi(phi);
    
    Vector_9 phi_voigt;
    double _phi[9];
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    for (int i=0; i<9; i++){_phi[i] = phi_voigt(i);}
    
    
    //The value of the gradient of the micro-displacement tensor in the current configuration.
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
    
    Matrix_3x3 F0;
    define_deformation_gradient(F0);
    
    //Define the value to compute the gradient w.r.t.
    Vector_9 voigt;
    deformation_measures::voigt_3x3_tensor(grad_u,voigt);
    
    //Initialize the finite difference gradient calculator
    std::vector<double> x0;
    x0.resize(nhat);
    
    for (int i=0; i<nhat; i++){ x0[i] = voigt[i];}
    
    std::cout << "\ndcauchydgrad_u\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_cauchy_stress_grad_u, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the gradient analytically
    
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Vector_27  grad_phi;
    Matrix_3x9 grad_phi_mat;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(_grad_u, F);
    
    deformation_measures::assemble_chi(_phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array,Matrix_3x3::Identity(),grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
                              
    _r = DcauchyDgrad_u;
    
    //std::cout << "r:\n"  << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_total_derivatives_DcauchyDgrad_u & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_total_derivatives_DcauchyDgrad_u & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_total_derivatives_DsDgrad_u(std::ofstream &results){
    /*!==================================================
    |    test_compute_total_derivatives_DsDgrad_u    |
    ==================================================
 
    Test the computation of the total derivatives of the stress 
    measures w.r.t. the gradients and magnitudes of the degrees of freedom.
 
    testing DsDgrad_u
 
    */
    
    const int mhat = 9;
    const int nhat = 9;
    
    Matrix_9x9  r; //The expected result
    Matrix_9x9 _r; //The function output
    
    //The gradient of the macro-deformation
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    
    double _grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }
    
    //The value of the micro-displacement
    Matrix_3x3 phi;
    define_phi(phi);
    
    Vector_9 phi_voigt;
    double _phi[9];
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    for (int i=0; i<9; i++){_phi[i] = phi_voigt(i);}
    
    
    //The value of the gradient of the micro-displacement tensor in the current configuration.
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
    
    //Define the value to compute the gradient w.r.t.
    Vector_9 voigt;
    deformation_measures::voigt_3x3_tensor(grad_u,voigt);
    
    //Initialize the finite difference gradient calculator
    std::vector<double> x0;
    x0.resize(nhat);
    
    for (int i=0; i<nhat; i++){ x0[i] = voigt[i];}
    
    std::cout << "\ndsdgrad_u\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_s_stress_grad_u, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the gradient analytically
    
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Vector_27  grad_phi;
    Matrix_3x9 grad_phi_mat;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(_grad_u, F);
    
    deformation_measures::assemble_chi(_phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array,Matrix_3x3::Identity(),grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
                              
    _r = DsDgrad_u;
    
    //std::cout << "r:\n"  << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_total_derivatives_DsDgrad_u & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_total_derivatives_DsDgrad_u & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_total_derivatives_DmDgrad_u(std::ofstream &results){
    /*!==================================================
    |    test_compute_total_derivatives_DmDgrad_u    |
    ==================================================
 
    Test the computation of the total derivatives of the stress 
    measures w.r.t. the gradients and magnitudes of the degrees of freedom.
 
    testing DmDgrad_u
 
    */
    
    const int mhat = 27;
    const int nhat = 9;
    
    Matrix_27x9  r; //The expected result
    Matrix_27x9 _r; //The function output
    
    //The gradient of the macro-deformation
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    
    double _grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }
    
    //The value of the micro-displacement
    Matrix_3x3 phi;
    define_phi(phi);
    
    Vector_9 phi_voigt;
    double _phi[9];
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    for (int i=0; i<9; i++){_phi[i] = phi_voigt(i);}
    
    
    //The value of the gradient of the micro-displacement tensor in the current configuration.
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
    
    //Define the value to compute the gradient w.r.t.
    Vector_9 voigt;
    deformation_measures::voigt_3x3_tensor(grad_u,voigt);
    
    //Initialize the finite difference gradient calculator
    std::vector<double> x0;
    x0.resize(nhat);
    
    for (int i=0; i<nhat; i++){ x0[i] = voigt[i];}
    
    std::cout << "\ndmdgrad_u\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_m_stress_grad_u, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the gradient analytically
    
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Vector_27  grad_phi;
    Matrix_3x9 grad_phi_mat;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(_grad_u, F);
    
    deformation_measures::assemble_chi(_phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array,Matrix_3x3::Identity(),grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
                              
    _r = DmDgrad_u;
    
    //std::cout << "r:\n"  << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_total_derivatives_DmDgrad_u & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_total_derivatives_DmDgrad_u & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_total_derivatives_DcauchyDgrad_phi(std::ofstream &results){
    /*!=========================================================
    |    test_compute_total_derivatives_DcauchyDgrad_phi    |
    =========================================================
 
    Test the computation of the total derivatives of the stress 
    measures w.r.t. the gradients and magnitudes of the degrees of freedom.
 
    testing DcauchyDgrad_phi
 
    */
    
    const int mhat = 9;
    const int nhat = 27;
    
    Matrix_9x27  r; //The expected result
    Matrix_9x27 _r; //The function output
    
    //The gradient of the macro-deformation
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    
    double _grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }
    
    //The value of the micro-displacement
    Matrix_3x3 phi;
    define_phi(phi);
    
    Vector_9 phi_voigt;
    double _phi[9];
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    for (int i=0; i<9; i++){_phi[i] = phi_voigt(i);}
    
    
    //The value of the gradient of the micro-displacement tensor in the current configuration.
    Vector_27 grad_phi;
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
    
    Matrix_3x9 grad_phi_mat;
    deformation_measures::assemble_grad_chi(grad_phi_data_array, Matrix_3x3::Identity(), grad_phi_mat);
    
    //Define the value to compute the gradient w.r.t.
    Vector_27 voigt;
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,voigt);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    //Initialize the finite difference gradient calculator
    std::vector<double> x0;
    x0.resize(nhat);
    
    for (int i=0; i<nhat; i++){ x0[i] = voigt[i];}
    
    std::cout << "\ndcauchydgrad_phi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_cauchy_stress_grad_phi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the gradient analytically
    
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(_grad_u, F);
    
    deformation_measures::assemble_chi(_phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
                              
    _r = DcauchyDgrad_phi;
    
    //std::cout << "r:\n"  << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_total_derivatives_DcauchyDgrad_phi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_total_derivatives_DcauchyDgrad_phi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_total_derivatives_DsDgrad_phi(std::ofstream &results){
    /*!=========================================================
    |    test_compute_total_derivatives_DsDgrad_phi    |
    =========================================================
 
    Test the computation of the total derivatives of the stress 
    measures w.r.t. the gradients and magnitudes of the degrees of freedom.
 
    testing DsDgrad_phi
 
    */
    
    const int mhat = 9;
    const int nhat = 27;
    
    Matrix_9x27  r; //The expected result
    Matrix_9x27 _r; //The function output
    
    //The gradient of the macro-deformation
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    
    double _grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }
    
    //The value of the micro-displacement
    Matrix_3x3 phi;
    define_phi(phi);
    
    Vector_9 phi_voigt;
    double _phi[9];
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    for (int i=0; i<9; i++){_phi[i] = phi_voigt(i);}
    
    
    //The value of the gradient of the micro-displacement tensor in the current configuration.
    Vector_27 grad_phi;
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
    
    Matrix_3x9 grad_phi_mat;
    deformation_measures::assemble_grad_chi(grad_phi_data_array, Matrix_3x3::Identity(), grad_phi_mat);
    
    //Define the value to compute the gradient w.r.t.
    Vector_27 voigt;
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,voigt);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    //Initialize the finite difference gradient calculator
    std::vector<double> x0;
    x0.resize(nhat);
    
    for (int i=0; i<nhat; i++){ x0[i] = voigt[i];}
    
    std::cout << "\ndsdgrad_phi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_s_stress_grad_phi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the gradient analytically
    
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(_grad_u, F);
    
    deformation_measures::assemble_chi(_phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
                                                    
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
                              
    _r = DsDgrad_phi;
    
    //std::cout << "r:\n"  << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_total_derivatives_DsDgrad_phi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_total_derivatives_DsDgrad_phi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_total_derivatives_DmDgrad_phi(std::ofstream &results){
    /*!=========================================================
    |    test_compute_total_derivatives_DsDgrad_phi    |
    =========================================================
 
    Test the computation of the total derivatives of the stress 
    measures w.r.t. the gradients and magnitudes of the degrees of freedom.
 
    testing DmDgrad_phi
 
    */
    
    const int mhat = 27;
    const int nhat = 27;
    
    Matrix_27x27  r; //The expected result
    Matrix_27x27 _r; //The function output
    
    //The gradient of the macro-deformation
    Matrix_3x3 grad_u;
    define_grad_u(grad_u);
    
    double _grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }
    
    //The value of the micro-displacement
    Matrix_3x3 phi;
    define_phi(phi);
    
    Vector_9 phi_voigt;
    double _phi[9];
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    for (int i=0; i<9; i++){_phi[i] = phi_voigt(i);}
    
    
    //The value of the gradient of the micro-displacement tensor in the current configuration.
    Vector_27 grad_phi;
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
    
    Matrix_3x9 grad_phi_mat;
    deformation_measures::assemble_grad_chi(grad_phi_data_array, Matrix_3x3::Identity(), grad_phi_mat);
    
    //Define the value to compute the gradient w.r.t.
    Vector_27 voigt;
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,voigt);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    //Initialize the finite difference gradient calculator
    std::vector<double> x0;
    x0.resize(nhat);
    
    for (int i=0; i<nhat; i++){ x0[i] = voigt[i];}
    
    std::cout << "\ndmdgrad_phi\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    auto t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_m_stress_grad_phi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    auto t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the gradient analytically
    
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(_grad_u, F);
    
    deformation_measures::assemble_chi(_phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data_array, F, grad_chi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
                                                    
    t1 = Clock::now();
    std::cout << "Analytic Jacobian (includes all other jacobian and stress calculations): " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
                              
    _r = DmDgrad_phi;
    
    //std::cout << "r:\n"  << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_total_derivatives_DmDgrad_phi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_total_derivatives_DmDgrad_phi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_internal_force(std::ofstream &results){
    /*!
    =====================================
    |    test_compute_internal_force    |
    =====================================

    Test the computation of the internal force. 

    */

    double  r[3]; //The expected result
    double _r[3]; //The function output

    //Define the test function
    double dNdx[3];
    define_dNdx(dNdx);

    //Define the stress tensor as an eigen matrix
    Vector_9 cauchy;
    define_cauchy(cauchy);

    //Define the stress tensor as a std::vector
    std::vector<double> cauchy_vec;
    balance_equations::map_eigen_to_vector(cauchy,cauchy_vec);

    //Define the stress tensor as a matrix
    Matrix_3x3 cauchy_mat;
    deformation_measures::undo_voigt_3x3_tensor(cauchy,cauchy_mat);

    //Compute the expected value of the internal force
    for (int i=0; i<3; i++){
        r[i] = 0;
        for (int j=0; j<3; j++){
            r[i] -= dNdx[j]*cauchy_mat(j,i);
        }
    }

    //Compute the value from the function
    balance_equations::compute_internal_force(dNdx, cauchy, _r);
    bool tot_result = true;

    for (int i=0; i<3; i++){tot_result *= 1e-6>fabs(r[i] - _r[i]);}

    balance_equations::compute_internal_force(dNdx, cauchy_vec, _r);
    for (int i=0; i<3; i++){tot_result *= 1e-6>fabs(r[i] - _r[i]);}

    for (int i=0; i<3; i++){balance_equations::compute_internal_force(i, dNdx, cauchy, _r[i]);}
    for (int i=0; i<3; i++){tot_result *= 1e-6>fabs(r[i] - _r[i]);}

    for (int i=0; i<3; i++){balance_equations::compute_internal_force(i, dNdx, cauchy_vec, _r[i]);}
    for (int i=0; i<3; i++){tot_result *= 1e-6>fabs(r[i] - _r[i]);}

    if (tot_result){
        results << "test_compute_internal_force & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_force & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_internal_couple(std::ofstream &results){
    /*!
    ======================================
    |    test_compute_internal_couple    |
    ======================================

    Test the computation of the internal force. 

    */

    Vector_9  r;    //The expected result
    double   _r[9]; //The function output

    //Define the test function
    double N;
    define_N(N);

    double dNdx[3];
    define_dNdx(dNdx);

    //Define the stress tensors as an eigen matrix
    Vector_9 cauchy;
    define_cauchy(cauchy);

    Vector_9 s;
    define_s(s);

    Vector_27 m;
    define_m(m);

    //Define the stress tensor as a std::vector
    std::vector<double> cauchy_vec;
    balance_equations::map_eigen_to_vector(cauchy,cauchy_vec);

    std::vector<double> s_vec;
    balance_equations::map_eigen_to_vector(s,s_vec);

    std::vector<double> m_vec;
    balance_equations::map_eigen_to_vector(m,m_vec);

    //Define the stress tensors as a matrix
    Matrix_3x3 cauchy_mat;
    deformation_measures::undo_voigt_3x3_tensor(cauchy,cauchy_mat);

    Matrix_3x3 s_mat;
    deformation_measures::undo_voigt_3x3_tensor(s,s_mat);

    Matrix_3x9 m_mat;
    deformation_measures::undo_voigt_3x9_tensor(m,m_mat);

    //Compute the expected value of the internal couple
    deformation_measures::voigt_3x3_tensor( N*(cauchy_mat - s_mat), r);

    Eigen::Vector3d dNdx_vec;
    for (int i=0; i<3; i++){dNdx_vec[i] = dNdx[i];}

    Vector_9 divm;
    divm = dNdx_vec.transpose()*m_mat;
    divm.row(3).swap(divm.row(6));
    divm.row(4).swap(divm.row(7));
    divm.row(5).swap(divm.row(8));

    r -= divm;

    //Compute the value of the function
    balance_equations::compute_internal_couple(N, dNdx, cauchy, s, m, _r);

    bool tot_result = true;
    for (int i=0; i<9; i++){tot_result *= 1e-6>fabs(r[i] - _r[i]);}

    balance_equations::compute_internal_couple(N, dNdx, cauchy_vec, s_vec, m_vec, _r);
    for (int i=0; i<9; i++){tot_result *= 1e-6>fabs(r[i] - _r[i]);}

    Matrix_3x3 _r_mat;
    double _r_ij;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            balance_equations::compute_internal_couple(i, j, N, dNdx, cauchy, s, m, _r_ij);
            _r_mat(i,j) = _r_ij;
        }
    }
    Vector_9 __r;
    deformation_measures::voigt_3x3_tensor(_r_mat, __r);
    tot_result *= r.isApprox(__r);

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            balance_equations::compute_internal_couple(i, j, N, dNdx, cauchy_vec, s_vec, m_vec, _r_ij);
            _r_mat(i,j) = _r_ij;
        }
    }
    deformation_measures::voigt_3x3_tensor(_r_mat, __r);
    tot_result *= r.isApprox(__r);

    if (tot_result){
        results << "test_compute_internal_couple & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_couple & False\\\\\n\\hline\n";
    }

    return 1;
}


int test_compute_internal_force_jacobian(std::ofstream &results){
    /*!==============================================
    |    test_compute_internal_force_jacobian    |
    ==============================================
    
    Test the computation of the jacobian associated with the internal force.
    
    */
    
    Matrix_3x12  r; //The expected result
    Matrix_3x12 _r; //The function output
    
    //Timers
    auto t0 = Clock::now();
    auto t1 = Clock::now();
    
    //The DOF vector
    std::vector<double> U = {1.22,2.1,4.1,-2.3,.124,7.2,-8.2,.28,7.21,2.1,-9.2,3.1};
    
    //The test function values
    double N;
    define_N(N);

    double dNdx[3];
    define_dNdx(dNdx);
    
    //The interpolation function values
    double eta;
    define_eta(eta);

    double detadx[3];
    define_detadx(detadx);

    //Compute the numeric values
    std::cout << "\nJacobian of Balance of Linear Momentum\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_balance_of_linear_momentum_U, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the analytic values
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(grad_u, F);
    
    deformation_measures::assemble_chi(phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data, F, grad_chi);
    
    Matrix_3x9 grad_phi_mat;
    Vector_27  grad_phi;
    deformation_measures::assemble_grad_chi(grad_phi_data, Matrix_3x3::Identity(), grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
    
    balance_equations::compute_internal_force_jacobian(N, dNdx, eta, detadx, DcauchyDgrad_u, dcauchydchi, DcauchyDgrad_phi, _r);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";

    //std::cout << " r:\n" <<  r << "\n";
    //std::cout << "_r:\n" << _r << "\n";

    bool tot_result = r.isApprox(_r,1e-6);

    //Compute the jacobian using a vector form    
    std::vector<std::vector<double>> _DcauchyDgrad_u;
    std::vector<std::vector<double>> _DcauchyDphi;
    std::vector<std::vector<double>> _DcauchyDgrad_phi;
    std::vector<std::vector<double>> __r;
    
    balance_equations::map_eigen_to_vector(DcauchyDgrad_u,   _DcauchyDgrad_u);
    balance_equations::map_eigen_to_vector(dcauchydchi,      _DcauchyDphi);
    balance_equations::map_eigen_to_vector(DcauchyDgrad_phi, _DcauchyDgrad_phi);
    balance_equations::compute_internal_force_jacobian(N, dNdx, eta, detadx, _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi, __r);
    balance_equations::map_vector_to_eigen(__r, _r);
    tot_result *= r.isApprox(_r,1e-6);

    if (tot_result){
        results << "test_compute_internal_force_jacobian & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_force_jacobian & False\\\\\n\\hline\n";
    }
}

int test_compute_internal_couple_jacobian(std::ofstream &results){
    /*!===============================================
    |    test_compute_internal_couple_jacobian    |
    ===============================================
    
    Test the computation of the jacobian associated with the internal couple.
    
    */
    
    Matrix_9x12  r; //The expected result
    Matrix_9x12 _r; //The function output
    
    //Timers
    auto t0 = Clock::now();
    auto t1 = Clock::now();
    
    //The DOF vector
    std::vector<double> U = {1.22,2.1,4.1,-2.3,.124,7.2,-8.2,.28,7.21,2.1,-9.2,3.1};
    
    //The shape function values
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);

    double eta;
    define_eta(eta);

    double detadx[3];
    define_detadx(detadx);
    
    //Compute the numeric values
    std::cout << "\nJacobian of Balance of First Moment of Momentum\n";
    std::cout << "Finite Difference vs. Analytic Jacobian\n";
    t0 = Clock::now();
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_balance_of_first_moment_of_momentum_U, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    t1 = Clock::now();
    std::cout << "Finite Difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    for (int i=0; i<r_vec.size(); i++){
        for (int j=0; j<r_vec[i].size(); j++){
            r(j,i) = r_vec[i][j];
        }
    }
    
    //Compute the analytic values
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    t0 = Clock::now();
    deformation_measures::get_deformation_gradient(grad_u, F);
    
    deformation_measures::assemble_chi(phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data, F, grad_chi);
    
    Matrix_3x9 grad_phi_mat;
    Vector_27  grad_phi;
    deformation_measures::assemble_grad_chi(grad_phi_data, Matrix_3x3::Identity(), grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
    
    balance_equations::compute_internal_couple_jacobian(N, dNdx, eta, detadx,
                                                        DcauchyDgrad_u, dcauchydchi, DcauchyDgrad_phi,
                                                        DsDgrad_u,      dsdchi,      DsDgrad_phi,
                                                        DmDgrad_u,      dmdchi,      DmDgrad_phi,
                                                        _r);
    t1 = Clock::now();
    std::cout << "Analytic Jacobian: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "\n";
    
    //std::cout << " r:\n" <<  r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    //Compute the jacobian using a vector form
    std::vector<double> _cauchy;
    std::vector<double> _s;
    std::vector<double> _m;
 
    std::vector<std::vector<double>> _DcauchyDgrad_u;
    std::vector<std::vector<double>> _DcauchyDphi;
    std::vector<std::vector<double>> _DcauchyDgrad_phi;

    std::vector<std::vector<double>> _DsDgrad_u;
    std::vector<std::vector<double>> _DsDphi;
    std::vector<std::vector<double>> _DsDgrad_phi;

    std::vector<std::vector<double>> _DmDgrad_u;
    std::vector<std::vector<double>> _DmDphi;
    std::vector<std::vector<double>> _DmDgrad_phi;

    std::vector<std::vector<double>> __r;
    
    balance_equations::map_eigen_to_vector(cauchy, _cauchy);
    balance_equations::map_eigen_to_vector(s, _s);
    balance_equations::map_eigen_to_vector(m, _m);

    balance_equations::map_eigen_to_vector(DcauchyDgrad_u,   _DcauchyDgrad_u);
    balance_equations::map_eigen_to_vector(dcauchydchi,      _DcauchyDphi);
    balance_equations::map_eigen_to_vector(DcauchyDgrad_phi, _DcauchyDgrad_phi);

    balance_equations::map_eigen_to_vector(DsDgrad_u,   _DsDgrad_u);
    balance_equations::map_eigen_to_vector(dsdchi,      _DsDphi);
    balance_equations::map_eigen_to_vector(DsDgrad_phi, _DsDgrad_phi);

    balance_equations::map_eigen_to_vector(DmDgrad_u,   _DmDgrad_u);
    balance_equations::map_eigen_to_vector(dmdchi,      _DmDphi);
    balance_equations::map_eigen_to_vector(DmDgrad_phi, _DmDgrad_phi);

    balance_equations::compute_internal_couple_jacobian(N, dNdx, eta, detadx,
                                                        _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                                        _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                                                        _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                                                        __r);
    balance_equations::map_vector_to_eigen(__r, _r);
    tot_result *= r.isApprox(_r,1e-6);

    int sot_to_voigt_map[3][3];
    deformation_measures::get_sot_to_voigt_map(sot_to_voigt_map);

    double __r_ijA;
    int Ihat;

    for (int i=0; i<3; i++){

        for (int j=0; j<3; j++){
            Ihat = sot_to_voigt_map[i][j];

            for (int dof_num=0; dof_num<12; dof_num++){

                balance_equations::compute_internal_couple_jacobian(i,j,dof_num, N, dNdx, eta, detadx,
                                                                    _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                                                                    _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                                                                    _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                                                                    __r_ijA);
                __r[Ihat][dof_num] = __r_ijA;
            }
        }
    }
    balance_equations::map_vector_to_eigen(__r, _r);
    tot_result *= r.isApprox(_r,1e-6);

    if (tot_result){
        results << "test_compute_internal_couple_jacobian & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_couple_jacobian & False\\\\\n\\hline\n";
    }
}

int test_micromorphic_material_library(std::ofstream &results){
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
    
    double dNdx[3];
    define_dNdx(dNdx);

    double eta;
    define_eta(eta);

    double detadx[3];
    define_detadx(detadx);
    
    //Compute the values explicitly
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
    
    Matrix_9x9   dPK2dF;
    Matrix_9x9   dPK2dchi;
    Matrix_9x27  dPK2dgrad_chi;
    Matrix_9x9   dSIGMAdF;
    Matrix_9x9   dSIGMAdchi;
    Matrix_9x27  dSIGMAdgrad_chi;
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;
    
    Matrix_9x9   dcauchydF;
    Matrix_9x9   dcauchydchi;
    Matrix_9x27  dcauchydgrad_chi;
    Matrix_9x9   dsdF;
    Matrix_9x9   dsdchi;
    Matrix_9x27  dsdgrad_chi;
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;
    
    Matrix_9x9   DcauchyDgrad_u;
    Matrix_9x27  DcauchyDgrad_phi;
    Matrix_9x9   DsDgrad_u;
    Matrix_9x27  DsDgrad_phi;
    Matrix_27x9  DmDgrad_u;
    Matrix_27x27 DmDgrad_phi;
    
    define_parameters(params);
    
    deformation_measures::get_deformation_gradient(grad_u, F);
    
    deformation_measures::assemble_chi(phi, chi);
    deformation_measures::assemble_grad_chi(grad_phi_data, F, grad_chi);
    
    Matrix_3x9 grad_phi_mat;
    Vector_27  grad_phi;
    deformation_measures::assemble_grad_chi(grad_phi_data, Matrix_3x3::Identity(), grad_phi_mat);
    deformation_measures::voigt_3x9_tensor(grad_phi_mat,grad_phi);
    
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                               dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                               dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                               dMdF,     dMdchi,     dMdgrad_chi);
    
    deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
                                               
    deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s, m,
                                                                 dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                 dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, dcauchydchi, dcauchydgrad_chi,
                                                                 dsdF,   dsdchi,   dsdgrad_chi,   dmdF,      dmdchi,      dmdgrad_chi);
    
    deformation_measures::compute_total_derivatives(F,               grad_phi,
                                                    dcauchydF,       dcauchydgrad_chi,
                                                    dsdF,            dsdgrad_chi,
                                                    dmdF,            dmdgrad_chi,
                                                    DcauchyDgrad_u,  DcauchyDgrad_phi,
                                                    DsDgrad_u,       DsDgrad_phi,
                                                    DmDgrad_u,       DmDgrad_phi);
                                                    
    //Compute the values using the version stored in the material library.
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();

    auto material = factory.GetMaterial("LinearElasticity");
    
    //Output
    
    Vector_9  _cauchy;
    Vector_9  _s;
    Vector_27 _m;
    
    Matrix_9x9   _DcauchyDgrad_u;
    Matrix_9x9   _DcauchyDphi;
    Matrix_9x27  _DcauchyDgrad_phi;
    Matrix_9x9   _DsDgrad_u;
    Matrix_9x9   _DsDphi;
    Matrix_9x27  _DsDgrad_phi;
    Matrix_27x9  _DmDgrad_u;
    Matrix_27x9  _DmDphi;
    Matrix_27x27 _DmDgrad_phi;

    //Vector output
    std::vector<double> _cauchy_v;
    std::vector<double> _s_v;
    std::vector<double> _m_v;

    std::vector<std::vector<double>> _DcauchyDgrad_u_v;
    std::vector<std::vector<double>> _DcauchyDphi_v;
    std::vector<std::vector<double>> _DcauchyDgrad_phi_v;
    std::vector<std::vector<double>> _DsDgrad_u_v;
    std::vector<std::vector<double>> _DsDphi_v;
    std::vector<std::vector<double>> _DsDgrad_phi_v;
    std::vector<std::vector<double>> _DmDgrad_u_v;
    std::vector<std::vector<double>> _DmDphi_v;
    std::vector<std::vector<double>> _DmDgrad_phi_v;

    
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
                             _cauchy, _s, _m, ADD_TERMS);
                             
    bool tot_result = true;
    
    tot_result *= cauchy.isApprox(_cauchy);
    tot_result *= s.isApprox(_s);
    tot_result *= m.isApprox(_m);
    
    //Compare the result from the vector form of the stress calculation alone
    material->evaluate_model(time, fparams, grad_u, phi, grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _cauchy_v, _s_v, _m_v, ADD_TERMS_v);

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(cauchy[i]-_cauchy_v[i]);
    }

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(s[i]-_s_v[i]);
    }

    for (int i=0; i<27; i++){
        tot_result *= 1e-6>fabs(m[i]-_m_v[i]);
    }

    //Compare the results from the jacobian calculation.
    material->evaluate_model(time, fparams, grad_u, phi, grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _cauchy, _s, _m,
                             _DcauchyDgrad_u, _DcauchyDphi,  _DcauchyDgrad_phi,
                             _DsDgrad_u,      _DsDphi,       _DsDgrad_phi,
                             _DmDgrad_u,      _DmDphi,       _DmDgrad_phi,
                             ADD_TERMS,       ADD_JACOBIANS);
    
    tot_result *= cauchy.isApprox(_cauchy);
    tot_result *= s.isApprox(_s);
    tot_result *= m.isApprox(_m);
    
    tot_result *= DcauchyDgrad_u.isApprox(_DcauchyDgrad_u);
    tot_result *= dcauchydchi.isApprox(_DcauchyDphi); //d(x)dchi = D(x)Dphi
    tot_result *= DcauchyDgrad_phi.isApprox(_DcauchyDgrad_phi);
    
    tot_result *= DsDgrad_u.isApprox(_DsDgrad_u);
    tot_result *= dsdchi.isApprox(_DsDphi); //d(x)dchi = D(x)Dphi
    tot_result *= DsDgrad_phi.isApprox(_DsDgrad_phi);
    
    tot_result *= DmDgrad_u.isApprox(_DmDgrad_u);
    tot_result *= dmdchi.isApprox(_DmDphi); //d(x)dchi = D(x)Dphi
    tot_result *= DmDgrad_phi.isApprox(_DmDgrad_phi);

    //Compare the result from the vector form of the jacobian calculation
    material->evaluate_model(time, fparams, grad_u, phi, grad_phi_data, SDVS, ADDDOF, ADD_grad_DOF,
                             _cauchy_v, _s_v, _m_v,
                             _DcauchyDgrad_u_v, _DcauchyDphi_v,  _DcauchyDgrad_phi_v,
                             _DsDgrad_u_v,      _DsDphi_v,       _DsDgrad_phi_v,
                             _DmDgrad_u_v,      _DmDphi_v,       _DmDgrad_phi_v,
                             ADD_TERMS_v,       ADD_JACOBIANS_v);

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(cauchy[i]-_cauchy_v[i]);
    }

    for (int i=0; i<9; i++){
        tot_result *= 1e-6>fabs(s[i]-_s_v[i]);
    }

    for (int i=0; i<27; i++){
        tot_result *= 1e-6>fabs(m[i]-_m_v[i]);
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DcauchyDgrad_u(i,j) - _DcauchyDgrad_u_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(dcauchydchi(i,j) - _DcauchyDphi_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<27; j++){
            tot_result *= 1e-6>fabs(DcauchyDgrad_phi(i,j) - _DcauchyDgrad_phi_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DsDgrad_u(i,j) - _DsDgrad_u_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(dsdchi(i,j) - _DsDphi_v[i][j]);
        }
    }

    for (int i=0; i<9; i++){
        for (int j=0; j<27; j++){
            tot_result *= 1e-6>fabs(DsDgrad_phi(i,j) - _DsDgrad_phi_v[i][j]);
        }
    }

    for (int i=0; i<27; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(DmDgrad_u(i,j) - _DmDgrad_u_v[i][j]);
        }
    }
    
    for (int i=0; i<27; i++){
        for (int j=0; j<9; j++){
            tot_result *= 1e-6>fabs(dmdchi(i,j) - _DmDphi_v[i][j]);
        }
    }
    
    for (int i=0; i<27; i++){
        for (int j=0; j<27; j++){
            tot_result *= 1e-6>fabs(DmDgrad_phi(i,j) - _DmDgrad_phi_v[i][j]);
        }
    }



    if (tot_result){
        results << "test_micromorphic_material_library & True\\\\\n\\hline\n";
    }
    else {
        results << "test_micromorphic_material_library & False\\\\\n\\hline\n";
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

    std::ofstream results;
    //Open the results file
    results.open ("results.tex");

    //!Run the test functions
    
    //!Test of the stiffness tensors
    test_compute_A_voigt(results);
    test_compute_B_voigt(results);
    test_compute_C_voigt(results);
    test_compute_D_voigt(results);
    
    //!Tests of the stress computations
    test_compute_PK2_stress(results);
    test_compute_symmetric_stress(results);
    test_compute_higher_order_stress(results);
    
    //!Tests of the non-zero gradients w.r.t. the derived deformation measures.
    test_compute_dPK2dRCG(results);
    test_compute_dSIGMAdRCG(results);
    test_compute_dPK2dPsi(results);
    test_compute_dSIGMAdPsi(results);
    test_compute_dPK2dGamma(results);
    test_compute_dSIGMAdGamma(results);
    test_compute_dMdGamma(results);
    
    //!Tests of the reference configuration jacobians w.r.t. the deformation gradient.
    test_compute_dPK2dF(results);
    test_compute_dSIGMAdF(results);
    test_compute_dMdF(results);
    
    //!Tests of the reference configuration jacobians w.r.t. the micro-displacement tensor.
    test_compute_dPK2dchi(results);
    test_compute_dSIGMAdchi(results);
    test_compute_dMdchi(results);
    
    //!Tests of the reference configuration jacobians w.r.t. the gradient of the micro-displacement tensor.
    test_compute_dPK2dgrad_chi(results);
    test_compute_dSIGMAdgrad_chi(results);
    test_compute_dMdgrad_chi(results);

    //!Tests of the current configuration jacobians w.r.t. the deformation gradient
    test_compute_dcauchydF(results);
    test_compute_dsdF(results);
    test_compute_dmdF(results);

    //!Tests of the current configuration jacobians w.r.t. the micro-displacement tensor
    test_compute_dcauchydchi(results);
    test_compute_dsdchi(results);
    test_compute_dmdchi(results);
    
    //!Tests of the current configuation jacobians w.r.t. the gradient of the micro-displacement tensor
    test_compute_dcauchydgrad_chi(results);
    test_compute_dsdgrad_chi(results);
    test_compute_dmdgrad_chi(results);
    
    //!Tests of the computation of the total derivatives
    test_compute_total_derivatives_DcauchyDgrad_u(results);
    test_compute_total_derivatives_DsDgrad_u(results);
    test_compute_total_derivatives_DmDgrad_u(results);
    
    test_compute_total_derivatives_DcauchyDgrad_phi(results);
    test_compute_total_derivatives_DsDgrad_phi(results);
    test_compute_total_derivatives_DmDgrad_phi(results);

    //!Test the computation of the internal force
    test_compute_internal_force(results);
    test_compute_internal_couple(results);
    
    //!Test of the computation of the balance equation derivatives
    test_compute_internal_force_jacobian(results);
    test_compute_internal_couple_jacobian(results);
    
    //!Test of the instance of the model in the material library
    test_micromorphic_material_library(results);

    //Close the results file
    results.close();
}

