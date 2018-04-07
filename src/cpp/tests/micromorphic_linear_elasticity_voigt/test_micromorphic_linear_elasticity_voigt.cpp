/*!============================================================================
   |                                                                          |
   |                      test_deformation_measures.cpp                       |
   |                                                                          |
   ----------------------------------------------------------------------------
   | The unit test file for deformation_measures.h/cpp. This file tests the   |
   | classes and functions defined in deformation_measures.h/cpp.             |
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
#include<finite_difference.h>
#include<deformation_measures.h>
#include<micromorphic_linear_elasticity_voigt.h>

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

std::vector<double> parse_inv_sot(std::vector<double> Avec){
    /*!=======================
    |    parse_inv_sot    |
    =======================
    
    Parse the inverse of a second order tensor to compute it's derivative 
    using finite differences.
    
    */
    
    Vector_9   V;
    Matrix_3x3 A;
    
    std::vector<double> Ainvvec;
    Ainvvec.resize(9);
    
    for (int i=0; i<9; i++){V(i) = Avec[i];}
    
    deformation_measures::undo_voigt_3x3_tensor(V,A);
    deformation_measures::voigt_3x3_tensor(A.inverse(),V);
    
    for (int i=0; i<9; i++){Ainvvec[i] = V(i);}
    
    return Ainvvec;
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

std::vector<double> parse_RCG_F(std::vector<double> Fvec){
    /*!=====================
    |    parse_RCG_F    |
    =====================
    
    Parse RCG with respect to the 
    deformation gradient.
    
    */
    
    Vector_9   F_voigt;
    Matrix_3x3 F;
    
    for (int i=0; i<9; i++){F_voigt(i) = Fvec[i];}
    
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    
    Matrix_3x3 RCG;
    Vector_9   RCG_voigt;
    std::vector<double> RCG_vec;
    RCG_vec.resize(9);
    deformation_measures::get_right_cauchy_green(F,RCG);
    deformation_measures::voigt_3x3_tensor(RCG,RCG_voigt);
    
    for (int i=0; i<9; i++){RCG_vec[i] = RCG_voigt(i);}
    return RCG_vec;
}

std::vector<double> parse_Psi_F(std::vector<double> Fvec){
    /*!=====================
    |    parse_Psi_F    |
    =====================
    
    Parse Psi with respect to the 
    deformation gradient.
    
    */
    
    Vector_9   F_voigt;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    define_chi(chi);
    
    for (int i=0; i<9; i++){F_voigt(i) = Fvec[i];}
    
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    
    Matrix_3x3 Psi;
    Vector_9   Psi_voigt;
    std::vector<double> Psi_vec;
    Psi_vec.resize(9);
    deformation_measures::get_psi(F,chi,Psi);
    deformation_measures::voigt_3x3_tensor(Psi,Psi_voigt);
    
    for (int i=0; i<9; i++){Psi_vec[i] = Psi_voigt(i);}
    return Psi_vec;
}

std::vector<double> parse_Psi_chi(std::vector<double> chivec){
    /*!=======================
    |    parse_Psi_chi    |
    =======================
    
    Parse Psi with respect to the 
    micro-displacement tensor.
    
    */
    
    Vector_9   chi_voigt;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    define_deformation_gradient(F);
    
    for (int i=0; i<9; i++){chi_voigt(i) = chivec[i];}
    
    deformation_measures::undo_voigt_3x3_tensor(chi_voigt,chi);
    
    Matrix_3x3 Psi;
    Vector_9   Psi_voigt;
    std::vector<double> Psi_vec;
    Psi_vec.resize(9);
    deformation_measures::get_psi(F,chi,Psi);
    deformation_measures::voigt_3x3_tensor(Psi,Psi_voigt);
    
    for (int i=0; i<9; i++){Psi_vec[i] = Psi_voigt(i);}
    return Psi_vec;
}

std::vector<double> parse_Gamma_F(std::vector<double> Fvec){
    /*!=======================
    |    parse_Gamma_F    |
    =======================
    
    Parse Gamma with respect to the 
    deformation gradient.
    
    */
    
    Vector_9   F_voigt;
    Matrix_3x3 F;
    Matrix_3x9 grad_chi;
    define_grad_phi(grad_chi);
    
    for (int i=0; i<9; i++){F_voigt(i) = Fvec[i];}
    
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    
    Matrix_3x9 Gamma;
    Vector_27  Gamma_voigt;
    std::vector<double> Gamma_vec;
    Gamma_vec.resize(27);
    deformation_measures::get_gamma(F,grad_chi,Gamma);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    for (int i=0; i<27; i++){Gamma_vec[i] = Gamma_voigt(i);}
    return Gamma_vec;
}

std::vector<double> parse_grad_chi_grad_phi(std::vector<double> grad_phivec){
    /*!=======================
    |    parse_Gamma_F    |
    =======================
    
    Parse Gamma with respect to the 
    deformation gradient.
    
    */
    
    Vector_27  grad_phi_voigt;
    Matrix_3x9 grad_phi;
    Matrix_3x3 F;
    define_deformation_gradient(F);
    
    for (int i=0; i<27; i++){grad_phi_voigt(i) = grad_phivec[i];}
    
    Matrix_3x9 grad_chi;
    Vector_27  grad_chi_voigt;
    
    deformation_measures::perform_right_positive_cyclic_permutation(grad_phi_voigt);
    deformation_measures::undo_voigt_3x9_tensor(grad_phi_voigt,grad_phi);
    grad_chi = F.transpose()*grad_phi;
    deformation_measures::voigt_3x9_tensor(grad_chi,grad_chi_voigt);
    deformation_measures::perform_left_positive_cyclic_permutation(grad_chi_voigt);
    
    std::vector<double> grad_chi_vec;
    grad_chi_vec.resize(27);
    
    for (int i=0; i<27; i++){grad_chi_vec[i] = grad_chi_voigt(i);}
    return grad_chi_vec;
}

std::vector<double> parse_grad_chi_F(std::vector<double> Fvec){
    /*!=======================
    |    parse_Gamma_F    |
    =======================
    
    Parse Gamma with respect to the 
    deformation gradient.
    
    */
    
    Matrix_3x9 grad_phi;
    Matrix_3x3 F;
    Vector_9   F_voigt;
    
    double grad_phi_tmp[9][3] = {{0.268677941492,  0.0956706826016, 0.325269776806},
                                 {0.638467644829,  0.997489522961,  0.173494630408},
                                 {0.81125862445,   0.224706569702,  0.514463248854},
                                 {0.748926548687,  0.1079833322,    0.797424154102},
                                 {0.0782798344482, 0.324585039367,  0.032983492952},
                                 {0.84359324251,   0.602202553902,  0.165660665757},
                                 {0.117005479066,  0.476022110021,  0.00730531427893},
                                 {0.792163708315,  0.968709140644,  0.0352176264845},
                                 {0.585823963851,  0.611515761521,  0.0102566326533}};
    
    //Put grad_phi into matrix form
    deformation_measures::assemble_grad_chi(grad_phi_tmp,Matrix_3x3::Identity(),grad_phi);
    
    for (int i=0; i<9; i++){F_voigt(i) = Fvec[i];}
    deformation_measures::undo_voigt_3x3_tensor(F_voigt,F);
    
    Matrix_3x9 grad_chi;
    Vector_27  grad_chi_voigt;
    Vector_27  grad_phi_voigt;
    
    deformation_measures::voigt_3x9_tensor(grad_phi,grad_phi_voigt);
    deformation_measures::perform_right_positive_cyclic_permutation(grad_phi_voigt);
    deformation_measures::undo_voigt_3x9_tensor(grad_phi_voigt,grad_phi);
    grad_chi = F.transpose()*grad_phi;
    deformation_measures::voigt_3x9_tensor(grad_chi,grad_chi_voigt);
    deformation_measures::perform_left_positive_cyclic_permutation(grad_chi_voigt);
    
    std::vector<double> grad_chi_vec;
    grad_chi_vec.resize(27);
    
    for (int i=0; i<27; i++){grad_chi_vec[i] = grad_chi_voigt(i);}
    return grad_chi_vec;
}

std::vector<double> parse_Gamma_grad_chi(std::vector<double> grad_chivec){
    /*!==============================
    |    parse_Gamma_grad_chi    |
    ==============================
    
    Parse Psi with respect to the 
    deformation gradient.
    
    */
    
    Vector_27  grad_chi_voigt;
    Matrix_3x3 F;
    define_deformation_gradient(F);
    Matrix_3x9 grad_chi;
    define_grad_phi(grad_chi);
    
    for (int i=0; i<27; i++){grad_chi_voigt(i) = grad_chivec[i];}
    
    deformation_measures::undo_voigt_3x9_tensor(grad_chi_voigt,grad_chi);
    
    Matrix_3x9 Gamma;
    Vector_27  Gamma_voigt;
    std::vector<double> Gamma_vec;
    Gamma_vec.resize(27);
    deformation_measures::get_gamma(F,grad_chi,Gamma);
    deformation_measures::voigt_3x9_tensor(Gamma,Gamma_voigt);
    
    for (int i=0; i<27; i++){Gamma_vec[i] = Gamma_voigt(i);}
    return Gamma_vec;
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

    define_A(A);     //Set the expected result
    micro_material::compute_A_voigt(0.191519450379, 0.62210877104, _A);

    bool tot_result = A.isApprox(_A);

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

    define_B(B);     //Set the expected result
    micro_material::compute_B_voigt(0.4377277390071145, 0.7799758081188035,
                                    0.2725926052826416, 0.2764642551430967,
                                    0.7853585837137692, _B);

    bool tot_result = B.isApprox(_B);

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

    define_C(C);     //Set the expected result
    micro_material::compute_C_voigt(0.8018721775350193, 0.9581393536837052,
                                    0.8759326347420947, 0.35781726995786667,
                                    0.5009951255234587, 0.6834629351721363,
                                    0.7127020269829002, 0.37025075479039493,
                                    0.5611961860656249, 0.5030831653078097,
                                    0.013768449590682241, _C);

    bool tot_result = C.isApprox(_C);

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

    define_D(D);     //Set the expected result
    micro_material::compute_D_voigt(0.2764642551430967, 0.7853585837137692, _D);

    bool tot_result = D.isApprox(_D);

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

int test_map_stresses_to_current_configuration(std::ofstream &results){

    Vector_9   cauchy; //The expected cauchy stress
    Vector_9  _cauchy; //The result of the function
    Vector_9   s; //The expected symmetric stress
    Vector_9  _s; //The result of the function
    Vector_27  m; //The expected higher order stress
    Vector_27 _m; //The result of the function

    Vector_9 PK2;   //The second piola kirchhoff stress
    Vector_9 SIGMA; //The symmetric stress
    Vector_27 M;    //The higher order stress

    Matrix_3x3 F;
    Matrix_3x3 chi;
    define_deformation_gradient(F);
    define_chi(chi);

    //Extract the stresses in the reference configuration
    define_PK2(PK2);
    define_SIGMA(SIGMA);
    define_M(M);

    //Extract the stresses in the current configuration
    define_cauchy(cauchy);
    define_s(s);
    define_m(m);

    micro_material::map_stresses_to_current_configuration(F,chi,PK2,SIGMA,M,_cauchy,_s,_m);

    bool tot_result = cauchy.isApprox(_cauchy,1e-6);
    tot_result *= s.isApprox(_s,1e-6);
    tot_result *= m.isApprox(_m,1e-6);

    if (tot_result){
        results << "test_map_stresses_to_current_configuration & True\\\\\n\\hline\n";
    }
    else {
        results << "test_map_stresses_to_current_configuration & False\\\\\n\\hline\n";
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
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_RCG, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dRCG_vec = fd.numeric_gradient();
    
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
    micro_material::compute_dPK2dRCG(RCG0, RCG0.inverse(), Gamma, Gamma_voigt,
                                     E, E_micro, E_voigt, E_micro_voigt,
                                     A, B, C, D, _dPK2dRCG);
    
    bool tot_result = dPK2dRCG.isApprox(_dPK2dRCG,1e-6);
    
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
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_RCG, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdRCG_vec = fd.numeric_gradient();
    
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
    micro_material::compute_dPK2dRCG(RCG0, RCG0.inverse(), Gamma, Gamma_voigt,
                                     E, E_micro, E_voigt, E_micro_voigt,
                                     A, B, C, D, terms, dPK2dRCG);
    micro_material::compute_dSIGMAdRCG(terms, _dSIGMAdRCG);
    
    bool tot_result = dSIGMAdRCG.isApprox(_dSIGMAdRCG,1e-6);
    
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
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_Psi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dPsi_vec = fd.numeric_gradient();
    
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
    micro_material::compute_dPK2dPsi(RCG.inverse(), E_micro, E_voigt, E_micro_voigt,
                                     B, D, _dPK2dPsi);
    
    bool tot_result = dPK2dPsi.isApprox(_dPK2dPsi,1e-6);
    
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
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_Psi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdPsi_vec = fd.numeric_gradient();
    
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
    micro_material::compute_dPK2dPsi(RCG.inverse(), E_micro, E_voigt, E_micro_voigt,
                                     B, D, terms, dPK2dPsi);
    micro_material::compute_dSIGMAdPsi(terms, _dSIGMAdPsi);
    
    bool tot_result = dSIGMAdPsi.isApprox(_dSIGMAdPsi,1e-6);
    
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
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_pk2_stress_Gamma, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPK2dGamma_vec = fd.numeric_gradient();
    
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
    micro_material::compute_dPK2dGamma(RCG.inverse(), Gamma0, Gamma_voigt, C, _dPK2dGamma);
    
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
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_symmetric_stress_Gamma, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dSIGMAdGamma_vec = fd.numeric_gradient();
    
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
    Matrix_9x27 dPK2dGamma;
    Matrix_9x27 terms[2];
    micro_material::compute_dPK2dGamma(RCG.inverse(), Gamma0, Gamma_voigt, C, terms, dPK2dGamma);
    micro_material::compute_dSIGMAdGamma(terms, _dSIGMAdGamma);
    
    bool tot_result = dSIGMAdGamma.isApprox(_dSIGMAdGamma,1e-6);
    
    if (tot_result){
        results << "test_compute_dSIGMAdGamma & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dSIGMAdGamma & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dAinvdA(std::ofstream &results){
    /*!==============================
    |    test_compute_dAinvdA    |
    ==============================
    
    Test the computation of the derivative of the inverse 
    of a second order tensor w.r.t. A.
    
    */
    
    Matrix_9x9  r; //The expected result
    Matrix_9x9 _r; //The function output
    
    Matrix_3x3 A; //The matrix to invert.
    A << 0.69646919,  0.41872705,  0.60380783,  0.41872705,  0.71946897,
         0.5539681 ,  0.60380783,  0.5539681 ,  0.4809319;
         
    Vector_9 A_voigt;
    deformation_measures::voigt_3x3_tensor(A,A_voigt);
    
    std::vector<double> Avec;
    Avec.resize(9);
    for (int i=0; i<9; i++){Avec[i] = A_voigt(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_inv_sot, 2, Avec , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    for (int i=0; i<r_vec.size(); i++){
        
        for (int j=0; j<r_vec[i].size(); j++){
            
            r(i,j) = r_vec[j][i];
            
        }
    }
    
    micro_material::compute_dAinvdA(A.inverse(),_r);
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_compute_dAinvdA & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dAinvdA & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dRCGdF(std::ofstream &results){
    /*!=============================
    |    test_compute_dRCGdF    |
    =============================
    
    Test the computation of the jacobian of the 
    right cauchy green deformation measure w.r.t. 
    the deformation gradient.
    */
    
    Matrix_9x9 dRCGdF; //The expected result
    std::vector< std::vector<double> > dRCGdF_vec; //Output from the finite_difference
    SpMat _dRCGdF(9,9); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F0;
    define_deformation_gradient(F0);
    
    Vector_9 F0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F0,F0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F0_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_RCG_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dRCGdF_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dRCGdF_vec.size(); i++){
        for (int j=0; j<dRCGdF_vec[i].size(); j++){
            dRCGdF(j,i) = dRCGdF_vec[i][j];
        }
    }
    
    micro_material::compute_dRCGdF(F0,_dRCGdF);
    
    Matrix_9x9 _dRCGdF_dense = _dRCGdF;
    bool tot_result = dRCGdF.isApprox(_dRCGdF_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dRCGdF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dRCGdF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPsidF(std::ofstream &results){
    /*!=============================
    |    test_compute_dPsidF    |
    =============================
    
    Test the computation of the jacobian of the 
    micro deformation measure Psi w.r.t. 
    the deformation gradient.
    */
    
    Matrix_9x9 dPsidF; //The expected result
    std::vector< std::vector<double> > dPsidF_vec; //Output from the finite_difference
    SpMat _dPsidF(9,9); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F0;
    Matrix_3x3 chi;
    define_deformation_gradient(F0);
    define_chi(chi);
    
    Vector_9 F0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F0,F0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F0_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_Psi_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPsidF_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPsidF_vec.size(); i++){
        for (int j=0; j<dPsidF_vec[i].size(); j++){
            dPsidF(j,i) = dPsidF_vec[i][j];
        }
    }
    
    micro_material::compute_dPsidF(chi,_dPsidF);
    
    Matrix_9x9 _dPsidF_dense = _dPsidF;
    bool tot_result = dPsidF.isApprox(_dPsidF_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dPsidF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPsidF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dPsidchi(std::ofstream &results){
    /*!===============================
    |    test_compute_dPsidchi    |
    ===============================
    
    Test the computation of the jacobian of the 
    micro deformation measure Psi w.r.t. 
    the micro-displacement.
    */
    
    Matrix_9x9 dPsidchi; //The expected result
    std::vector< std::vector<double> > dPsidchi_vec; //Output from the finite_difference
    SpMat _dPsidchi(9,9); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F;
    Matrix_3x3 chi0;
    define_deformation_gradient(F);
    define_chi(chi0);
    
    Vector_9 chi0_vec; //The vector form of chi
    deformation_measures::voigt_3x3_tensor(chi0,chi0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = chi0_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_Psi_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dPsidchi_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dPsidchi_vec.size(); i++){
        for (int j=0; j<dPsidchi_vec[i].size(); j++){
            dPsidchi(j,i) = dPsidchi_vec[i][j];
        }
    }
    
    micro_material::compute_dPsidchi(F,_dPsidchi);
    
    Matrix_9x9 _dPsidchi_dense = _dPsidchi;
    bool tot_result = dPsidchi.isApprox(_dPsidchi_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dPsidchi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dPsidchi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dGammadF(std::ofstream &results){
    /*!===============================
    |    test_compute_dGammadF    |
    ===============================
    
    Test the computation of the jacobian of the 
    micro deformation measure Gamma w.r.t. 
    the deformation gradient.
    */
    
    Matrix_27x9 dGammadF; //The expected result
    std::vector< std::vector<double> > dGammadF_vec; //Output from the finite_difference
    SpMat _dGammadF(27,9); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F0;
    Matrix_3x9 grad_chi;
    Vector_27  grad_chi_voigt;
    define_deformation_gradient(F0);
    define_grad_phi(grad_chi); //Note: grad_chi == grad_phi
    deformation_measures::voigt_3x9_tensor(grad_chi,grad_chi_voigt);
    
    Vector_9 F0_vec; //The vector form of F
    deformation_measures::voigt_3x3_tensor(F0,F0_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F0_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_Gamma_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dGammadF_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dGammadF_vec.size(); i++){
        for (int j=0; j<dGammadF_vec[i].size(); j++){
            dGammadF(j,i) = dGammadF_vec[i][j];
        }
    }
    
    micro_material::compute_dGammadF(grad_chi_voigt,_dGammadF);
    
    Matrix_27x9 _dGammadF_dense = _dGammadF;
    bool tot_result = dGammadF.isApprox(_dGammadF_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dGammadF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dGammadF & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dGammadgrad_chi(std::ofstream &results){
    /*!======================================
    |    test_compute_dGammadgrad_chi    |
    ======================================
    
    Test the computation of the jacobian of the 
    micro deformation measure Gamma w.r.t. 
    the micro-deformation measure Gamma.
    */
    
    Matrix_27x27 dGammadgrad_chi; //The expected result
    std::vector< std::vector<double> > dGammadgrad_chi_vec; //Output from the finite_difference
    SpMat _dGammadgrad_chi(27,27); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F;
    Matrix_3x9 grad_chi0;
    Vector_27  grad_chi0_voigt;
    define_deformation_gradient(F);
    define_grad_phi(grad_chi0); //Note: grad_chi == grad_phi
    deformation_measures::voigt_3x9_tensor(grad_chi0,grad_chi0_voigt);
    
    Vector_27 grad_chi0_vec; //The vector form of grad_chi
    deformation_measures::voigt_3x9_tensor(grad_chi0,grad_chi0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = grad_chi0_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_Gamma_grad_chi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dGammadgrad_chi_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dGammadgrad_chi_vec.size(); i++){
        for (int j=0; j<dGammadgrad_chi_vec[i].size(); j++){
            dGammadgrad_chi(j,i) = dGammadgrad_chi_vec[i][j];
        }
    }
    
    micro_material::compute_dGammadgrad_chi(F,_dGammadgrad_chi);
    
    Matrix_27x27 _dGammadgrad_chi_dense = _dGammadgrad_chi;
    bool tot_result = dGammadgrad_chi.isApprox(_dGammadgrad_chi_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dGammadgrad_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dGammadgrad_chi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dgrad_chidgrad_phi(std::ofstream &results){
    /*!=========================================
    |    test_compute_dgrad_chidgrad_phi    |
    =========================================
    
    Test the computation of the gradient of chi w.r.t. 
    the gradient of phi.
    
    Note: We assume that the phi is the gradient in the 
          current configuration.
    */
    
    Matrix_27x27 dgrad_chidgrad_phi; //The expected result
    std::vector< std::vector<double> > dgrad_chidgrad_phi_vec; //Output from the finite_difference
    SpMat _dgrad_chidgrad_phi(27,27); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F;
    define_deformation_gradient(F);
    Matrix_3x9 grad_phi0;
    double grad_phi0_tmp[9][3] = {{0.268677941492,  0.0956706826016, 0.325269776806},
                                  {0.638467644829,  0.997489522961,  0.173494630408},
                                  {0.81125862445,   0.224706569702,  0.514463248854},
                                  {0.748926548687,  0.1079833322,    0.797424154102},
                                  {0.0782798344482, 0.324585039367,  0.032983492952},
                                  {0.84359324251,   0.602202553902,  0.165660665757},
                                  {0.117005479066,  0.476022110021,  0.00730531427893},
                                  {0.792163708315,  0.968709140644,  0.0352176264845},
                                  {0.585823963851,  0.611515761521,  0.0102566326533}};
    
    //Put grad_phi into matrix form
    deformation_measures::assemble_grad_chi(grad_phi0_tmp,Matrix_3x3::Identity(),grad_phi0);
    
    Vector_27 grad_phi0_vec; //The vector form of grad_chi
    deformation_measures::voigt_3x9_tensor(grad_phi0,grad_phi0_vec);
    
    std::vector<double> x0;
    x0.resize(27);
    for (int i=0; i<27; i++){x0[i] = grad_phi0_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_grad_chi_grad_phi, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dgrad_chidgrad_phi_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dgrad_chidgrad_phi_vec.size(); i++){
        for (int j=0; j<dgrad_chidgrad_phi_vec[i].size(); j++){
            dgrad_chidgrad_phi(j,i) = dgrad_chidgrad_phi_vec[i][j];
        }
    }
    
    micro_material::compute_dgrad_chidgrad_phi(F,_dgrad_chidgrad_phi);
    
    Matrix_27x27 _dgrad_chidgrad_phi_dense = _dgrad_chidgrad_phi;
    bool tot_result = dgrad_chidgrad_phi.isApprox(_dgrad_chidgrad_phi_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dgrad_chidgrad_phi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dgrad_chidgrad_phi & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_compute_dgrad_chidF(std::ofstream &results){
    /*!==================================
    |    test_compute_dgrad_chidF    |
    ==================================
    
    Test the computation of the gradient of chi w.r.t. 
    the deformation gradient.
    
    Note: We assume that the phi is the gradient in the 
          current configuration.
    */
    
    Matrix_27x9 dgrad_chidF; //The expected result
    std::vector< std::vector<double> > dgrad_chidF_vec; //Output from the finite_difference
    SpMat _dgrad_chidF(27,9); //The function result
    
    
    //Define the required fundamental measures
    Matrix_3x3 F0;
    define_deformation_gradient(F0);
    Matrix_3x9 grad_phi;
    define_grad_phi(grad_phi);
    Vector_27 grad_phi_voigt;
    double grad_phi_tmp[9][3] = {{0.268677941492,  0.0956706826016, 0.325269776806},
                                 {0.638467644829,  0.997489522961,  0.173494630408},
                                 {0.81125862445,   0.224706569702,  0.514463248854},
                                 {0.748926548687,  0.1079833322,    0.797424154102},
                                 {0.0782798344482, 0.324585039367,  0.032983492952},
                                 {0.84359324251,   0.602202553902,  0.165660665757},
                                 {0.117005479066,  0.476022110021,  0.00730531427893},
                                 {0.792163708315,  0.968709140644,  0.0352176264845},
                                 {0.585823963851,  0.611515761521,  0.0102566326533}};
    
    //Put grad_phi into matrix form
    deformation_measures::assemble_grad_chi(grad_phi_tmp,Matrix_3x3::Identity(),grad_phi);
    deformation_measures::voigt_3x9_tensor(grad_phi,grad_phi_voigt);
    
    Vector_9 F_vec; //The vector form of grad_chi
    deformation_measures::voigt_3x3_tensor(F0,F_vec);
    
    std::vector<double> x0;
    x0.resize(9);
    for (int i=0; i<9; i++){x0[i] = F_vec(i);}
    
    //Initialize the finite difference operator
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_grad_chi_F, 2, x0 , 1e-6);
    
    //Compute the numeric gradient
    dgrad_chidF_vec = fd.numeric_gradient();
    
    //Populate the expected result for easy comparison.
    for (int i=0; i<dgrad_chidF_vec.size(); i++){
        for (int j=0; j<dgrad_chidF_vec[i].size(); j++){
            dgrad_chidF(j,i) = dgrad_chidF_vec[i][j];
        }
    }
    
    micro_material::compute_dgrad_chidF(grad_phi_voigt,_dgrad_chidF);
    
    Matrix_27x9 _dgrad_chidF_dense = _dgrad_chidF;
    bool tot_result = dgrad_chidF.isApprox(_dgrad_chidF_dense,1e-6);
    
    if (tot_result){
        results << "test_compute_dgrad_chidF & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_dgrad_chidF & False\\\\\n\\hline\n";
    }

    return 1;
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
    test_compute_A_voigt(results);
    test_compute_B_voigt(results);
    test_compute_C_voigt(results);
    test_compute_D_voigt(results);
    test_compute_PK2_stress(results);
    test_compute_symmetric_stress(results);
    test_compute_higher_order_stress(results);
    test_map_stresses_to_current_configuration(results);
    test_compute_dPK2dRCG(results);
    test_compute_dSIGMAdRCG(results);
    test_compute_dPK2dPsi(results);
    test_compute_dSIGMAdPsi(results);
    test_compute_dPK2dGamma(results);
    test_compute_dAinvdA(results);
    test_compute_dRCGdF(results);
    test_compute_dPsidF(results);
    test_compute_dPsidchi(results);
    test_compute_dGammadF(results);
    test_compute_dGammadgrad_chi(results);
    test_compute_dgrad_chidgrad_phi(results);
    test_compute_dgrad_chidF(results);

    //Close the results file
    results.close();
}

