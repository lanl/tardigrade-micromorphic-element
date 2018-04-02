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
    tripletList.push_back(T(3,3,1.24421754207966));
    tripletList.push_back(T(3,6,1.24421754207966));
    tripletList.push_back(T(4,4,1.24421754207966));
    tripletList.push_back(T(4,7,1.24421754207966));
    tripletList.push_back(T(5,5,1.24421754207966));
    tripletList.push_back(T(5,8,1.24421754207966));
    tripletList.push_back(T(6,3,1.24421754207966));
    tripletList.push_back(T(6,6,1.24421754207966));
    tripletList.push_back(T(7,4,1.24421754207966));
    tripletList.push_back(T(7,7,1.24421754207966));
    tripletList.push_back(T(8,5,1.24421754207966));
    tripletList.push_back(T(8,8,1.24421754207966));

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
    tripletList.push_back(T(3,3,0.552928510286193));
    tripletList.push_back(T(3,6,0.552928510286193));
    tripletList.push_back(T(4,4,0.552928510286193));
    tripletList.push_back(T(4,7,0.552928510286193));
    tripletList.push_back(T(5,5,0.552928510286193));
    tripletList.push_back(T(5,8,0.552928510286193));
    tripletList.push_back(T(6,3,0.552928510286193));
    tripletList.push_back(T(6,6,0.552928510286193));
    tripletList.push_back(T(7,4,0.552928510286193));
    tripletList.push_back(T(7,7,0.552928510286193));
    tripletList.push_back(T(8,5,0.552928510286193));
    tripletList.push_back(T(8,8,0.552928510286193));

    D.setFromTriplets(tripletList.begin(), tripletList.end());
    return;
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

    define_parameters(params);
    define_deformation_gradient(F);
    define_chi(chi);
    define_grad_phi(grad_chi);  //Note: grad_chi = grad_phi

//    define_PK2(PK2);
    micro_material::get_stress(t, dt, params, F, chi, grad_chi, SDVS, _PK2);

    //std::cout << PK2;

    bool tot_result = PK2.isApprox(_PK2);

    if (tot_result){
        results << "test_compute_PK2_voigt & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_PK2_voigt & False\\\\\n\\hline\n";
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

    //Close the results file
    results.close();
}

