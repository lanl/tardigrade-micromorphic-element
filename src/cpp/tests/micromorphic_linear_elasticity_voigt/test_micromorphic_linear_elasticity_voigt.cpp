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
    /*====================
      |    define_chi    |
      ====================

      Define the micro-deformation tensor to be used.

    */

    Matrix_3x3 phi;
    define_phi(phi);

    chi = phi + Matrix_3x3::Identity();
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

    //Close the results file
    results.close();
}

