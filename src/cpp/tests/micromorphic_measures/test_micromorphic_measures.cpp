/*!============================================================================
   |                                                                          |
   |                      test_micromorphic_measures.cpp                      |
   |                                                                          |
   ----------------------------------------------------------------------------
   | The unit test file for micromorphic_measures.h/cpp. This file tests the  |
   | classes and functions defined in micromorphic_measures.h/cpp.            |
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
#include<micromorphic_measures.h>

void voigt_3x3(const Matrix_3x3 &A, Vector_9 &v){
    /*===================
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
    /*=======================
      |    define_grad_u    |
      =======================

      Define the gradient of u to be used

    */

    grad_u << 0.69646919,  0.28613933,  0.22685145,  0.55131477,  0.71946897,
              0.42310646,  0.9807642 ,  0.68482974,  0.4809319;

    return;
}

void define_phi(Matrix_3x3 &phi){
    /*====================
      |    define_phi    |
      ====================

      Define the values of phi to be used.

    */

    phi << 0.39211752,  0.34317802,  0.72904971,  0.43857224,  0.0596779,
           0.39804426,  0.73799541,  0.18249173,  0.17545176;
}

void define_grad_phi(Matrix_3x9 &grad_phi){
    /*=========================
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
    /*=====================================
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

int test_get_deformation_gradient(std::ofstream &results){
    /*=======================================
      |    test_get_deformation_gradient    |
      =======================================

      Test the get_deformation_gradient function in 
      micromorphic_measures.h/cpp

    */

    #
    Matrix_3x3 grad_u; //gradient of displacements
    Matrix_3x3 _F;     //F returned by function being tested
    Matrix_3x3 F;      //Baseline F

    define_grad_u(grad_u);          //Populate grad_u
    define_deformation_gradient(F); //Populate F

    //Prepare grad_u for function
    double _grad_u[3][3];

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _grad_u[i][j] = grad_u(i,j);
        }
    }

    //Execute the function
    micromorphic_measures::get_deformation_gradient(_grad_u,_F);

    //Compare the results
    bool tot_result = F.isApprox(_F);

    if (tot_result){
        results << "test_get_deformation_gradient & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_deformation_gradient & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_assemble_chi(std::ofstream &results){
    /*===========================
      |    test_assemble_chi    |
      ===========================

      Test the assembly of the micro-deformation 
      tensor chi.

    */

    Matrix_3x3 chi;  //The reference value of chi
    Matrix_3x3 _chi; //The value of chi from the function being tested
    Matrix_3x3 phi;  //The reference value of phi
    define_phi(phi);
    Vector_9 phi_voigt;     //phi in voigt notation
    voigt_3x3(phi,phi_voigt); //Put phi in voigt notation

    double v[9]; //Form of phi required for the function being tested
    for(int i=0; i<9; i++){v[i] = phi_voigt[i];}

    define_chi(chi); //Set the reference value of chi

    micromorphic_measures::assemble_chi(v,_chi);

    bool tot_result = chi.isApprox(_chi);

    if (tot_result){
        results << "test_assemble_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_assemble_chi & False\\\\\n\\hline\n";
    }

    return 1;

}

int test_assemble_grad_chi(std::ofstream &results){
    /*================================
      |    test_assemble_grad_chi    |
      ================================

      Test the assembly of the gradient of the micro-deformation 
      tensor chi.

    */

    Matrix_3x9 grad_chi;  //The expected result
    Matrix_3x9 _grad_chi; //The result from the function

    Matrix_3x3 F;         //The deformation gradient
    define_deformation_gradient(F); //Populate the deformation gradient

    define_grad_phi(grad_chi); //The gradient of phi=grad_chi

    Eigen::Matrix<double,9,3> grad_chi_data; //The incoming data
    grad_chi_data << 0.53155137,  0.53182759,  0.63440096,  0.09210494,  0.43370117,
                     0.43086276,  0.31728548,  0.41482621,  0.86630916,  0.4936851 ,
                     0.42583029,  0.31226122,  0.72244338,  0.32295891,  0.36178866,
                     0.84943179,  0.72445532,  0.61102351,  0.50183668,  0.62395295,
                     0.1156184 ,  0.42635131,  0.89338916,  0.94416002,  0.22826323,
                     0.29371405,  0.63097612;

    double grad_chi_data_array[9][3]; //The format expected by the function

    for (int i=0; i<9; i++){

        for (int j=0; j<3; j++){

            grad_chi_data_array[i][j] = grad_chi_data(i,j);

        }

    }

    micromorphic_measures::assemble_grad_chi(grad_chi_data_array, F, _grad_chi);

    bool tot_result = grad_chi.isApprox(_grad_chi,1e-7);

    if (tot_result){
        results << "test_assemble_chi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_assemble_chi & False\\\\\n\\hline\n";
    }

    return 1;

}

int test_get_right_cauchy_green(std::ofstream &results){
    /*=====================================
      |    test_get_right_cauchy_green    |
      =====================================

      Test the computation of the right cauchy green 
      deformation tensor.

    */

    Matrix_3x3 RCG;  //The expected result
    Matrix_3x3 _RCG; //The result from the function

    Matrix_3x3 F;    //The deformation gradient
    define_deformation_gradient(F); //Set the value of the deformation gradient.

    RCG = F.transpose()*F;

    micromorphic_measures::get_right_cauchy_green(F,_RCG);

    bool tot_result = RCG.isApprox(_RCG);

    if (tot_result){
        results << "test_get_right_cauchy_green & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_right_cauchy_green & False\\\\\n\\hline\n";
    }

    return 1;

}


int test_get_left_cauchy_green(std::ofstream &results){
    /*====================================
      |    test_get_left_cauchy_green    |
      ====================================

      Test the computation of the left cauchy green 
      deformation tensor.

    */

    Matrix_3x3 LCG;  //The expected result
    Matrix_3x3 _LCG; //The result from the function

    Matrix_3x3 F;    //The deformation gradient
    define_deformation_gradient(F); //Set the value of the deformation gradient.

    LCG = F*F.transpose();

    micromorphic_measures::get_left_cauchy_green(F,_LCG);

    bool tot_result = LCG.isApprox(_LCG);

    if (tot_result){
        results << "test_get_left_cauchy_green & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_left_cauchy_green & False\\\\\n\\hline\n";
    }

    return 1;

}

int test_get_lagrange_strain(std::ofstream &results){
    /*==================================
      |    test_get_lagrange_strain    |
      ==================================

      Test the computation of the lagrange strain.

    */

    Matrix_3x3 E;  //The expected result
    Matrix_3x3 _E; //The result from the function

    Matrix_3x3 F;  //The deformation gradient
    define_deformation_gradient(F); //Set the value of the deformation gradient

    E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());

    micromorphic_measures::get_lagrange_strain(F,_E);

    bool tot_result = E.isApprox(_E);

    if (tot_result){
        results << "test_get_lagrange_strain & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_lagrange_strain & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_get_almansi_strain(std::ofstream &results){
    /*=================================
      |    test_get_almansi_strain    |
      =================================

      Test the computation of the almansi strain.

    */

    Matrix_3x3 e;  //The expected result
    Matrix_3x3 _e; //The result from the function

    Matrix_3x3 F;  //The deformation gradient
    define_deformation_gradient(F); //Set the value of the deformation gradient

    e = 0.5*(Matrix_3x3::Identity() - (F*F.transpose()).inverse());

    micromorphic_measures::get_almansi_strain(F,_e);

    bool tot_result = e.isApprox(_e);

    if (tot_result){
        results << "test_get_almansi_strain & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_almansi_strain & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_get_small_strain(std::ofstream &results){
    /*===============================
      |    test_get_small_strain    |
      ===============================

      Test the computation of the small strain.

    */

    Matrix_3x3 epsilon;  //The expected result
    Matrix_3x3 _epsilon; //The result from the function

    Matrix_3x3 grad_u;   //The gradient of the displacement
    define_grad_u(grad_u); //Define the gradient of the displacement

    epsilon = 0.5*(grad_u + grad_u.transpose());

    double _grad_u[3][3]; //Convert grad_u to the form expected by the function

    for (int i=0; i<3; i++){

        for (int j=0; j<3; j++){

            _grad_u[i][j] = grad_u(i,j);

        }

    }

    micromorphic_measures::get_small_strain(_grad_u,_epsilon);

    bool tot_result = epsilon.isApprox(_epsilon);

    if (tot_result){
        results << "test_get_small_strain & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_small_strain & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_get_psi(std::ofstream &results){
    /*======================
      |    test_get_psi    |
      ======================

      Test the computation of psi.

    */

    Matrix_3x3 psi;  //The expected result
    Matrix_3x3 _psi; //The result from the function

    Matrix_3x3 chi; //The micro-deformation tensor
    Matrix_3x3 F;   //The deformation gradient

    define_deformation_gradient(F); //Populate the deformation gradient
    define_chi(chi);                //Populate chi

    psi = F.transpose()*chi;

    micromorphic_measures::get_psi(F,chi,_psi);

    bool tot_result = psi.isApprox(_psi);

    if (tot_result){
        results << "test_get_psi & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_psi & False\\\\\n\\hline\n";
    }

    return 1;

}

int test_get_gamma(std::ofstream &results){
    /*========================
      |    test_get_gamma    |
      ========================

      Test the computation of gamma.

    */

    Matrix_3x9 Gamma;  //The expected result
    Matrix_3x9 _Gamma; //The result from the function

    Matrix_3x9 grad_chi; //The gradient of chi
    Matrix_3x3 F;        //The deformation gradient

    define_deformation_gradient(F); //Populate the deformation gradient
    define_grad_phi(grad_chi);      //Populate the gradient of chi (grad_phi = grad_chi)

    Gamma = F.transpose()*grad_chi;

    micromorphic_measures::get_gamma(F,grad_chi,_Gamma);

    bool tot_result = Gamma.isApprox(_Gamma);

    if (tot_result){
        results << "test_get_gamma & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_gamma & False\\\\\n\\hline\n";
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
    test_get_deformation_gradient(results);
    test_assemble_chi(results);
    test_assemble_grad_chi(results);
    test_get_right_cauchy_green(results);
    test_get_left_cauchy_green(results);
    test_get_lagrange_strain(results);
    test_get_almansi_strain(results);
    test_get_small_strain(results);
    test_get_psi(results);
    test_get_gamma(results);

    //Close the results file
    results.close();
}

