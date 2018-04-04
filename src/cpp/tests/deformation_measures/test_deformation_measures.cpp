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
#include<vector>
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
    /*!====================
      |    define_chi    |
      ====================

      Define the micro-deformation tensor to be used.

    */

    Matrix_3x3 phi;
    define_phi(phi);

    chi = phi + Matrix_3x3::Identity();
}

int test_get_deformation_gradient(std::ofstream &results){
    /*!=======================================
       |    test_get_deformation_gradient    |
       =======================================

       Test the get_deformation_gradient function in 
       deformation_measures.h/cpp

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
    deformation_measures::get_deformation_gradient(_grad_u,_F);

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
    /*!===========================
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

    deformation_measures::assemble_chi(v,_chi);

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
    /*!================================
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

    deformation_measures::assemble_grad_chi(grad_chi_data_array, F, _grad_chi);

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
    /*!=====================================
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

    deformation_measures::get_right_cauchy_green(F,_RCG);

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
    /*!====================================
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

    deformation_measures::get_left_cauchy_green(F,_LCG);

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
    /*!==================================
       |    test_get_lagrange_strain    |
       ==================================

       Test the computation of the lagrange strain.

    */

    Matrix_3x3 E;  //The expected result
    Matrix_3x3 _E; //The result from the function

    Matrix_3x3 F;  //The deformation gradient
    define_deformation_gradient(F); //Set the value of the deformation gradient

    E = 0.5*(F.transpose()*F - Matrix_3x3::Identity());

    deformation_measures::get_lagrange_strain(F,_E);

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
    /*!=================================
       |    test_get_almansi_strain    |
       =================================

       Test the computation of the almansi strain.

    */

    Matrix_3x3 e;  //The expected result
    Matrix_3x3 _e; //The result from the function

    Matrix_3x3 F;  //The deformation gradient
    define_deformation_gradient(F); //Set the value of the deformation gradient

    e = 0.5*(Matrix_3x3::Identity() - (F*F.transpose()).inverse());

    deformation_measures::get_almansi_strain(F,_e);

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
    /*!===============================
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

    deformation_measures::get_small_strain(_grad_u,_epsilon);

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
    /*!======================
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

    deformation_measures::get_psi(F,chi,_psi);

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
    /*!========================
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

    deformation_measures::get_gamma(F,grad_chi,_Gamma);

    bool tot_result = Gamma.isApprox(_Gamma);

    if (tot_result){
        results << "test_get_gamma & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_gamma & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_voigt_3x3_symm_tensor(std::ofstream &results){
    /*!====================================
       |    test_voigt_3x3_symm_tensor    |
       ====================================

       Test the mapping of a symmetric tensor to voigt
       notation.

    */

    Matrix_3x3 A; //The tensor to map
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9; //Populate the tensor

    Vector_6 v;  //The expected result
    Vector_6 _v; //The function result

    v(0) = A(0,0);
    v(1) = A(1,1);
    v(2) = A(2,2);
    v(3) = A(1,2);
    v(4) = A(0,2);
    v(5) = A(0,1);

    deformation_measures::voigt_3x3_symm_tensor(A,_v);

    bool tot_result = v.isApprox(_v);

    if (tot_result){
        results << "test_voigt_3x3_symm_tensor & True\\\\\n\\hline\n";
    }
    else {
        results << "test_voigt_3x3_symm_tensor & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_voigt_3x3_tensor(std::ofstream &results){
    /*!===============================
       |    test_voigt_3x3_tensor    |
       ===============================

       Test the mapping of a general tensor to voigt
       notation.

    */

    Matrix_3x3 A; //The tensor to map
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9; //Populate the tensor

    Vector_9 v;  //The expected result
    Vector_9 _v; //The function result

    v(0) = A(0,0);
    v(1) = A(1,1);
    v(2) = A(2,2);
    v(3) = A(1,2);
    v(4) = A(0,2);
    v(5) = A(0,1);
    v(6) = A(2,1);
    v(7) = A(2,0);
    v(8) = A(1,0);

    deformation_measures::voigt_3x3_tensor(A,_v);

    bool tot_result = v.isApprox(_v);

    if (tot_result){
        results << "test_voigt_3x3_tensor & True\\\\\n\\hline\n";
    }
    else {
        results << "test_voigt_3x3_tensor & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_undo_voigt_3x3_tensor(std::ofstream &results){
    /*!====================================
       |    test_undo_voigt_3x3_tensor    |
       ====================================

       Test the mapping of a general 3x3 tensor in voigt 
       notation back to the matrix form.

    */

    Matrix_3x3 A;  //The expected result
    Matrix_3x3 _A; //The function result
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9; //Populate the tensor

    Vector_9 v;  //The matrix in voigt form

    deformation_measures::voigt_3x3_tensor(A,v);

    deformation_measures::undo_voigt_3x3_tensor(v,_A);
    
    bool tot_result = A.isApprox(_A);

    if (tot_result){
        results << "test_undo_voigt_3x3_tensor & True\\\\\n\\hline\n";
    }
    else {
        results << "test_undo_voigt_3x3_tensor & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_voigt_3x9_tensor(std::ofstream &results){
    /*!===============================
       |    test_voigt_3x9_tensor    |
       ===============================

       Test the computation of the voigt notation 
       for a third order tensor.

    */

    Matrix_3x9 A; //The tensor to map
    Vector_27 v;  //The expected result
    Vector_27 _v; //The function result
    int a = 1;
    for (int i=0; i<3; i++){

        for (int j=0; j<9; j++){

            A(i,j) = a;
            v(i*9+j) = a;
            a++;

        }

    }

    deformation_measures::voigt_3x9_tensor(A,_v);

    bool tot_result = v.isApprox(_v);

    if (tot_result){
        results << "test_voigt_3x9_tensor & True\\\\\n\\hline\n";
    }
    else {
        results << "test_voigt_3x9_tensor & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_get_micro_strain(std::ofstream &results){
    /*!===============================
       |    test_get_micro_strain    |
       ===============================

       Test the computation of the deformation strain.

    */

    Matrix_3x3 E_micro;  //The expected result
    Matrix_3x3 _E_micro; //The function result

    Matrix_3x3 chi; //The micro-deformation tensor
    Matrix_3x3 F;   //The deformation gradient

    define_deformation_gradient(F); //Populate the deformation gradient
    define_chi(chi);                //Populate chi

    E_micro = F.transpose()*chi - Matrix_3x3::Identity();

    deformation_measures::get_micro_strain(F.transpose()*chi,_E_micro);

    bool tot_result = E_micro.isApprox(_E_micro);

    if (tot_result){
        results << "test_get_micro_strain & True\\\\\n\\hline\n";
    }
    else {
        results << "test_get_micro_strain & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_undo_voigt_3x9_tensor(std::ofstream &results){
    /*!====================================
       |    test_undo_voigt_3x9_tensor    |
       ====================================

       Test undoing the computation of the voigt notation 
       for a third order tensor.

    */

    Matrix_3x9 A;  //The expected result
    Matrix_3x9 _A; //The function result
    Vector_27 v;  //The expected result

    deformation_measures::voigt_3x9_tensor(A,v);
    deformation_measures::undo_voigt_3x9_tensor(v,_A);
    
    bool tot_result = A.isApprox(_A);

    if (tot_result){
        results << "test_undo_voigt_3x9_tensor & True\\\\\n\\hline\n";
    }
    else {
        results << "test_undo_voigt_3x9_tensor & False\\\\\n\\hline\n";
    }

    return 1;
}

int test_perform_left_positive_cyclic_permutation(std::ofstream &results){
    /*!=======================================================
       |    test_perform_left_positive_cyclic_permutation    |
       =======================================================
    
        Test performing a left positive cyclic permutation on a 3x9 (third order)
        tensor.
    
    */
    
    Vector_27 A;  //The expected result
    Vector_27 _A; //The function result
    
    A  << 0, 10, 20, 19, 18, 9, 11, 2, 1, 3, 13, 23, 22, 21, 12, 14, 5, 4, 6, 16, 26, 25, 24, 15, 17, 8, 7;
    _A << 0, 4, 8, 5, 2, 1, 7, 6, 3, 9, 13, 17, 14, 11, 10, 16, 15, 12, 18, 22, 26, 23, 20, 19, 25, 24, 21;
    
    deformation_measures::perform_left_positive_cyclic_permutation(_A);
    
    bool tot_result = A.isApprox(_A);
    
    if (tot_result){
        results << "test_perform_left_positive_cyclic_permutation & True\\\\\n\\hline\n";
    }
    else {
        results << "test_perform_left_positive_cyclic_permutation & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_perform_right_positive_cyclic_permutation(std::ofstream &results){
    /*!========================================================
       |    test_perform_right_positive_cyclic_permutation    |
       ========================================================
    
        Test performing a right positive cyclic permutation on a 3x9 (third order)
        tensor.
    
    */
    
    Vector_27 A;  //The expected result
    Vector_27 _A; //The function result
    
    A  << 0, 12, 24, 15, 6, 3, 21, 18, 9, 1, 13, 25, 16, 7, 4, 22, 19, 10, 2, 14, 26, 17, 8, 5, 23, 20, 11;
    _A << 0, 4, 8, 5, 2, 1, 7, 6, 3, 9, 13, 17, 14, 11, 10, 16, 15, 12, 18, 22, 26, 23, 20, 19, 25, 24, 21;
    
    deformation_measures::perform_right_positive_cyclic_permutation(_A);
    
    bool tot_result = A.isApprox(_A);
    
    if (tot_result){
        results << "test_perform_right_positive_cyclic_permutation & True\\\\\n\\hline\n";
    }
    else {
        results << "test_perform_right_positive_cyclic_permutation & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_dot_2ot_4ot(std::ofstream &results){
    /*!==========================
    |    test_dot_2ot_4ot    |
    ==========================
    
    Test the dot product of a second order tensor with an index of a fourth order tensor.
    
    */
    
    Matrix_9x9  r[4]; //The expected results
    Matrix_9x9 _r[4]; //The output of the function
    
    //Load the expected results
    r[0] << 135,  147,  159,  150,  141,  138,  156,  153,  144,  486,  534,
            582,  546,  510,  498,  570,  558,  522,  999, 1083, 1167, 1104,
           1041, 1020, 1146, 1125, 1062,  594,  642,  690,  654,  618,  606,
            678,  666,  630,  189,  201,  213,  204,  195,  192,  210,  207,
            198,  162,  174,  186,  177,  168,  165,  183,  180,  171,  810,
            894,  978,  915,  852,  831,  957,  936,  873,  621,  705,  789,
            726,  663,  642,  768,  747,  684,  378,  426,  474,  438,  402,
            390,  462,  450,  414;
            
    r[1] << 45,   57,   69,   60,   51,   48,   66,   63,   54,  450,  498,
           546,  510,  474,  462,  534,  522,  486, 1341, 1425, 1509, 1446,
          1383, 1362, 1488, 1467, 1404,  774,  822,  870,  834,  798,  786,
           858,  846,  810,  207,  219,  231,  222,  213,  210,  228,  225,
           216,  126,  138,  150,  141,  132,  129,  147,  144,  135,  774,
           858,  942,  879,  816,  795,  921,  900,  837,  207,  291,  375,
           312,  249,  228,  354,  333,  270,  126,  174,  222,  186,  150,
           138,  210,  198,  162;
           
    r[2] << 15,   45,   75,   48,   21,   18,   72,   69,   42,  366,  486,
           606,  498,  390,  378,  594,  582,  474, 1203, 1413, 1623, 1434,
          1245, 1224, 1602, 1581, 1392,  690,  810,  930,  822,  714,  702,
           918,  906,  798,  177,  207,  237,  210,  183,  180,  234,  231,
           204,   96,  126,  156,  129,  102,   99,  153,  150,  123,  636,
           846, 1056,  867,  678,  657, 1035, 1014,  825,   69,  279,  489,
           300,  111,   90,  468,  447,  258,   42,  162,  282,  174,   66,
            54,  270,  258,  150;
            
    r[3] << 5,   41,   77,   50,   23,   14,   68,   59,   32,  338,  482,
          626,  518,  410,  374,  590,  554,  446, 1157, 1409, 1661, 1472,
         1283, 1220, 1598, 1535, 1346,  662,  806,  950,  842,  734,  698,
          914,  878,  770,  167,  203,  239,  212,  185,  176,  230,  221,
          194,   86,  122,  158,  131,  104,   95,  149,  140,  113,  590,
          842, 1094,  905,  716,  653, 1031,  968,  779,   23,  275,  527,
          338,  149,   86,  464,  401,  212,   14,  158,  302,  194,   86,
           50,  266,  230,  122;
           
    //Load the generating matrices
    Matrix_3x3 SOT; //The generating second order tensor
    Matrix_9x9 FOT; //The generating fourth order tensor
    
    SOT << 0, 1, 2, 3, 4, 5, 6, 7, 8;
    FOT << 0,  4,  8,  5,  2,  1,  7,  6,  3, 36, 40, 44, 41, 38, 37, 43, 42,
          39, 72, 76, 80, 77, 74, 73, 79, 78, 75, 45, 49, 53, 50, 47, 46, 52,
          51, 48, 18, 22, 26, 23, 20, 19, 25, 24, 21,  9, 13, 17, 14, 11, 10,
          16, 15, 12, 63, 67, 71, 68, 65, 64, 70, 69, 66, 54, 58, 62, 59, 56,
          55, 61, 60, 57, 27, 31, 35, 32, 29, 28, 34, 33, 30;
    
    for (int i=0; i<4; i++){deformation_measures::dot_2ot_4ot(i,SOT,FOT,_r[i]);}
    
    bool tot_result = true;
    for (int i=0; i<4; i++){tot_result *= r[i].isApprox(_r[i]);}
    
    if (tot_result){
        results << "test_dot_2ot_4ot & True\\\\\n\\hline\n";
    }
    else {
        results << "test_dot_2ot_4ot & False\\\\\n\\hline\n";
    }
    
    return 1;
    
}

int test_dot_2ot_5ot(std::ofstream &results){
    /*!=======================
    |    test_dot_2ot_5ot    |
    ==========================
    
    Compute the dot product of a second order tensor 
    and a fifth order tensor.
    
    */
    
    Matrix_9x27  r_m1[5]; //The expected results for mode 1
    Matrix_9x27 _r_m1[5]; //The function output for mode 1
    
    Matrix_27x9  r_m2[5]; //The expected results for mode 2
    Matrix_27x9 _r_m2[5]; //The function output for mode 2
    
    r_m1[0] << 405,  417,  429,  420,  411,  408,  426,  423,  414,  432,  444,
               456,  447,  438,  435,  453,  450,  441,  459,  471,  483,  474,
               465,  462,  480,  477,  468, 1458, 1506, 1554, 1518, 1482, 1470,
              1542, 1530, 1494, 1566, 1614, 1662, 1626, 1590, 1578, 1650, 1638,
              1602, 1674, 1722, 1770, 1734, 1698, 1686, 1758, 1746, 1710, 2997,
              3081, 3165, 3102, 3039, 3018, 3144, 3123, 3060, 3186, 3270, 3354,
              3291, 3228, 3207, 3333, 3312, 3249, 3375, 3459, 3543, 3480, 3417,
              3396, 3522, 3501, 3438, 1782, 1830, 1878, 1842, 1806, 1794, 1866,
              1854, 1818, 1890, 1938, 1986, 1950, 1914, 1902, 1974, 1962, 1926,
              1998, 2046, 2094, 2058, 2022, 2010, 2082, 2070, 2034,  567,  579,
               591,  582,  573,  570,  588,  585,  576,  594,  606,  618,  609,
               600,  597,  615,  612,  603,  621,  633,  645,  636,  627,  624,
               642,  639,  630,  486,  498,  510,  501,  492,  489,  507,  504,
               495,  513,  525,  537,  528,  519,  516,  534,  531,  522,  540,
               552,  564,  555,  546,  543,  561,  558,  549, 2430, 2514, 2598,
              2535, 2472, 2451, 2577, 2556, 2493, 2619, 2703, 2787, 2724, 2661,
              2640, 2766, 2745, 2682, 2808, 2892, 2976, 2913, 2850, 2829, 2955,
              2934, 2871, 1863, 1947, 2031, 1968, 1905, 1884, 2010, 1989, 1926,
              2052, 2136, 2220, 2157, 2094, 2073, 2199, 2178, 2115, 2241, 2325,
              2409, 2346, 2283, 2262, 2388, 2367, 2304, 1134, 1182, 1230, 1194,
              1158, 1146, 1218, 1206, 1170, 1242, 1290, 1338, 1302, 1266, 1254,
              1326, 1314, 1278, 1350, 1398, 1446, 1410, 1374, 1362, 1434, 1422,
              1386;
              
    r_m1[1] << 135,  147,  159,  150,  141,  138,  156,  153,  144,  162,  174,
               186,  177,  168,  165,  183,  180,  171,  189,  201,  213,  204,
               195,  192,  210,  207,  198, 1350, 1398, 1446, 1410, 1374, 1362,
              1434, 1422, 1386, 1458, 1506, 1554, 1518, 1482, 1470, 1542, 1530,
              1494, 1566, 1614, 1662, 1626, 1590, 1578, 1650, 1638, 1602, 4023,
              4107, 4191, 4128, 4065, 4044, 4170, 4149, 4086, 4212, 4296, 4380,
              4317, 4254, 4233, 4359, 4338, 4275, 4401, 4485, 4569, 4506, 4443,
              4422, 4548, 4527, 4464, 2322, 2370, 2418, 2382, 2346, 2334, 2406,
              2394, 2358, 2430, 2478, 2526, 2490, 2454, 2442, 2514, 2502, 2466,
              2538, 2586, 2634, 2598, 2562, 2550, 2622, 2610, 2574,  621,  633,
               645,  636,  627,  624,  642,  639,  630,  648,  660,  672,  663,
               654,  651,  669,  666,  657,  675,  687,  699,  690,  681,  678,
               696,  693,  684,  378,  390,  402,  393,  384,  381,  399,  396,
               387,  405,  417,  429,  420,  411,  408,  426,  423,  414,  432,
               444,  456,  447,  438,  435,  453,  450,  441, 2322, 2406, 2490,
              2427, 2364, 2343, 2469, 2448, 2385, 2511, 2595, 2679, 2616, 2553,
              2532, 2658, 2637, 2574, 2700, 2784, 2868, 2805, 2742, 2721, 2847,
              2826, 2763,  621,  705,  789,  726,  663,  642,  768,  747,  684,
               810,  894,  978,  915,  852,  831,  957,  936,  873,  999, 1083,
              1167, 1104, 1041, 1020, 1146, 1125, 1062,  378,  426,  474,  438,
               402,  390,  462,  450,  414,  486,  534,  582,  546,  510,  498,
               570,  558,  522,  594,  642,  690,  654,  618,  606,  678,  666,
               630;
               
    r_m1[2] << 45,   57,   69,   60,   51,   48,   66,   63,   54,  126,  138,
              150,  141,  132,  129,  147,  144,  135,  207,  219,  231,  222,
              213,  210,  228,  225,  216, 1098, 1146, 1194, 1158, 1122, 1110,
             1182, 1170, 1134, 1422, 1470, 1518, 1482, 1446, 1434, 1506, 1494,
             1458, 1746, 1794, 1842, 1806, 1770, 1758, 1830, 1818, 1782, 3609,
             3693, 3777, 3714, 3651, 3630, 3756, 3735, 3672, 4176, 4260, 4344,
             4281, 4218, 4197, 4323, 4302, 4239, 4743, 4827, 4911, 4848, 4785,
             4764, 4890, 4869, 4806, 2070, 2118, 2166, 2130, 2094, 2082, 2154,
             2142, 2106, 2394, 2442, 2490, 2454, 2418, 2406, 2478, 2466, 2430,
             2718, 2766, 2814, 2778, 2742, 2730, 2802, 2790, 2754,  531,  543,
              555,  546,  537,  534,  552,  549,  540,  612,  624,  636,  627,
              618,  615,  633,  630,  621,  693,  705,  717,  708,  699,  696,
              714,  711,  702,  288,  300,  312,  303,  294,  291,  309,  306,
              297,  369,  381,  393,  384,  375,  372,  390,  387,  378,  450,
              462,  474,  465,  456,  453,  471,  468,  459, 1908, 1992, 2076,
             2013, 1950, 1929, 2055, 2034, 1971, 2475, 2559, 2643, 2580, 2517,
             2496, 2622, 2601, 2538, 3042, 3126, 3210, 3147, 3084, 3063, 3189,
             3168, 3105,  207,  291,  375,  312,  249,  228,  354,  333,  270,
              774,  858,  942,  879,  816,  795,  921,  900,  837, 1341, 1425,
             1509, 1446, 1383, 1362, 1488, 1467, 1404,  126,  174,  222,  186,
              150,  138,  210,  198,  162,  450,  498,  546,  510,  474,  462,
              534,  522,  486,  774,  822,  870,  834,  798,  786,  858,  846,
              810;
              
    r_m1[3] << 15,   45,   75,   48,   21,   18,   72,   69,   42,   96,  126,
              156,  129,  102,   99,  153,  150,  123,  177,  207,  237,  210,
              183,  180,  234,  231,  204, 1014, 1134, 1254, 1146, 1038, 1026,
             1242, 1230, 1122, 1338, 1458, 1578, 1470, 1362, 1350, 1566, 1554,
             1446, 1662, 1782, 1902, 1794, 1686, 1674, 1890, 1878, 1770, 3471,
             3681, 3891, 3702, 3513, 3492, 3870, 3849, 3660, 4038, 4248, 4458,
             4269, 4080, 4059, 4437, 4416, 4227, 4605, 4815, 5025, 4836, 4647,
             4626, 5004, 4983, 4794, 1986, 2106, 2226, 2118, 2010, 1998, 2214,
             2202, 2094, 2310, 2430, 2550, 2442, 2334, 2322, 2538, 2526, 2418,
             2634, 2754, 2874, 2766, 2658, 2646, 2862, 2850, 2742,  501,  531,
              561,  534,  507,  504,  558,  555,  528,  582,  612,  642,  615,
              588,  585,  639,  636,  609,  663,  693,  723,  696,  669,  666,
              720,  717,  690,  258,  288,  318,  291,  264,  261,  315,  312,
              285,  339,  369,  399,  372,  345,  342,  396,  393,  366,  420,
              450,  480,  453,  426,  423,  477,  474,  447, 1770, 1980, 2190,
             2001, 1812, 1791, 2169, 2148, 1959, 2337, 2547, 2757, 2568, 2379,
             2358, 2736, 2715, 2526, 2904, 3114, 3324, 3135, 2946, 2925, 3303,
             3282, 3093,   69,  279,  489,  300,  111,   90,  468,  447,  258,
              636,  846, 1056,  867,  678,  657, 1035, 1014,  825, 1203, 1413,
             1623, 1434, 1245, 1224, 1602, 1581, 1392,   42,  162,  282,  174,
               66,   54,  270,  258,  150,  366,  486,  606,  498,  390,  378,
              594,  582,  474,  690,  810,  930,  822,  714,  702,  918,  906,
              798;
              
    r_m1[4] << 5,   41,   77,   50,   23,   14,   68,   59,   32,   86,  122,
             158,  131,  104,   95,  149,  140,  113,  167,  203,  239,  212,
             185,  176,  230,  221,  194,  986, 1130, 1274, 1166, 1058, 1022,
            1238, 1202, 1094, 1310, 1454, 1598, 1490, 1382, 1346, 1562, 1526,
            1418, 1634, 1778, 1922, 1814, 1706, 1670, 1886, 1850, 1742, 3425,
            3677, 3929, 3740, 3551, 3488, 3866, 3803, 3614, 3992, 4244, 4496,
            4307, 4118, 4055, 4433, 4370, 4181, 4559, 4811, 5063, 4874, 4685,
            4622, 5000, 4937, 4748, 1958, 2102, 2246, 2138, 2030, 1994, 2210,
            2174, 2066, 2282, 2426, 2570, 2462, 2354, 2318, 2534, 2498, 2390,
            2606, 2750, 2894, 2786, 2678, 2642, 2858, 2822, 2714,  491,  527,
             563,  536,  509,  500,  554,  545,  518,  572,  608,  644,  617,
             590,  581,  635,  626,  599,  653,  689,  725,  698,  671,  662,
             716,  707,  680,  248,  284,  320,  293,  266,  257,  311,  302,
             275,  329,  365,  401,  374,  347,  338,  392,  383,  356,  410,
             446,  482,  455,  428,  419,  473,  464,  437, 1724, 1976, 2228,
            2039, 1850, 1787, 2165, 2102, 1913, 2291, 2543, 2795, 2606, 2417,
            2354, 2732, 2669, 2480, 2858, 3110, 3362, 3173, 2984, 2921, 3299,
            3236, 3047,   23,  275,  527,  338,  149,   86,  464,  401,  212,
             590,  842, 1094,  905,  716,  653, 1031,  968,  779, 1157, 1409,
            1661, 1472, 1283, 1220, 1598, 1535, 1346,   14,  158,  302,  194,
              86,   50,  266,  230,  122,  338,  482,  626,  518,  410,  374,
             590,  554,  446,  662,  806,  950,  842,  734,  698,  914,  878,
             770;
             
    r_m2[0] << 405,  417,  429,  420,  411,  408,  426,  423,  414,  513,  525,
               537,  528,  519,  516,  534,  531,  522,  621,  633,  645,  636,
               627,  624,  642,  639,  630,  540,  552,  564,  555,  546,  543,
               561,  558,  549,  459,  471,  483,  474,  465,  462,  480,  477,
               468,  432,  444,  456,  447,  438,  435,  453,  450,  441,  594,
               606,  618,  609,  600,  597,  615,  612,  603,  567,  579,  591,
               582,  573,  570,  588,  585,  576,  486,  498,  510,  501,  492,
               489,  507,  504,  495, 1134, 1182, 1230, 1194, 1158, 1146, 1218,
              1206, 1170, 1566, 1614, 1662, 1626, 1590, 1578, 1650, 1638, 1602,
              1998, 2046, 2094, 2058, 2022, 2010, 2082, 2070, 2034, 1674, 1722,
              1770, 1734, 1698, 1686, 1758, 1746, 1710, 1350, 1398, 1446, 1410,
              1374, 1362, 1434, 1422, 1386, 1242, 1290, 1338, 1302, 1266, 1254,
              1326, 1314, 1278, 1890, 1938, 1986, 1950, 1914, 1902, 1974, 1962,
              1926, 1782, 1830, 1878, 1842, 1806, 1794, 1866, 1854, 1818, 1458,
              1506, 1554, 1518, 1482, 1470, 1542, 1530, 1494, 1863, 1947, 2031,
              1968, 1905, 1884, 2010, 1989, 1926, 2619, 2703, 2787, 2724, 2661,
              2640, 2766, 2745, 2682, 3375, 3459, 3543, 3480, 3417, 3396, 3522,
              3501, 3438, 2808, 2892, 2976, 2913, 2850, 2829, 2955, 2934, 2871,
              2241, 2325, 2409, 2346, 2283, 2262, 2388, 2367, 2304, 2052, 2136,
              2220, 2157, 2094, 2073, 2199, 2178, 2115, 3186, 3270, 3354, 3291,
              3228, 3207, 3333, 3312, 3249, 2997, 3081, 3165, 3102, 3039, 3018,
              3144, 3123, 3060, 2430, 2514, 2598, 2535, 2472, 2451, 2577, 2556,
              2493;
              
    r_m2[1] << 135,  147,  159,  150,  141,  138,  156,  153,  144,  405,  417,
               429,  420,  411,  408,  426,  423,  414,  675,  687,  699,  690,
               681,  678,  696,  693,  684,  432,  444,  456,  447,  438,  435,
               453,  450,  441,  189,  201,  213,  204,  195,  192,  210,  207,
               198,  162,  174,  186,  177,  168,  165,  183,  180,  171,  648,
               660,  672,  663,  654,  651,  669,  666,  657,  621,  633,  645,
               636,  627,  624,  642,  639,  630,  378,  390,  402,  393,  384,
               381,  399,  396,  387,  378,  426,  474,  438,  402,  390,  462,
               450,  414, 1458, 1506, 1554, 1518, 1482, 1470, 1542, 1530, 1494,
              2538, 2586, 2634, 2598, 2562, 2550, 2622, 2610, 2574, 1566, 1614,
              1662, 1626, 1590, 1578, 1650, 1638, 1602,  594,  642,  690,  654,
               618,  606,  678,  666,  630,  486,  534,  582,  546,  510,  498,
               570,  558,  522, 2430, 2478, 2526, 2490, 2454, 2442, 2514, 2502,
              2466, 2322, 2370, 2418, 2382, 2346, 2334, 2406, 2394, 2358, 1350,
              1398, 1446, 1410, 1374, 1362, 1434, 1422, 1386,  621,  705,  789,
               726,  663,  642,  768,  747,  684, 2511, 2595, 2679, 2616, 2553,
              2532, 2658, 2637, 2574, 4401, 4485, 4569, 4506, 4443, 4422, 4548,
              4527, 4464, 2700, 2784, 2868, 2805, 2742, 2721, 2847, 2826, 2763,
               999, 1083, 1167, 1104, 1041, 1020, 1146, 1125, 1062,  810,  894,
               978,  915,  852,  831,  957,  936,  873, 4212, 4296, 4380, 4317,
              4254, 4233, 4359, 4338, 4275, 4023, 4107, 4191, 4128, 4065, 4044,
              4170, 4149, 4086, 2322, 2406, 2490, 2427, 2364, 2343, 2469, 2448,
              2385;
              
    r_m2[2] << 45,   57,   69,   60,   51,   48,   66,   63,   54,  369,  381,
              393,  384,  375,  372,  390,  387,  378,  693,  705,  717,  708,
              699,  696,  714,  711,  702,  450,  462,  474,  465,  456,  453,
              471,  468,  459,  207,  219,  231,  222,  213,  210,  228,  225,
              216,  126,  138,  150,  141,  132,  129,  147,  144,  135,  612,
              624,  636,  627,  618,  615,  633,  630,  621,  531,  543,  555,
              546,  537,  534,  552,  549,  540,  288,  300,  312,  303,  294,
              291,  309,  306,  297,  126,  174,  222,  186,  150,  138,  210,
              198,  162, 1422, 1470, 1518, 1482, 1446, 1434, 1506, 1494, 1458,
             2718, 2766, 2814, 2778, 2742, 2730, 2802, 2790, 2754, 1746, 1794,
             1842, 1806, 1770, 1758, 1830, 1818, 1782,  774,  822,  870,  834,
              798,  786,  858,  846,  810,  450,  498,  546,  510,  474,  462,
              534,  522,  486, 2394, 2442, 2490, 2454, 2418, 2406, 2478, 2466,
             2430, 2070, 2118, 2166, 2130, 2094, 2082, 2154, 2142, 2106, 1098,
             1146, 1194, 1158, 1122, 1110, 1182, 1170, 1134,  207,  291,  375,
              312,  249,  228,  354,  333,  270, 2475, 2559, 2643, 2580, 2517,
             2496, 2622, 2601, 2538, 4743, 4827, 4911, 4848, 4785, 4764, 4890,
             4869, 4806, 3042, 3126, 3210, 3147, 3084, 3063, 3189, 3168, 3105,
             1341, 1425, 1509, 1446, 1383, 1362, 1488, 1467, 1404,  774,  858,
              942,  879,  816,  795,  921,  900,  837, 4176, 4260, 4344, 4281,
             4218, 4197, 4323, 4302, 4239, 3609, 3693, 3777, 3714, 3651, 3630,
             3756, 3735, 3672, 1908, 1992, 2076, 2013, 1950, 1929, 2055, 2034,
             1971;
             
    r_m2[3] << 15,   45,   75,   48,   21,   18,   72,   69,   42,  339,  369,
              399,  372,  345,  342,  396,  393,  366,  663,  693,  723,  696,
              669,  666,  720,  717,  690,  420,  450,  480,  453,  426,  423,
              477,  474,  447,  177,  207,  237,  210,  183,  180,  234,  231,
              204,   96,  126,  156,  129,  102,   99,  153,  150,  123,  582,
              612,  642,  615,  588,  585,  639,  636,  609,  501,  531,  561,
              534,  507,  504,  558,  555,  528,  258,  288,  318,  291,  264,
              261,  315,  312,  285,   42,  162,  282,  174,   66,   54,  270,
              258,  150, 1338, 1458, 1578, 1470, 1362, 1350, 1566, 1554, 1446,
             2634, 2754, 2874, 2766, 2658, 2646, 2862, 2850, 2742, 1662, 1782,
             1902, 1794, 1686, 1674, 1890, 1878, 1770,  690,  810,  930,  822,
              714,  702,  918,  906,  798,  366,  486,  606,  498,  390,  378,
              594,  582,  474, 2310, 2430, 2550, 2442, 2334, 2322, 2538, 2526,
             2418, 1986, 2106, 2226, 2118, 2010, 1998, 2214, 2202, 2094, 1014,
             1134, 1254, 1146, 1038, 1026, 1242, 1230, 1122,   69,  279,  489,
              300,  111,   90,  468,  447,  258, 2337, 2547, 2757, 2568, 2379,
             2358, 2736, 2715, 2526, 4605, 4815, 5025, 4836, 4647, 4626, 5004,
             4983, 4794, 2904, 3114, 3324, 3135, 2946, 2925, 3303, 3282, 3093,
             1203, 1413, 1623, 1434, 1245, 1224, 1602, 1581, 1392,  636,  846,
             1056,  867,  678,  657, 1035, 1014,  825, 4038, 4248, 4458, 4269,
             4080, 4059, 4437, 4416, 4227, 3471, 3681, 3891, 3702, 3513, 3492,
             3870, 3849, 3660, 1770, 1980, 2190, 2001, 1812, 1791, 2169, 2148,
             1959;
             
    r_m2[4] << 5,   41,   77,   50,   23,   14,   68,   59,   32,  329,  365,
             401,  374,  347,  338,  392,  383,  356,  653,  689,  725,  698,
             671,  662,  716,  707,  680,  410,  446,  482,  455,  428,  419,
             473,  464,  437,  167,  203,  239,  212,  185,  176,  230,  221,
             194,   86,  122,  158,  131,  104,   95,  149,  140,  113,  572,
             608,  644,  617,  590,  581,  635,  626,  599,  491,  527,  563,
             536,  509,  500,  554,  545,  518,  248,  284,  320,  293,  266,
             257,  311,  302,  275,   14,  158,  302,  194,   86,   50,  266,
             230,  122, 1310, 1454, 1598, 1490, 1382, 1346, 1562, 1526, 1418,
            2606, 2750, 2894, 2786, 2678, 2642, 2858, 2822, 2714, 1634, 1778,
            1922, 1814, 1706, 1670, 1886, 1850, 1742,  662,  806,  950,  842,
             734,  698,  914,  878,  770,  338,  482,  626,  518,  410,  374,
             590,  554,  446, 2282, 2426, 2570, 2462, 2354, 2318, 2534, 2498,
            2390, 1958, 2102, 2246, 2138, 2030, 1994, 2210, 2174, 2066,  986,
            1130, 1274, 1166, 1058, 1022, 1238, 1202, 1094,   23,  275,  527,
             338,  149,   86,  464,  401,  212, 2291, 2543, 2795, 2606, 2417,
            2354, 2732, 2669, 2480, 4559, 4811, 5063, 4874, 4685, 4622, 5000,
            4937, 4748, 2858, 3110, 3362, 3173, 2984, 2921, 3299, 3236, 3047,
            1157, 1409, 1661, 1472, 1283, 1220, 1598, 1535, 1346,  590,  842,
            1094,  905,  716,  653, 1031,  968,  779, 3992, 4244, 4496, 4307,
            4118, 4055, 4433, 4370, 4181, 3425, 3677, 3929, 3740, 3551, 3488,
            3866, 3803, 3614, 1724, 1976, 2228, 2039, 1850, 1787, 2165, 2102,
            1913;
             
    //Load the generating matrices
    Matrix_3x3  SOT;
    Matrix_9x27 FOT_m1;
    Matrix_27x9 FOT_m2;
    
    SOT << 0, 1, 2, 3, 4, 5, 6, 7, 8;
    FOT_m1 << 0,   4,   8,   5,   2,   1,   7,   6,   3,   9,  13,  17,  14,
             11,  10,  16,  15,  12,  18,  22,  26,  23,  20,  19,  25,  24,
             21, 108, 112, 116, 113, 110, 109, 115, 114, 111, 117, 121, 125,
            122, 119, 118, 124, 123, 120, 126, 130, 134, 131, 128, 127, 133,
            132, 129, 216, 220, 224, 221, 218, 217, 223, 222, 219, 225, 229,
            233, 230, 227, 226, 232, 231, 228, 234, 238, 242, 239, 236, 235,
            241, 240, 237, 135, 139, 143, 140, 137, 136, 142, 141, 138, 144,
            148, 152, 149, 146, 145, 151, 150, 147, 153, 157, 161, 158, 155,
            154, 160, 159, 156,  54,  58,  62,  59,  56,  55,  61,  60,  57,
             63,  67,  71,  68,  65,  64,  70,  69,  66,  72,  76,  80,  77,
             74,  73,  79,  78,  75,  27,  31,  35,  32,  29,  28,  34,  33,
             30,  36,  40,  44,  41,  38,  37,  43,  42,  39,  45,  49,  53,
             50,  47,  46,  52,  51,  48, 189, 193, 197, 194, 191, 190, 196,
            195, 192, 198, 202, 206, 203, 200, 199, 205, 204, 201, 207, 211,
            215, 212, 209, 208, 214, 213, 210, 162, 166, 170, 167, 164, 163,
            169, 168, 165, 171, 175, 179, 176, 173, 172, 178, 177, 174, 180,
            184, 188, 185, 182, 181, 187, 186, 183,  81,  85,  89,  86,  83,
             82,  88,  87,  84,  90,  94,  98,  95,  92,  91,  97,  96,  93,
             99, 103, 107, 104, 101, 100, 106, 105, 102;
             
    FOT_m2 << 0,   4,   8,   5,   2,   1,   7,   6,   3,  36,  40,  44,  41,
             38,  37,  43,  42,  39,  72,  76,  80,  77,  74,  73,  79,  78,
             75,  45,  49,  53,  50,  47,  46,  52,  51,  48,  18,  22,  26,
             23,  20,  19,  25,  24,  21,   9,  13,  17,  14,  11,  10,  16,
             15,  12,  63,  67,  71,  68,  65,  64,  70,  69,  66,  54,  58,
             62,  59,  56,  55,  61,  60,  57,  27,  31,  35,  32,  29,  28,
             34,  33,  30,  81,  85,  89,  86,  83,  82,  88,  87,  84, 117,
            121, 125, 122, 119, 118, 124, 123, 120, 153, 157, 161, 158, 155,
            154, 160, 159, 156, 126, 130, 134, 131, 128, 127, 133, 132, 129,
             99, 103, 107, 104, 101, 100, 106, 105, 102,  90,  94,  98,  95,
             92,  91,  97,  96,  93, 144, 148, 152, 149, 146, 145, 151, 150,
            147, 135, 139, 143, 140, 137, 136, 142, 141, 138, 108, 112, 116,
            113, 110, 109, 115, 114, 111, 162, 166, 170, 167, 164, 163, 169,
            168, 165, 198, 202, 206, 203, 200, 199, 205, 204, 201, 234, 238,
            242, 239, 236, 235, 241, 240, 237, 207, 211, 215, 212, 209, 208,
            214, 213, 210, 180, 184, 188, 185, 182, 181, 187, 186, 183, 171,
            175, 179, 176, 173, 172, 178, 177, 174, 225, 229, 233, 230, 227,
            226, 232, 231, 228, 216, 220, 224, 221, 218, 217, 223, 222, 219,
            189, 193, 197, 194, 191, 190, 196, 195, 192;
    
    for (int i=0; i<5; i++){deformation_measures::dot_2ot_5ot(i,SOT,FOT_m1,_r_m1[i]);}
    for (int i=0; i<5; i++){deformation_measures::dot_2ot_5ot(i,SOT,FOT_m2,_r_m2[i]);}
    
    bool tot_result = true;
    for (int i=0; i<5; i++){tot_result *= r_m1[i].isApprox(_r_m1[i]);}
    for (int i=0; i<5; i++){tot_result *= r_m2[i].isApprox(_r_m2[i]);}
    
    if (tot_result){
        results << "test_dot_2ot_5ot & True\\\\\n\\hline\n";
    }
    else {
        results << "test_dot_2ot_5ot & False\\\\\n\\hline\n";
    }
    
    return 1;
}

int test_two_sot_to_fot(std::ostream &results){
    /*!=============================
    |    test_two_sot_to_fot    |
    =============================
    
    Test mapping two second order tensors to a 
    fourth order tensor.
    */
    
    Matrix_9x9 C;  //The expected result
    Matrix_9x9 _C; //The function value
    
    //Initialize the generating tensors
    Matrix_3x3 A;
    Matrix_3x3 B;
    
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    B << 10, 11, 12, 13, 14, 15, 16, 17, 18;
    
    //Load the expected result
    C << 10,  14,  18,  15,  12,  11,  17,  16,  13,  50,  70,  90,  75,
         60,  55,  85,  80,  65,  90, 126, 162, 135, 108,  99, 153, 144,
        117,  60,  84, 108,  90,  72,  66, 102,  96,  78,  30,  42,  54,
         45,  36,  33,  51,  48,  39,  20,  28,  36,  30,  24,  22,  34,
         32,  26,  80, 112, 144, 120,  96,  88, 136, 128, 104,  70,  98,
        126, 105,  84,  77, 119, 112,  91,  40,  56,  72,  60,  48,  44,
         68,  64,  52;
    
    deformation_measures::two_sot_to_fot(A,B,_C);
    
    bool tot_result = C.isApprox(_C);
    
    if (tot_result){
        results << "test_two_sot_to_fot & True\\\\\n\\hline\n";
    }
    else {
        results << "test_two_sot_to_fot & False\\\\\n\\hline\n";
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
    test_voigt_3x3_symm_tensor(results);
    test_voigt_3x3_tensor(results);
    test_voigt_3x9_tensor(results);
    test_get_micro_strain(results);
    test_undo_voigt_3x3_tensor(results);
    test_undo_voigt_3x9_tensor(results);
    test_perform_left_positive_cyclic_permutation(results);
    test_perform_right_positive_cyclic_permutation(results);
    test_dot_2ot_4ot(results);
    test_dot_2ot_5ot(results);
    test_two_sot_to_fot(results);

    //Close the results file
    results.close();
}

