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
#include<vector>
#include<fstream>
#include<finite_difference.h>
#include<deformation_measures.h>
#include<balance_equations.h>

void get_map_sot_voigt(double (&sot_voigt)[9][2]){
    /*!===========================
    |    get_map_sot_voigt    |
    ===========================
    
    Get the map of the second order tensor to voigt notation.
    
    */
    
    Matrix_3x3 A;
    A << 0,1,2,3,4,5,6,7,8;
    
    double indices[9][2] = {{0,0},{0,1},{0,2},
                            {1,0},{1,1},{1,2},
                            {2,0},{2,1},{2,2}};
    
    Vector_9 Avoigt;
    deformation_measures::voigt_3x3_tensor(A,Avoigt);
    
    for (int i=0; i<9; i++){
        sot_voigt[i][0] = indices[(int)(Avoigt[i]+0.5)][0];
        sot_voigt[i][1] = indices[(int)(Avoigt[i]+0.5)][1];
    }
    
    return;
}

void get_map_tot_voigt(double (&tot_voigt)[27][3]){
    /*!===========================
    |    get_map_tot_voigt    |
    ===========================
    
    Get the map of a third order tensor in the form A_ijk -> A_{i}{jk} to voigt notation.
    
    */
    
    double sot_voigt[9][2];
    get_map_sot_voigt(sot_voigt);
    
    for (int i=0; i<3; i++){
        
        for (int j=0; j<9; j++){
            
            tot_voigt[i*9+j][0] = i;
            tot_voigt[i*9+j][1] = sot_voigt[j][0];
            tot_voigt[i*9+j][2] = sot_voigt[j][1];            
        }
        
    }
    
    return;
}

void find_sot_index(const int &i, const int &j, int &Ihat){
    /*!========================
    |    find_sot_index    |
    ========================
    
    Map a pair of indices for a second order tensor to 
    the voigt index.
    
    */
    
    double sot_voigt[9][2];
    get_map_sot_voigt(sot_voigt);
    
    for (int indx=0; indx<9; indx++){
        if((sot_voigt[indx][0]==i) && (sot_voigt[indx][1]==j)){Ihat = indx; break;}
    }
    
    return;
}

void find_tot_index(const int &i, const int &j, const int &k, int &Ihat){
    /*!========================
    |    find_tot_index    |
    ========================
    
    Map a set of indices for a third order tensor to 
    the voigt index.
    
    */
    
    double tot_voigt[27][3];
    get_map_tot_voigt(tot_voigt);
    
    for (int indx=0; indx<27; indx++){
        if((tot_voigt[indx][0]==i) && (tot_voigt[indx][1]==j) && (tot_voigt[indx][2]==k)){Ihat = indx; break;}
    }
    
    return;
}

void define_N(double &N){
    /*!==================
    |    define_N    |
    ==================
    
    Define the value of the test function to be used.
    
    */
    
    N = 0.261;
    return;
}

void define_dNdx(double (&dNdx)[3]){
    /*!=====================
    |    define_dNdx    |
    =====================
    
    Define the gradient of the test function 
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

void define_cauchy(Vector_9 &cauchy){
    /*!=======================
    |    define_cauchy    |
    =======================
    
    Define the cauchy stress to be used.
    
    */
    
    cauchy[0] =  4.8247383;
    cauchy[1] =  0.69369043;
    cauchy[2] =  1.87181647;
    cauchy[3] =  3.98965695;
    cauchy[4] = -3.60138975;
    cauchy[5] = -4.47294891;
    cauchy[6] =  3.25726929;
    cauchy[7] = -0.81253737;
    cauchy[8] = -1.37661323;
    
    return;
}

void define_s(Vector_9 &s){
    /*!==================
    |    define_s    |
    ==================
    
    Define the symmetric stress to be used.
    
    */
    
    s[0] =  0.94496726;
    s[5] =  0.73268522;
    s[6] =  0.36190426;
    s[8] =  0.73268522;
    s[1] =  0.17534639;
    s[3] =  0.35118444;
    s[7] =  0.36190426;
    s[6] =  0.35118444;
    s[2] =  0.59818451;
    
    return;
}

void define_m(Vector_27 &m){
    /*!==================
    |    define_m    |
    ==================
    
    Define the higher order couple stress.
    
    */
    
    m[ 0] = 0.49291984245;
    m[ 1] = 0.463682849187;
    m[ 2] = 0.410503481618;
    m[ 3] = 0.790191457966;
    m[ 4] = 0.747691366301;
    m[ 5] = 0.415153945938;
    m[ 6] = 0.494002422078;
    m[ 7] = 0.917159660011;
    m[ 8] = 0.202756449146;
    m[ 9] = 0.842117576964;
    m[10] = 0.105880945737;
    m[11] = 0.345590957821;
    m[12] = 0.786366853468;
    m[13] = 0.343248849409;
    m[14] = 0.481767549881;
    m[15] = 0.898090172949;
    m[16] = 0.319073694596;
    m[17] = 0.414916901482;
    m[18] = 0.731387402298;
    m[19] = 0.388900418949;
    m[20] = 0.064710443124;
    m[21] = 0.327182261845;
    m[22] = 0.337523957493;
    m[23] = 0.950942676106;
    m[24] = 0.419581796783;
    m[25] = 0.164816718675;
    m[26] = 0.521393759903;
    
    return;
}

void define_body_force(double (&b)[3]){
    /*!===========================
    |    define_body_force    |
    ===========================
    
    Define the body force.
    
    */
    
    b[0] =  0.251;
    b[1] = -2.312;
    b[2] = 10.832;
    
    return;
}

void define_body_couple(double (&l)[9]){
    /*!============================
    |    define_body_couple    |
    ============================
    
    Define the body couple.
    
    */
    
    l[0] =  0.128;
    l[1] = -2.123;
    l[2] =  8.172;
    l[3] =  0.271;
    l[4] =  0.781;
    l[5] = -6.271;
    l[6] =  4.721;
    l[7] =  0.927;
    l[8] =  1.271;
    
    return;
}

void define_acceleration(double (&a)[3]){
    /*!=============================
    |    define_acceleration    |
    =============================
    
    Define the acceleration.
    
    */
    
    a[0] =  0.821;
    a[1] = -5.261;
    a[2] =  2.312;
    
    return;
}

void define_micro_gyration(double (&omega)[9]){
    /*!===============================
    |    define_micro_gyration    |
    ===============================
    
    Define the micro-gyration tensor.
    
    */
    
    omega[0] = 0.773449687115;
    omega[1] = 0.0971455366725;
    omega[2] = 0.0972789699315;
    omega[3] = 0.62718100678;
    omega[4] = 0.740238358734;
    omega[5] = 0.341171867764;
    omega[6] = 0.436016627415;
    omega[7] = 0.0326331393607;
    omega[8] = 0.229705402608;
    
    return;
}

void define_DcauchyDgrad_u(Matrix_9x9 &DcauchyDgrad_u){
    /*!===============================
    |    define_DcauchyDgrad_u    |
    ===============================
    
    Define the gradient of the cauchy stress w.r.t. the 
    gradient of u.
    
    */
    
    int initial = 1;
    
    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            DcauchyDgrad_u(i,j) = initial;
            initial++;
        }
    }
}

void define_DcauchyDphi(Matrix_9x9 &DcauchyDphi){
    /*!============================
    |    define_DcauchyDphi    |
    ============================
    
    Define the gradient of the cauchy stress w.r.t. the 
    micro-displacement tensor.
    
    */
    
    int initial = 27;
    
    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            DcauchyDphi(i,j) = initial;
            initial++;
        }
    }
}

void define_DcauchyDgrad_phi(Matrix_9x27 &DcauchyDgrad_phi){
    /*!=================================
    |    define_DcauchyDgrad_phi    |
    =================================
    
    Define the gradient of the cauchy stress w.r.t. the 
    micro-displacement tensor.
    
    */
    
    int initial = 2*27;
    
    for (int i=0; i<9; i++){
        for (int j=0; j<27; j++){
            DcauchyDgrad_phi(i,j) = initial;
            initial++;
        }
    }
}

void define_DsDgrad_u(Matrix_9x9 &DsDgrad_u){
    /*!==========================
    |    define_DsDgrad_u    |
    ==========================
    
    Define the gradient of the symmetric stress w.r.t. the 
    gradient of u.
    
    */
    
    int initial = 87;
    
    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            DsDgrad_u(i,j) = initial;
            initial++;
        }
    }
}

void define_DsDphi(Matrix_9x9 &DsDphi){
    /*!=======================
    |    define_DsDphi    |
    =======================
    
    Define the gradient of the symmetric stress w.r.t. the 
    micro-displacement tensor.
    
    */
    
    int initial = 13;
    
    for (int i=0; i<9; i++){
        for (int j=0; j<9; j++){
            DsDphi(i,j) = initial;
            initial++;
        }
    }
}

void define_DsDgrad_phi(Matrix_9x27 &DsDgrad_phi){
    /*!============================
    |    define_DsDgrad_phi    |
    ============================
    
    Define the gradient of the symmetric stress w.r.t. the 
    micro-displacement tensor.
    
    */
    
    int initial = 2*82;
    
    for (int i=0; i<9; i++){
        for (int j=0; j<27; j++){
            DsDgrad_phi(i,j) = initial;
            initial++;
        }
    }
}

void define_DmDgrad_u(Matrix_27x9 &DmDgrad_u){
    /*!==========================
    |    define_DsDgrad_u    |
    ==========================
    
    Define the gradient of the higher-order stress w.r.t. the 
    gradient of u.
    
    */
    
    int initial = 91;
    
    for (int i=0; i<27; i++){
        for (int j=0; j<9; j++){
            DmDgrad_u(i,j) = initial;
            initial++;
        }
    }
}

void define_DmDphi(Matrix_27x9 &DmDphi){
    /*!=======================
    |    define_DsDphi    |
    =======================
    
    Define the gradient of the higher-order stress w.r.t. the 
    micro-displacement tensor.
    
    */
    
    int initial = 32;
    
    for (int i=0; i<27; i++){
        for (int j=0; j<9; j++){
            DmDphi(i,j) = initial;
            initial++;
        }
    }
}

void define_DmDgrad_phi(Matrix_27x27 &DmDgrad_phi){
    /*!============================
    |    define_DmDgrad_phi    |
    ============================
    
    Define the gradient of the higher-order stress w.r.t. the 
    micro-displacement tensor.
    
    */
    
    int initial = 2*17;
    
    for (int i=0; i<27; i++){
        for (int j=0; j<27; j++){
            DmDgrad_phi(i,j) = initial;
            initial++;
        }
    }
}

std::vector<double> parse_grad_u_U(std::vector<double> U){
    /*!========================
    |    parse_grad_u_U    |
    ========================
    
    Parse the construction of the gradient of u w.r.t. U.
    in the reference configuration.
    */
    
    double detadx[3];
    define_detadx(detadx);
    
    Matrix_3x3 grad_u;
    
    grad_u(0,0) = U[0]*detadx[0];
    grad_u(1,1) = U[1]*detadx[1];
    grad_u(2,2) = U[2]*detadx[2];
    grad_u(1,2) = U[1]*detadx[2];
    grad_u(0,2) = U[0]*detadx[2];
    grad_u(0,1) = U[0]*detadx[1];
    grad_u(2,1) = U[2]*detadx[1];
    grad_u(2,0) = U[2]*detadx[0];
    grad_u(1,0) = U[1]*detadx[0];
    
    Vector_9 grad_u_voigt;
    deformation_measures::voigt_3x3_tensor(grad_u,grad_u_voigt);
    
    std::vector<double> result;
    result.resize(9);
    for (int i=0; i<9; i++){result[i] = grad_u_voigt(i);}
    return result;
}

std::vector<double> parse_grad_u_U_current(std::vector<double> U){
    /*!
    ================================
    |    parse_grad_u_U_current    |
    ================================
    
    Parse the construction of the gradient of u w.r.t. U.
    in the current configuration.
    */
    
    //Use detadx to define detadX
    double detadX[3];
    define_detadx(detadX);
    
    Matrix_3x3 dudX;
    
    dudX(0,0) = U[0]*detadX[0];
    dudX(1,1) = U[1]*detadX[1];
    dudX(2,2) = U[2]*detadX[2];
    dudX(1,2) = U[1]*detadX[2];
    dudX(0,2) = U[0]*detadX[2];
    dudX(0,1) = U[0]*detadX[1];
    dudX(2,1) = U[2]*detadX[1];
    dudX(2,0) = U[2]*detadX[0];
    dudX(1,0) = U[1]*detadX[0];
    
    Matrix_3x3 Finv;
    Finv = (Matrix_3x3::Identity() + dudX).inverse();
    
    Matrix_3x3 grad_u;
    grad_u = dudX*Finv;
    
    Vector_9 grad_u_voigt;
    deformation_measures::voigt_3x3_tensor(grad_u,grad_u_voigt);
    
    std::vector<double> result;
    result.resize(9);
    for (int i=0; i<9; i++){result[i] = grad_u_voigt(i);}
    return result;
}

std::vector<double> parse_phi_U(std::vector<double> U){
    /*!=====================
    |    parse_phi_U    |
    =====================
    
    Parse the construction of phi w.r.t. U.
    
    */
    
    double eta;
    define_eta(eta);
    
    Matrix_3x3 phi;
    
    phi(0,0) = U[ 3]*eta;
    phi(1,1) = U[ 4]*eta;
    phi(2,2) = U[ 5]*eta;
    phi(1,2) = U[ 6]*eta;
    phi(0,2) = U[ 7]*eta;
    phi(0,1) = U[ 8]*eta;
    phi(2,1) = U[ 9]*eta;
    phi(2,0) = U[10]*eta;
    phi(1,0) = U[11]*eta;
    
    Vector_9 phi_voigt;
    deformation_measures::voigt_3x3_tensor(phi,phi_voigt);
    
    std::vector<double> result;
    result.resize(9);
    for (int i=0; i<9; i++){result[i] = phi_voigt(i);}
    return result;
}

std::vector<double> parse_grad_phi_U(std::vector<double> U){
    /*!==========================
    |    parse_grad_phi_U    |
    ==========================
    
    Parse the construction of the gradient of phi w.r.t. U.
    
    */
    
    double detadx[3];
    define_detadx(detadx);
    
    int Jhat;
    
    Matrix_3x3 phi;
    Matrix_3x9 grad_phi;
    
    phi(0,0) = U[ 3];
    phi(1,1) = U[ 4];
    phi(2,2) = U[ 5];
    phi(1,2) = U[ 6];
    phi(0,2) = U[ 7];
    phi(0,1) = U[ 8];
    phi(2,1) = U[ 9];
    phi(2,0) = U[10];
    phi(1,0) = U[11];
    
    for (int i=0; i<3; i++){
        for (int I=0; I<3; I++){
            for (int j=0; j<3; j++){
                find_sot_index(I,j,Jhat);
                
                grad_phi(i,Jhat) = phi(i,I)*detadx[j];
                
            }
        }
    }
    
    Vector_27 grad_phi_voigt;
    deformation_measures::voigt_3x9_tensor(grad_phi,grad_phi_voigt);
    
    std::vector<double> result;
    result.resize(27);
    for (int i=0; i<27; i++){result[i] = grad_phi_voigt(i);}
    return result;
}

std::vector<double> parse_grad_phi_U_current(std::vector<double> U){
    /*!
    ==================================
    |    parse_grad_phi_U_current    |
    ==================================
    
    Parse the construction of the gradient of phi w.r.t. U in the 
    current configuration.
    
    */
    
    //Note we use define_detadX to define detadX
    double detadX[3];
    define_detadx(detadX);
    
    int Jhat;
    
    Matrix_3x3 phi;
    Matrix_3x9 dphidX;
    Matrix_3x9 grad_phi;
    
    phi(0,0) = U[ 3];
    phi(1,1) = U[ 4];
    phi(2,2) = U[ 5];
    phi(1,2) = U[ 6];
    phi(0,2) = U[ 7];
    phi(0,1) = U[ 8];
    phi(2,1) = U[ 9];
    phi(2,0) = U[10];
    phi(1,0) = U[11];
    
    for (int i=0; i<3; i++){
        for (int I=0; I<3; I++){
            for (int J=0; J<3; J++){
                find_sot_index(I,J,Jhat);
                
                dphidX(i,Jhat) = phi(i,I)*detadX[J];
                
            }
        }
    }
    
    //Construct the inverse deformation gradient
    Matrix_3x3 dudX;
    
    dudX(0,0) = U[0]*detadX[0];
    dudX(1,1) = U[1]*detadX[1];
    dudX(2,2) = U[2]*detadX[2];
    dudX(1,2) = U[1]*detadX[2];
    dudX(0,2) = U[0]*detadX[2];
    dudX(0,1) = U[0]*detadX[1];
    dudX(2,1) = U[2]*detadX[1];
    dudX(2,0) = U[2]*detadX[0];
    dudX(1,0) = U[1]*detadX[0];
    
    Matrix_3x3 Finv;
    Finv = (Matrix_3x3::Identity() + dudX).inverse();
    
    int Khat;
    
    for (int i=0; i<3; i++){
        for (int I=0; I<3; I++){
            for (int j=0; j<3; j++){
                find_sot_index(I,j,Jhat);
                
                grad_phi(i,Jhat) = 0.;
                
                for (int J=0; J<3; J++){
                    find_sot_index(I,J,Khat);
                    grad_phi(i,Jhat) += dphidX(i,Khat)*Finv(J,j);
                }
                
            }
        }
    }
    
    Vector_27 grad_phi_voigt;
    deformation_measures::voigt_3x9_tensor(grad_phi,grad_phi_voigt);
    
    std::vector<double> result;
    result.resize(27);
    for (int i=0; i<27; i++){result[i] = grad_phi_voigt(i);}
    return result;
}

int test_compute_internal_force(std::ofstream &results){
    /*!=====================================
    |    test_compute_internal_force    |
    =====================================
    
    Test the computation of the internal force.
    
    This test makes sure that the indices are lining up with 
    what is expected. To make sure that your model is also 
    correct you should test this in your code as well.
    
    */
    
    double  r[3]; //The expected result.
    double _r[3]; //The function output.
    
    double dNdx[3];
    define_dNdx(dNdx);
    
    Vector_9 cauchy;
    define_cauchy(cauchy);
    
    Matrix_3x3 cauchy_mat;
    deformation_measures::undo_voigt_3x3_tensor(cauchy,cauchy_mat);
    
    //Compute the expected value
    for (int i=0; i<3; i++){
        r[i] = 0.;
        for (int j=0; j<3; j++){
            r[i] -= dNdx[j]*cauchy_mat(j,i);
        }
    }
    
    bool tot_result = true;
    
    //Compute the vector result.
    balance_equations::compute_internal_force(dNdx, cauchy, _r);
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    //Zero out the result.
    for (int i=0; i<3; i++){_r[i] = 0.;}
    
    //Compute the component result.
    for (int i=0; i<3; i++){balance_equations::compute_internal_force(i, dNdx, cauchy, _r[i]);}
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    if (tot_result){
        results << "test_compute_internal_force & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_force & False\\\\\n\\hline\n";
    }
}

int test_compute_body_force(std::ofstream &results){
    /*!=================================
    |    test_compute_body_force    |
    =================================
    
    Test the computation of the body force.
    
    This test makes sure that the indices are lining up with 
    what is expected. To make sure that your model is also 
    correct you should test this in your code as well.
    
    */
    
    double  r[3]; //The expected result.
    double _r[3]; //The function output.
    
    double N;
    define_N(N);
    
    double density;
    define_density(density);
    
    double b[3];
    define_body_force(b);
    
    for (int i=0; i<3; i++){
        r[i] = N*density*b[i];
    }

    bool tot_result = true;
    
    //Test the vector form
    balance_equations::compute_body_force(N, density, b, _r);
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    //Zero out the result.
    for (int i=0; i<3; i++){_r[i] = 0.;}
    
    //Test the component form
    for (int i=0; i<3; i++){balance_equations::compute_body_force(i, N, density, b, _r[i]);}
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    if (tot_result){
        results << "test_compute_body_force & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_body_force & False\\\\\n\\hline\n";
    }
}

int test_compute_kinematic_force(std::ofstream &results){
    /*!======================================
    |    test_compute_kinematic_force    |
    ======================================
    
    Test the computation of the kinematic force.
    
    This test makes sure that the indices are lining up with 
    what is expected. To make sure that your model is also 
    correct you should test this in your code as well.
    
    */
    
    double  r[3]; //The expected result.
    double _r[3]; //The function output.
    
    double N;
    define_N(N);
    
    double density;
    define_density(density);
    
    double a[3];
    define_acceleration(a);
    
    for (int i=0; i<3; i++){
        r[i] = -N*density*a[i];
    }
    
    bool tot_result = true;
    
    //Test the vector form
    balance_equations::compute_kinematic_force(N, density, a, _r);
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}

    //Zero out the result.
    for (int i=0; i<3; i++){_r[i] = 0.;}
    
    //Test the component form
    for (int i=0; i<3; i++){balance_equations::compute_kinematic_force(i, N, density, a, _r[i]);}
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    
    if (tot_result){
        results << "test_compute_acceleration & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_acceleration & False\\\\\n\\hline\n";
    }
}

int test_compute_internal_couple(std::ofstream &results){
    /*!======================================
    |    test_compute_internal_couple    |
    ======================================
    
    Test the computation of the internal couple.
    
    This test makes sure that the indices are lining up with 
    what is expected. To make sure that your model is also 
    correct you should test this in your code as well.
    
    */
    
    double  r[9]; //The expected result.
    double _r[9]; //The function output.
    
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);
    
    Vector_9 cauchy;
    define_cauchy(cauchy);
    
    Vector_9 s;
    define_s(s);
    
    Vector_27 m;
    define_m(m);
    
    Matrix_3x3 cauchy_mat;
    deformation_measures::undo_voigt_3x3_tensor(cauchy,cauchy_mat);
    
    Matrix_3x3 s_mat;
    deformation_measures::undo_voigt_3x3_tensor(s,s_mat);
    
    Matrix_3x9 m_mat;
    deformation_measures::undo_voigt_3x9_tensor(m,m_mat);
    
    //Compute the divergence of the higher-order couple stress
    Matrix_3x3 divm_mat;
    Eigen::Matrix<double,1,3> dNdx_mat;
    for (int i=0; i<3; i++){dNdx_mat(0,i) = dNdx[i];}
    
    Vector_9 divm = dNdx_mat*m_mat;
    deformation_measures::undo_voigt_3x3_tensor(divm,divm_mat);
    
    //Compute the expected value
    Matrix_3x3 r_mat = -(N*(cauchy_mat - s_mat) - divm_mat.transpose());
    Vector_9 r_vec;
    deformation_measures::voigt_3x3_tensor(r_mat,r_vec);
    
    for (int i=0; i<9; i++){r[i] = r_vec[i];}
    
    bool tot_result = true;
    
    //std::cout << "r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << r[i] << " ";
    //}
    //std::cout << "\n";
    
    //Compute the vector result.
    balance_equations::compute_internal_couple(N, dNdx, cauchy, s, m, _r);
    for (int i=0; i<9; i++){tot_result *= 1e-6>fabs(r[i]-_r[i]);}
    
    //std::cout << "Vector _r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << _r[i] << " ";
    //}
    //std::cout << "\n";
    
    //Zero out the result.
    for (int i=0; i<3; i++){_r[i] = 0.;}
    
    //Compute the component result.
    int Ihat;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            find_sot_index(i,j,Ihat);
            balance_equations::compute_internal_couple(i, j, N, dNdx, cauchy, s, m, _r[Ihat]);
        }
    }
    
    //std::cout << "Component _r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << _r[i] << " ";
    //}
    //std::cout << "\n";
    
    for (int i=0; i<9; i++){tot_result *= 1e-6>fabs(r[i]-_r[i]);}
    
    if (tot_result){
        results << "test_compute_internal_couple & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_couple & False\\\\\n\\hline\n";
    }
}

int test_compute_body_couple(std::ofstream &results){
    /*!==================================
    |    test_compute_body_couple    |
    ==================================
    
    Test the computation of the body couple.
    
    This test makes sure that the indices are lining up with 
    what is expected. To make sure that your model is also 
    correct you should test this in your code as well.
    
    */
    
    double  r[9]; //The expected result.
    double _r[9]; //The function output.
    
    double N;
    define_N(N);
    
    double density;
    define_density(density);
    
    double l[9];
    define_body_couple(l);
    
    int Ihat;
    int Jhat;
    
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            find_sot_index(i,j,Ihat);
            find_sot_index(j,i,Jhat);
            r[Ihat] = N*density*l[Jhat];
        }
    }
    
    //std::cout << "Expected r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << r[i] << " ";
    //}
    //std::cout << "\n";
    

    bool tot_result = true;
    
    //Test the vector form
    balance_equations::compute_body_couple(N, density, l, _r);
    for (int i=0; i<9; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    //std::cout << "Vector _r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << _r[i] << " ";
    //}
    //std::cout << "\n";
    
    
    //Zero out the result.
    for (int i=0; i<9; i++){_r[i] = 0.;}
    
    //Test the component form
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            find_sot_index(i,j,Ihat);
            balance_equations::compute_body_couple(i, j, N, density, l, _r[Ihat]);
        }
    }
    for (int i=0; i<9; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    //std::cout << "Component _r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << _r[i] << " ";
    //}
    //std::cout << "\n";
    
    
    if (tot_result){
        results << "test_compute_body_couple & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_body_couple & False\\\\\n\\hline\n";
    }
}

int test_compute_kinematic_couple(std::ofstream &results){
    /*!=======================================
    |    test_compute_kinematic_couple    |
    =======================================
    
    Test the computation of the kinematic couple.
    
    This test makes sure that the indices are lining up with 
    what is expected. To make sure that your model is also 
    correct you should test this in your code as well.
    
    */
    
    double  r[9]; //The expected result.
    double _r[9]; //The function output.
    
    double N;
    define_N(N);
    
    double density;
    define_density(density);
    
    double omega[9];
    define_micro_gyration(omega);
    
    int Ihat;
    int Jhat;
    
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            find_sot_index(i,j,Ihat);
            find_sot_index(j,i,Jhat);
            r[Ihat] = -N*density*omega[Jhat];
        }
    }
    
    //std::cout << "Expected r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << r[i] << " ";
    //}
    //std::cout << "\n";
    

    bool tot_result = true;
    
    //Test the vector form
    balance_equations::compute_kinematic_couple(N, density, omega, _r);
    for (int i=0; i<9; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    //std::cout << "Vector _r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << _r[i] << " ";
    //}
    //std::cout << "\n";
    
    
    //Zero out the result.
    for (int i=0; i<9; i++){_r[i] = 0.;}
    
    //Test the component form
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            find_sot_index(i,j,Ihat);
            balance_equations::compute_kinematic_couple(i, j, N, density, omega, _r[Ihat]);
        }
    }
    for (int i=0; i<9; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    //std::cout << "Component _r: ";
    //for (int i=0; i<9; i++){
    //    std::cout << _r[i] << " ";
    //}
    //std::cout << "\n";
    
    
    if (tot_result){
        results << "test_compute_kinematic_couple & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_kinematic_couple & False\\\\\n\\hline\n";
    }
}

int test_compute_internal_force_jacobian(std::ofstream &results){
    /*!==============================================
    |    test_compute_internal_force_jacobian    |
    ==============================================
    
    Test the computation of the jacobian of the internal force.
    
    */
    
    double r[3][12]; //The expected result.
    Matrix_3x12 _r;  //The function output.
    
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);

    double eta;
    define_eta(eta);

    double detadx[3];
    define_detadx(detadx);
    
    double sot_map[9][2];
    get_map_sot_voigt(sot_map);
    
    double tot_map[27][3];
    get_map_tot_voigt(tot_map);
    
    Matrix_9x9  DcauchyDgrad_u;
    define_DcauchyDgrad_u(DcauchyDgrad_u);
    
    Matrix_9x9  DcauchyDphi;
    define_DcauchyDphi(DcauchyDphi);
    
    Matrix_9x27 DcauchyDgrad_phi;
    define_DcauchyDgrad_phi(DcauchyDgrad_phi);
    
    //Define the required Kronecker delta-like term for the first term
    Eigen::Matrix<double, 3, 12> K_eye1;
    K_eye1 = Eigen::Matrix<double, 3, 12>::Zero();
    K_eye1(0,0) = 1.;
    K_eye1(1,1) = 1.;
    K_eye1(2,2) = 1.;
    
    int Ihat;
    int Jhat;
    
    for (int j=0; j<3; j++){
        for (int K=0; K<12; K++){
            
            r[j][K] = 0;
            
            for (int i=0; i<3; i++){
                
                find_sot_index(i,j,Ihat);
                
                for (int I=0; I<3; I++){
                    
                    for (int l=0; l<3; l++){
                        
                        find_sot_index(I,l,Jhat);
                        
                        r[j][K] -= dNdx[i]*DcauchyDgrad_u(Ihat,Jhat)*K_eye1(I,K)*detadx[l];
                        
                    }
                    
                }
                
            }
            
        }
    }
    
    Matrix_3x3 eye;
    eye = Matrix_3x3::Identity();
    
    //Define the required Kronecker delta-like term for the second term
    Eigen::Matrix<double, 12, 12> K_eye2;
    K_eye2 = Eigen::Matrix<double, 12, 12>::Identity();
    K_eye2(0,0) = 0;
    K_eye2(1,1) = 0;
    K_eye2(2,2) = 0;
    
    for (int j=0; j<3; j++){
        for (int K=0; K<12; K++){
        
            for (int m=0; m<3; m++){
                for (int M=0; M<3; M++){
                
                    find_sot_index(m,M,Jhat);
                
                    for (int i=0; i<3; i++){
                    
                        find_sot_index(i,j,Ihat);
                    
                        if(K>=3){
                            r[j][K] -= dNdx[i]*DcauchyDphi(Ihat,Jhat)*eta*K_eye2(Jhat+3,K);
                        }
                        
                    }
                }
            }
        }
    }
    
    int idx;
    
    //Compute the third term
    for (int j=0; j<3; j++){
        for (int K=0; K<12; K++){
            for (int m=0; m<3; m++){
                for (int M=0; M<3; M++){
                    
                    find_sot_index(m,M,idx);
                    
                    for (int i=0; i<3; i++){
                        
                        find_sot_index(i,j,Ihat);
                        
                        for (int l=0; l<3; l++){
                            
                            find_tot_index(m,M,l,Jhat);
                    
                            if(K>=3){
                                r[j][K] -= dNdx[i]*DcauchyDgrad_phi(Ihat,Jhat)*detadx[l]*K_eye2(idx+3,K);
                            }
                        }
                    }
                }
            }
        }
    }
    
//    std::cout << "r:\n";
//    for (int i=0; i<3; i++){
//        for (int j=0; j<12; j++){
//            std::cout << r[i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//    
//    std::cout << "_r:\n";
//    for (int i=0; i<3; i++){
//        for (int j=0; j<12; j++){
//            std::cout << _r(i,j) << " ";
//        }
//        std::cout << "\n";
//    }
//    
//    std::cout << "error:\n";
//    for (int i=0; i<3; i++){
//        for (int j=0; j<12; j++){
//            std::cout << r[i][j]-_r(i,j) << " ";
//        }
//        std::cout << "\n";
//    }
    
    bool tot_result = true;
    
    //Compute the value in the vector function
    balance_equations::compute_internal_force_jacobian(N, dNdx, eta, detadx, DcauchyDgrad_u, DcauchyDphi, DcauchyDgrad_phi, _r);
    for (int i=0; i<3; i++){
        for (int j=0; j<12; j++){
            tot_result *= 1e-9>fabs(r[i][j]-_r(i,j));
        }
    }
    
    //std::cout << "Vector _r:" << _r << "\n";
    
    //Zero out the result matrix.
    _r = Matrix_3x12::Zero();
    
    //Compute the value for the component function
    for (int i=0; i<3; i++){
        
        for (int A=0; A<12; A++){
            balance_equations::compute_internal_force_jacobian(i, A, N, dNdx, eta, detadx, DcauchyDgrad_u, DcauchyDphi, DcauchyDgrad_phi, _r(i,A));
        }
    }
    
    for (int i=0; i<3; i++){
        for (int j=0; j<12; j++){
            tot_result *= 1e-9>fabs(r[i][j]-_r(i,j));
        }
    }
    
    //std::cout << "Component _r:" << _r << "\n";
    
    
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
    
    Test the computation of the jacobian associated with the 
    internal couple.
    
    */
    
    Matrix_9x12  r; //The expected result
    Matrix_9x12 _r; //The function output
    
    //Define the test function and it's gradient
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);

    //Define the interpolation function and it's gradient
    double eta;
    define_eta(N);

    double detadx[3];
    define_detadx(detadx);
    
    //Compute the mapping from the PDE degrees of freedom to the FE degrees of freedom
    SpMat  dgrad_udU(9,12);
    SpMat  dphidU(9,12);
    SpMat dgrad_phidU(27,12);
    
    balance_equations::construct_dgrad_udU(detadx,dgrad_udU);
    balance_equations::construct_dphidU(eta,dphidU);
    balance_equations::construct_dgrad_phidU(detadx,dgrad_phidU);
    
    //Get the stress jacobians w.r.t. the degrees of freedom
    Matrix_9x9  DcauchyDgrad_u;
    Matrix_9x9  DcauchyDphi;
    Matrix_9x27 DcauchyDgrad_phi;
    
    Matrix_9x9  DsDgrad_u;
    Matrix_9x9  DsDphi;
    Matrix_9x27 DsDgrad_phi;
    
    Matrix_27x9  DmDgrad_u;
    Matrix_27x9  DmDphi;
    Matrix_27x27 DmDgrad_phi;
    
    define_DcauchyDgrad_u(DcauchyDgrad_u);
    define_DcauchyDphi(DcauchyDphi);
    define_DcauchyDgrad_phi(DcauchyDgrad_phi);
    
    define_DsDgrad_u(DsDgrad_u);
    define_DsDphi(DsDphi);
    define_DsDgrad_phi(DsDgrad_phi);
    
    define_DmDgrad_u(DmDgrad_u);
    define_DmDphi(DmDphi);
    define_DmDgrad_phi(DmDgrad_phi);
    
    //Construct the jacobians of the stress tensors
    Matrix_9x12  DcauchyDU;
    Matrix_9x12  DsDU;
    Matrix_27x12 DmDU;
    
    DcauchyDU = DcauchyDgrad_u*dgrad_udU + DcauchyDphi*dphidU + DcauchyDgrad_phi*dgrad_phidU;
    DsDU      = DsDgrad_u*dgrad_udU      + DsDphi*dphidU      + DsDgrad_phi*dgrad_phidU;
    DmDU      = DmDgrad_u*dgrad_udU      + DmDphi*dphidU      + DmDgrad_phi*dgrad_phidU;
    
    r = -N*(DcauchyDU - DsDU);
    
    int Ihat;
    int Jhat;
    
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            
            find_sot_index(i,j,Ihat);
            
            for (int alpha=0; alpha<12; alpha++){
                
                for (int k=0; k<3; k++){
                    find_tot_index(k,j,i,Jhat);
                    r(Ihat,alpha) += dNdx[k]*DmDU(Jhat,alpha);
                }
                
            }
        }
    }
    
    //Test vector computation
    balance_equations::compute_internal_couple_jacobian(N, dNdx, eta, detadx,
                                                        DcauchyDgrad_u, DcauchyDphi, DcauchyDgrad_phi,
                                                        DsDgrad_u,      DsDphi,      DsDgrad_phi,
                                                        DmDgrad_u,      DmDphi,      DmDgrad_phi,
                                                        _r);
    bool tot_result = r.isApprox(_r,1e-6);
    
    //std::cout << "Expected  r:\n" <<  r << "\n";
    //std::cout << "Vector _r:\n" << _r << "\n";
    
    //Test component computation
    _r = Matrix_9x12::Zero();
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            for (int A=0; A<12; A++){
                find_sot_index(i,j,Ihat);
                balance_equations::compute_internal_couple_jacobian(i, j, A, N, dNdx, eta, detadx,
                                                                    DcauchyDgrad_u, DcauchyDphi, DcauchyDgrad_phi,
                                                                    DsDgrad_u,      DsDphi,      DsDgrad_phi,
                                                                    DmDgrad_u,      DmDphi,      DmDgrad_phi,
                                                                    _r(Ihat,A));
            }
        }
    }
    
    tot_result *= r.isApprox(_r,1e-6);
    
    //std::cout << "Component _r:\n" << _r << "\n";
    
    
    if (tot_result){
        results << "test_compute_internal_couple_jacobian & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_couple_jacobian & False\\\\\n\\hline\n";
    }
}

int test_construct_dgrad_udU(std::ofstream &results){
    /*!==================================
    |    test_construct_dgrad_udU    |
    ==================================
    
    Test the construction of the derivative of the gradient of u w.r.t. 
    the DOF vector in the reference configuration.
    
    */
    
    const int m = 9;
    const int n = 12;
    
    Matrix_9x12  r; //The expected output
    Matrix_9x12 _r; //The function result
    SpMat       _rsparse(9,12);
    
    double detadx[3];
    define_detadx(detadx);
    
    std::vector<double> U = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    //Compute the finite difference
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_grad_u_U, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            r(i,j) = r_vec[j][i];
        }
    }
    
    balance_equations::construct_dgrad_udU(detadx,_rsparse);
    _r = _rsparse;
    
    //std::cout << " r:\n" << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);

    balance_equations::construct_dgrad_udU(detadx, _r);
    tot_result *= r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_construct_dgrad_udU & True\\\\\n\\hline\n";
    }
    else {
        results << "test_construct_dgrad_udU & False\\\\\n\\hline\n";
    }
}

int test_construct_dgrad_udU_current(std::ofstream &results){
    /*!
    ==========================================
    |    test_construct_dgrad_udU_current    |
    ==========================================
    
    Test the construction of the derivative of the gradient of u w.r.t. 
    the DOF vector in the current configuration.
    
    */
    
    const int m = 9;
    const int n = 12;
    
    Matrix_9x12  r; //The expected output
    Matrix_9x12 _r; //The function result
    SpMat       _rsparse(9,12);
    
    //Use define_detadx to define detadX
    double detadX[3];
    define_detadx(detadX);
    
    std::vector<double> U = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    //Compute the finite difference
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_grad_u_U_current, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            r(i,j) = r_vec[j][i];
        }
    }
    
    //Form grad_u
    Matrix_3x3 dudX;
    for (int i=0; i<3; i++){
        for (int I=0; I<3; I++){
            dudX(i,I) = detadX[I]*U[i];
        }
    }
    
    //Form the inverse of the deformation gradient
    Matrix_3x3 Finv;
    Finv = (Matrix_3x3::Identity() + dudX).inverse();
    
    //Extract u
    double u[3];
    for (int i=0; i<3; i++){u[i] = U[i];}
    
    //Compute grad_u
    double eye[3][3] = {{1,0,0},
                        {0,1,0},
                        {0,0,1}};
    double grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            grad_u[i][j] = eye[i][j] - Finv(i,j);
        }
    }
    
    balance_equations::construct_dgrad_udU(u,grad_u,detadX,_r);
    
    //std::cout << " r:\n" << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_construct_dgrad_udU_current & True\\\\\n\\hline\n";
    }
    else {
        results << "test_construct_dgrad_udU_current & False\\\\\n\\hline\n";
    }
} 
    
int test_construct_dphidU(std::ofstream &results){
    /*!===============================
    |    test_construct_dphidU    |
    ===============================
    
    Test the construction of the derivative of phi w.r.t. 
    the DOF vector.
    
    */
    
    const int m = 9;
    const int n = 12;
    
    Matrix_9x12  r; //The expected output
    Matrix_9x12 _r; //The function result
    SpMat       _rsparse(9,12);
    
    double eta;
    define_eta(eta);
    
    std::vector<double> U = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    //Compute the finite difference
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_phi_U, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            r(i,j) = r_vec[j][i];
        }
    }
    
    balance_equations::construct_dphidU(eta,_rsparse);
    _r = _rsparse;
    
    //std::cout << " r:\n" << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);

    balance_equations::construct_dphidU(eta,_r);
    tot_result *= r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_construct_dphidU & True\\\\\n\\hline\n";
    }
    else {
        results << "test_construct_dphidU & False\\\\\n\\hline\n";
    }
}
    
int test_construct_dgrad_phidU(std::ofstream &results){
    /*!====================================
    |    test_construct_dgrad_phidU    |
    ====================================
    
    Test the construction of the derivative of the gradient of phi w.r.t. 
    the DOF vector in the reference configuration.
    
    */
    
    const int m = 27;
    const int n = 12;
    
    Matrix_27x12  r; //The expected output
    Matrix_27x12 _r; //The function result
    SpMat       _rsparse(27,12);
    
    double detadx[3];
    define_detadx(detadx);
    
    std::vector<double> U = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    //Compute the finite difference
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_grad_phi_U, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            r(i,j) = r_vec[j][i];
        }
    }
    
    balance_equations::construct_dgrad_phidU(detadx,_rsparse);
    _r = _rsparse;
    
    //std::cout << " r:\n" << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);

    balance_equations::construct_dgrad_phidU(detadx,_r);

    tot_result *= r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_construct_dgrad_phidU & True\\\\\n\\hline\n";
    }
    else {
        results << "test_construct_dgrad_phidU & False\\\\\n\\hline\n";
    }
} 

int test_construct_dgrad_phidU_current(std::ofstream &results){
    /*!
    ============================================
    |    test_construct_dgrad_phidU_current    |
    ============================================
    
    Test the construction of the derivative of the gradient of phi w.r.t. 
    the DOF vector in the current configuration.
    
    */
    
    const int m = 27;
    const int n = 12;
    
    Matrix_27x12  r; //The expected output
    Matrix_27x12 _r; //The function result
    SpMat       _rsparse(27,12);
    
    //Use define_detadx to define detadX
    double detadX[3];
    define_detadx(detadX);
    
    std::vector<double> U = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    //Compute the finite difference
    finite_difference::FiniteDifference fd;
    fd = finite_difference::FiniteDifference(parse_grad_phi_U_current, 2, U , 1e-6);
    
    //Compute the numeric gradient
    std::vector<std::vector<double>> r_vec = fd.numeric_gradient();
    
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            r(i,j) = r_vec[j][i];
        }
    }
    
    //Form grad_u
    Matrix_3x3 dudX;
    for (int i=0; i<3; i++){
        for (int I=0; I<3; I++){
            dudX(i,I) = detadX[I]*U[i];
        }
    }
    
    //Form the inverse of the deformation gradient
    Matrix_3x3 Finv;
    Finv = (Matrix_3x3::Identity() + dudX).inverse();
    
    //Extract u
    double u[3];
    for (int i=0; i<3; i++){u[i] = U[i];}
    
    double phi[9];
    for (int i=0; i<9; i++){phi[i] = U[i+3];}
    
    //Compute grad_u
    double eye[3][3] = {{1,0,0},
                        {0,1,0},
                        {0,0,1}};
    double grad_u[3][3];
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            grad_u[i][j] = eye[i][j] - Finv(i,j);
        }
    }
    
    Matrix_9x12 dgrad_udU;
    balance_equations::construct_dgrad_udU(u,grad_u,detadX,dgrad_udU);
    balance_equations::construct_dgrad_phidU(phi,Finv,detadX,dgrad_udU,_r);
    
    //std::cout << " r:\n" << r << "\n";
    //std::cout << "_r:\n" << _r << "\n";
    
    bool tot_result = r.isApprox(_r,1e-6);
    
    if (tot_result){
        results << "test_construct_dgrad_phidU_current & True\\\\\n\\hline\n";
    }
    else {
        results << "test_construct_dgrad_phidU_current & False\\\\\n\\hline\n";
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
    
    //!Test the force and couple computations
    test_compute_internal_force(results);
    test_compute_body_force(results);
    test_compute_kinematic_force(results);
    
    test_compute_internal_couple(results);
    test_compute_body_couple(results);
    test_compute_kinematic_couple(results);
    
    test_compute_internal_force_jacobian(results);
    test_compute_internal_couple_jacobian(results);
    
    test_construct_dgrad_udU(results);
    test_construct_dgrad_udU_current(results);
    test_construct_dphidU(results);
    test_construct_dgrad_phidU(results);
    test_construct_dgrad_phidU_current(results);

    //Close the results file
    results.close();
}
