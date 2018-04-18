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
            
            tot_voigt[i*3+j][0] = i;
            tot_voigt[i*3+j][1] = sot_voigt[j][0];
            tot_voigt[i*3+j][2] = sot_voigt[j][1];            
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
        if((tot_voigt[indx][0]==i) && (tot_voigt[indx][1]==j) && (tot_voigt[indx][2])){Ihat = indx; break;}
    }
    
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
    /*!============================
    |    define_DcauchyDphi    |
    ============================
    
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
    
    //Compute the function output.
    for (int i=0; i<3; i++){balance_equations::compute_internal_force(i, dNdx, cauchy, _r[i]);}
    
    bool tot_result = true;
    
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
    
    for (int i=0; i<3; i++){balance_equations::compute_body_force(i, N, density, b, _r[i]);}
    
    bool tot_result = true;
    
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
    
    for (int i=0; i<3; i++){balance_equations::compute_kinematic_force(i, N, density, a, _r[i]);}
    
    bool tot_result = true;
    
    for (int i=0; i<3; i++){tot_result *= 1e-9>fabs(r[i]-_r[i]);}
    
    if (tot_result){
        results << "test_compute_acceleration & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_acceleration & False\\\\\n\\hline\n";
    }
}

int test_compute_internal_force_jacobian(std::ofstream &results){
    /*!==============================================
    |    test_compute_internal_force_jacobian    |
    ==============================================
    
    Test the computation of the jacobian of the internal force.
    
    */
    
    double  r[3][12]; //The expected result.
    double _r[3][12]; //The function output.
    
    double N;
    define_N(N);
    
    double dNdx[3];
    define_dNdx(dNdx);
    
    double sot_map[9][2];
    get_map_sot_voigt(sot_map);
    
    Matrix_9x9  DcauchyDgrad_u;
    define_DcauchyDgrad_u(DcauchyDgrad_u);
    
    Matrix_9x9  DcauchyDphi;
    define_DcauchyDphi(DcauchyDphi);
    
    Matrix_9x27 DcauchyDgrad_phi;
    define_DcauchyDgrad_phi(DcauchyDgrad_phi);
    
    //Define the required Kronecker delta-like term
    Eigen::Matrix<double, 3, 12> K_eye;
    K_eye = Eigen::Matrix<double, 3, 12>::Zero();
    K_eye(0,0) = 1.;
    K_eye(1,1) = 1.;
    K_eye(2,2) = 1.;
    
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
                        
                        r[j][K] -= dNdx[i]*DcauchyDgrad_u(Ihat,Jhat)*K_eye(I,K)*dNdx[l];
                        
                    }
                    
                }
                
            }
            
        }
    }
    
    for (int i=0; i<3; i++){
        
        for (int K=0; K<12; K++){
    
            balance_equations::compute_internal_force_jacobian(i, K, N, dNdx, DcauchyDgrad_u, DcauchyDphi, DcauchyDgrad_phi, _r[i][K]);
            
        }
    }
    
    std::cout << "r:\n";
    for (int i=0; i<3; i++){
        for (int j=0; j<12; j++){
            std::cout << r[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "_r:\n";
    for (int i=0; i<3; i++){
        for (int j=0; j<12; j++){
            std::cout << _r[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    bool tot_result = true;
    
    for (int i=0; i<3; i++){
        for (int j=0; j<12; j++){
            tot_result *= 1e-9>fabs(r[i][j]-_r[i][j]);
        }
    }
    
    if (tot_result){
        results << "test_compute_internal_force_jacobian & True\\\\\n\\hline\n";
    }
    else {
        results << "test_compute_internal_force_jacobian & False\\\\\n\\hline\n";
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
    test_compute_internal_force(results);
    test_compute_body_force(results);
    test_compute_kinematic_force(results);
    
    test_compute_internal_force_jacobian(results);

    //Close the results file
    results.close();
}