/*Newton-Krylov test program*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <newton_krylov.h>

std::vector<double> residual_function1(const std::vector<double> &x){
    /*Residual function for test nonlinear problem*/
    //std::cout << "\nIn residual calculation\n";
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    
    std::vector< double > r;
    r.resize(3);
    
    r[0] = 3.*x1-cos(x2*x3) - 1.5;
    r[1] = 4.*pow(x1,2.)-625.*pow(x2,2.)+2.*x3-1;
    r[2] = 20.*x3+exp(-x1*x2)+9.;
    
    //std::cout << "Residual:" << "\nr1: " << r[0] << "\nr2: " << r[1] << "\nr3: " << r[2] << "\n";
    
    return r;
}

std::vector<double> residual_function2(const std::vector<double> &x){
    /*Residual function for test nonlinear problem*/
    //std::cout << "\nIn residual calculation\n";
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    
    std::vector< double > r;
    r.resize(3);
    
    r[0] = pow(x1,2.)-2.*x1+pow(x2,2.)-x3+1.;
    r[1] = x1*pow(x2,2.)-x1-3.*x2+x2*x3+2.;
    r[2] = x1*pow(x3,2.)-3.*x3+x2*pow(x3,2.)+x1*x2;
    
    //std::cout << "Residual:" << "\nr1: " << r[0] << "\nr2: " << r[1] << "\nr3: " << r[2] << "\n";
    
    return r;
}

std::vector<double> residual_function3(const std::vector<double> &x){
    /*Residual function for test nonlinear problem*/
    double x1 = x[0];
    std::vector<double> r;
    r.resize(1);
    
    r[0] = x1*(x1-2.)*(x1+7.);
    return r;
}

void run_solver1(std::ofstream &results){
    /*Set initial vector*/
    std::vector< double> u;
    std::vector< double> r;
    u.resize(3);
    r.resize(3);
    bool test_result = true;
    u[0] = 1.;
    u[1] = 1.;
    u[2] = 1.;
    r[0] =  0.83328161;
    r[1] =  0.03533462;
    r[2] = -0.49854928;
    
    krylov::Solver solver1 = krylov::Solver(residual_function1, u, 10,10);
    solver1.solve();
    
    for(int i=0; i<3; i++){test_result *= 1e-6>fabs(solver1.u[i]-r[i]);}
    
    //std::cout << "Newton-Krylov Solution fxn 1:\n" << "x1: " << solver1.u[0] << "\nx2: " << solver1.u[1] << "\nx3: " << solver1.u[2];
    //std::cout << "\ntest result: " << test_result << "\n";;
    
    if(test_result){results << "test_newton_krylov_1 & True\\\\\n\\hline\n";}
    else{results << "\ntest_newton_krylov_1 & False\\\\\n\\hline\n";}
}

void run_solver2(std::ofstream &results,int mode=1){
    /*Set initial vector*/
    std::vector< double > u;
    std::vector< double > r;
    u.resize(3);
    r.resize(3);
    bool test_result = true;
    if(mode==1){
        u[0] = 1.;
        u[1] = 2.;
        u[2] = 3.;
        
        r[0] = 1.;
        r[1] = 1.;
        r[2] = 1.;
    }
    else if(mode==2){
        u[0] = 0.;
        u[1] = 0.;
        u[2] = 0.;
        
        r[0] = 1.09894258;
        r[1] = 0.36761668;
        r[2] = 0.14493166;
    }
    krylov::Solver solver1 = krylov::Solver(residual_function2, u, 50, 10);
    solver1.solve();
    
    std::cout << "Newton-Krylov Solution fxn 2:\n" << "x1: " << solver1.u[0] << "\nx2: " << solver1.u[1] << "\nx3: " << solver1.u[2];
    for(int i=0; i<3; i++){test_result *= 1e-6>fabs(solver1.u[i]-r[i]);}
    std::cout << "\ntest_result: " << test_result << "\n";
    
    if(mode==1){
        if(test_result){results << "test_newton_krylov_2_mode_1 & True\\\\\n\\hline\n";}
        else{results << "\ntest_newton_krylov_2_mode_1 & False\\\\\n\\hline\n";}
    }
    else if(mode==2){
        if(test_result){results << "test_newton_krylov_2_mode_2 & True\\\\\n\\hline\n";}
        else{results << "\ntest_newton_krylov_2_mode_2 & False\\\\\n\\hline\n";}
    }
}

void run_solver3(std::ofstream &results){
    /*Run the third nonlinear problem*/
    std::vector< double > u;
    double answer = 0.;
    u.resize(1);
    u[0] = 0.5;
    krylov::Solver solver1 = krylov::Solver(residual_function3, u, 50, 10, false);
    solver1.solve();
    
    if(1e-7>fabs(solver1.u[0]-answer)){results << "test_newton_krylov_3 & True\\\\\n\\hline\n";}
    else{results << "\ntest_newton_krylov_3 & False\\\\\n\\hline\n";}
}

int main(){
    /*Main function for Newton-Krylov test*/
    
    std::ofstream results;
    //Open the results file
    results.open ("results.tex");
    
    run_solver1(results);
    run_solver2(results,1);
    run_solver2(results,2);
    run_solver3(results);
    
    //Close the results file
    results.close();
}