/*Newton-Krylov test program*/

#include <iostream>
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

void run_solver1(){
    /*Set initial vector*/
    std::vector< double> u;
    u.resize(3);
    u[0] = 1.;
    u[1] = 1.;
    u[2] = 1.;
    
    krylov::Solver solver1 = krylov::Solver(residual_function1, u, 10,10);
    solver1.solve();
    
    std::cout << "Newton-Krylov Solution fxn 1:\n" << "x1: " << solver1.u[0] << "\nx2: " << solver1.u[1] << "\nx3: " << solver1.u[2];
}

void run_solver2(int mode=1){
    /*Set initial vector*/
    std::vector< double> u;
    u.resize(3);
    if(mode==1){
        u[0] = 1.;
        u[1] = 2.;
        u[2] = 3.;
    }
    else if(mode==2){
        u[0] = 0.;
        u[1] = 0.;
        u[2] = 0.;
    }
    krylov::Solver solver1 = krylov::Solver(residual_function2, u, 50, 10);
    solver1.solve();
    
    std::cout << "Newton-Krylov Solution fxn 2:\n" << "x1: " << solver1.u[0] << "\nx2: " << solver1.u[1] << "\nx3: " << solver1.u[2];
}

void run_solver3(){
    /*Run the third nonlinear problem*/
    std::vector< double > u;
    u.resize(1);
    u[0] = 0.5;
    krylov::Solver solver1 = krylov::Solver(residual_function3, u, 50, 10, true);
    solver1.solve();
    std::cout << "Newton-Krylov Solution fxn 3:\n" << "x: " << solver1.u[0];
}

int main(){
    /*Main function for Newton-Krylov test*/
    
    run_solver1();
    run_solver2(1);
    run_solver2(2);
    run_solver3();
}