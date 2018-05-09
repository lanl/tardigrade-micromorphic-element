#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include<deformation_measures.h>

typedef std::chrono::high_resolution_clock Clock;

int main(){
    /*
    ==============
    |    main    |
    ==============
    
    Tests different methods for computing tensor multiplications and reports the 
    relative speeds.
    */

    Matrix_27x27 dMdgrad_chi;
    double J = 1.2;
    Matrix_3x3 F;
    Matrix_3x3 chi;
    Matrix_27x27 dmdgrad_chi;

    dMdgrad_chi = Matrix_27x27::Random();
    F           = Matrix_3x3::Random();
    chi         = Matrix_3x3::Random();

    double _dMdgrad_chi[729];
    double _F[9];
    double _chi[9];
    double _dmdgrad_chi[729];

    std::vector<double> __dMdgrad_chi(729);
    std::vector<double> __F(9);
    std::vector<double> __chi(9);
    std::vector<double> __dmdgrad_chi(729);

    std::array<double,729> ___dMdgrad_chi;
    std::array<double,9>   ___F;
    std::array<double,9>   ___chi;
    std::array<double,729> ___dmdgrad_chi;

    for (int i=0; i<27; i++){
        for (int j=0; j<27; j++){
            _dMdgrad_chi[27*i + j]   = dMdgrad_chi(i,j);
            __dMdgrad_chi[27*i + j]  = dMdgrad_chi(i,j);
            ___dMdgrad_chi[27*i + j] = dMdgrad_chi(i,j);
        }
    }

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            _F[3*i + j] = F(i,j);
            _chi[3*i + j] = chi(i,j);

            __F[3*i + j] = F(i,j);
            __chi[3*i + j] = chi(i,j);

            ___F[3*i + j] = F(i,j);
            ___chi[3*i + j] = chi(i,j);
        }
    }

    int n = 1;

    auto t0 = Clock::now();
    auto t1 = Clock::now();
    
    t0 = Clock::now();
    for (int _n=0; _n<n; _n++){
        deformation_measures::map_dAdgrad_chi_to_dadgrad_chi(dMdgrad_chi, J, F, chi, dmdgrad_chi);
    }
    t1 = Clock::now();
    std::cout << "Eigen Matrices: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/((double)n) << "\n";

    t0 = Clock::now();
    for (int _n=0; _n<n; _n++){
        deformation_measures::map_dAdgrad_chi_to_dadgrad_chi(_dMdgrad_chi, J, _F, _chi, _dmdgrad_chi);
    }
    t1 = Clock::now();
    std::cout << "1D Arrays: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/((double)n) << "\n";
    
    t0 = Clock::now();
    for (int _n=0; _n<n; _n++){
        deformation_measures::map_dAdgrad_chi_to_dadgrad_chi(__dMdgrad_chi, J, __F, __chi, __dmdgrad_chi);
    }
    t1 = Clock::now();
    std::cout << "1D std::vectors: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/((double)n) << "\n";
    
    t0 = Clock::now();
    for (int _n=0; _n<n; _n++){
        deformation_measures::map_dAdgrad_chi_to_dadgrad_chi(___dMdgrad_chi, J, ___F, ___chi, ___dmdgrad_chi);
    }
    t1 = Clock::now();
    std::cout << "1D std::arrays: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/((double)n) << "\n";

    return 1;
}
