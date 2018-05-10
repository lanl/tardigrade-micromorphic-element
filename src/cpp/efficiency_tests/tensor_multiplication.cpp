#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include<deformation_measures.h>

typedef std::chrono::high_resolution_clock Clock;
typedef std::vector<double> stdv;

int main(){
    /*
    ==============
    |    main    |
    ==============
    
    Tests different methods for computing tensor multiplications and reports the 
    relative speeds.
    */

    //Deformation measures
    double J = 1.2;
    Matrix_3x3 F;
    Matrix_3x3 chi;

    //Stress measures
    Vector_9  PK2_voigt;
    Vector_9  SIGMA_voigt;
    Vector_27 M_voigt;

    Vector_9  cauchy_voigt;
    Vector_9  s_voigt;
    Vector_27 m_voigt;

    //PK2 Jacobians
    Matrix_9x9  dPK2dF;
    Matrix_9x9  dPK2dchi;
    Matrix_9x27 dPK2dgrad_chi;

    //SIGMA Jacobians
    Matrix_9x9  dSIGMAdF;
    Matrix_9x9  dSIGMAdchi;
    Matrix_9x27 dSIGMAdgrad_chi;

    //Higher order stress Jacobians
    Matrix_27x9  dMdF;
    Matrix_27x9  dMdchi;
    Matrix_27x27 dMdgrad_chi;

    //Cauchy jacobians
    Matrix_9x9  dcauchydF;
    Matrix_9x9  dcauchydchi;
    Matrix_9x27 dcauchydgrad_chi;

    //s Jacobians
    Matrix_9x9  dsdF;
    Matrix_9x9  dsdchi;
    Matrix_9x27 dsdgrad_chi;

    //m Jacobians
    Matrix_27x9  dmdF;
    Matrix_27x9  dmdchi;
    Matrix_27x27 dmdgrad_chi;

    //Set random values
    F            = Matrix_3x3::Random();
    chi          = Matrix_3x3::Random();

    PK2_voigt    = Vector_9::Random();
    SIGMA_voigt  = Vector_9::Random();
    M_voigt      = Vector_27::Random();

    cauchy_voigt = Vector_9::Random();
    s_voigt      = Vector_9::Random();
    m_voigt      = Vector_27::Random();

    dPK2dF        = Matrix_9x9::Random();
    dPK2dchi      = Matrix_9x9::Random();
    dPK2dgrad_chi = Matrix_9x27::Random();

    dSIGMAdF        = Matrix_9x9::Random();
    dSIGMAdchi      = Matrix_9x9::Random();
    dSIGMAdgrad_chi = Matrix_9x27::Random();

    dMdF        = Matrix_27x9::Random();
    dMdchi      = Matrix_27x9::Random();
    dMdgrad_chi = Matrix_27x27::Random();

    int n = 10000;

    auto t0 = Clock::now();
    auto t1 = Clock::now();
    
    t0 = Clock::now();
    for (int _n=0; _n<n; _n++){
        deformation_measures::map_jacobians_to_current_configuration(F, chi,
            PK2_voigt,    SIGMA_voigt, M_voigt,
            cauchy_voigt, s_voigt,     m_voigt,
            dPK2dF,       dPK2dchi,    dPK2dgrad_chi,
            dSIGMAdF,     dSIGMAdchi,  dSIGMAdgrad_chi,
            dMdF,         dMdchi,      dMdgrad_chi,
            dcauchydF,    dcauchydchi, dcauchydgrad_chi,
            dsdF,         dsdchi,      dsdgrad_chi,
            dmdF,         dmdchi,      dmdgrad_chi);
    }
    t1 = Clock::now();
    std::cout << "Eigen Matrices: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/((double)n) << "\n";

    return 1;
}
