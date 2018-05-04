#include <Eigen/Dense>
#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

//Global variables
const double sot_to_voigt[3][3] = {{0, 5, 4},
                                   {8, 1, 3},
                                   {7, 6, 2}};

void explicit_loop(const int n, const Eigen::Matrix3d &A, const Eigen::Matrix3d &B, const Eigen::Matrix<double, 9, 9> C, Eigen::Matrix<double, 9, 9> &D){
    /*
    =======================
    |    explicit_loop    |
    =======================
    
    Compute the contraction explicitly.
    
    */
    
    int I;
    int J;
    int M;

    for (int _n=0; _n<n; _n++){
    
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                I = sot_to_voigt[i][j];
            
                for (int k=0; k<3; k++){
                    for (int l=0; l<3; l++){
                        J = sot_to_voigt[k][l];
                    
                        D(I,J) = 0.;
                    
                        for (int m=0; m<3; m++){
                            for (int n=0; n<3; n++){
                                M = sot_to_voigt[m][n];
                            
                                D(I,J) += A(i,m)*C(M,J)*B(j,n);
                            }
                        }
                    }
                }
            }
        }
    }
}

int main(){
    /*
    ==============
    |    main    |
    ==============
    
    Tests different methods for computing tensor multiplications and reports the 
    relative speeds.
    */

    int n = 10000;
    
    Eigen::Matrix3d A;
    A = Eigen::Matrix3d::Random();
    
    Eigen::Matrix3d B;
    B = Eigen::Matrix3d::Random();
    
    Eigen::Matrix<double, 9, 9> C;
    C = Eigen::Matrix<double, 9, 9>::Random();
    
    Eigen::Matrix<double, 9, 9> D;
    
    auto t0 = Clock::now();
    auto t1 = Clock::now();
    
    t0 = Clock::now();
    explicit_loop(n,A,B,C,D);
    t1 = Clock::now();
    std::cout << "Explicit Loop: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/((double)n) << "\n";
    return 1;
}
