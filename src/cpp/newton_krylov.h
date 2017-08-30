/*Collection of routines to perform a Newton-Krylov iteration*/

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>

namespace krylov{
    
    typedef std::vector< double > (*residual_function)(const std::vector<double>&);
    typedef std::vector< double > Vector;
    typedef std::vector< std::vector< double > > Matrix;

    class Solver{
        /*Class which defines a krylov solver*/
        
        public:
            unsigned int NKi  = 0;     //Overall Iteration number
            unsigned int imax = 40;    //Maximum number of Newton iterations
            unsigned int kmax = 10;    //Maximum number of GMRES iterations
            bool verbose=false;        //Output root finding messages
            bool inconv_flg=false;     //Inconvergence flag
            
            double NKtol = 1e-7;       //The Newton-Raphson tolerance
            double tol  = 1e-7;        //The iteration tolerance
            
            residual_function R_fxn;   //The residual function
            
            std::vector< double > R;   //The residual vector
            
            std::vector< double > u;    //Initial value
            double norm_u;              //Norm of u
            
            /*Constructors*/
            Solver();
            Solver(residual_function,std::vector<double>,
                   unsigned int im=10,unsigned int km=10,
                   bool V=false);
            
            /*Newton Solver*/
            void solve();
            void perform_iteration();
            
            /*Residual computation*/
            virtual std::vector<double> get_residual(std::vector<double> u){
                /*Compute the residual*/
                //std::cout << "nothing happened\n";
                return R_fxn(u);
            }
            
            /*Derivative Computation*/
            std::vector< double > residual_derivative(const std::vector< double >&);
            std::vector< double > get_utrial(const double&, const std::vector< double >&);
            double get_h(std::vector< double >);
            
            
            /*gmres functions*/
            std::vector< double > gmres();
            
            std::vector< std::vector< double > > setup_q();
            std::vector< std::vector< double > > setup_H();
            std::vector< double > set_b();         
            std::vector< double > solve_triangular(const unsigned int&, const std::vector< std::vector< double >>&, const std::vector< double >&);
            
            /*Vector functions*/
            double vector_norm(const std::vector< double >&);
            void normalize_vector(std::vector< double >&);
            void normalize_vector(std::vector< double >&,const double&);
            double inner_product(const std::vector< double >&,const std::vector< double >&);
            //Matrix ATA(const Matrix&);
            
            /*Debugging functions*/
            void print_matrix(const std::vector< std::vector< double > >&);
            void print_vector(const std::vector< double >&);
    };
}