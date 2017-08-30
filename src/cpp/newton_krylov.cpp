/*Collection of routines to perform a Newton-Krylov iteration*/

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <typeinfo>
#include <newton_krylov.h>

namespace krylov{
    
    Solver::Solver(){
        /*Blank initilizer*/
    }
    
    Solver::Solver(residual_function rin, std::vector< double > u_in,
                   unsigned int im, unsigned int km,
                   bool V){
        /*Initializer*/
        R_fxn = rin;
        imax = im;
        kmax = km;
        verbose = V;
        
        u      = u_in;
        norm_u = vector_norm(u);
        
        if(u.size()<kmax){
            std::cout << "\nWarning: Reducing gkmax to size of u\n";
            kmax = u.size();
        }
        
    }
    
    std::vector< double > Solver::get_utrial(const double &h, const std::vector<double>& s){
        /*Get the u value used for derivative computation*/
        std::vector< double > u_np1;
        u_np1.resize(u.size());
        for(int i=0; i<u.size(); i++){
            u_np1[i] = u[i]+h*s[i];
        }
        return u_np1;
    }
    
    std::vector< double > Solver::residual_derivative(const std::vector< double >& s){
        /*Compute the derivative of the residual*/
        double h = get_h(s);
        
        std::vector< double > u_np1 = get_utrial(h,s);
        
        std::vector< double > Dn = get_residual(u_np1);
        
        for(int i=0; i<u.size(); i++){
            Dn[i] = (Dn[i]-R[i])/h;
        }
        
        return Dn;
    }
    
    double Solver::get_h(std::vector< double > s){
        /*Get the step size*/
        double norm_s = vector_norm(s);
        double h = 0;
        
        if((norm_s>0) && (norm_u>0)){
            h = sqrt(tol)*norm_u/norm_s;
        }
        else if((norm_s>0) && (norm_u<tol)){
            h = sqrt(tol)*norm_s;
        }
        return h;
    }
    
    void Solver::solve(){
        /*Attempt to solve the problem*/
        
        /*Compute the initial residual*/
        
        R = get_residual(u);

        double Rnorm_0 = vector_norm(R);
        
        double Rnorm = Rnorm_0;
        
        double rel_tol = Rnorm/Rnorm_0;
        
        while((fabs(Rnorm/Rnorm_0)>NKtol) && (fabs(Rnorm)>1e-9) && (NKi<imax) && (!inconv_flg)){
            
            if(verbose){
                std::cout << "\n***************************\n";
                std::cout << "NK Iteration: "<<NKi+1<<"\n";
                std::cout << "u:\n";
                print_vector(u);
                std::cout << "R:\n";
                print_vector(R);
                std::cout << "Rnorm:          " << Rnorm << "\n";
                std::cout << "Relative Rnorm: " << Rnorm/Rnorm_0 << "\n";
                std::cout << "\n***************************\n";
            }
            
            //std::cout << "R:\n";
            //print_vector(R);
            
            perform_iteration();
            
            R = get_residual(u);
            
            //std::cout << "R:\n";
            //print_vector(R);
            
            Rnorm = vector_norm(R);
            
            NKi++;
        }
        
        if((NKi>=imax)||(inconv_flg)){
            std::cout << "Error: Newton-Krylov failed to converge\n";
            inconv_flg = true;
        }
        else{
            std::cout << "***********************\n";
            std::cout << "Newton-Krylov Converged\nRelative Norm: "<< Rnorm/Rnorm_0 << "\nNorm: " << Rnorm << "\n";
            std::cout << "***********************\n";
        }
        
        
    }
    
    void Solver::perform_iteration(){
        /*Perform an iteration of the krylov solver*/
        //std::cout << "R:\n";
        //print_vector(R);
        std::vector< double > s = gmres();
        //std::cout << "u:\n";
        //print_vector(u);
        for(int i=0; i<u.size(); i++){
            u[i] = u[i]+s[i];
        }
        //std::cout << "uip1:\n";
        //print_vector(u);
    }
    
    double Solver::vector_norm(const std::vector<double> &vec){
        /*Compute the norm of a vector*/
        double norm = 0;
        for(int i=0; i<vec.size(); i++){
            norm += pow(vec[i],2);
        }
        norm = sqrt(norm);
        return norm;
    }
    
    void Solver::normalize_vector(std::vector<double> &vec){
        /*Compute a normalized vector*/
        double norm = vector_norm(vec);
        normalize_vector(vec,norm);
    }
    
    void Solver::normalize_vector(std::vector<double> &vec, const double &norm){
        /*Compute a normalized vector*/
        for(int i=0; i<vec.size(); i++){
            vec[i] = vec[i]/norm;
        }
    }
    
    double Solver::inner_product(const std::vector<double> &vec1, const std::vector<double> &vec2){
        /*Compute the inner product*/
        double ip = 0.;
        
        //std::cout << "vec1:\n";
        //print_vector(vec1);
        
        //std::cout << "vec2:\n";
        //print_vector(vec2);
        
        for(int i=0; i<vec1.size(); i++){
            //std::cout << "vec[i]*vec2[i]: " << vec1[i]*vec2[i] << "\n";
            
            ip += vec1[i]*vec2[i];
        }
        
        //std::cout << "ip: " << ip << "\n";
        return ip;
    }
    
    std::vector< std::vector< double > > Solver::setup_q(){
        /*Setup the q data structure*/
        std::vector< std::vector< double > > q; //Normalized residual structure
        q.resize(kmax+1);
        for(int i=0; i<=kmax; i++){
            q[i].resize(u.size());
        }
        return q;
    }
    
    std::vector< std::vector< double > > Solver::setup_H(){
        /*Setup the H matrix*/
        std::vector< std::vector< double > > H;
        H.resize(kmax+1);
        for(int i=0; i<kmax+1; i++){
            H[i].resize(kmax);
        }
        //std::cout << "H0:\n";
        //print_matrix(H);
        return H;
    }
    
    std::vector< double > Solver::set_b(){
        /*Set the value of b*/
        
        std::vector< double > b;
        b.resize(R.size());
        for(int i=0; i<R.size(); i++){
            b[i] = -R[i];
        }
        return b;
    }
    
    std::vector<double> Solver::gmres(){
        /*Perform gmres iteration*/
        //std::cout << "In GMRES\n";
        
        /*Set-up required values*/
        std::vector< std::vector< double > > q;// = setup_q();
        q.reserve(kmax+1);
        
        std::vector<double> beta;
        beta.resize(kmax+1);
        
        //std::cout << "initial beta:\n";
        //print_vector(beta);
        
        std::vector< double > cj;
        cj.resize(kmax);
        
        std::vector< double > sj;
        sj.resize(kmax);
        
        //std::cout << "R:\n";
        //print_vector(R);
        
        //std::vector<double> b = set_b();
        q.push_back(set_b());
        
        //std::cout << "b:\n";
        //print_vector(b);
        
        double bnorm = vector_norm(q[0]);
        //std::cout << "bnorm: " << bnorm << "\n";
        
        normalize_vector(q[0],bnorm);
        
        //std::cout << "q.size(): " << q.size() << "\n";
        
        //std::cout << "q:\n";
        //print_matrix(q);
        
        beta[0] = bnorm;
        
        std::vector< std::vector< double > > H = setup_H();

        //std::cout << "beta:\n";
        //print_vector(beta);
        
        //std::cout << "kmax: " << kmax << "\n";
        
        unsigned int k = 0;
        double temp1 = 0;
        double temp2 = 0;
        //Begin loop
        while((k<kmax)&&(fabs(beta[k])>tol)){//&&(k<u.size())){
            
            q.push_back(residual_derivative(q[k]));
            
            //std::cout << "q[k+1]:";
            //print_vector(q[k+1]);
            
            for(int j=0; j<=k; j++){
                H[j][k] = inner_product(q[j],q[k+1]);
                for(int i=0; i<u.size(); i++){
                    q[k+1][i] -= H[j][k]*q[j][i];
                }
            }
            
            H[k+1][k] = vector_norm(q[k+1]);
            
            //std::cout << "H[k+1][k]: " << H[k+1][k] << "\n";
            
            //std::cout << "q[k+1]:\n";
            //print_vector(q[k+1]);
            
            normalize_vector(q[k+1],H[k+1][k]);
            
            //std::cout << "q[k+1]:\n";
            //print_vector(q[k+1]);
            
//            std::cout << "Pre rotation H:\n";
//            print_matrix(H);
            
                      
//            for(int j=0; j<=k; j++){
//                temp1 = sqrt(pow(H[j][j],2)+pow(H[j+1][j],2));
//                //std::cout << "temp1: " << temp1 << "\n";
//                cj[j] =   H[j][j]/temp1;
//                sj[j] = H[j+1][j]/temp1;
//            }

            /*Apply rotations to the current column*/
//            std::cout << "old rotations cj:\n";
//            print_vector(cj);
//            std::cout << "old rotations sj:\n";
//            print_vector(sj);
            
            for(int j=0; j<k; j++){
                temp1     =  cj[j]*H[j][k]+sj[j]*H[j+1][k];
                temp2     = -sj[j]*H[j][k]+cj[j]*H[j+1][k];
                H[j][k]   = temp1;
                H[j+1][k] = temp2;
            }
            
            /*Get new rotation*/
            temp1 = sqrt(pow(H[k][k],2)+pow(H[k+1][k],2));
            cj[k] =   H[k][k]/temp1;
            sj[k] = H[k+1][k]/temp1;
            
//            std::cout << "Pre new rotation H:\n";
//            print_matrix(H);
            
            for(int j=k; j<=k; j++){
                temp1     =  cj[j]*H[j][k]+sj[j]*H[j+1][k];
                temp2     = -sj[j]*H[j][k]+cj[j]*H[j+1][k];
                H[j][k]   = temp1;
                H[j+1][k] = temp2;
            }

  
//            std::cout << "Post new rotation H:\n";
//            print_matrix(H);

            temp1     =  cj[k]*beta[k]+sj[k]*beta[k+1];
            temp2     = -sj[k]*beta[k]+cj[k]*beta[k+1];
            beta[k]   = temp1;
            beta[k+1] = temp2;
            
//            for(int j=0; j<=k; j++){
//                temp1     =  cj[j]*beta[k]+sj[j]*beta[k+1];
//                temp2     = -sj[j]*beta[k]+cj[j]*beta[k+1];
//                beta[k]   = temp1;
//                beta[k+1] = temp2;
//            }
            
//            std::cout << "beta:\n";
//            print_vector(beta);
            
            k++;
            
        }
        
//        std::cout << "H:\n";
//        print_matrix(H);
//        
//        std::cout << "q:\n";
//        print_matrix(q);
//        
//        std::cout << "beta:\n";
//        print_vector(beta);
//        
//        std::cout << "k: " << k << "\n";
        
//        std::cout << "y:\n";
//        print_vector(y);
        
        std::vector< double > s;
        s.resize(u.size());
        
        if(k>0){
        
            std::vector< double > y = solve_triangular(k,H,beta);
        
            for(int j=0; j<k; j++){
                for(int i=0; i<u.size(); i++){
                    s[i] += y[j]*q[j][i];
                }
            }
        }
        else{
            std::cout << "************************************************\n";
            std::cout << "Error: No significant orthonormal vectors found.\n"
                      << "       Try increasing NKtol if appropriate.\n"
                      << "       Decreasing tol may also be helpful.\n"
                      << "       Ending iteration.\n";
            std::cout << "************************************************\n";
            inconv_flg = true;
        }
        
        //std::cout << "s:\n";
        //print_vector(s);
        
        //assert(1==0);
        
        return s;
  
    }  
    
    std::vector< double > Solver::solve_triangular(const unsigned int& kub, const std::vector< std::vector< double > > &H, const std::vector< double > &beta){
        /*Solve the triangular matrix*/
/*        std::vector <double> y;
        y.resize(kmax);
        
        
        int kmj = 0;
        int k = kmax-1;
        
        //std::cout << "k: " << k << "\n";
        
        y[k] = beta[k]/H[k][k];
        
        //std::cout << "y[k]: " << y[k] << "\n";
        
        for(int j=1; j<=k; j++){
            //std::cout << "j: " << j << "\n";
            y[k-j] = beta[k-j];
            for(int i=(k-j+1); i<=k; i++){
                //std::cout << "i: " << i << "\n";
                //std::cout << "H[k-j][i]: " << H[k-j][i] << "\n";
                //std::cout << "y[i]:      " << y[i] << "\n";
                y[k-j] -= H[k-j][i]*y[i];
            }
            y[k-j] = y[k-j]/H[k-j][k-j];
        }

        return y;*/
        std::vector <double> y;
        y.resize(kub);
        
        
        int kmj = 0;
        int k = kub-1;
        
        //std::cout << "k: " << k << "\n";
        
        y[k] = beta[k]/H[k][k];
        
        //std::cout << "y[k]: " << y[k] << "\n";
        
        for(int j=1; j<=k; j++){
            //std::cout << "j: " << j << "\n";
            y[k-j] = beta[k-j];
            for(int i=(k-j+1); i<=k; i++){
                //std::cout << "i: " << i << "\n";
                //std::cout << "H[k-j][i]: " << H[k-j][i] << "\n";
                //std::cout << "y[i]:      " << y[i] << "\n";
                y[k-j] -= H[k-j][i]*y[i];
            }
            y[k-j] = y[k-j]/H[k-j][k-j];
        }

        return y;
    }
    
    void Solver::print_matrix(const std::vector< std::vector< double > > &M){
        /*Print out the values of the matrix*/
        
        for(int i=0; i<M.size(); i++){
            for(int j=0; j<M[0].size(); j++){
                std::cout << "\t" << M[i][j];
            }
            std::cout << "\n";
            //assert(1==0);
        }
    }
    
    void Solver::print_vector(const std::vector< double > &V){
        /*Print out all of the values of a vector*/
        for(int i=0; i<V.size(); i++){
            std::cout << "\t" << V[i];
        }
        std::cout << "\n";
        //assert(1==0);
    }
}