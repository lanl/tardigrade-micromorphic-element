/*!=======================================================
  |                                                     |
  |                finite_difference.h                  |
  |                                                     |
  -------------------------------------------------------
  | The header file for the definition of a finite      |
  | difference gradient calculator                      |
  -------------------------------------------------------
  | Notes: The functions defined here compute the       |
  |        finite difference derivatives and gradients  |
  |        of scalar and vector functions.              |
  =======================================================*/
  
#include <vector>
  
namespace finite_difference{
  
    typedef std::vector< double > (*incoming_function)(std::vector< double >);
  
    class FiniteDifference{
        /*!=================================
        |    class: FiniteDifference    |
        =================================
        
        Class which defines finite difference operations to 
        compute derivatives and gradients numerically.
        
        */
        
        public:
            int accuracy_order                      = 2;              //!The accuracy of the finite difference calculation
            std::vector< double > deltah_multiplier = {-1,0,1};       //!The multipliers used to give the step size in the 
                                                                      //!finite difference calculation.
            std::vector< double > coefficients      = {-0.5, 0, 0.5}; //!The coefficients for each term in the finite difference 
                                                                      //!calculation.
            incoming_function fxn;                                    //!The function on which the finite difference will be calculated
            std::vector< double > x0;                                 //!The vector about which to compute the derivative
            double h=1e-6;                                            //!The step size for the gradient
            
            //!=
            //!| Constructors
            //!=
            
            FiniteDifference(){}
            FiniteDifference(incoming_function _fxn, int _accuracy_order, const std::vector< double > &_x0 ,double _h){
                /*!==========================
                |    FiniteDifference    |
                ==========================
                
                The constructor for the finite difference 
                operator where the accuracy order is chosen.
                
                */
                
                accuracy_order    = _accuracy_order;
                fxn               = _fxn;
                x0                = _x0;
                h                 = _h;
                deltah_multiplier = deltah_multipliers[accuracy_order/2-1];
                coefficients      = first_derivative_coeffs[accuracy_order/2-1];
            }
            
            //!=
            //!| Methods
            //!=
            
            std::vector< double > finite_difference(std::vector< double > h_vec){
                /*!===================================
                |        finite_difference        |
                ===================================
                
                Comptute the finite difference derivative 
                of the vector valued function fxn in the 
                direction of h_vec;
                
                */
                
                double h_den = norm(h_vec);
                std::vector< double > result;   //!The returned value from the function
                std::vector< double > solution; //!The overall solution
                solution.resize(0);
                
                std::vector< double > xi(h_vec.size(),0.); //!The perturbed value of x0
                double Ci     = 0;                         //!The corresponding coefficient value
                double deltah = 0;                         //!The corresponding deltah multiplier value
                
                for(int i=0; i<deltah_multiplier.size(); i++){
                    Ci     = coefficients[i];            //Update the current coefficient value
                    
                    if(fabs(Ci)>C_tol){
                        
                        //Update xi
                        for(int j=0; j<h_vec.size(); j++){xi[j] = x0[j]+deltah_multiplier[i]*h_vec[j];}
                        
                        //Compute the value of the result
                        result = fxn(xi);
                                                
                        //Change the size of the solution to match the result if required
                        if(solution.size()==0){solution.resize(result.size()); solution = std::vector<double>(result.size(),0.);}
                        
                        //Add the contribution of the current perturbation to the solution
                        for(int j=0; j<result.size(); j++){solution[j] += Ci*result[j];}
                    }
                }
                
                //Divide the solution vector by the denomenator factor
                for(int i=0; i<solution.size(); i++){solution[i] = solution[i]/h_den;}
                
                return solution;
            }
            
            std::vector< std::vector< double > > numeric_gradient(){
                /*!==========================================
                |            numeric_gradient            |
                ==========================================
                
                Compute the numeric gradient of a function
                
                */
                
                //Initialize the gradient output
                std::vector< std::vector< double > > gradient;
                gradient.resize(x0.size());
                
                std::vector< double > h_vec(x0.size(),0);
                int indx = 0;
                
                //Compute the gradient
                for(int i=0; i<x0.size(); i++){
                    h_vec[i] = h; //Update the perturbation vector
                    if(i>0){
                        h_vec[i-1] = 0;
                    }
                    indx++;
                    
                    gradient[i] = finite_difference(h_vec);
                }
                
                return gradient;
            }
            
        private:
            std::vector< std::vector< double > > deltah_multipliers      = {{-1.,0.,1.},{-2.,-1.,0.,1.,2.},{-3.,-2.,-1.,0.,1.,2.,3.},{-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.}}; //!All deltah multipliers
            std::vector< std::vector< double > > first_derivative_coeffs = {{-0.5,0,0.5},{1./12,-2./3,0.,2./3,-1./12},
                                                                            {-1./60,3./20,-3./4,0.,3./4,-3./20,1./60},
                                                                            {-1./280,-4./105,-1./5,-4./5,0,4./5,-1./5,4./105,-1./280}}; //!All coefficients for first order derivatives
            double C_tol = 1e-10; //!The tolerance on the coefficients
            double norm(std::vector< double > V){
                /*!==============================
                |            norm            |
                ==============================
                
                Compute the norm of a vector
                
                */
                
                double val = 0;
                for(int i=0; i<V.size(); i++){val += V[i]*V[i];}
                return sqrt(val);
            }
    };
    
}