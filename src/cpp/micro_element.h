/*!=======================================================
  |                                                     |
  |                  micro_element.h                    |
  |                                                     |
  -------------------------------------------------------
  | The header file for the definition of a             |
  | micromorphic continuum element.                     |
  =======================================================
  | Dependencies:                                       |
  | tensor:      The class which defines tensor access  |
  |              to an underlying Eigen matrix. This    |
  |              may result in a somewhat slower result |
  |              however it should allow for a faster   |
  |              implementation.                        |
  =======================================================*/
  
#include <iostream>
#include <vector>
#include <tensor.h>
  
namespace micro_element
{
    class Hex8{
        /*!===
         |
         | H E X 8
         |
        ===
        
        An eight noded hexehedral element.
        
        The nodes have all degrees of freedom i.e. the 
        deformation u and the micro-deformation phi are 
        defined at all nodes.
        
        The element is defined such that different constitutive 
        equations can be implemented without extreme difficulty.
        
        */
        
        public:
        
            //!==
            //!|
            //!| Attribute Definitions
            //!|
            //!==
            
            //!=
            //!| Element definitions
            //!=
            
            std::vector< std::vector< double > > reference_coords = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1}
                                                                     {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //!The reference coordinates of the element's nodes
            std::vector< std::vector< double > > current_coords   = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1}
                                                                     {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //!The current coordinates of the element's nodes
            std::vector< std::vector< double > > local_coords     = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1}
                                                                     {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //!The local coordinates of the element's nodes
                                                          
            //!=
            //!| Gauss Quadrature
            //!=
            int number_gauss_points = 8; //!The number of gauss points
            
            std::vector< std::vector< double > > points  = {{-0.57735026918962573, -0.57735026918962573, -0.57735026918962573},
                                                            {-0.57735026918962573, -0.57735026918962573,  0.57735026918962573},
                                                            {-0.57735026918962573,  0.57735026918962573, -0.57735026918962573},
                                                            {-0.57735026918962573,  0.57735026918962573,  0.57735026918962573},
                                                            { 0.57735026918962573, -0.57735026918962573, -0.57735026918962573},
                                                            { 0.57735026918962573, -0.57735026918962573,  0.57735026918962573},
                                                            { 0.57735026918962573,  0.57735026918962573, -0.57735026918962573},
                                                            { 0.57735026918962573,  0.57735026918962573,  0.57735026918962573}}; //!The local coordinates of the gauss points
                                                           
            std::vector< double >                weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; //!The weights for the gauss quadrature
            
            //!=
            //!| Nodal degree of freedom variables
            //!=
            
            std::vector< double > dof;       //!The current value of the degree of freedom vector
            std::vector< double > Delta_dof; //!The change in the degree of freedom vector between the current increment and the last
            
            //!=
            //!| Residual and Element Stiffness
            //!=
            
            std::vector< double > RHS;      //!The, ``right hand side,'' of the linearized equation (i.e. the residual vector)
            RHS.resize(96); //Allocate the required memory for the right hand side vector
            std::vector< double > AMATRX;   //!The negative of the element stiffness matrix (i.e. -dRHSddof)
            
            //Allocate the required memory for the AMATRX
            AMATRX.resize(96);
            for(int i=0; i<96; i++){
                AMATRX[i].resize(96);
            }
            
            //!=
            //!| Stresses
            //!=
            
            std::vector<tensor::Tensor(std::vector<int> shape = {3,3})>   PK2;   //!The second Piola-Kirchhoff stresses at the gauss points
            std::vector<tensor::Tensor(std::vector<int> shape = {3,3})>   SIGMA; //!The symmetric stress in the reference configuration at the gauss points
            std::vector<tensor::Tensor(std::vector<int> shape = {3,3,3})> M;     //!The higher order stress at the gauss points
            
            //Resize the stress measure vectors
            PK2.resize(number_gauss_points);
            SIGMA.resize(number_gauss_points);
            M.resize(number_gauss_points);
            
            //!==
            //!|
            //!| Constructors
            //!|
            //!==
        
            //Constructors
            Hex8();
            Hex8(std::vector< std::vector< double > >);
            Hex8(std::vector< std::vector< double > >, std::vector< std::vector< double > >);
            
            //!==
            //!|
            //!| Operators
            //!|
            //!==
            
            Hex8& operator=(const Hex8&);
            
            //!==
            //!|
            //!| Methods
            //!|
            //!==
            
            //!=
            //!| Shape Functions
            //!=
            
            double shape_function(int,const std::vector< double >&);
            
            std::vector< double > local_gradient_shape_function(int,const std::vector< double >&);
            
            std::vector< double > global_gradient_shape_function(bool, int, const std::vector< double >&);
            
            tensor::Tensor & compute_jacobian(bool, const std::vector< double >&);
            
            //!=
            //!| Fundamental Deformation Measures
            //!=
            
            tensor::Tensor & compute_deformation_gradient();
            tensor::Tensor & compute_microdisplacement();
            tensor::Tensor & compute_gradient_microdisplacement();
            
            //!=
            //!| Micromorphic Deformation Measures
            //!=
            
            void set_deformation_measures();
            void compute_right_cauchy_green();
            void compute_Psi();            
            void compute_Gamma();
            
            //!=
            //!| Constitutive Model Interface
            //!=
            
            void set_stresses(int);
            
            //!=
            //!| Residuals
            //!=
            
            //!|=> Forces
            
            void add_internal_nodal_forces();
            void add_extermal_nodal_forces();
            void add_kinematic_nodal_forces();
            
            //!|=> Moments
            
            void add_internal_moments();
            void add_external_moments();
            void add_kinematic_moments();
            
            //!=
            //!| Tangents
            //!=
            
            //!|=> Balance of Linear Momentum
            
            void add_dFint_ddof();
            void add_dFext_ddof();
            void add_dFkin_ddof();
            
            //!|=> Balance of First Moment of Momentum
            
            void add_dMint_ddof();
            void add_dMext_ddof();
            void add_dMkin_ddof();
            
        private:
            //!=
            //!| Common solution variables
            //!=
            //!  These are kept private because they are to be utilized 
            //!  during computations over the nodes. They should not be
            //!  used generally since it may not be known outside of the 
            //!  computational routines exactly which gauss point they 
            //!  are defined for.
            
            tensor::Tensor F        = tensor::eye();                                      //!The deformation gradient
            tensor::Tensor chi      = tensor::eye();                                      //!The microdisplacement tensor
            tensor::Tensor grad_chi = tensor::Tensor(std::vector<int> shape = {3,3,3});   //!The gradient of the microdisplacement tensor
            
            tensor::Tensor C        = tensor::Tensor(std::vector<int> shape = {3,3});     //!The Right Cauchy-Green deformation tensor
            tensor::Tensor Cinv     = tensor::Tensor(std::vector<int> shape = {3,3});     //!The inverse of teh Right Cauchy-Green deformation tensor
            tensor::Tensor Psi      = tensor::Tensor(std::vector<int> shape = {3,3});     //!The micro deformation tensor
            tensor::Tensor Gamma    = tensor::Tensor(std::vector<int> shape = {3,3,3});   //!The higher order deformation tensor
            
            double detF         = 1.                                                      //!The determinant of the jacobian of the deformation gradient
    };
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
    
    vector_dyadic_product(const std::vector< double > &, const std::vector< double >&);
}