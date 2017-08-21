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
            
            std::vector< std::vector< double > > reference_coords = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
                                                                     {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //!The reference coordinates of the element's nodes
            std::vector< std::vector< double > > current_coords   = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
                                                                     {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //!The current coordinates of the element's nodes
            std::vector< std::vector< double > > local_coords     = {{-1,-1,-1},{ 1,-1,-1},{ 1, 1,-1},{-1, 1,-1},
                                                                     {-1,-1, 1},{ 1,-1, 1},{ 1, 1, 1},{-1, 1, 1}}; //!The local coordinates of the element's nodes
            std::vector< tensor::Tensor > node_phis;
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
            
            std::vector< std::vector< double > > dof_at_nodes;       //!The current value of the degree of freedom vector
                                                                     //!at each node. 
            std::vector< std::vector< double > > Delta_dof_at_nodes; //!The change in the degree of freedom vector between the current increment and the last
                                                                     //!at each node.
            
            //!=
            //!| Residual and Element Stiffness
            //!=
            
            std::vector< double > RHS;                     //!The, ``right hand side,'' of the linearized equation (i.e. the residual vector)
            std::vector< std::vector< double > > AMATRX;   //!The negative of the element stiffness matrix (i.e. -dRHSddof)
            
            //!=
            //!| Stresses
            //!=
            
            std::vector< tensor::Tensor >   PK2;   //!The second Piola-Kirchhoff stresses at the gauss points
            std::vector< tensor::Tensor >   SIGMA; //!The symmetric stress in the reference configuration at the gauss points
            std::vector< tensor::Tensor >   M;     //!The higher order stress at the gauss points
            
            //!==
            //!|
            //!| Constructors
            //!|
            //!==
        
            //Constructors
            Hex8();
            Hex8(std::vector< double >);
            Hex8(std::vector< double >, std::vector< double >, std::vector< double >);
            
            //!==
            //!|
            //!| Operators
            //!|
            //!==
            
            //Hex8& operator=(const Hex8&);
            
            //!==
            //!|
            //!| Methods
            //!|
            //!==
            
            //!=
            //!| Shape Functions
            //!=
            
            double                shape_function(int,const std::vector< double >&);
            std::vector< double > local_gradient_shape_function(int,const std::vector< double >&);
            std::vector< double > global_gradient_shape_function(bool, int, const std::vector< double >&);
            tensor::Tensor        compute_jacobian(bool, const std::vector< double >&);
            
            void set_shape_function_values();
            void set_shape_functions();
            void set_local_gradient_shape_functions();
            void set_global_gradient_shape_functions(bool);
            
            //!=
            //!| Fundamental Deformation Measures
            //!=
            
            void set_fundamental_measures();
            void compute_deformation_gradient();
            void compute_microdisplacement();
            void compute_gradient_microdisplacement();
            
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
            
            void set_stresses();
            
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
            
            //!=
            //!| Element Integration 
            //!=
            
            void update_gauss_point();
            void integrate_element();
            
        private:
            //!=
            //!| Common solution variables
            //!=
            //!  These are kept private because they are to be utilized 
            //!  during computations over the nodes. They should not be
            //!  used generally since it may not be known outside of the 
            //!  computational routines exactly which gauss point they 
            //!  are defined for.
            
            std::vector< int > tot_shape = {3,3,3};
            std::vector< int > sot_shape = {3,3};
            
            int gpt_num             = -1;                          //!The current gauss point number
            std::vector< double >                Ns;               //!The shape function values
            std::vector< std::vector< double > > dNdxis;           //!The derivatives of the shape function with respect
                                                                   //!to the local coordinates.
            std::vector< std::vector< double > > dNdxs;            //!The derivatives of the shape function with respect
                                                                   //!to the current configuration.
            std::vector< std::vector< double > > dNdXs;            //!The derivatives of the shape function with respect 
                                                                   //!to the reference configuration.
            tensor::Tensor J        = tensor::eye();               //!The jacobian of transformation e.g. dxi dX
            
            std::vector< double > Jhatdet  = {1,1,1,1,1,1,1,1};    //!The determinant of the jacobian of transformation 
                                                                   //!to the reference configuration. e.g. det(dxi dX)
            
            tensor::Tensor F        = tensor::eye();               //!The deformation gradient
            tensor::Tensor chi      = tensor::eye();               //!The microdisplacement tensor
            tensor::Tensor grad_chi = tensor::Tensor(tot_shape);   //!The gradient of the microdisplacement tensor
            
            tensor::Tensor C        = tensor::Tensor(sot_shape);   //!The Right Cauchy-Green deformation tensor
            tensor::Tensor Cinv     = tensor::Tensor(sot_shape);   //!The inverse of teh Right Cauchy-Green deformation tensor
            tensor::Tensor Psi      = tensor::Tensor(sot_shape);   //!The micro deformation tensor
            tensor::Tensor Gamma    = tensor::Tensor(tot_shape);   //!The higher order deformation tensor
            
            double Fdet             = 1.;                          //!The determinant of the jacobian of the deformation gradient
            
            //!=
            //!| Parse incoming vectors
            //!=
            
            std::vector< std::vector< double > > parse_incoming_vectors(int,const std::vector< double > &);
    };
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
    
    tensor::Tensor vector_dyadic_product(const std::vector< double > &, const std::vector< double >&);
}