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
            std::vector< tensor::Tensor23 > node_phis;                                                             //!The nodal microdisplacements
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
            //!| Constitutive model parameters
            //!=
            
            std::vector< int > iparams;
            std::vector< double > fparams;
            
            //!=
            //!| Stresses
            //!=
            
            std::vector< tensor::Tensor23 >   PK2;   //!The second Piola-Kirchhoff stresses at the gauss points
            std::vector< tensor::Tensor23 >   SIGMA; //!The symmetric stress in the reference configuration at the gauss points
            std::vector< tensor::Tensor33 >   M;     //!The higher order stress at the gauss points
            
            //!==
            //!|
            //!| Constructors
            //!|
            //!==
        
            //Constructors
            Hex8();
            Hex8(std::vector< double >);
            Hex8(std::vector< double >, std::vector< double >, std::vector< double >, std::vector< double > _fparams = {}, std::vector< int > _iparams = {});
            
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
            
            void update_shape_function_values();
            
            void set_shape_functions();
            void set_shape_function(int);
            
            void set_local_gradient_shape_functions();
            void set_local_gradient_shape_function(int);
            
            void set_jacobian(bool);
            
            void set_global_gradient_shape_functions(bool);
            void set_global_gradient_shape_function(bool, int);
            
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
            
            void set_stresses(bool set_tangents = false);
            
            //!=
            //!| Residuals
            //!=
            
            //!|=> Forces
            
            void add_all_forces();
            void add_internal_nodal_forces();
            void add_extermal_nodal_forces();
            void add_kinematic_nodal_forces();
            
            //!|=> Moments
            
            void add_all_moments();
            void add_internal_moments();
            void add_external_moments();
            void add_kinematic_moments();
            
            //!=
            //!| Tangents
            //!=
            
            //!|=> Compute derivatives of fundamental deformation measures
            
            void set_fundamental_tangents();
            void set_dFdU();
            void set_dchidU();
            void set_dgradchi_dU();
            
            //!|=> Compute derivatives of micromorphic deformation measures
            
            void set_deformation_tangents();
            void set_dCdU();
            void set_dPsidU();
            void set_dGammadU();
            
            //!|=> Compute derivatives of stress measures
            
            void set_stress_tangents();
            void set_dPK2dU();
            void set_dSIGMAdU();
            void set_dMdU();
            
            //!|=> Balance of Linear Momentum
            
            void set_force_tangent();
            void add_dFintdU();
            void add_dFextdU();
            void add_dFkindU();
            
            //!|=> Balance of First Moment of Momentum
            
            void set_moment_tangent();
            void add_dMintdU();
            void add_dMextdU();
            void add_dMkindU();
            
            //!|=> Tangent utilities
            void reset_tangents();
            
            //!=
            //!| Element Integration 
            //!=
            
            void update_gauss_point(bool set_tangents = false);
            void integrate_element(bool set_tangents = false);
            
            //!=
            //!| Test functions
            //!=
            
            void set_gpt_num(int);
            double get_N(int);
            std::vector< double > get_dNdxi(int);
            tensor::Tensor23 get_jacobian(int);
            std::vector< double > get_dNdx(bool,int);
            double get_Jhatdet(bool);
            
            tensor::Tensor23 get_F();
            tensor::Tensor23 get_chi();
            tensor::Tensor33 get_grad_chi();
            
            tensor::Tensor23 get_C();
            tensor::Tensor23 get_Psi();
            tensor::Tensor33 get_Gamma();
            
            tensor::BaseTensor<3,288> get_dFdU();
            tensor::BaseTensor<3,288> get_dchidU();
            tensor::BaseTensor<9,288> get_dgrad_chidU();
            
            tensor::BaseTensor<3,288> get_dCdU();
            tensor::BaseTensor<3,288> get_dPsidU();
            tensor::BaseTensor<9,288> get_dGammadU();
            
            tensor::BaseTensor<3,288> get_dPK2dU();
            tensor::BaseTensor<3,288> get_dSIGMAdU();
            tensor::BaseTensor<9,288> get_dMdU();
            
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
            
            int gpt_num             = -1;                               //!The current gauss point number
            std::vector< double >                Ns;                    //!The shape function values
            std::vector< std::vector< double > > dNdxis;                //!The derivatives of the shape function with respect
                                                                        //!to the local coordinates.
            std::vector< std::vector< double > > dNdxs;                 //!The derivatives of the shape function with respect
                                                                        //!to the current configuration.
            std::vector< std::vector< double > > dNdXs;                 //!The derivatives of the shape function with respect 
                                                                        //!to the reference configuration.
            tensor::Tensor23 J = tensor::Tensor23({3,3});               //!The jacobian of transformation e.g. dxi dX
            tensor::Tensor23 Jinv = tensor::Tensor23({3,3});            //!The inverse of the jacobian of transformation e.g. (dxi dX)^(-1)
            double Jhatdet  = 0;                                        //!The determinant of the jacobian of transformation 
                                                                        //!to the reference configuration. e.g. det(dxi dX)
            
            tensor::Tensor23 F = tensor::Tensor23({3,3});               //!The deformation gradient
            tensor::Tensor23 chi = tensor::Tensor23({3,3});             //!The microdisplacement tensor
            tensor::Tensor33 grad_chi = tensor::Tensor33({3,3,3});      //!The gradient of the microdisplacement tensor
            
            tensor::Tensor23 C = tensor::Tensor23({3,3});               //!The Right Cauchy-Green deformation tensor
            tensor::Tensor23 Cinv = tensor::Tensor23({3,3});            //!The inverse of teh Right Cauchy-Green deformation tensor
            tensor::Tensor23 Psi = tensor::Tensor23({3,3});             //!The micro deformation tensor
            tensor::Tensor33 Gamma = tensor::Tensor33({3,3,3});         //!The higher order deformation tensor
            
            double Fdet = 0;                                            //!The determinant of the jacobian of the deformation gradient
            
            
            //!|=> Tangents
            
            //!Fundamental deformation measures
            
            bool                      dFdU_flag        = false;
            bool                      dchidU_flag      = false;
            bool                      dgrad_chidU_flag = false;
            
            tensor::BaseTensor<3,288> dFdU        = tensor::BaseTensor<3,288>({3,3,96});   //!The derivative of F w.r.t. the degree of freedom vector
            tensor::BaseTensor<3,288> dchidU      = tensor::BaseTensor<3,288>({3,3,96});   //!The derivative of chi w.r.t. the degree of freedom vector
            tensor::BaseTensor<9,288> dgrad_chidU = tensor::BaseTensor<9,288>({3,3,3,96}); //!The derivative of grad_chi w.r.t. the degree of freedom vector
            
            //!Derived deformation measures
            
            bool                      dCdU_flag     = false;
            bool                      dPsidU_flag   = false;
            bool                      dGammadU_flag = false;
            
            tensor::BaseTensor<3,288> dCdU        = tensor::BaseTensor<3,288>({3,3,96});   //!The derivative of C w.r.t. the degree of freedom vector
            tensor::BaseTensor<3,288> dPsidU      = tensor::BaseTensor<3,288>({3,3,96});   //!The derivative of Psi w.r.t. the degree of freedom vector
            tensor::BaseTensor<9,288> dGammadU    = tensor::BaseTensor<9,288>({3,3,3,96}); //!The derivative of Gamma w.r.t. the degree of freedom vector
            
            //!Stress derivatives
            
            bool                      dPK2dU_flag   = false;
            bool                      dSIGMAdU_flag = false;
            bool                      dMdU_flag     = false;
            
            tensor::BaseTensor<3,288> dPK2dU      = tensor::BaseTensor<3,288>({3,3,96});   //!The derivative of the second Piola-Kirchhoff stress w.r.t the degree of freedom vector
            tensor::BaseTensor<3,288> dSIGMAdU    = tensor::BaseTensor<3,288>({3,3,96});   //!The derivative of the symmetric stress tensor in the reference configuration w.r.t. the degree of freedom vector
            tensor::BaseTensor<9,288> dMdU        = tensor::BaseTensor<9,288>({3,3,3,96}); //!The derivative of the higher order stress tensor in the reference configuration w.r.t. the degree of freedom vector
            
            //!Stress tangents from constitutive model
            
            bool             stress_tangent_flag = false;
            
            tensor::Tensor43 dPK2dC       = tensor::Tensor43({3,3,3,3});     //!The derivative of the second Piola-Kirchhoff stress w.r.t. the right Cauchy-Green deformation tensor
            tensor::Tensor43 dSIGMAdC     = tensor::Tensor43({3,3,3,3});     //!The derivative of the symmetric stress w.r.t. the right Cauchy-Green deformation tensor
            tensor::Tensor53 dMdC         = tensor::Tensor53({3,3,3,3,3});   //!The derivative of the higher order stress w.r.t. the right Cauchy-Green deformation tensor
            
            tensor::Tensor43 dPK2dPsi     = tensor::Tensor43({3,3,3,3});     //!The derivative of the second Piola-Kirchhoff stress w.r.t. the micro-deformation tensor
            tensor::Tensor43 dSIGMAdPsi   = tensor::Tensor43({3,3,3,3});     //!The derivative of the symmetric stress w.r.t. the micro-deformation tensor
            tensor::Tensor53 dMdPsi       = tensor::Tensor53({3,3,3,3,3});   //!The derivative of the higher order stress w.r.t. the micro-deformation tensor
            
            tensor::Tensor53 dPK2dGamma   = tensor::Tensor53({3,3,3,3,3});   //!The derivative of the second Piola-Kirchhoff stress w.r.t. micro-gradient deformation tensor
            tensor::Tensor53 dSIGMAdGamma = tensor::Tensor53({3,3,3,3,3});   //!The derivative of the symmetric stress w.r.t. micro-gradient deformation tensor
            tensor::Tensor63 dMdGamma     = tensor::Tensor63({3,3,3,3,3,3}); //!The derivative of the higher order stress w.r.t. micro-gradient deformation tensor
            
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
    
    tensor::Tensor23 vector_dyadic_product(const std::vector< double > &, const std::vector< double >&); //!Dyadic product for second order tensors
}