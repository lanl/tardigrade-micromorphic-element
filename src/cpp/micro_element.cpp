/*!=======================================================
  |                                                     |
  |                 micro_element.cpp                   |
  |                                                     |
  -------------------------------------------------------
  | The source file for the definition of a             |
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
#include <micro_element.h>
  
namespace micro_element
{
    
    /*!==
    |
    | Constructors
    |
    ==*/
    
    Hex8::Hex8(){
        //Resize the RHS and AMATRX containers
        RHS.resize(96); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX.resize(96);
        for(int i=0; i<96; i++){
            AMATRX[i].resize(96);
        }
        
        //Resize the stress measure vectors
        PK2.resize(number_gauss_points);
        SIGMA.resize(number_gauss_points);
        M.resize(number_gauss_points);
    }
    
    Hex8::Hex8(std::vector< double > rcs){
        /*!====================
        |        Hex8       |
        =====================
        
        The constructor for a hexehedral element when 
        given the nodal reference coordinates.
        
        Input:
           rcs: A vector of doubles which are the
                coordinates of the nodes.
                
                The nodes are ordered in a counter clockwise 
                manner i.e.
                
               4,8      3,7
                o--------o
                |        |
                |        |
                |        |
                o--------o
               1,5      2,6
               
               where the comma indicates the ``upper layer'' of the
               hexehedral element.
        */
        
        //Resize the RHS and AMATRX containers
        RHS.resize(96); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX.resize(96);
        for(int i=0; i<96; i++){
            AMATRX[i].resize(96);
        }
        
        //Resize the phi vector
        std::vector< int >> shape = {3,3};
        node_phis.reshape(reference_coords.size());
        for(int n=0; n<reference_coords.size(), n++){
            node_phis[n] = tensor::Tensor(shape);
        }
        
        //Resize the stress measure vectors
        PK2.resize(number_gauss_points);
        SIGMA.resize(number_gauss_points);
        M.resize(number_gauss_points);
        
        
        
        reference_coords = parse_incoming_vector(1,rcs);
        current_coords   = parse_incoming_vector(1,rcs);
    }
    
    Hex8::Hex8(std::vector< double > rcs, std::vector< double > U, std::vector< double > dU){
        /*!====================
        |        Hex8       |
        =====================
        
        The constructor for a hexehedral element when 
        given the nodal reference coordinates.
        
        Input:
            rcs: A vector doubles which are the coordinates of 
                 the nodes.
                
                 The nodes are ordered in a counter clockwise 
                 manner i.e.
                
                4,8      3,7
                 o--------o
                 |        |
                 |        |
                 |        |
                 o--------o
                1,5      2,6
               
                where the comma indicates the ``upper layer'' of the
                hexehedral element.
               
            U:  The vector of changes of the degree of freedom 
                vector. It is organized such that the degrees of 
                freedom are in the order of the nodes.
        */
        
        //Resize the RHS and AMATRX containers
        RHS.resize(96); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX.resize(96);
        for(int i=0; i<96; i++){
            AMATRX[i].resize(96);
        }
        
        //Resize the phi vector
        std::vector< int >> shape = {3,3};
        node_phis.reshape(reference_coords.size());
        for(int n=0; n<reference_coords.size(), n++){
            node_phis[n] = tensor::Tensor(shape);
        }
        
        //Resize the stress measure vectors
        PK2.resize(number_gauss_points);
        SIGMA.resize(number_gauss_points);
        M.resize(number_gauss_points);
        
        //Break the rcs, U, and, dU vectors into a vector of vectors
        //where each subvector are the values associated with a given 
        //node.
        reference_coords   = parse_incoming_vector(1,rcs);
        dof_at_nodes       = parse_incoming_vector(2,U);
        Delta_dof_at_nodes = parse_incoming_vector(2,dU);
        
        //Assign the current value of the nodal coordinates
        for(int n=0; n<reference_coords.size(), n++){
            current_coords[n].resize(3);
            for(int i=0; i<3; i++){
                current_coords[n][i] = reference_coords[n][i]+dof_at_nodes[n][i];
            }
        }
        
        //Assign the values of phi at the nodes
        for(int n=0; n<reference_coords.size(), n++){
            //!NOTE: Assumes dof vector is set up as phi_11, phi_22, phi_33,
            //!                                      phi_23, phi_13, phi_12,
            //!                                      phi_32, phi_31, phi_21
            node_phis[n](0,0) = dof_at_nodes[n][0];
            node_phis[n](1,1) = dof_at_nodes[n][1];
            node_phis[n](2,2) = dof_at_nodes[n][2];
            node_phis[n](1,2) = dof_at_nodes[n][3];
            node_phis[n](0,2) = dof_at_nodes[n][4];
            node_phis[n](0,1) = dof_at_nodes[n][5];
            node_phis[n](2,1) = dof_at_nodes[n][6];
            node_phis[n](2,0) = dof_at_nodes[n][7];
            node_phis[n](1,0) = dof_at_nodes[n][8];
        }
    }
    
    //!==
    //!|
    //!| Operators
    //!|
    //!==
    
    Hex8& Hex8::operator=(const Hex8& hex8_in){
        /*!=======================
        |      operator=       |
        ========================
        
        Copy operator to allow for copying 
        hexahedral elements*/
        
        reference_coords = hex8_in.reference_coords;
        current_coords   = hex8_in.current_coords;
    }
    
    //!==
    //!|
    //!| Methods
    //!|
    //!==
    
    //!=
    //!| Shape Functions
    //!=
    
    Hex8::set_shape_function_values(){
        /*!===================================
        |    set_shape_function_values    |
        ===================================
        
        Set all of the shape function values to
        a value consistent with the current gauss 
        point (number gpt_num).
        
        */
        
        set_shape_functions();
        
    }
    
    double Hex8::shape_function(int n, const std::vector< double > &xi){
        /*!=========================
        |      shape_function    |
        ==========================
        
        Compute the value of the shape function
        at a given node and xi (local coordinate) 
        locations
        
        Input:
            n  : Node number
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        return 0.125*(1+xi[0]*local_coords[n][0])*(1+xi[1]*local_coords[n][1])*(1+xi[2]*local_coords[n][2]);
    }
    
    Hex8::set_shape_functions(){
        /*!===============================
        |      set_shape_functions    |
        ===============================
        
        Compute the value of the shape function
        at the current gauss point. Set the private 
        variable Ns with the results.
        
        Input:
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        for(n=0; n<reference_coords.size(); n++){
            Ns[n] = shape_function(n,xi);
        }
        return;
    }
    
    std::vector< double > Hex8::local_gradient_shape_function(int n, const std::vector< double > &xi){
        /*!=======================================
        |    local_gradient_shape_function    |
        =======================================
        
        Compute the value of the gradient of the 
        shape function at a given node and xi 
        (local coordinate) locations locally.
        
        Input:
            n  : Node number
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        //Initialize the shape of the output tensor
        std::vector< double > dNdxi;
        dNdxi.resize(3);
        
        dNdxi[0] = 0.125*local_coords[n][0]*(1+xi[1]*local_coords[n][1])*(1+xi[2]*local_coords[n][2]);
        dNdxi[1] = 0.125*(1+xi[0]*local_coords[n][0])*local_coords[n][1]*(1+xi[2]*local_coords[n][2]);
        dNdxi[2] = 0.125*(1+xi[0]*local_coords[n][0])*(1+xi[1]*local_coords[n][1])*local_coords[n][2];
        
        return dNdxi;
    }
    
    std::vector< double > Hex8::global_gradient_shape_function(bool mode, int n, const std::vector< double > &xi){
        /*!=======================================
        |    global_gradient_shape_function    |
        ========================================
        
        Compute the value of the gradient of the 
        shape function at a given node and xi 
        (local coordinate) locations globally.
        
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
            n   : Node number
            xi  : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        //Initialize the output
        std::vector< double > dNdx;
        dNdx.resize(3);
        
        tensor::Tensor J          = compute_jacobian(mode,xi);           //Compute the jacobian
        tensor::Tensor Jinv       = J.inverse();                         //Invert the jacobian
        double         Jdet       = J.det();                             //Get the determinant of the jacobian
        std::vector<double> dNdxi = local_gradient_shape_function(n,xi); //Compute the local gradient of the shape function
        
        //Compute the global gradient w.r.t. either the reference or global x
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                dNdx[i] += dNdxi[j]*Jinv(j,i);
            }
        }                       
        return dNdx;
    }
    
    tensor::Tensor Hex8::compute_jacobian(bool mode, const std::vector< double > & xi){
        /*!=========================
        |    compute_jacobian    |
        ==========================
        
        Compute the value of the jacobian of 
        transformation between the given 
        coordinates and the local coordinates.
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
            xi  : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        //Initialize the output
        std::vector< int > shape = {3,3};
        tensor::Tensor J = tensor::Tensor(shape);
        
        //Initialize local gradient value
        std::vector< double > local_gradient;
        local_gradient.resize(3);

        //Set up a temporary renaming variable depending on the mode
        std::vector< std::vector< double > > coordinates;
        
        if(mode){
            coordinates = current_coords;
        }
        else{
            coordinates = reference_coords;
        }
        
        //Add the contributions of each node in the element
        for(int n = 0; n<coordinates.size(); n++){
            local_gradient = local_gradient_shape_function(n,xi); //Get the local gradient of the given nodes shape function
            J += vector_dyadic_product(local_gradient,coordinates[n]); //Compute the vector dyadic product and add it to the jacobian
        }
        
        Jhatdet[gpt_num] = J.det();
        
        return J;
    }
    
    void set_local_gradient_shape_functions(){
        /*!================================================
        |      set_local_gradient_shape_functions      |
        ================================================
        
        Set all of the gradients of the shape functions.
        
        This sets the value of the gradients of the 
        shape functions at a given gauss point.
        
        */
        
        //Populate dNdxis
        for(int n=0; i<reference_coords.size(); n++){
            dNdxis[n] = local_gradient_shape_function(n, xi);
        }
        
        return;
    }
    
    void get_global_gradient_shape_functions(bool mode){
        /*!================================================
        |      set_global_gradient_shape_functions      |
        ================================================
        
        Get all of the gradients of the shape functions 
        with respect to global coordinates.
        
        This sets the value of the gradients of the 
        shape functions at the current gauss point.
        
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
        
        */
        
        //Populate dNdxs
        
        if(mode==0){
            for(int n=0; i<reference_coords.size(); n++){
                dNdXs[n] = global_gradient_shape_function(mode, n, xi);
            }
        }
        else if(mode==1){
            for(int n=0; i<reference_coords.size(); n++){
                dNdxs[n] = global_gradient_shape_function(mode, n, xi);
            }
        }
        else{
            std::cout "Error: Options are 0 and 1";
            assert(1==0);
        }
        return;
    }
    
    //!=
    //!| Fundamental Deformation Measures
    //!=
    
    void Hex8::set_fundamental_measures(){
        /*!========================================
        |       set_fundamental_measures       |
        ========================================
        
        Compute the fundamental deformation measures
        at the current gauss point. These are stored 
        in the private attributes.
        
        This function allows all of the measures 
        to be updated simultaneously and ensures 
        they are evaluated at a consistent location.
        
        Note that the current gauss point is defined 
        by gpt_num.
            
        */
        
        compute_deformation_gradient();
        compute_microdisplacement();
        compute_gradient_microdisplacement();
    }
    
    void Hex8::compute_deformation_gradient(){
        /*!======================================
        |    compute_deformation_gradient    |
        ======================================
        
        Compute the deformation gradient from the 
        reference and current coordinates at the 
        current gauss point defined by gpt_num.
        
        */
        
        //Initialize vectors
        std::vector< double > shape = {3,3};
        tensor::Tensor dxdxi = tensor::Tensor(shape); //!The derivative of the current coordinates w.r.t. the local coordinates
        tensor::Tensor dXdxi = tensor::Tensor(shape); //!The derivative of the reference coordinates w.r.t. the local coordinates
        
        //Compute the derivatives of the reference and current coordinates w.r.t. xi
        for(int n=0; n<reference_coords.size(); n++){
            dxdxi = vector_dyadic_product(current_coordinates[n],  dNdxis[n]);
            dXdxi = vector_dyadic_product(reference_coordinates[n],dNdxis[n]);
        }
        
        tensor::Tensor dxidX = dXdxi.inverse() //!The derivative of the local coordinates w.r.t. the reference coordinates
        
        //Reset F to zero
        F = tensor::Tensor(shape); //!The deformation gradient
        
        for(int i = 0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int k=0; k<3; k++){
                    F(i,j) += dxdxi(i,k)*dxidX(k,j);
                }
            }
        }
        
        //Compute the determinant of F to ensure they are consistent
        Fdet = F.det();
        
        return;
    }
    
    void Hex8::compute_microdisplacement(){
        /*!===================================
        |    compute_microdisplacement    |
        ===================================
        
        Compute the microdisplacement measure 
        chi at the current gauss point.
        
        */
        
        //Initialize vectors
        std::vector< double > shape = {3,3};
        
        //Reset chi to zero
        chi = tensor::Tensor(shape);
        
        //Interpolate the nodal phis to xi
        for(int n=0; n<reference_coords.size(); n++){
            chi = Ns[n]*node_phis[n];
        }
        chi += tensor::eye();
        return;
    }
    
    void Hex8::compute_gradient_microdisplacement(){
        /*!=============================================
        |    compute_gradient_microdisplacement    |
        ============================================
        
        Compute the gradient of the microdisplacement 
        tensor with respect to the reference coordinates
        at the gauss point.
        */
        
        //Initialize vectors
        std::vector< int > shape = {3,3};
        std::vector< int > tot_shape = {3,3,3};
        tensor::Tensor chi_n = tensor::Tensor(shape); //!The value of chi at a node
        tensor::Tensor I     = tensor::eye(); //!The second order identity tensor
        
        //Set grad_chi to zero
        grad_chi = tensor::Tensor(tot_shape);
        
        for(int n=0; n<reference_coords.size(); n++){
            chi_n = node_phi[n]+I:
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        grad_chi += chi_n(i,j)*dNdXs[n][k];
                    }
                }
            }
        }
        
        return;
    }
    
    //!=
    //!| Micromorphic Deformation Measures
    //!=
    
    void compute_right_cauchy_green(){
        /*!====================================
        |    compute_right_cauchy_green    |
        ====================================
        
        Compute the right Cauchy-Green deformation 
        tensor using the deformation gradient 
        stored in the private attributes.
        
        */
        
        //Zero the contents of C
        std::vector< int > shape = {3,3};
        C = tensor::Tensor(shape);
        
        //Form the right Cauchy-Green deformation tensor
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int i=0; i<3; i++){
                    C(I,J) += F(i,I)*F(i,J);
                }
            }
        }
        
        //Set the inverse of C so they are consistent
        Cinv = C.inverse();
    }
    
    void compute_Psi(){
        /*!=====================
        |    compute_Psi    |
        =====================
        
        Compute micro deformation measure Psi 
        using the deformation gradient and chi 
        stored in the private attributes.
        
        */
        
        //Zero the contents of Psi
        std::vector< int > shape = {3,3};
        Psi = tensor::Tensor(shape);
        
        //Form Psi
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int i=0; i<3; i++){
                    Psi(I,J) += F(i,I)*chi(i,J);
                }
            }
        }
    }
    
    void compute_Gamma(){
        /*!=======================
        |    compute_Gamma    |
        =======================
        
        Compute micro deformation measure Gamma 
        using the deformation gradient and the 
        gradient of chi stored in the private 
        attributes.
        
        */
        
        //Zero the contents of Gamma
        std::vector< int > shape = {3,3,3};
        Gamma = tensor::Tensor(shape);
        
        //Form Gamma
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    for(int i=0; i<3; i++){
                        Gamma(I,J,K) += F(i,I)*grad_chi(i,J,K);
                    }
                }
            }
        }
    }
    
    //!=
    //!| Constitutive Model Interface
    //!=
    
    void set_stresses(){
        /*!========================
        |     set_stresses     |
        ========================
        
        Set the stresses computed from the 
        constitutive model at the current 
        gauss point.
        
        */
        
        micro_material::get_stress(C,Psi,Gamma,PK2[gpt_num],SIGMA[gpt_num],M[gpt_num]);
    }
    
    //!=
    //!| Residuals
    //!=
    
    //!|=> Forces
    
    void Hex8::add_internal_nodal_forces(){
        /*!=====================================
        |     add_internal_nodal_forces     |
        =====================================
        
        Add the contributions of the internal 
        nodal forces to the right hand side vector.
        
        */
        std::vector< double > integral_value;
        //Put the force residuals in the RHS vector
        for(int n = 0; n<reference_coords.size; n++){
            
            for(int j=(n*3); j<((n+1)*3); j++){
                
                for(int I=0; I<3; I++){
                    for(int J=0; J<3; J++){
                        RHS[j] += -dNdXs[n][I]*PK2[n](I,J)*F(j,J)*Jhatdet[n]*weight[n];
                    }
                }
                
            }
        }
    }
    
    //!|=> Moments
    
    void Hex8::add_internal_moments(){
        /*!=====================================
        |       add_internal_moments        |
        =====================================
        
        Add the contributions of the internal 
        moments to the right hand side vector.
        
        */
        
        std::vector< double > shape = {3,3};
        tensor::Tensor mu_int = tensor::Tensor(shape);
        
        int initial_index = 3*reference_coord.size();
        //Create the moment residual at the different nodes
        for(int n=0; n<reference_coords.size(); n++){
            
            //Create the current value of mu_int
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            mu_int(i,j) += -Ns[n]*F(iI)*(SIGMA[n](IJ) - PK2[n](I,J))*F(j,J)*Jhatdet[n]*weight[n];
                            
                        }
                    }
                    
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            for(int K=0; K<3; K++){
                            mu_int(i,j) += -dNdXs[n][K]*F(jJ)*chi(iI)*M(KJI)*Jhatdet[n]*weight[n];
                            }
                        }
                    }
                }
            }
            
            //Populate the right hand side tensor
            RHS[initial_index+n*9+0] = mu_int(0,0);
            RHS[initial_index+n*9+1] = mu_int(1,1);
            RHS[initial_index+n*9+2] = mu_int(2,2);
            RHS[initial_index+n*9+3] = mu_int(1,2);
            RHS[initial_index+n*9+4] = mu_int(0,2);
            RHS[initial_index+n*9+5] = mu_int(0,1);
            RHS[initial_index+n*9+6] = mu_int(2,1);
            RHS[initial_index+n*9+7] = mu_int(2,0);
            RHS[initial_index+n*9+8] = mu_int(1,0);
        }
    }
    
    //!==
    //!|
    //!| Private Methods
    //!|
    //!=
            
    std::vector< std::vector< double > > parse_incoming_vectors(int mode; std::vector< double > incoming){
        /*!================================
        |    parse_incoming_vectors    |
        ================================
        
        Takes incoming vectors in different formats and parses them 
        into a vector of vectors form. This allows communication 
        between the driver program (Abaqus or otherwise) and the 
        element.
        
        */
        
        //!Variable defintions
        std::vector< double >::const_iterator first;        //! The iterator which identifies the start of a subvector
        std::vector< double >::const_iterator last;         //! The iterator which identifies the end of a subvector
        std::vector< std::vector< double > > parsed_vector; //! The parsed vector which is returned
        parsed_vector.resize(reference_coords.size());
        
        if(mode==1){//Parse an incoming coordinate vector
            factor = 3;
        if(mode==2){//Parse an incoming dof vector
            factor = 8;
        }
        
        //Extract the subvectors
        for(int n=0; n<reference_coords.size(); n++){
            first = incoming.begin()+n*factor;
            last  = incoming.begin()+(n+1)*factor;
            parsed_vector[n] = std::vector< double > new_vec(first,last);
        }
        
        return parsed_vector;
    }
    
    //!=
    //!| Parse DOF vectors
    //!=
            
            
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
    
    tensor::Tensor vector_dyadic_product(const std::vector< double > & V1, const std::vector< double >& V2){
        /*================================
        |    vector_dyadic_product    |
        ===============================
        
        Compute the dyadic product of two vectors.
        
        i.e. A_{ij} = V1_i V2_j
        
        Input:
            V1: Vector 1
            V2: Vector 2
        
        */
        
        //Compute the dimension of the resulting tensor
        int rows = V1.size();
        int cols = V2.size();
        
        //Initialize the tensor data and shape variables
        Eigen::MatrixXd data = Eigen::MatrixXd::Zero(rows,cols);
        std::vector< int > shape = {rows, cols};
        
        //Compute the dyadic product of the two vectors
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                data(i,j) = V1[i]*V2[j];
            }
        }
        
        tensor::Tensor T = tensor::Tensor(shape,data);
        return T;
    }
}