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
    
    std::vector< double > Hex8::get_shape_functions(const std::vector< double > &xi){
        /*!===============================
        |      get_shape_functions    |
        ===============================
        
        Compute the value of the shape function
        at a given node and xi (local coordinate) 
        locations
        
        Input:
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        std::vector< double > Ns;
        Ns.resize(reference_coords.size());
        
        for(n=0; n<reference_coords.size(); n++){
            Ns[n] = shape_function(n,xi);
        }
        return Ns;
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
        
        return J;
    }
    
    std::vector< std::vector< double > > get_local_gradient_shape_functions(const std::vector< double > &xi){
        /*!================================================
        |      get_local_gradient_shape_functions      |
        ================================================
        
        Get all of the gradients of the shape functions.
        
        This returns the value of the gradients of the 
        shape functions at a given value of xi.
        
        */
        
        std::vector< std::vector< double > > dNdxis; //!The derivative of the shape functions with respect to the local coordinates.
        dNdxis.resize(reference_coords.size());      //Resize the vector
        
        //Populate dNdxis
        for(int n=0; i<reference_coords.size(); n++){
            dNdxis[n] = local_gradient_shape_function(n, xi);
        }
        
        return dNdxis;
    }
    
    std::vector< std::vector< double > > get_global_gradient_shape_functions(bool mode, const std::vector< double > &xi){
        /*!================================================
        |      get_global_gradient_shape_functions      |
        ================================================
        
        Get all of the gradients of the shape functions 
        with respect to global coordinates.
        
        This returns the value of the gradients of the 
        shape functions at a given value of xi.
        
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
            xi  : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        std::vector< std::vector< double > > dNdXs; //!The derivative of the shape functions with respect to the local coordinates.
        dNdxis.resize(reference_coords.size());      //Resize the vector
        
        //Populate dNdxis
        for(int n=0; i<reference_coords.size(); n++){
            dNdxis[n] = global_gradient_shape_function(mode, n, xi);
        }
        
        return dNdxis;
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
        std::vector< std::vector< double > > dNdxis = get_local_gradient_shape_functions(local_coordinates[gpt_num]); //!The derivative of the shape functions with respect to the local coordinates.
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
        std::vector< double > Ns = get_shape_functions(local_coordinates[gpt_num]); //!The value of the shape functions at the gauss point.
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
        std::vector< std::vector< double > > dNdXs = get_gloal_gradient_shape_functions(local_coordinates[gpt_num]); //!The derivative of the shape functions with respect to the global coordinates.
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