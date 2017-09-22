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
#include <micromorphic_linear_elasticity.h>
#include <ctime>
  
namespace micro_element
{
    
    /*!==
    |
    | Constructors
    |
    ==*/
    
    Hex8::Hex8(){
        //Resize the RHS and AMATRX containers
        RHS.resize(96,0.); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX.resize(96);
        for(int i=0; i<96; i++){
            AMATRX[i].resize(96,0.);
        }
        
        //Resize the phi vector
        node_phis.resize(reference_coords.size());
        for(int n=0; n<reference_coords.size(); n++){
            node_phis[n] = tensor::Tensor23({3,3});
        }
        
        //Resize the stress measure vectors
        PK2.resize(number_gauss_points);
        SIGMA.resize(number_gauss_points);
        M.resize(number_gauss_points);
        
        //Initialize stress measures to zero
        for(int i=0; i<number_gauss_points; i++){
            PK2[i]   = tensor::Tensor23({3,3});
            SIGMA[i] = tensor::Tensor23({3,3});
            M[i]     = tensor::Tensor33({3,3,3});
        }
        
        //Resize the private vector attributes
        Ns.resize(local_coords.size());
        dNdxis.resize(local_coords.size());
        dNdxs.resize(local_coords.size());
        dNdXs.resize(local_coords.size());
        
        //Set the dimension of the gradient vectors to 3
        for(int n=0; n<8; n++){dNdxis[n].resize(3);}
        for(int n=0; n<8; n++){dNdxs[n].resize(3);}
        for(int n=0; n<8; n++){dNdXs[n].resize(3);}
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
        RHS.resize(96,0.); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX.resize(96);
        for(int i=0; i<96; i++){
            AMATRX[i].resize(96,0.);
        }
        
        //Resize the phi vector
        node_phis.resize(reference_coords.size());
        for(int n=0; n<reference_coords.size(); n++){
            node_phis[n] = tensor::Tensor23({3,3});
        }
        
        //Resize the stress measure vectors
        PK2.resize(number_gauss_points);
        SIGMA.resize(number_gauss_points);
        M.resize(number_gauss_points);
        
        //Initialize stress measures to zero
        for(int i=0; i<number_gauss_points; i++){
            PK2[i]   = tensor::Tensor23({3,3});
            SIGMA[i] = tensor::Tensor23({3,3});
            M[i]     = tensor::Tensor33({3,3,3});
        }
        
        reference_coords = parse_incoming_vectors(1,rcs);
        current_coords   = parse_incoming_vectors(1,rcs);
        
        //Resize the private vector attributes
        Ns.resize(local_coords.size());
        dNdxis.resize(local_coords.size());
        dNdxs.resize(local_coords.size());
        dNdXs.resize(local_coords.size());
        
        //Set the dimension of the gradient vectors to 3
        for(int n=0; n<8; n++){dNdxis[n].resize(3);}
        for(int n=0; n<8; n++){dNdxs[n].resize(3);}
        for(int n=0; n<8; n++){dNdXs[n].resize(3);}
    }
    
    Hex8::Hex8(std::vector< double > rcs, std::vector< double > U, std::vector< double > dU,
               std::vector< double > _fparams, std::vector< int > _iparams){
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
            
            dU: The change in the degree of freedom vector over the 
                last increment. The organization is the same as for 
                U.
                
            fparams: Parameters for the constitutive model which are 
                     floating point numbers.
                     
            iparams: Parameters for the constitutive model which are 
                     integer numbers.
        */
        
        //Resize the RHS and AMATRX containers
        RHS.resize(96,0.); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX.resize(96);
        for(int i=0; i<96; i++){
            AMATRX[i].resize(96,0.);
        }
        
        //Resize the phi vector
        node_phis.resize(reference_coords.size());
        for(int n=0; n<reference_coords.size(); n++){
            node_phis[n] = tensor::Tensor23({3,3});
        }
        
        //Resize the stress measure vectors
        PK2.resize(number_gauss_points);
        SIGMA.resize(number_gauss_points);
        M.resize(number_gauss_points);
        
        //Initialize stress measures to zero
        for(int i=0; i<number_gauss_points; i++){
            PK2[i]   = tensor::Tensor23({3,3});
            SIGMA[i] = tensor::Tensor23({3,3});
            M[i]     = tensor::Tensor33({3,3,3});
        }
        
        //Resize the private vector attributes
        Ns.resize(local_coords.size());
        dNdxis.resize(local_coords.size());
        dNdxs.resize(local_coords.size());
        dNdXs.resize(local_coords.size());
        
        //Set the dimension of the gradient vectors to 3
        for(int n=0; n<8; n++){dNdxis[n].resize(3);}
        for(int n=0; n<8; n++){dNdxs[n].resize(3);}
        for(int n=0; n<8; n++){dNdXs[n].resize(3);}
        
        //Break the rcs, U, and, dU vectors into a vector of vectors
        //where each subvector are the values associated with a given 
        //node.
        reference_coords   = parse_incoming_vectors(1,rcs);
        dof_at_nodes       = parse_incoming_vectors(2,U);
        Delta_dof_at_nodes = parse_incoming_vectors(2,dU);
        
        //Assign the current value of the nodal coordinates
        for(int n=0; n<reference_coords.size(); n++){
            current_coords[n].resize(3);
            for(int i=0; i<3; i++){
                current_coords[n][i] = reference_coords[n][i]+dof_at_nodes[n][i];
            }
        }
        
        //Assign the values of phi at the nodes
        for(int n=0; n<reference_coords.size(); n++){
            //!NOTE: Assumes dof vector is set up as phi_11, phi_22, phi_33,
            //!                                      phi_23, phi_13, phi_12,
            //!                                      phi_32, phi_31, phi_21
            node_phis[n](0,0) = dof_at_nodes[n][ 3];
            node_phis[n](1,1) = dof_at_nodes[n][ 4];
            node_phis[n](2,2) = dof_at_nodes[n][ 5];
            node_phis[n](1,2) = dof_at_nodes[n][ 6];
            node_phis[n](0,2) = dof_at_nodes[n][ 7];
            node_phis[n](0,1) = dof_at_nodes[n][ 8];
            node_phis[n](2,1) = dof_at_nodes[n][ 9];
            node_phis[n](2,0) = dof_at_nodes[n][10];
            node_phis[n](1,0) = dof_at_nodes[n][11];
        }
        
        //Set the material parameters
        fparams = _fparams;
        iparams = _iparams;
    }
    
    //!==
    //!|
    //!| Operators
    //!|
    //!==
    
//    Hex8& Hex8::operator=(const Hex8& hex8_in){
//        /*!=======================
//        |      operator=       |
//        ========================
//        
//        Copy operator to allow for copying 
//        hexahedral elements*/
//        
//        reference_coords = hex8_in.reference_coords;
//        current_coords   = hex8_in.current_coords;
//    }
    
    //!==
    //!|
    //!| Methods
    //!|
    //!==
    
    //!=
    //!| Shape Functions
    //!=
    
    void Hex8::update_shape_function_values(){
        /*!======================================
        |    update_shape_function_values    |
        ======================================
        
        Set all of the shape function values to
        a value consistent with the current gauss 
        point (number gpt_num).
        
        */
        
        set_shape_functions();                  //!Set N for each node at the current gauss point in private attributes
        set_local_gradient_shape_functions();   //!Set dNdxi for each node at the current gauss point in private attributes
        set_global_gradient_shape_functions(0); //!Set dNdX for each node at the current gauss point in private attributes
        set_global_gradient_shape_functions(1); //!Set dNdx for each node at the current gauss point in private attributes
        
        
    }
    
    void Hex8::set_shape_function(int n){
        /*!==============================
        |      set_shape_function    |
        ==============================
        
        Compute the value of the shape function at 
        node n at the gauss point.
        
        Input:
            n  : Node number        
        */
        
        Ns[n] = 0.125*(1+points[gpt_num][0]*local_coords[n][0])*(1+points[gpt_num][1]*local_coords[n][1])*(1+points[gpt_num][2]*local_coords[n][2]);
        return;
    }
    
    void Hex8::set_shape_functions(){
        /*!===============================
        |      set_shape_functions    |
        ===============================
        
        Compute the value of the shape function
        at the current gauss point. Set the private 
        variable Ns with the results.
        
        Input:
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        for(int n=0; n<reference_coords.size(); n++){
            set_shape_function(n);
        }
        return;
    }
    
    void Hex8::set_local_gradient_shape_function(int n){
        /*!===========================================
        |    set_local_gradient_shape_function    |
        ===========================================
        
        Compute the value of the gradient of the 
        shape function at a given node and xi 
        (local coordinate) locations locally.
        
        Input:
            n  : Node number
        */
        
        //Initialize the output vector
        dNdxis[n][0] = 0;
        dNdxis[n][1] = 0;
        dNdxis[n][2] = 0;
        
        dNdxis[n][0] = 0.125*local_coords[n][0]*(1+points[gpt_num][1]*local_coords[n][1])*(1+points[gpt_num][2]*local_coords[n][2]);
        dNdxis[n][1] = 0.125*(1+points[gpt_num][0]*local_coords[n][0])*local_coords[n][1]*(1+points[gpt_num][2]*local_coords[n][2]);
        dNdxis[n][2] = 0.125*(1+points[gpt_num][0]*local_coords[n][0])*(1+points[gpt_num][1]*local_coords[n][1])*local_coords[n][2];
        
        return;
    }
    
    void Hex8::set_global_gradient_shape_function(bool mode, int n){
        /*!============================================
        |    set_global_gradient_shape_function    |
        ============================================
        
        Compute the value of the gradient of the 
        shape function at a given node and xi 
        (local coordinate) locations globally.
        
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
            n   : Node number
        
        */
        //Initialize the gradient
        if(mode==0){     dNdXs[n] = {0,0,0};}
        else if(mode==1){dNdxs[n] = {0,0,0};}
        
        //Compute the global gradient w.r.t. either the reference or global x
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                if(mode==0)     {dNdXs[n][i] += Jinv(j,i)*dNdxis[n][j];} //Compute the gradient of the shape function associated with node n w.r.t. the current coordinates
                else if(mode==1){dNdxs[n][i] += Jinv(j,i)*dNdxis[n][j];} //Compute the gradient of the shape function associated with node n w.r.t. the reference coordinates
            }
        }
        return;
    }
    
    void Hex8::set_jacobian(bool mode){
        /*!======================
        |    set_jacobian    |
        ======================
        
        Set the value of the jacobian of 
        transformation between the given 
        coordinates and the local coordinates.
        
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)        
        */
        
        //Initialize the jacobian value i.e. set to zero
        J.data.setZero();

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
            //J += vector_dyadic_product(dNdxis[n],coordinates[n]); //!Compute the vector dyadic product and add it to the jacobian
            //                                                        //!Jacobian computed consistent with Felippa AFEM Ch11 11.8
            J += vector_dyadic_product(coordinates[n],dNdxis[n]);   //!Compute the vector dyadic product and add it to the jacobian
                                                                    //!Jacobian computed consistent with Belytchko E4.3.8
        }
        
        Jinv    = J.inverse();
        if(mode){
            Jhatdet = J.det();
        }
        
        return;
    }
    
    void Hex8::set_local_gradient_shape_functions(){
        /*!================================================
        |      set_local_gradient_shape_functions      |
        ================================================
        
        Set all of the local gradients of the shape functions.
        
        This sets the value of the gradients of the 
        shape functions at a given gauss point.
        
        */
        
        //Populate dNdxis
        for(int n=0; n<reference_coords.size(); n++){
            set_local_gradient_shape_function(n);
        }
        
        return;
    }
    
    void Hex8::set_global_gradient_shape_functions(bool mode){
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
        
        set_jacobian(mode); //Set the jacobian and its inverse to the correct values
        
        for(int n=0; n<reference_coords.size(); n++){
            set_global_gradient_shape_function(mode, n);
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
        tensor::Tensor23 dxdxi({3,3}); //!The derivative of the current coordinates w.r.t. the local coordinates
        tensor::Tensor23 dXdxi({3,3}); //!The derivative of the reference coordinates w.r.t. the local coordinates
        
        //Compute the derivatives of the reference and current coordinates w.r.t. xi
        for(int n=0; n<reference_coords.size(); n++){
            dxdxi += vector_dyadic_product(current_coords[n],  dNdxis[n]);
            dXdxi += vector_dyadic_product(reference_coords[n],dNdxis[n]);
        }
        
        tensor::Tensor23 dxidX = dXdxi.inverse(); //!The derivative of the local coordinates w.r.t. the reference coordinates
        
        //Reset F to zero
        F.data.setZero(); //!The deformation gradient
        
        for(int i = 0; i<3; i++){
            for(int J=0; J<3; J++){
                for(int k=0; k<3; k++){
                    F(i,J) += dxdxi(i,k)*dxidX(k,J);
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
        //Reset chi to zero
        chi.data.setZero();
        
        //Interpolate the nodal phis to xi
        for(int n=0; n<reference_coords.size(); n++){
            chi += Ns[n]*node_phis[n];
        }
        chi(0,0) += 1;
        chi(1,1) += 1;
        chi(2,2) += 1;
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
        tensor::Tensor23 chi_n({3,3});          //!The value of chi at a node
        tensor::Tensor23 I     = tensor::eye(); //!The second order identity tensor
        
        //Set grad_chi to zero
        grad_chi.data.setZero();
        
        for(int n=0; n<reference_coords.size(); n++){
            chi_n = node_phis[n]+I;
            for(int i=0; i<3; i++){
                for(int J=0; J<3; J++){
                    for(int K=0; K<3; K++){
                        grad_chi(i,J,K) += chi_n(i,J)*dNdXs[n][K];
                    }
                }
            }
        }
        
        return;
    }
    
    //!=
    //!| Micromorphic Deformation Measures
    //!=
    
    void Hex8::set_deformation_measures(){
        /*!==================================
        |    set_deformation_measures    |
        ==================================
        
        Compute the deformation measures which 
        will be passed to the constitutive equation 
        subroutine to produce the stress. These measures 
        are the Right Cauchy-Green deformation measure, 
        the microdeformation tensor Psi, and the micro-
        gradient deformation measures Gamma.
        
        */
        
        compute_right_cauchy_green();
        compute_Psi();
        compute_Gamma();
    }
    
    void Hex8::compute_right_cauchy_green(){
        /*!====================================
        |    compute_right_cauchy_green    |
        ====================================
        
        Compute the right Cauchy-Green deformation 
        tensor using the deformation gradient 
        stored in the private attributes.
        
        */
        
        //Zero the contents of C
        C.data.setZero();
        
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
    
    void Hex8::compute_Psi(){
        /*!=====================
        |    compute_Psi    |
        =====================
        
        Compute micro deformation measure Psi 
        using the deformation gradient and chi 
        stored in the private attributes.
        
        */
        
        //Zero the contents of Psi
        Psi.data.setZero();
        
        //Form Psi
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int i=0; i<3; i++){
                    Psi(I,J) += F(i,I)*chi(i,J);
                }
            }
        }
    }
    
    void Hex8::compute_Gamma(){
        /*!=======================
        |    compute_Gamma    |
        =======================
        
        Compute micro deformation measure Gamma 
        using the deformation gradient and the 
        gradient of chi stored in the private 
        attributes.
        
        */
        
        //Zero the contents of Gamma
        Gamma.data.setZero();
        
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
    
    void Hex8::set_stresses(){
        /*!========================
        |     set_stresses     |
        ========================
        
        Set the stresses computed from the 
        constitutive model at the current 
        gauss point.
        
        */
        
        micro_material::get_stress(fparams,iparams,C,Psi,Gamma,PK2[gpt_num],SIGMA[gpt_num],M[gpt_num]);
        
        return;
    }
    
    //!=
    //!| Residuals
    //!=
    
    //!|=> Forces
    
    void Hex8::add_all_forces(){
        /*!========================
        |    add_all_forces    |
        ========================
        
        Add all of the forces to the residual 
        vector for a given gauss point.
        
        */
        
        add_internal_nodal_forces();
    }
    
    void Hex8::add_internal_nodal_forces(){
        /*!=====================================
        |     add_internal_nodal_forces     |
        =====================================
        
        Add the contributions of the internal 
        nodal forces to the right hand side vector.
        
        */
        
        //Put the force residuals in the RHS vector
        for(int n = 0; n<reference_coords.size(); n++){
            for(int j=0; j<3; j++){
                for(int I=0; I<3; I++){
                    for(int J=0; J<3; J++){
                        RHS[j+n*12] += -dNdXs[n][I]*PK2[gpt_num](I,J)*F(j,J)*Jhatdet*weights[gpt_num];
                    }
                }
            }
        }
        
        return;
    }
    
    //!|=> Moments
    
    void Hex8::add_all_moments(){
        /*!=========================
        |    add_all_moments    |
        =========================
        
        Add all of the moments to the 
        residual vector.
        
        */
        
        add_internal_moments();
    }
    
    void Hex8::add_internal_moments(){
        /*!=====================================
        |       add_internal_moments        |
        =====================================
        
        Add the contributions of the internal 
        moments to the right hand side vector.
        
        */
        
        //Initialize the internal moment tensor
        tensor::Tensor23 mu_int({3,3}); //!The internal moment tensor
        
        int initial_index = 3;//*reference_coords.size();
        //Create the moment residual at the different nodes
        for(int n=0; n<reference_coords.size(); n++){
            
            mu_int.data.setZero();
            
            //Create the current value of mu_int
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            mu_int(i,j) += -Ns[n]*F(i,I)*(SIGMA[gpt_num](I,J) - PK2[gpt_num](I,J))*F(j,J)*Jhatdet*weights[gpt_num];
                        }
                    }
                    
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            for(int K=0; K<3; K++){
                                mu_int(i,j) += -dNdXs[n][K]*F(j,J)*chi(i,I)*M[gpt_num](K,J,I)*Jhatdet*weights[gpt_num];
                            }
                        }
                    }
                }
            }
            
            RHS[initial_index+n*12+0] += mu_int(0,0);
            RHS[initial_index+n*12+1] += mu_int(1,1);
            RHS[initial_index+n*12+2] += mu_int(2,2);
            RHS[initial_index+n*12+3] += mu_int(1,2);
            RHS[initial_index+n*12+4] += mu_int(0,2);
            RHS[initial_index+n*12+5] += mu_int(0,1);
            RHS[initial_index+n*12+6] += mu_int(2,1);
            RHS[initial_index+n*12+7] += mu_int(2,0);
            RHS[initial_index+n*12+8] += mu_int(1,0);
        }
    }
    
    //!=
    //!| Tangents
    //!=
            
    //!|=> Compute derivatives of fundamental deformation measures
    
    void Hex8::set_fundamental_tangents(){
        /*!==================================
        |    set_fundamental_tangents    |
        ==================================
        
        Set the tangents of the fundamental 
        deformation measures.
        
        */
        
        set_dFdU();
        set_dchidU();
        set_dgradchi_dU();
        
        return;
    }
    
    void Hex8::set_dFdU(){
        /*!==================
        |    set_dFdU    |
        ==================
        
        Set the derivative of the deformation gradient w.r.t.
        the degree of freedom vector.
        
        Note that currently a large part of this tangent is 
        populated with zeros. Future updates should investigate 
        how to reduce this unused space.
        
        */
        
        for(int n=0; n<8; n++){//Iterate through the nodes
            //Set the derivatives of the deformation gradient w.r.t. the 
            //degree of freedom vector
            for(int i=0; i<3; i++){
                for(int J=0; J<3; J++){
                    for(int K=0; K<3; K++){
                        if(K==i){dFdU(i,J,K+n*12) = dNdXs[n][J];}
                    }
                }
            }
        }
        return;
    }
    
    void Hex8::set_dchidU(){
        /*!====================
        |    set_dchidU    |
        ====================
        
        Set the derivative of the micro-displacement 
        tensor with respect to the degree of freedom 
        vector.
        
        */
        
        for(int n=0; n<8; n++){//Iterate through the nodes
            //Set the derivatives of the deformation gradient w.r.t. the 
            //degree of freedom vector
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    dchidU(i,j,3+12*n+3*i+j) = Ns[n];
                }
            }
        }
        return;
    }
    
    void Hex8::set_dgradchi_dU(){
        /*!=========================
        |    set_dgradchi_dU    |
        =========================
        
        Set the derivative of the derivative 
        of the gradient of the micro-deformation 
        tensor with respect to the degree of 
        freedom vector.
        
        */
        for(int n=0; n<8; n++){//Iterate through the nodes
            //Set the derivatives of the deformation gradient w.r.t. the 
            //degree of freedom vector
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        dgrad_chidU(i,j,k,3+12*n+i*3+j) = dNdXs[n][k];
                    }
                }
            }
        }
        
        
    }
    
     //!|=> Compute derivatives of derived deformation measures
    
    void Hex8::set_dCdU(){
        /*!==================
        |    set_dCdU    |
        ==================
        
        Set the derivative of the right 
        Cauchy-Green deformation tensor 
        with respect to the degree of 
        freedom vector.
        
        Note: Could potentially remove some 
              multipliciations by intelegently 
              iterating through K.
        
        */
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<96; K++){
                    
                    for(int i=0; i<3; i++){
                        dCdU(I,J,K) += F(i,J)*dFdU(i,I,K) + F(i,I)*dFdU(i,J,K);
                    }
                }
            }
        }
        return;
    }
    
    void Hex8::set_dPsidU(){
        /*!====================
        |    set_dPsidU    |
        ====================
        
        Set the derivative of the 
        micro-deformation tensor 
        with respect to the degree 
        of freedom vector.
        
        Note: Could potentially remove some 
              multipliciations by intelegently 
              iterating through K.
        
        */
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<96; K++){
                    
                    for(int i=0; i<3; i++){
                        dPsidU(I,J,K) += chi(i,J)*dFdU(i,I,K) + F(i,I)*dchidU(i,J,K);
                    }
                    
                }
            }
        }
        return;
    }
    
    void Hex8::set_dGammadU(){
        /*!======================
        |    set_dGammadU    |
        ======================
        
        Set the derivative of Gamma 
        with respect to the degree 
        of freedom vector.
        
        Note: Could potentially remove some 
              multipliciations by intelegently 
              iterating through K.
        
        */
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int L=0; L<3; L++){
                    for(int K=0; K<96; K++){
                        
                        for(int i=0; i<3; i++){
                            dGammadU(I,J,L,K) += grad_chi(i,J,L)*dFdU(i,I,K) + F(i,I)*dgrad_chidU(i,J,L,K);
                        }
                    }
                }
            }
        }
        return;
    }
    
    //!|=> Stress tangents
    
    void Hex8::set_dPK2dU(){
        /*!====================
        |    set_dPK2dU    |
        ====================
        
        Set the derivative of the second 
        Piola-Kirchhoff stress with respect 
        to the degree of freedom vector.
        
        Note: Efficiencies may be gained by 
              more inteligently iterating through 
              the degrees of freedom.
        
        */
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<96; K++){
                    
                    //Add the contributions of the lower order stress tensors
                    for(int L=0; L<3; L++){
                        for(int M=0; M<3; M++){
                            dPK2dU(I,J,K) += dPK2dC(I,J,L,M)*dCdU(L,M,K) + dPK2dPsi(I,J,L,M)*dPsidU(L,M,K);
                        }
                    }
                    
                    //Add the contributions of the higher order stress tensor
                    for(int L=0; L<3; L++){
                        for(int M=0; M<3; M++){
                            for(int N=0; N<3; N++){
                                dPK2dU(I,J,K) += dPK2dGamma(I,J,L,M,N)*dGammadU(L,M,N,K);
                            }
                        }
                    }
                }
            }
        }
        return;
    }
    
    void Hex8::set_dSIGMAdU(){
        /*!======================
        |    set_dSIGMAdU    |
        ======================
        
        Set the derivative of the symmetric 
        stress tensor in the reference configuration 
        w.r.t. the degree of freedom vector.
        
        Note: Efficiencies may be gained by 
              more inteligently iterating through 
              the degrees of freedom.
        
        */
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<96; K++){
                    
                    //Add the contributions of the lower order deformation tensors
                    for(int L=0; L<3; L++){
                        for(int M=0; M<3; M++){
                            dSIGMAdU(I,J,K) += dSIGMAdC(I,J,L,M)*dCdU(L,M,K) + dSIGMAdPsi(I,J,L,M)*dPsidU(L,M,K);
                        }
                    }
                    
                    //Add the contributions of the higher order deformation tensor
                    for(int L=0; L<3; L++){
                        for(int M=0; M<3; M++){
                            for(int N=0; N<3; N++){
                                dSIGMAdU(I,J,K) += dSIGMAdGamma(I,J,L,M,N)*dGammadU(L,M,N,K);
                            }
                        }
                    }
                }
            }
        }
        return;
    }
    
    void Hex8::set_dMdU(){
        /*!====================
        |    set_dMdU    |
        ====================
        
        Set the derivative of the higher
        order stress with respect to the 
        degree of freedom vector.
        
        Note: Efficiencies may be gained by 
              more inteligently iterating through 
              the degrees of freedom.
        
        */
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    for(int L=0; L<96; L++){
                        
                        //Add the contributions of the lower order deformation tensors
                        for(int O=0; O<3; O++){
                            for(int P=0; P<3; P++){
                                dMdU(I,J,K,L) += dMdC(I,J,K,O,P)*dCdU(O,P,L) + dMdPsi(I,J,K,O,P)*dPsidU(O,P,L);
                            }
                        }
                        
                        //Add the contributions of the higher order deformation tensor
                        for(int O=0; O<3; O++){
                            for(int P=0; P<3; P++){
                                for(int Q=0; Q<3; Q++){
                                    dMdU(I,J,K,L) += dMdGamma(I,J,K,O,P,Q)*dGammadU(O,P,Q,L);
                                }
                            }
                        }
                        
                    }
                }
            }
        }
        return;
    }
    
    //!|=> Add contribution of gauss point to the element stiffness matrix
    
    void Hex8::add_dFintdU(){
        /*!=====================
        |    add_dFintdU    |
        =====================
        
        Add the contribution of the gauss point 
        to the derivative of the balance of linear 
        momentum to the overall element stiffness.
        
        */
        
        for(int n=0; n<8; n++){

            for(int j=0; j<3; j++){
                
                for(int K=0; K<96; K++){
                    
                    for(int I=0; I<3; I++){
                        for(int J=0; J<3; J++){
                            AMATRX[j+n*12][K] -= dNdXs[n][I]*(dPK2dU(I,J,K)*F(j,J) + PK2[gpt_num](I,J)*dFdU(j,J,K))*weights[gpt_num]*Jhatdet;
                        }
                    }
                }
            }
        }
        return;
    }
    
    void Hex8::add_dMintdU(){
        /*!=====================
        |    add_dMintdU    |
        =====================
        
        Add the derivative of the internal 
        stress term in the balance of first 
        moment of momentum with respect to 
        the degree of freedom vector to the 
        element stiffness matrix.
        
        */
        
        for(int n=0; n<8; n++){
            
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int L=0; L<96; L++){
                        
                        for(int I=0; I<3; I++){
                            for(int J=0; J<3; J++){
                                
                                AMATRX[3+n*12+3*i+j][L] -= ( Ns[n]*(dFdU(i,I,L)*(SIGMA[gpt_num](I,J) - PK2[gpt_num](I,J))*F(j,J)
                                                            +F(i,I)*(dSIGMAdU(I,J,L) - dPK2dU(I,J,L))*F(j,J) + F(i,I)*(SIGMA[gpt_num](I,J) - PK2[gpt_num](I,J))*dFdU(j,J,L)))*weights[gpt_num]*Jhatdet;
                                
                            }
                        }
                        
                        for(int I=0; I<3; I++){
                            for(int J=0; J<3; J++){
                                for(int K=0; K<3; K++){
                                    AMATRX[3+n*12+3*i+j][L] -= dNdXs[n][K]*(dFdU(j,J,L)*chi(i,I)*M[gpt_num](K,J,I) + F(j,J)*dchidU(i,I,L)*M[gpt_num](K,J,I) + F(j,J)*chi(i,I)*dMdU(K,J,I,L))*weights[gpt_num]*Jhatdet;
                                }
                            }
                        }
                        
                    }
                }
            }
        }
    }
    
    //!=
    //!| Element Integration 
    //!=
    
    void Hex8::update_gauss_point(){
        /*!============================
        |    update_gauss_point    |
        ============================
        
        Update the deformation and stress measures 
        at a gauss point
        
        */
        
        //!Compute the shape function values
        update_shape_function_values();
    
        //!Set the fundamental deformation measures
        set_fundamental_measures();
    
        //!Set the deformation measures
        set_deformation_measures();
    
        //!Set the stresses at the gauss point
        set_stresses();
        return;
    }
    
    void Hex8::integrate_element(){
        /*!===========================
        |    integrate_element    |
        ===========================
        
        Integrate the element
        
        */
        
        for(int i=0; i<number_gauss_points; i++){
            gpt_num = i;          //Update the gauss point number
            update_gauss_point(); //Update all of the gauss points
            add_all_forces();     //Update the RHS vector from all of the forces
            add_all_moments();    //Update the RHS vector from all of the stresses
        }
        
        return;
    }
    
    //!=
    //!| Test functions
    //!=
    
    void Hex8::set_gpt_num(int _gpt_num){
        /*!=====================
        |    set_gpt_num    |
        =====================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Set the gauss point number.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        
        gpt_num = _gpt_num;
    }
    
    double Hex8::get_N(int _n){
        /*!===============
        |    get_N    |
        ===============
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the shape function 
        from node n.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        
        return Ns[_n];
    }
    
    std::vector< double > Hex8::get_dNdxi(int _n){
        /*!===================
        |    get_dNdxi    |
        ===================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the shape function 
        from node n.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        
        return dNdxis[_n];
    }
    
    tensor::Tensor23 Hex8::get_jacobian(int _mode){
        /*!======================
        |    get_jacobian    |
        ======================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the jacobian for 
        either the reference or current 
        coordinates.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        
        set_jacobian(_mode);
        
        return J;
    }
    
    double Hex8::get_Jhatdet(bool _mode){
        /*!=====================
        |    get_Jdethat    |
        =====================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the determinant of the 
        jacobian.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        
        set_jacobian(_mode);
        
        return Jhatdet;
    }
    
    std::vector< double > Hex8::get_dNdx(bool _mode, int _node){
        /*!==================
        |    get_dNdx    |
        ==================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the gradient of the 
        shape function with respect to the 
        global coordinates.
        
        Input:
            mode: 0 for reference coordinates
                  1 for current coordinates
            node: Which node's shape function 
                  to return.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        
        if(_mode){return dNdxs[_node];}
        else{     return dNdXs[_node];}
    }
    
    tensor::Tensor23 Hex8::get_F(){
        /*!===============
        |    get_F    |
        ===============
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the deformation 
        gradient.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        return F;
    }
    
    tensor::Tensor23 Hex8::get_chi(){
        /*!=================
        |    get_chi    |
        =================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the microdisplacement 
        mapping.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        return chi;
    }
    
    tensor::Tensor33 Hex8::get_grad_chi(){
        /*!======================
        |    get_grad_chi    |
        ======================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the gradient of the 
        microdisplacement mapping.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        return grad_chi;
    }
    
    tensor::Tensor23 Hex8::get_C(){
        /*!===============
        |    get_C    |
        ===============
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the right Cauchy-Green
        deformation tensor.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        return C;
    }
    
    tensor::Tensor23 Hex8::get_Psi(){
        /*!=================
        |    get_Psi    |
        =================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the microdeformation 
        tensor Psi.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        return Psi;
    }
    
    tensor::Tensor33 Hex8::get_Gamma(){
        /*!===================
        |    get_Gamma    |
        ===================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the microdeformation 
        tensor Gamma.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        */
        return Gamma;
    }
    
    tensor::BaseTensor<3,288> Hex8::get_dFdU(){
        /*!==================
        |    get_dFdU    |
        ==================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        deformation gradient.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dFdU;
    }
    
    //!==
    //!|
    //!| Private Methods
    //!|
    //!=
    
    //!=
    //!| Parse DOF vectors
    //!=
            
    std::vector< std::vector< double > > Hex8::parse_incoming_vectors(int mode, const std::vector< double > & incoming){
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
        int factor;                                         //! The number of dof associated with a given mode
        
        if(mode==1){//Parse an incoming coordinate vector
            factor = 3;
        }
        else if(mode==2){//Parse an incoming dof vector
            factor = 12;
        }
        else{//Unrecognized mode
            std::cout << "\nError: The mode value of " << mode << " is not recognized.\n";
            assert(1==0);
        }
        
        //Extract the subvectors
        for(int n=0; n<reference_coords.size(); n++){
            first = incoming.begin()+n*factor;
            last  = incoming.begin()+(n+1)*factor;
            std::vector< double > new_vec(first,last);
            parsed_vector[n] = new_vec;
        }
        
        return parsed_vector;
    }
    
            
            
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
    
    tensor::Tensor23 vector_dyadic_product(const std::vector< double > & V1, const std::vector< double >& V2){
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
        Eigen::Matrix<double,3,3> data = Eigen::MatrixXd::Zero(rows,cols);
        std::vector< int > shape = {rows, cols};
        
        //Compute the dyadic product of the two vectors
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                data(i,j) = V1[i]*V2[j];
            }
        }
        
        tensor::Tensor23 T = tensor::Tensor23(shape,data);
        return T;
    }
}