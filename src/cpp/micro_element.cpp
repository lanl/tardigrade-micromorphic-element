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
#include <fstream>
#include <vector>
#include <tensor.h>
#include <micro_element.h>
#include <tardigrade_micromorphic_linear_elasticity.h>
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
        RHS    = Eigen::MatrixXd::Zero(96,1);//Matrix_Xd_Map(RHS_PTR,96,1); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX = Eigen::MatrixXd::Zero(96,96);//Matrix_Xd_Map(AMATRX_PTR,96,96);
        
        //Initialize the mini_mass matrix
        mini_mass = Matrix_8d::Zero(8,8);
        
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
        RHS    = Eigen::MatrixXd::Zero(96,1); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX = Eigen::MatrixXd::Zero(96,96);
        
        //Initialize the mini_mass matrix
        mini_mass = Matrix_8d::Zero(8,8);
        
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
               std::vector<double> _fparams, std::vector<int> _iparams){
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
                freedom are in the order of the nodes. Each node's 
                dof is organized:
                
                u1, u2, u3, phi11, phi22, phi33, phi23, phi13, phi12, phi32, phi31, phi21
            
            dU: The change in the degree of freedom vector over the 
                last increment. The organization is the same as for 
                U.
                
            fparams: Parameters for the constitutive model which are 
                     floating point numbers.
                     
            iparams: Parameters for the constitutive model which are 
                     integer numbers.
        */
        
        //Resize the RHS and AMATRX containers
        RHS    = Eigen::MatrixXd::Zero(96,1); //Allocate the required memory for the right hand side vector
        //Allocate the required memory for the AMATRX
        AMATRX = Eigen::MatrixXd::Zero(96,96);
        
        //Initialize the mini_mass matrix
        mini_mass = Matrix_8d::Zero(8,8);
        
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
        fparams.resize(_fparams.size(),1);
        iparams.resize(_iparams.size(),1);
        for(int i=0; i<_fparams.size(); i++){fparams(i) = _fparams[i];}
        for(int i=0; i<_iparams.size(); i++){iparams(i) = _iparams[i];}
    }
    
    Hex8::Hex8(double *_RHS,          double *_AMATRX,      Vector &_SVARS, energy_vector &ENERGY,
               Vector &PROPS,         Matrix_RM &COORDS,    Vector &U,      Vector &DU,
               Vector &V,             Vector &A,            double TIME[2], double DTIME, 
               int KSTEP,             int KINC,             int JELEM,      params_vector &PARAMS,
               Matrixi_RM &JDLTYP,    Vector &ADLMAG,       double *PREDEF, int NPREDF,
               lflags_vector &LFLAGS, Matrix_RM &DDLMAG,    double PNEWDT,  Vectori &JPROPS,
               double PERIOD,         std::string output_fn){
        /*!====================
        |        Hex8       |
        =====================
        
        The constructor for a hexehedral element when 
        used for an Abaqus UMAT
        
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
            
            RHS:    A pointer to the RHS array
            AMATRX: A pointer to the AMATRX array
            SVARS:  A pointer to the SVARS array
            ENERGY: A Eigen::Map which contians a pointer to the ENERGY array
            PROPS:  A Eigen::Map which contains a pointer to the PROPS array
            COORDS: A Eigen::Map which contains a pointer to the COORDS array
            U:      A Eigen::Map which contains a pointer to the U array
            DU:     A Eigen::Map which contains a pointer to the DU array
            V:      A Eigen::Map which contains a pointer to the V array
            A:      A Eigen::Map which contains a pointer to the A array
            TIME:   The current step and total time
            DTIME:  The change in time
            KSTEP:  The step count
            KINC:   The increment count
            JELEM:  The user-defined element number
            PARAMS: Integration parameters
            JDLTYP: A Eigen::Map which contains a pointer to the JDLTYP array
            ADLMAG: A Eigen::Map which contains a pointer to the ADLMAG array
            PREDEF: The pointer to the pre-defined field array
            NPREDF: The number of pre-defined fields
            LFLAGS: A Eigen::Map which contains a pointer to the LFLAGS array
            DDLMAG: A Eigen::Map which containts a pointer to the DDLMAG array
            PNEWDT: The new time-step ratio compared to the current DTIME
            JPROPS: A Eigen::Map which contains a pointer to the JPROPS array
            PERIOD: Time period of the current step.
        */
        
        //Assign the RHS and AMATRX arrays
        RHS = Matrix_Xd_Map(_RHS,96,1);
        //Allocate the required memory for the AMATRX
        AMATRX = Matrix_Xd_Map(_AMATRX,96,96);
        
        //Initialize the mini_mass matrix
        mini_mass = Matrix_8d::Zero(8,8);

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
        reference_coords   = parse_incoming_vectors(1,COORDS);
        dof_at_nodes       = parse_incoming_vectors(2,U);
        Delta_dof_at_nodes = parse_incoming_vectors(2,DU);
        
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
        
        //Set the state variables
        SVARS   = _SVARS;

        //Set the material parameters
        fparams = PROPS;
        iparams = JPROPS;

        //Set the output filename
        output_name = output_fn;
        step_num    = KSTEP;
        inc_num     = KINC;
        el_num      = JELEM;
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
    
    void Hex8::update_shape_function_values(bool compute_mass){
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
        if(compute_mass){
            update_mini_mass();
        }
        
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
        if(!mode){
            Jhatdet = J.det(); //!Jhatdet is the jacobian of transformation for volumes in the reference coordinates to volumes in 
                               //!the local (or master) coordinate system.
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
    
    void Hex8::update_mini_mass(){
        /*!==========================
        |    update_mini_mass    |
        ==========================
        
        Update the condensed mass matrix. 
        It should be noted that because the 
        density is assumed to be constant 
        over the element, the mini mass 
        matrix will be populated by the 
        density times Jdethat. This may not 
        always be the case if future 
        enhancements are performed.
        
        */
        
        for(int i=0; i<number_gauss_points; i++){
            for(int j=0; j<number_gauss_points; j++){
                mini_mass(i,j) += Ns[i]*Ns[j]*fparams(0)*Jhatdet*weights[gpt_num];
            }
        }
        
    }
    
    void Hex8::set_mass_matrix(){
        /*!=========================
        |    set_mass_matrix    |
        =========================
        
        Set the correct values in AMATRX to the 
        expected values of the mass matrix. This 
        is the approach used in Abaqus and may not 
        be totally general.
        
        */
        
        for(int i=0; i<number_gauss_points; i++){
            for(int j=0; j<number_gauss_points; j++){
                for(int k=0; k<12; k++){
                    AMATRX(k+12*i,k+12*j) = mini_mass(i,j);
                }
            }
        }
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
    
    void Hex8::set_stresses(bool set_tangents){
        /*!========================
        |     set_stresses     |
        ========================
        
        Set the stresses computed from the 
        constitutive model at the current 
        gauss point.
        
        */
        
        if(set_tangents){
            micro_material::get_stress(     fparams,       iparams,
                                                  C,           Psi,        Gamma,
                                       PK2[gpt_num],SIGMA[gpt_num],   M[gpt_num],
                                             dPK2dC,      dPK2dPsi,   dPK2dGamma,
                                           dSIGMAdC,    dSIGMAdPsi, dSIGMAdGamma,
                                               dMdC,        dMdPsi,     dMdGamma);
                                           
            stress_tangent_flag = true; //Alert the element that the stress tangents have been set
        
        }
        else{
            micro_material::get_stress(fparams,iparams,C,Psi,Gamma,PK2[gpt_num],SIGMA[gpt_num],M[gpt_num]);
        }

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
                        RHS(j+n*12) += -dNdXs[n][I]*PK2[gpt_num](I,J)*F(j,J)*Jhatdet*weights[gpt_num];
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
            
            RHS(initial_index+n*12+0) += mu_int(0,0);
            RHS(initial_index+n*12+1) += mu_int(1,1);
            RHS(initial_index+n*12+2) += mu_int(2,2);
            RHS(initial_index+n*12+3) += mu_int(1,2);
            RHS(initial_index+n*12+4) += mu_int(0,2);
            RHS(initial_index+n*12+5) += mu_int(0,1);
            RHS(initial_index+n*12+6) += mu_int(2,1);
            RHS(initial_index+n*12+7) += mu_int(2,0);
            RHS(initial_index+n*12+8) += mu_int(1,0);
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
        
        dFdU_flag = true; //Alert the element that dFdU has been set
        
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
            
            dchidU(0,0, 3+12*n) = Ns[n];
            dchidU(1,1, 4+12*n) = Ns[n];
            dchidU(2,2, 5+12*n) = Ns[n];
            dchidU(1,2, 6+12*n) = Ns[n];
            dchidU(0,2, 7+12*n) = Ns[n];
            dchidU(0,1, 8+12*n) = Ns[n];
            dchidU(2,1, 9+12*n) = Ns[n];
            dchidU(2,0,10+12*n) = Ns[n];
            dchidU(1,0,11+12*n) = Ns[n];
            
        }
        
        dchidU_flag = true; //Alert the element that dchidU has been set
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
            
            for(int K=0; K<3; K++){
                
                dgrad_chidU(0,0,K, 3+12*n) = dNdXs[n][K];
                dgrad_chidU(1,1,K, 4+12*n) = dNdXs[n][K];
                dgrad_chidU(2,2,K, 5+12*n) = dNdXs[n][K];
                dgrad_chidU(1,2,K, 6+12*n) = dNdXs[n][K];
                dgrad_chidU(0,2,K, 7+12*n) = dNdXs[n][K];
                dgrad_chidU(0,1,K, 8+12*n) = dNdXs[n][K];
                dgrad_chidU(2,1,K, 9+12*n) = dNdXs[n][K];
                dgrad_chidU(2,0,K,10+12*n) = dNdXs[n][K];
                dgrad_chidU(1,0,K,11+12*n) = dNdXs[n][K];
                
            }
        }
        
        bool dgrad_chidU_flag = true; //Alert the element that dgrad_chidU has been set
        
        return;
    }
        
    //!|=> Compute the derivatives of the deformation measures
        
    void Hex8::set_deformation_tangents(){
        /*!==================================
        |    set_deformation_tangents    |
        ==================================
        
        Set the tangents of the deformation 
        measures.
        
        */
        
        set_dCdU();
        set_dPsidU();
        set_dGammadU();
        
        return;
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
        
        dCdU_flag = true; //Alert the element that dCdU has been set
        
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
        
        dPsidU_flag = true; //Alert the element that dPsidU has been set
        
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
        
        dGammadU_flag = true; //Alert the element that dGammadU has been set
        
        return;
    }
    
    //!|=> Stress tangents
    
    void Hex8::set_stress_tangents(){
        /*!=============================
        |    set_stress_tangents    |
        =============================
        
        Set the tangents of the stress 
        measures with respect to the 
        degree of freedom vector.
        
        */
        
        set_dPK2dU();
        set_dSIGMAdU();
        set_dMdU();
        return;
    }
    
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
        
        dPK2dU_flag = true; //Alert the element that dPK2dU has been set
        
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
        
        dSIGMAdU_flag = true; //Alert the element that dSigmadU has been set
        
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
        
        dMdU_flag = true; //Alert the element that dMdU has been set
        
        return;
    }
    
    //!|=> Add contribution of gauss point to the element stiffness matrix
    
    void Hex8::set_force_tangent(){
        /*!===========================
        |    set_force_tangent    |
        ===========================
        
        Set the tangent of the linear momentum 
        with respect to the degree of freedom
        vector.
        
        */
        
        add_dFintdU();
        return;
    }
    
    void Hex8::add_dFintdU(){
        /*!=====================
        |    add_dFintdU    |
        =====================
        
        Add the contribution of the gauss point 
        to the derivative of the balance of linear 
        momentum to the overall element stiffness.
        
        */
        
        for(int n=0; n<reference_coords.size(); n++){

            for(int j=0; j<3; j++){
                
                for(int I=0; I<3; I++){
                    
                    for(int J=0; J<3; J++){
                        
                        for(int K=0; K<96; K++){
                            
                            AMATRX(j+n*12,K) -= dNdXs[n][I]*(dPK2dU(I,J,K)*F(j,J) + PK2[gpt_num](I,J)*dFdU(j,J,K))*weights[gpt_num]*Jhatdet;
                            
                        }
                    }
                }
            }
        }
        
        return;
    }
    
    void Hex8::set_moment_tangent(){
        /*!============================
        |    set_moment_tangent    |
        ============================
        
        Set the tangents for the first 
        moment of momentum with respect 
        to the degree of freedom vector.
        
        */
        
        add_dMintdU();
        
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
        
        double T; //!Internal variable to make code more readable
        
        int value_map[3][3] = { {0,5,4},
                                {8,1,3},
                                {7,6,2}}; //!Map of i,j indicices to the AMATRX index
        
        for(int n=0; n<8; n++){
            
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int L=0; L<96; L++){
                        
                        for(int I=0; I<3; I++){
                            for(int J=0; J<3; J++){
                                
                                AMATRX(3+n*12+value_map[i][j],L) -= ( Ns[n]*(dFdU(i,I,L)*(SIGMA[gpt_num](I,J) - PK2[gpt_num](I,J))*F(j,J)
                                                                                  +F(i,I)*(    dSIGMAdU(I,J,L) -     dPK2dU(I,J,L))*F(j,J)
                                                                                  +F(i,I)*(SIGMA[gpt_num](I,J) - PK2[gpt_num](I,J))*dFdU(j,J,L)))*weights[gpt_num]*Jhatdet;
                                
                            }
                        }
                        
                        for(int K=0; K<3; K++){
                            T = 0;
                            for(int I=0; I<3; I++){
                                for(int J=0; J<3; J++){
                                    T += dFdU(j,J,L)*chi(i,I)*M[gpt_num](K,J,I) + F(j,J)*dchidU(i,I,L)*M[gpt_num](K,J,I) + F(j,J)*chi(i,I)*dMdU(K,J,I,L);
                                }
                            }
                            AMATRX(3+n*12+value_map[i][j],L) -= dNdXs[n][K]*T*weights[gpt_num]*Jhatdet;
                        }
                        
                    }
                }
            }
        }
    }
    
    void Hex8::reset_tangents(){
        /*!========================
        |    reset_tangents    |
        ========================
        
        Reset the tangents which have been set
        
        */
        
        //!Reset fundamental deformation tangents
        if(dFdU_flag)        {dFdU.zero();                dFdU_flag        = false;}
        if(dchidU_flag)      {dchidU.zero();              dchidU_flag      = false;}
        if(dgrad_chidU_flag) {dgrad_chidU.zero();         dgrad_chidU_flag = false;}
        
        //!Reset deformation measure tangents
        if(dCdU_flag)        {dCdU.zero();                dCdU_flag        = false;}
        if(dPsidU_flag)      {dPsidU.zero();              dPsidU_flag      = false;}
        if(dGammadU_flag)    {dGammadU.zero();            dGammadU_flag    = false;}
        
        //!Reset stress tangents from material model
        if(stress_tangent_flag){
            
            //!Set derivatives with respect to C to zero
            dPK2dC.zero();
            dSIGMAdC.zero();
            dMdC.zero();
            
            //!Set derivatives with respect to Psi to zero
            dPK2dPsi.zero();
            dSIGMAdPsi.zero();
            dMdPsi.zero();
            
            //!Set derivatives with respect to Gamma to zero
            dPK2dGamma.zero();
            dSIGMAdGamma.zero();
            dMdGamma.zero();
            
            stress_tangent_flag = false;
        }
        
        //!Reset stress tangents
        
        if(dPK2dU_flag)        {dPK2dU.zero();                dPK2dU_flag        = false;}
        if(dSIGMAdU_flag)      {dSIGMAdU.zero();              dSIGMAdU_flag      = false;}
        if(dMdU_flag)          {dMdU.zero();                  dMdU_flag          = false;}
        
        return;
    }
    
    //!=
    //!| Element Integration 
    //!=
    
    void Hex8::update_gauss_point(bool set_tangents, bool compute_mass){
        /*!============================
        |    update_gauss_point    |
        ============================
        
        Update the deformation and stress measures 
        at a gauss point
        
        */
        
        //!Reset tangents if required
        reset_tangents();
        
        //!Compute the shape function values
        update_shape_function_values(compute_mass);
    
        //!Set the fundamental deformation measures
        set_fundamental_measures();
    
        //!Set the deformation measures
        set_deformation_measures();
        
        //!Set the stresses at the gauss point
        
        //!Set the tangents if required
        if(set_tangents){
            set_fundamental_tangents();
            set_deformation_tangents();
            set_stresses(set_tangents);
            set_stress_tangents();
        }
        else{
            set_stresses(set_tangents);
        }
        
        return;
    }
    
    void Hex8::integrate_element(bool set_tangents, bool ignore_RHS, bool compute_mass, bool output_stress){
        /*!===========================
        |    integrate_element    |
        ===========================
        
        Integrate the element
        
        */
        
        for(int i=0; i<number_gauss_points; i++){
            gpt_num = i;                                   //Update the gauss point number
            update_gauss_point(set_tangents,compute_mass); //Update all of the gauss points
            if(!ignore_RHS){
                add_all_forces();                          //Update the RHS vector from all of the forces
                add_all_moments();                         //Update the RHS vector from all of the stresses
            }
            if(set_tangents){                              //Update AMATRX from the forces and stresses
                set_force_tangent();
                set_moment_tangent();
            }
        }

        if(output_stress){write_output();}
        
        if(compute_mass){
            if(set_tangents){
                std::cout << "Error: Tangent and mass matrix cannot be simultaneously requested\n";
                assert(1==0);
            }
            set_mass_matrix();
        }

        return;
    }
    
    //!=
    //!| Output functions
    //!=

    void Hex8::write_output(){
        /*!======================
        |    write_output    |
        ======================

        Write the data to the output file. This 
        may be used as an error reporting method 
        in the future.

        */
        

        std::ofstream f;
        f.open(output_name + '_' + std::to_string(el_num) + "_" + std::to_string(step_num) + "_" + std::to_string(inc_num));

        for(int i=0; i<8; i++){
            //Output the PK2 stress
            f << "PK2\n";
            f << PK2[i](0,0) << ", " << PK2[i](1,1) << ", " << PK2[i](2,2) << ", " <<
                 PK2[i](1,2) << ", " << PK2[i](0,2) << ", " << PK2[i](0,1) << ", " <<
                 PK2[i](2,1) << ", " << PK2[i](2,0) << ", " << PK2[i](1,0) << "\n";

            //Output the symmetric stress
            f << "SIGMA\n";
            f << SIGMA[i](0,0) << ", " << SIGMA[i](1,1) << ", " << SIGMA[i](2,2) << ", " <<
                 SIGMA[i](1,2) << ", " << SIGMA[i](0,2) << ", " << SIGMA[i](0,1) << "\n";

            //Output the higher order couple stress
            f << "M\n";
            f << M[i](0,0,0) << ", " << M[i](1,1,0) << ", " << M[i](2,2,0) << ", " <<
                 M[i](1,2,0) << ", " << M[i](0,2,0) << ", " << M[i](0,1,0) << ", " <<
                 M[i](2,1,0) << ", " << M[i](2,0,0) << ", " << M[i](1,0,0) << "\n" <<
                 M[i](0,0,1) << ", " << M[i](1,1,1) << ", " << M[i](2,2,1) << ", " <<
                 M[i](1,2,1) << ", " << M[i](0,2,1) << ", " << M[i](0,1,1) << ", " <<
                 M[i](2,1,1) << ", " << M[i](2,0,1) << ", " << M[i](1,0,1) << "\n" <<
                 M[i](0,0,2) << ", " << M[i](1,1,2) << ", " << M[i](2,2,2) << ", " <<
                 M[i](1,2,2) << ", " << M[i](0,2,2) << ", " << M[i](0,1,2) << ", " <<
                 M[i](2,1,2) << ", " << M[i](2,0,2) << ", " << M[i](1,0,2) << "\n";

        }
        f.close();


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
    
    tensor::BaseTensor<3,288> Hex8::get_dchidU(){
        /*!====================
        |    get_dchidU    |
        ====================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        micro-displacement tensor.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dchidU;
    }
    
    tensor::BaseTensor<9,288> Hex8::get_dgrad_chidU(){
        /*!=========================
        |    get_dgrad_chidU    |
        =========================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        gradient of the micro-displacement tensor.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dgrad_chidU;
    }
    
    tensor::BaseTensor<3,288> Hex8::get_dCdU(){
        /*!==================
        |    get_dCdU    |
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
        right Cauchy-Green deformation 
        tensor.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dCdU;
    }
    
    tensor::BaseTensor<3,288> Hex8::get_dPsidU(){
        /*!====================
        |    get_dPsidU    |
        ====================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        micro-deformation gradient.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dPsidU;
    }
    
    tensor::BaseTensor<9,288> Hex8::get_dGammadU(){
        /*!======================
        |    get_dGammadU    |
        ======================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        gradient of the micro-deformation 
        gradient.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dGammadU;
    }
    
    tensor::BaseTensor<3,288> Hex8::get_dPK2dU(){
        /*!====================
        |    get_dPK2dU    |
        ====================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        second Piola-Kirchhoff stress tensor.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dPK2dU;
    }
    
    tensor::BaseTensor<3,288> Hex8::get_dSIGMAdU(){
        /*!======================
        |    get_dSIGMAdU    |
        ======================
        
        !!!!!!!!!!! WARNING !!!!!!!!!!!!!!!
        ! DO NOT USE THIS FUNCTION EXCEPT !
        ! TO TEST THE CODE!               !
        !                                 !
        ! ELEMENT INTEGRATION SHOULD BE   !
        ! PERFORMED USING EXISTING        !
        ! METHODS!                        !
        !!!!!!!!!!!!! WARNING !!!!!!!!!!!!!
        
        Get the value of the tangent of the 
        symmetric micro-stress.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dSIGMAdU;
    }
    
    tensor::BaseTensor<9,288> Hex8::get_dMdU(){
        /*!==================
        |    get_dMdU    |
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
        gradient of the higher order stress 
        tensor.
        
        Used to access the private variable 
        from outside the class. This should 
        not be done in general but is allowed 
        here for testing purposes.
        
        */
        
        return dMdU;
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
    
    std::vector< std::vector< double > > Hex8::parse_incoming_vectors(int mode, const Matrix_RM &incoming){
        /*!================================
        |    parse_incoming_vectors    |
        ================================
        
        Takes incoming Eigen::Matrix objects and maps them to 
        vectors.
        
        */
        
        //!Variable defintions
        std::vector< double >::const_iterator first;        //! The iterator which identifies the start of a subvector
        std::vector< double >::const_iterator last;         //! The iterator which identifies the end of a subvector
        std::vector< std::vector< double > > parsed_vector; //! The parsed vector which is returned
        parsed_vector.resize(reference_coords.size());
        int factor;                                         //! The number of dof associated with a given mode
        
        if(mode==1){//Parse an incoming coordinate Eigen::Map
            factor = 3;
        }
        else if(mode==2){//Parse an incoming dof Eigen::Map
            factor = 12;
        }
        else{//Unrecognized mode
            std::cout << "\nError: The mode value of " << mode << " is not recognized.\n";
            assert(1==0);
        }

        //Extract the subvectors
        for(int n=0; n<reference_coords.size(); n++){
            parsed_vector[n] = std::vector<double>(factor,0.);

            if(mode==1){
                for(int i=0; i<factor; i++){parsed_vector[n][i] = incoming(n,i);}
            }
            else if(mode==2){
                for(int i=0; i<factor; i++){parsed_vector[n][i] = incoming(i+n*factor);}
            }
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
