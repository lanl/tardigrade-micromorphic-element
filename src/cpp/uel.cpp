#include<Eigen/Dense>
#include <micro_element.h>

extern "C" void UEL(double *RHS,double *AMATRX,double *SVARS,double *ENERGY,
                    int NDOFEL,int NRHS,int NSVARS,double *PROPS,int NPROPS,
                    double *COORDS,int MCRD,int NNODE,double *U,double *DU,
                    double *V,double *A,int JTYPE,double TIME[2],double DTIME,
                    int KSTEP,int KINC,int JELEM,double *PARAMS,int NDLOAD,
                    int *JDLTYP,double *ADLMAG,double *PREDEF,int NPREDF,
                    int *LFLAGS,int MLVARX,double *DDLMAG,int MDLOAD,
                    double PNEWDT,int *JPROPS,int NJPROP,double PERIOD){
                        
                        /*!===========================================================================================
                          |                                           UEL                                            |
                          --------------------------------------------------------------------------------------------
                          |                                                                                          |
                          | The interface between Abaqus and the micromorphic finite element.                        |
                          ============================================================================================
                          |                                                                                          |
                          |RHS:         An array containing the contributions of the element                         |
                          |             to the right-hand-side vectors of the overall system                         |
                          |             of equations.                                                                |
                          |                                                                                          |
                          |             The right hand side is organized sequentially for each                       |
                          |             node in the following order for a total of 96 terms                          |
                          |                                                                                          |
                          |             RHS = [F1,F2,F3,M11,M22,M33,M23,M13,M12,M32,M31,M21,...]                     |
                          |                                                                                          |
                          |AMATRX:      An array containing the contribution of this element                         |
                          |             to the Jacobian (stiffness) or other matrix of the                           |
                          |             overall system of equations.                                                 |
                          |                                                                                          |
                          |             with F being the force balance equation and M being the                      |
                          |             balance of the first moment of momentum. The matrix is                       |
                          |             organized as follows for each node in the following order                    |
                          |             for a total of 96x96 terms                                                   |
                          |                                                                                          |
                          |                                                                                          |
                          |             AMATRX = -[[dF1du1, dF1du2, dF1du3, dF1dphi11, dF1dphi22, ...],              |
                          |                        [dF2du1, dF2du2, dF2du3, dF2dphi11, dF2dphi22, ...],              |
                          |                        [dF2du1, dF2du2, dF2du3, dF2dphi11, dF2dphi22, ...],              |
                          |                                                                                          |
                          |SVARS:       An array containing the values of the solution-dependent                     |
                          |             state variables associated with this element. Should be                      |
                          |             updated to the values at the end of the increment unless                     |
                          |             indicated in LFLAGS.                                                         |
                          |                                                                                          |
                          |ENERGY:      ENERGY[1] = Kinetic energy                                                   |
                          |             ENERGY[2] = Elastic strain energy                                            |
                          |             ENERGY[3] = Creep dissipation                                                |
                          |             ENERGY[4] = Plastic dissipation                                              |
                          |             ENERGY[5] = Viscous dissipation                                              |
                          |             ENERGY[6] = ``Artificial strain energy''                                     |
                          |             ENERGY[7] = Electrostatic energy                                             |
                          |             ENERGY[8] = Incremental work done by load applied within                     |
                          |                         the user element                                                 |
                          |                                                                                          |
                          |PNEWDT:      Ratio of suggested new time increment to the time increment                  |
                          |             currently being used.                                                        |
                          |                                                                                          |
                          |PROPS:       A floating point array containing the NPROPS real property                   |
                          |             values defined for use with this element.                                    |
                          |                                                                                          |
                          |JPROPS:      An integer array containing the NJPROP integer property                      |
                          |                                                                                          |
                          |COORDS:      An array containing the original coordinates of the nodes of                 |
                          |             the element                                                                  |
                          |                                                                                          |
                          |U, DU, V, A: Arrays containing the current estimates of the basic solution                |
                          |             variables at the nodes of the element.                                       |
                          |             U  = Total values of the variables.                                          |
                          |             DU = Incremental values of the variables                                     |
                          |             V  = Velocity (defined only for implicit dynamics)                           |
                          |             A  = Acceleration (defined only for implicit dynamics)                       |
                          |                                                                                          |
                          |             U = [u1,u2,u3,phi11,phi22,phi33,phi23,phi13,phi12,phi32,phi31,phi21,...]     |
                          |                                                                                          |
                          |JDLTYP:      An array containing the integers used to define distributed                  |
                          |             load types for the element. Loads of type Un are identified                  |
                          |             by the integer value n in JDLTYP; loads of type UnNU are                     |
                          |             identified by the negative integer value.                                    |
                          |                                                                                          |
                          |ADLMAG:      For genreal nonlinear steps ADLMAG is the total load magnitude               |
                          |             of the K1th distributed load at the end of the current increment             |
                          |             for distributed loads of type Un. For distributed loads of type              |
                          |             UnNU, the load magnitude is defined in UEL.                                  |
                          |                                                                                          |
                          |DDLMAG:      For general nonlinear steps DDLMAG contains the increments in the            |
                          |             magnitudes of the distributed loads that are currently active on             |
                          |             this element for distributed loads of type Un.                               |
                          |                                                                                          |
                          |PREDEF:      An array containing the values of predefined field variables,                |
                          |             such as temperature in an uncoupled stress/displacement analysis,            |
                          |             at the nodes of the element.                                                 |
                          |                                                                                          |
                          |             The first index K1 is either 1 or 2, with 1 indicating the value             |
                          |             of the field variable at the end of the increment and 2 indicating           |
                          |             the increment in the field variable. The second index K2, indicates          |
                          |             the variable: the temperature coresponds to index 1, and the                 |
                          |             predefined field variables correspond to indices 2 and above. In             |
                          |             cases where temperature is not defined, the predefined field                 |
                          |             variables begin with index 1. The third index, K3, indicates the             |
                          |             local node number on the element.                                            |
                          |                                                                                          |
                          |             PREDEF(K1, 1,K3) = Temperature                                               |
                          |             PREDEF(K1, 2,K3) = First predefined field variable                           |
                          |             PREDEF(K1, 3,K3) = Second predefined field variable                          |
                          |             Etc.               Any other predefined field variable                       |
                          |             PREDEF(K1,K2,K3) = Total or incremental value of the K2th predefined         |
                          |                                field variable at the K3th node of the element.           |
                          |PREDEF(1,K2,K3)  = Values of the variables at the end of the current                      |
                          |                   increment                                                              |
                          |PREDEF(2,K2,K3)  = Incremental values corresponding to the current                        |
                          |                   time increment                                                         |
                          |                                                                                          |
                          |PARAMS:      An array containing the parameters associated with the solution              |
                          |             procedure. The entries in this array depend on the solution procedure        |
                          |             currently being used when UEL is called, as indicated by the entries         |
                          |             in the LFLAGS array.                                                         |
                          |                                                                                          |
                          |             For implicit dynamics (LFLAGS(1) = 11 or 12) PARAMS contains the             |
                          |             integration operator values as                                               |
                          |                                                                                          |
                          |             PARAMS(1) = alpha                                                            |
                          |             PARAMS(2) = beta                                                             |
                          |             PARAMS(3) = gamma                                                            |
                          |                                                                                          |
                          |LFLAGS:      An array containing the flags that define the current solution               |
                          |             procedure and requirements for element calculations.                         |
                          |                                                                                          |
                          |             LFLAGS(1) = Defines the procedure type                                       |
                          |             LFLAGS(2) = 0 Small-displacement analysis                                    |
                          |             LFLAGS(2) = 1 Large-displacement analysis                                    |
                          |             LFLAGS(3) = 1 Normal implicit time incrementation procedure. User            |
                          |                           subroutine UEL must define the residual vector in              |
                          |                           RHS and the Jacobian matrix in AMATRX                          |
                          |             LFLAGS(3) = 2 Define the current stiffness matrix (AMATRX = K^{NM}           |
                          |                           = -dF^N/du^M or -dG^N/du^M) only.                              |
                          |             LFLAGS(3) = 3 Define the current damping matrix (AMATRX = C^{NM}             |
                          |                           = -dF^N/ddotu^M or -dG^N/ddotu^M) only.                        |
                          |             LFLAGS(3) = 4 Definen the current mass matrix (AMATRX = -dF^N/dddotu^M) only |
                          |             LFLAGS(3) = 5 Define the current residual or load vector (RHS = F^N) only.   |
                          |             LFLAGS(3) = 6 Define the current mass matrix and the residual vector for     |
                          |                           the initial acceleration calculation (or the calculation of    |
                          |                           accelerations after impact).
                          |             LFLAGS(3) = 100 Define peturbation quantities for output.
                          |             LFLAGS(4) = 0 The step is a general step
                          |             LFLAGS(4) = 1 The step is a linear perturbation step
                          |             LFLAGS(5) = 0 The current approximations to u^M, etc. were based on Newton
                          |                           corrections.
                          |             LFLAGS(5) = 1 The current approximations were found by extrapolation from
                          |                           the previous increment.
                          |TIME(1):     Current value of step time
                          |TIME(2):     Current value of total time
                          |DTIME:       Time increment
                          |PERIOD:      Time period of the current step
                          |NDOFEL:      Number of degrees of freedom in the element
                          |MLVARX:      Dimensioning parameter used when several displacement or 
                          |             right-hand-side vectors are used.
                          |NRHS:        Number of load vetors. NRHS is 1 in most nonlinear problems:
                          |             it is 2 for the modified Riks static procedure, and it is 
                          |             greater than 1 in some linear analysis procedures and during 
                          |             substructure generation.
                          |NSVARS:      User-defined number of solution-dependent state variables 
                          |             associated with the element.
                          |NPROPS:      User-defined number of real property values associated with 
                          |             the element.
                          |NJPROP:      User-defined number of integer property values assocated with 
                          |             the element.
                          |MCRD:        MCRD is defined as the maximum of the user-defined maximum 
                          |             number of coordinates needed at any node point and the 
                          |             value of the largest active degree of freedom of the user
                          |             element that is less than or equal to 3.
                          |NNODE:       User-defined number of nodes on the element
                          |JTYPE:       Integer defining the element type. This is the user-defined
                          |             integer value n in element type Un.
                          |KSTEP:       Current step number.
                          |KINC:        Current increment number.
                          |JELEM:       User-assigned element number.
                          |NDLOAD:      Identification number of the distributed load or flux
                          |             currently active on this element.
                          |MDLOAD:      Total number of distributed loads and/or flues defined 
                          |             on this element.
                          |NPREDF:      Number of predefined field variables, including temperature.
                          ===========================================================================================
                          
                        */
                        
                        //!Form Eigen::Map representations of the incoming vectors
                        
                        Matrix_Xd_Map RHS_mat(RHS,NDOFEL,NRHS);
                        Matrix_Xd_Map AMATRX_mat(AMATRX,NDOFEL,NDOFEL);
                        Vector_Xd_Map SVARS_vec(SVARS,NSVARS,1);
                        Vector_8d_Map ENERGY_vec(ENERGY,8,1);
                        Vector_Xd_Map PROPS_vec(PROPS,NPROPS,1);
                        Vector_Xd_Map JPROPS_vec(JPROPS,NJPROPS,1);
                        
                        Matrix_Xd_Map COORDS_mat(COORDS,NNODE,MCRD);
                        Vector_Xd_Map U_vec(U,NDOFEL,1);
                        Vector_Xd_Map DU_vec(DU,NDOFEL,1);
                        Vector_Xd_Map V_vec(V,NDOFEL,1);
                        Vector_Xd_Map A_vec(A,NDOFEL,1);
                        Vector_3d_Map PARAMS_vec(PARAMS,3,1);
                        
                        Matrix_Xd_Map JDLTYP_mat(JDLTYP,MDLOAD,1); //May not be totally general.
                        Vector_Xd_Map ADLMAG_vec(ADLMAG,MDLOAD,1);
                        
                        Vector_5i_Map LFLAGS_vec(LFLAGS,5,1);
                        
                        Matrix_Xd_Map DDLMAG_mat(DDLMAG,MDLOAD,1); //May not be totally general.
                        
                        Vector_Xd_Map JPROPS_vec(JPROPS,NJPROP,1);
                        
                        //!Parse the incoming element to the correct user subroutine
                        if(JTYPE==1){
                            compute_hex8(RHS_mat,AMATRX_mat,SVARS_vec,ENERGY_vec,
                                         PROPS_vec,COORDS_mat,U_vec,DU_vec,V_vec,A_vec,
                                         TIME,DTIME,KSTEP,KINC,JELEM,PARAMS_vec,
                                         JDLTYP_mat,ADLMAG_vec,PREDEF,NPREDEF,LFLAGS,
                                         DDLMAG,PNEWDT,JPROPS_vec,PERIOD);
                        }
                        else{
                            std::cout << "\nError: Element type not recognized";
                        }
                    }
                    
void compute_hex8(Matrix_Xd_Map &RHS, Matrix_Xd_Map &AMATRX, Vector_Xd_Map &SVARS_vec,
                  Vector_8d_Map &ENERGY, Vector_3d_Map &PROPS, Matrix_Xd_Map &COORDS, Vector_Xd_Map &U, 
                  Vector_Xd_Map &DU, Vector_Xd_Map &V, Vector_Xd_Map &A, double TIME[2], double DTIME, 
                  int KSTEP,int KINC,int JELEM, Vector_3d_Map &PARAMS, Matrix_Xd_Map &JDLTYP, 
                  Vector_Xd_Map &ADLMAG, double *PREDEF, int NPREDF, Vector_5i_Map &LFLAGS, 
                  Matrix_Xd_Map &DDLMAG, double PNEWDT, Vector_Xd_Map &JPROPS_vec, double PERIOD):
    /*!====================
    |   compute_hex8   |
    ====================
    
    Compute the response for a hex8 
    micromorphic element.
    
    */
    
    //!Reform the incoming coordinates to a form which the element can read
                        
    std::vector<double> flat_coords;
    for(int i=0; i<NNODE; i++){
        for(int j=0; j<dimension; j++){
            flat_coords[j+i*8] = COORDS_mat(i,j);
        }
    }
    
    //!Create the element
    micro_element::hex8 element(COORDS,RHS,AMATRX,SVARS,ENERGY,PROPS,U,DU,V,A,TIME,DTIME,KSTEP,KINC,JELEM,
                                PARAMS,JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,DDLMAG,PNEWDT,JPROPS,PERIOD);
    
    //!Compute the required values