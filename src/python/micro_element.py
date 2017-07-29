import numpy as np
import hex8
import unittest
import os
import finite_difference as fd
import micromorphic_linear_elasticity as micro_LE
from hex8 import T_to_V_mapping as T2V
from hex8 import V_to_T_mapping as V2T

"""===========================================
|                                         |
| Definition of a 3D Micromorphic Element |
|                                         |
===========================================
| Intended to mimic an Abaqus UEL to      |
| allow easy transfer between codes. This |
| code has not been written for speed but |
| for more clarity while minimizing the   |
| dependence upon specialized python      |
| functions in the non-testing areas.     |
==========================================="""

def UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,\
        PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,\
        KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,\
        LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD):
    """The main subroutine for the element. Intentionally mimics
    an Abaqus UEL for easy porting between versions.
    
    RHS:         An array containing the contributions of the element 
                 to the right-hand-side vectors of the overall system
                 of equations.
                 
                 The right hand side is organized sequentially for each 
                 node in the following order for a total of 96 terms
                 
                 RHS = [F1,F2,F3,M11,M22,M33,M23,M13,M12,M32,M31,M21,...]
              
    AMATRX:      An array containing the contribution of this element
                 to the Jacobian (stiffness) or other matrix of the 
                 overall system of equations.
                
                 with F being the force balance equation and M being the 
                 balance of the first moment of momentum. The matrix is 
                 organied as follows for each node in the following order
                 for a total of 96x96 terms
                 
                 
                 AMATRX = -[[dF1du1, dF1du2, dF1du3, dF1dphi11, dF1dphi22, ...],
                            [dF2du1, dF2du2, dF2du3, dF2dphi11, dF2dphi22, ...],
                            [dF2du1, dF2du2, dF2du3, dF2dphi11, dF2dphi22, ...],
              
    SVARS:       An array containing the values of the solution-dependent
                 state variables associated with this element. Should be 
                 updated to the values at the end of the increment unless
                 indicated in LFLAGS.
              
    ENERGY:      ENERGY[1] = Kinetic energy
                 ENERGY[2] = Elastic strain energy
                 ENERGY[3] = Creep dissipation
                 ENERGY[4] = Plastic dissipation
                 ENERGY[5] = Viscous dissipation
                 ENERGY[6] = ``Artificial strain energy''
                 ENERGY[7] = Electrostatic energy
                 ENERGY[8] = Incremental work done by load applied within 
                             the user element
              
    PNEWDT:      Ratio of suggested new time increment to the time increment 
                 currently being used.
              
    PROPS:       A floating point array containing the NPROPS real property
                 values defined for use with this element.
              
    JPROPS:      An integer array containing the NJPROP integer property
               
    COORDS:      An array containing the original coordinates of teh nodes of
                 the element
    
    U, DU, V, A: Arrays containing the current estimates of the basic solution
                 variables at the nodes of the element.
                 U  = Total values of the variables.
                 DU = Incremental values of the variables
                 V  = Velocity (defined only for implicit dynamics)
                 A  = Acceleration (defined only for implicit dynamics)
                 
                 U = [u1,u2,u3,phi11,phi22,phi33,phi23,phi13,phi12,phi32,phi31,phi21,...]
                 
    JDLTYP:      An array containing the integers used to define distributed 
                 load types for the element. Loads of type Un are identified 
                 by the integer value n in JDLTYP; loads of type UnNU are 
                 identified by the negative integer value.
                 
    ADLMAG:      For genreal nonlinear steps ADLMAG is the total load magnitude 
                 of the K1th distributed load at the end of the current increment
                 for distributed loads of type Un. For distributed loads of type 
                 UnNU, the load magnitude is defined in UEL.
                 
    DDLMAG:      For general nonlinear steps DDLMAG contains the increments in the 
                 magnitudes of the distributed loads that are currently active on 
                 this element for distributed loads of type Un.
                 
    PREDEF:      An array containing the values of predefined field variables, 
                 such as temperature in an uncoupled stress/displacement analysis, 
                 at the nodes of the element.
                 
                 The first index K1 is either 1 or 2, with 1 indicating the value 
                 of the field variable at the end of the increment and 2 indicating
                 the increment in the field variable. The second index K2, indicates 
                 the variable: the temperature coresponds to index 1, and the 
                 predefined field variables correspond to indices 2 and above. In 
                 cases where temperature is not defined, the predefined field 
                 variables begin with index 1. The third index, K3, indicates the 
                 local node number on the element.
                 
                 PREDEF(K1, 1,K3) = Temperature
                 PREDEF(K1, 2,K3) = First predefined field variable
                 PREDEF(K1, 3,K3) = Second predefined field variable
                 Etc.               Any other predefined field variable
                 PREDEF(K1,K2,K3) = Total or incremental value of the K2th predefined
                                    field variable at the K3th node of the element.
                 PREDEF(1,K2,K3)  = Values of the variables at the end of the current
                                    increment
                 PREDEF(2,K2,K3)  = Incremental values corresponding to the current 
                                    time increment
                                    
    PARAMS:      An array containing the parameters associated with the solution 
                 procedure. The entries in this array depend on the solution procedure 
                 currently being used when UEL is called, as indicated by the entries 
                 in the LFLAGS array.
                 
                 For implicit dynamics (LFLAGS(1) = 11 or 12) PARAMS contains the 
                 integration operator values as
                 
                 PARAMS(1) = alpha
                 PARAMS(2) = beta
                 PARAMS(3) = gamma
                 
    LFLAGS:      An array containing the flags that define the current solution 
                 procedure and requirements for element calculations.
                 
                 LFLAGS(1) = Defines the procedure type
                 LFLAGS(2) = 0 Small-displacement analysis
                 LFLAGS(2) = 1 Large-displacement analysis
                 LFLAGS(3) = 1 Normal implicit time incrementation procedure. User
                               subroutine UEL must define the residual vector in 
                               RHS and the Jacobian matrix in AMATRIX
                 LFLAGS(3) = 2 Define the current stiffness matrix (AMATRX = K^{NM}
                               = -dF^N/du^M or -dG^N/du^M) only.
                 LFLAGS(3) = 3 Define the current damping matrix (AMATRX = C^{NM}
                               = -dF^N/ddotu^M or -dG^N/ddotu^M) only.
                 LFLAGS(3) = 4 Definen the current mass matrix (AMATRX = -dF^N/dddotu^M) only
                 LFLAGS(3) = 5 Define the current residual or load vector (RHS = F^N) only.
                 LFLAGS(3) = 6 Define the current mass matrix and the residual vector for 
                               the initial acceleration calculation (or the calculation of 
                               accelerations after impact).
                 LFLAGS(3) = 100 Define peturbation quantities for output.
                 LFLAGS(4) = 0 The step is a general step
                 LFLAGS(4) = 1 The step is a linear perturbation step
                 LFLAGS(5) = 0 The current approximations to u^M, etc. were based on Newton
                               corrections.
                 LFLAGS(5) = 1 The current approximations were found by extrapolation from
                               the previous increment.
    TIME(1):     Current value of step time
    TIME(2):     Current value of total time
    DTIME:       Time increment
    PERIOD:      Time period of the current step
    NDOFEL:      Number of degrees of freedom in the element
    MLVARX:      Dimensioning parameter used when several displacement or 
                 right-hand-side vectors are used.
    NRHS:        Number of load vetors. NRHS is 1 in most nonlinear problems:
                 it is 2 for the modified Riks static procedure, and it is 
                 greater than 1 in some linear analysis procedures and during 
                 substructure generation.
    NSVARS:      User-defined number of solution-dependent state variables 
                 associated with the element.
    NPROPS:      User-defined number of real property values associated with 
                 the element.
    NJPROP:      User-defined number of integer property values assocated with 
                 the element.
    MCRD:        MCRD is defined as the maximum of the user-defined maximum 
                 number of coordinates needed at any node point and the 
                 value of the largest active degree of freedom of the user
                 element that is less than or equal to 3.
    NNODE:       User-defined number of nodes on the element
    JTYPE:       Integer defining the element type. This is the user-defined
                 integer value n in element type Un.
    KSTEP:       Current step number.
    KINC:        Current increment number.
    JELEM:       User-assigned element number.
    NDLOAD:      Identification number of the distributed load or flux
                 currently active on this element.
    MDLOAD:      Total number of distributed loads and/or flues defined 
                 on this element.
    NPREDF:      Number of predefined field variables, including temperature.
    """
        
    #Set the order of the gauss quadrature and get the gauss points and weights
    ORDER_QUAD = 2
    PQ,WQ      = hex8.get_gpw(ORDER_QUAD)
    UN,PHIN    = parse_dof_vector(U)
    DUN,DPHIN  = parse_dof_vector(DU)
    
    #Parse the LFLAGS array
    if(LFLAGS[2]==1):
        print "Derp"
        
###### Extract degrees of freedom ######
def parse_dof_vector(U): #Test function written
    """Parse the degree of freedom vector"""
    #Extract degrees of freedom
    node_dofs = [U[i:(i+12)] for i in range(0,8*12,12)]
    #Parse the nodal degrees of freedom
    node_us   = [node_dofs[i][:3]   for i in range(8)]
    node_phis = [node_dofs[i][3:12] for i in range(8)]
    return node_us,node_phis
    
###### Compute Fundamental Quantities ######

#TODO: Replace matrix form with vector form

def interpolate_dof(xi_vec,phi_vectors,nodal_global_coords_current,nodal_global_coords_reference): #Test function written
    """Interpolate the degree of freedom vectors to form the required gradients"""
    F        = compute_F(xi_vec,nodal_global_coords_current,nodal_global_coords_reference)
    chi      = compute_chi(xi_vec,phi_vectors)
    grad_chi = compute_grad_chi(xi_vec,phi_vectors,nodal_global_coords_reference)
    return F,chi,grad_chi

def compute_F(xi_vec,nodal_global_coords_current,nodal_global_coords_reference): #Test function written
    """Compute the deformation gradient at xi_loc following Belytschko Section 4, Example 4.3"""
    
    #Compute the gradients of the shape functions w.r.t. the local coordinates
    dNdxis =  [hex8.Hex8_local_grad_shape_function(n,xi_vec) for n in range(8)]

    #Compute the gradients of the global coordinates w.r.t. the local coordinates
    dxdxi    = sum([hex8.vector_dyadic_product(ngcc,dNdxi) for ngcc,dNdxi in zip(nodal_global_coords_current,dNdxis)])
    dXdxi    = sum([hex8.vector_dyadic_product(ngcr,dNdxi) for ngcr,dNdxi in zip(nodal_global_coords_reference,dNdxis)])
    
    #Invert the partial of the global referencec coordinates w.r.t. the local coordindates
    dxidX    = hex8.invert_3x3_matrix(dXdxi)[1]
    
    #Compute and return F
    return hex8.reduce_tensor_to_vector_form(hex8.matrix_dot(dxdxi,dxidX)) #TODO: make this a vector by nature
    
def compute_chi(xi_vec,phi_vectors): #Test function written
    """Compute the value of chi"""
    #Initialize phi
    phi    = np.zeros([9,])
    #Iterate through the nodes
    for i in range(8):
        #Compute the shape function
        N = hex8.Hex8_shape_function(i,xi_vec)
        for j in range(9):
            #Update phi
            phi[j] += N*phi_vectors[i][j]
    #Return the formed value of chi
    return form_chi(phi)
    
def form_chi(phi_vector): #Test function written (contained in compute_chi test)
    """Using the interpolated values of phi form the chi tensor"""
    #Create an empty array
    chi = np.empty([9,1])
    #Populate array
    chi[0] = 1.+phi_vector[0]
    chi[1] = 1.+phi_vector[1]
    chi[2] = 1.+phi_vector[2]
    chi[3] =    phi_vector[3]
    chi[4] =    phi_vector[4]
    chi[5] =    phi_vector[5]
    chi[6] =    phi_vector[6]
    chi[7] =    phi_vector[7]
    chi[8] =    phi_vector[8]
    return hex8.convert_M_to_V(chi,[3,3]) #TODO: Make this a vector by nature
    
def compute_grad_chi(xi_vec,phi_vectors,nodal_global_coords_reference): #Test function written
    """Compute the gradient of chi from the vector of phi values"""
    #Populate an empty matrix
    grad_chi = np.zeros([3,3,3])
    
    #Compute the gradients of the shape functions w.r.t. the local coordinates
    dNdxis =  [hex8.Hex8_local_grad_shape_function(n,xi_vec) for n in range(8)]
    
    #Compute the gradients of the referencec global coordinates w.r.t. the local coordinates
    dXdxi    = sum([hex8.vector_dyadic_product(ngcr,dNdxi) for ngcr,dNdxi in zip(nodal_global_coords_reference,dNdxis)])
    
    #Invert dXdxi to get dxidX
    dxidX  = hex8.invert_3x3_matrix(dXdxi)[1]
    
    #Compute the dNdXes
    dNdxs = [hex8.vector_dot_matrix(dNdxi,dxidX) for dNdxi in dNdxis]
    
    #Compute the gradients
    for n in range(8):
        chin = form_chi(phi_vectors[n])
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    grad_chi[i,j,k] += chin[T2V([i,j],[3,3])]*dNdxs[n][k]
    return hex8.reduce_tensor_to_vector_form(grad_chi) #TODO: Make this a vector by nature
    
###### Compute the derivatives of the fundamental quantities ######
    
def compute_fundamental_derivatives(xi_vec,nodal_global_coords_reference): #Test function written
    """Compute the derivative of the fundamental expressions with respect to the degree of freedom vector"""
    dFdU        = compute_dFdU(xi_vec,nodal_global_coords_reference)
    dchidU      = compute_dchidU(xi_vec)
    dgrad_chidU = compute_dgrad_chidU(xi_vec,nodal_global_coords_reference)
    return dFdU,dchidU,dgrad_chidU
    
def compute_dFdU(xi_vec,nodal_global_coords_reference): #Test function written
    """Compute the derivative of the deformation gradient with respect to the degree of freedom vector"""
    dFdU = np.zeros([9*96])
    for n in range(8):
        dFdUn = compute_dFndUn(n,xi_vec,nodal_global_coords_reference)
        for i in range(3):
            for I in range(3):
                for J in range(12):
                    dFdU[T2V([i,I,J+12*n],[3,3,96])] = dFdUn[T2V([i,I,J],[3,3,12])]
    return dFdU
    
def compute_dFndUn(n,xi_vec,nodal_global_coords_reference): #Test function written (part of test_compute_dFdU)
    """Compute the nth matrix in dFdU"""
    dFdUn   = np.zeros([9*12])
    dNdX    = hex8.Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords_reference)[0]
    Itmp    = np.reshape([1,0,0,0,0,0,0,0,0,0,0,0,\
                          0,1,0,0,0,0,0,0,0,0,0,0,\
                          0,0,1,0,0,0,0,0,0,0,0,0,],[3,12]).astype(float)
    
    for i in range(3):
        for I in range(3):
            for J in range(12):
                dFdUn[T2V((i,I,J),[3,3,12])] = Itmp[i,J]*dNdX[I]
    return dFdUn
    
def compute_dchidU(xi_vec): #Test function written
    """Compute the derivative of the micro-displacement tensor with respect to the degree of freedom vector"""
    X = np.zeros([9*96])
    for n in range(8):
        dchidUn = compute_dchidUn(n,xi_vec)
        for i in range(3):
            for I in range(3):
                for J in range(12):
                    X[T2V([i,I,J+12*n],[3,3,96])] = dchidUn[T2V([i,I,J],[3,3,12])]
    return X
    
def compute_dchidUn(n,xi_vec): #Test function written (part of compute_dchidU)
    """Compute the nth matrix in dchidU"""
    dchidUn = np.zeros([9*12])
    N = hex8.Hex8_shape_function(n,xi_vec)
    Itmp    = np.reshape([0,0,0,1,0,0,0,0,0,0,0,0,\
                          0,0,0,0,0,0,0,0,1,0,0,0,\
                          0,0,0,0,0,0,0,1,0,0,0,0,\
                          0,0,0,0,0,0,0,0,0,0,0,1,\
                          0,0,0,0,1,0,0,0,0,0,0,0,\
                          0,0,0,0,0,0,1,0,0,0,0,0,\
                          0,0,0,0,0,0,0,0,0,0,1,0,\
                          0,0,0,0,0,0,0,0,0,1,0,0,\
                          0,0,0,0,0,1,0,0,0,0,0,0],[3*3*12]).astype(float)
    for i in range(3):
        for I in range(3):
            for J in range(12):
                iIJ = T2V([i,I,J],[3,3,12])
                dchidUn[iIJ] = N*Itmp[iIJ]
    return dchidUn
    
def compute_dgrad_chidU(xi_vec,nodal_global_coords_reference): #Test function written
    """Compute the derivative of the gradient of the micro-displacement tensor with respect to the degree of freedom vector"""
    dgrad_chidU = np.zeros([27*96])
    for n in range(8):
        dgrad_chidUn = compute_dgrad_chidUn(n,xi_vec,nodal_global_coords_reference)
        for i in range(3):
            for I in range(3):
                for J in range(3):
                    for K in range(12):
                        dgrad_chidU[T2V([i,I,J,K+n*12],[3,3,3,96])] = dgrad_chidUn[T2V([i,I,J,K],[3,3,3,12])]
    return dgrad_chidU
    
def compute_dgrad_chidUn(n,xi_vec,nodal_global_coords_reference): #Test function written (part of compute_dgrad_chidU)
    """Compute the nth matrix in dgrad_chindU"""
    dgrad_chidUn = np.zeros([27*12])
    dNdX = hex8.Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords_reference)[0]
    Itmp    = np.reshape([0,0,0,1,0,0,0,0,0,0,0,0,\
                          0,0,0,0,0,0,0,0,1,0,0,0,\
                          0,0,0,0,0,0,0,1,0,0,0,0,\
                          0,0,0,0,0,0,0,0,0,0,0,1,\
                          0,0,0,0,1,0,0,0,0,0,0,0,\
                          0,0,0,0,0,0,1,0,0,0,0,0,\
                          0,0,0,0,0,0,0,0,0,0,1,0,\
                          0,0,0,0,0,0,0,0,0,1,0,0,\
                          0,0,0,0,0,1,0,0,0,0,0,0],[3*3*12]).astype(float)
    for i in range(3):
        for I in range(3):
            for J in range(3):
                for K in range(12):
                    iIJK = T2V([i,I,J,K],[3,3,3,12])
                    iIK  = T2V([i,I,K],[3,3,12])
                    dgrad_chidUn[iIJK] = Itmp[iIK]*dNdX[J]
    return dgrad_chidUn
    
###### Compute Deformation Measures ######
                    
def get_deformation_measures(F,chi,grad_chi): #Test function written
    """Compute and return the deformation measures"""
    C   = hex8.matrix_Tdot_V(F,F)
    Phi = hex8.matrix_Tdot_V(F,chi)
    
    Gamma = np.zeros([27,])
    for I in range(3):
        for J in range(3):
            for K in range(3):
                index = T2V([I,J,K],[3,3,3]) #Identify the index of the Gamma vector
                for i in range(3):
                    Findx = T2V([i,I],[3,3])      #Identify the index of the F vector
                    Gindx = T2V([i,J,K],[3,3,3]) #Identify the index of the grad_chi vector
                    Gamma[index] += F[Findx]*grad_chi[Gindx]
    #Gamma = hex8.matrix_Tdot_TOT(F,grad_chi)
    return C,Phi,Gamma
    
###### Compute Derivatives of Deformation Measures ######
def compute_DM_derivatives(F,chi,grad_chi,dFdU,dchidU,dgrad_chidU): #Test function written
    """Compute the derivatives of the deformation measures"""
    dCdU     = compute_dCdU(F,dFdU)
    dPsidU   = compute_dPsidU(F,chi,dFdU,dchidU)
    dGammadU = compute_dGammadU(F,grad_chi,dFdU,dgrad_chidU)
    return dCdU,dPsidU,dGammadU
    
def compute_dCdU(F,dFdU): #Test function written
    """Compute the derivative of the right Cauchy-Green Deformation Tensor"""
    dCdU = np.zeros([9*96])
    for n in range(8):
        dCdUn = compute_dCdUn(n,F,dFdU)
        for I in range(3):
            for J in range(3):
                for K in range(12):
                    dCdU[T2V([I,J,K+12*n],[3,3,96])] = dCdUn[T2V([I,J,K],[3,3,12])]
    return dCdU

def compute_dCdUn(n,F,dFdU): #Test function written (part of test_compute_dCdU)
    """Compute a submatrix of the derivative of the right Cauchy-Green Deformation Tensor
    n goes from 0 to 7
    """
    
    dCdUn = np.zeros([9*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(12):
                index = T2V([I,J,K],[3,3,12])
                for i in range(3):
                    
                    iJ  = T2V([i,J],[3,3])
                    iI  = T2V([i,I],[3,3])
                    iIK = T2V([i,I,K+12*n],[3,3,96])
                    iJK = T2V([i,J,K+12*n],[3,3,96])
                    
                    dCdUn[index] += F[iJ]*dFdU[iIK] + F[iI]*dFdU[iJK]
    return dCdUn
    
def compute_dCinvdU(Cinv,dCdU): #Test function written
    """Compute the derivative of the inverse of the 
    right Cauchy-Green deformation tensor with respect 
    to the degree of freedom vector"""
    dCinvdU = np.zeros([9*96])
    
    for n in range(8):
        dCinvdUn = compute_dCinvdUn(n,Cinv,dCdU)
        for I in range(3):
            for J in range(3):
                for K in range(12):
                    dCinvdU[T2V([I,J,K+12*n],[3,3,96])] = dCinvdUn[T2V([I,J,K],[3,3,12])]
    return dCinvdU
    
def compute_dCinvdUn(n,Cinv,dCdU): #Test function written (part of compute_dCinvdU)
    """Compute a submatrix of the derivative of the inverse 
    of the right Cauchy-Green deformation tensor with respect 
    to the degree of freedom vector. n goes from 0 to 7"""
    #Initialize matrix
    dCinvdUn = np.zeros([9*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(12):
                IJK = T2V([I,J,K],[3,3,12])
                
                for M in range(3):
                    for L in range(3):
                        IL  = T2V([I,L],[3,3])
                        MJ  = T2V([M,J],[3,3])
                        LMK = T2V([L,M,K+12*n],[3,3,96])
                        
                        dCinvdUn[IJK] += -Cinv[IL]*Cinv[MJ]*dCdU[LMK]
    return dCinvdUn
    
def compute_dPsidU(F,chi,dFdU,dchidU): #Test function written
    """Compute the derivative of the deformation measure
    Psi with respect to the deformation gradient"""
    dPsidU = np.zeros([9*96])
    
    for n in range(8):
        dPsidUn = compute_dPsidUn(n,F,chi,dFdU,dchidU)
        for I in range(3):
            for J in range(3):
                for K in range(12):
                    dPsidU[T2V([I,J,K+n*12],[3,3,96])] = dPsidUn[T2V([I,J,K],[3,3,12])]
    return dPsidU
    
def compute_dPsidUn(n,F,chi,dFdU,dchidU): #Test function written (part of test_compute_dPsidU)
    """Compute a submatrix of the derivative of the deformation measure Phi with 
    respect to the dof vector n goes from 0 to 7
    """
    dPsidUn = np.zeros([9*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(12):
                IJK = T2V([I,J,K],[3,3,12])
                
                for i in range(3):
                    iJ  = T2V([i,J],[3,3])
                    iIK = T2V([i,I,K+n*12],[3,3,96])
                    iI  = T2V([i,I],[3,3])
                    iJK = T2V([i,J,K+n*12],[3,3,96])
                    
                    dPsidUn[IJK] += chi[iJ]*dFdU[iIK] + F[iI]*dchidU[iJK]
    return dPsidUn
    
def compute_dGammadU(F,grad_chi,dFdU,dgrad_chidU): #Test function written
    """Compute the derivative of the deformation measure
    Gamma with respect to the deformation gradient"""
    dGammadU = np.zeros([27*96])
    
    for n in range(8):
        dGammadUn = compute_dGammadUn(n,F,grad_chi,dFdU,dgrad_chidU)
        for I in range(3):
            for J in range(3):
                for L in range(3):
                    for K in range(12):
                        dGammadU[T2V([I,J,L,K+n*12],[3,3,3,96])] = dGammadUn[T2V([I,J,L,K],[3,3,3,12])]
    return dGammadU
    
def compute_dGammadUn(n,F,grad_chi,dFdU,dgrad_chidU): #Test function written (contained in test_compute_dGammadU)
    """Compute a submatrix of the derivative of the
    deformation measure Gamma with respect to the
    degree of freedom vector. n goes from 0 to 7"""
 
    #Initialize the array
    dGammadUn = np.zeros([27*12])
    
    for I in range(3):
        for J in range(3):
            for L in range(3):
                for K in range(12):
                    IJLK = T2V([I,J,L,K],[3,3,3,12])
                    
                    for i in range(3):
                        iJL  = T2V([i,J,L],          [3,3,3])
                        iIK  = T2V([i,I,K+12*n],    [3,3,96])
                        iI   = T2V([i,I],              [3,3])
                        iJLK = T2V([i,J,L,K+12*n],[3,3,3,96])
                    
                        dGammadUn[IJLK] += grad_chi[iJL]*dFdU[iIK]+F[iI]*dgrad_chidU[iJLK]
    return dGammadUn    
    
###### Residual and Residual Gradient Calculations ######
    
def compute_residuals_gpt(xi_vec,node_us,node_phis,nodal_global_coords_current,nodal_global_coords_reference,params,state_variables):
    """Compute the residuals at a given gauss point"""
    F,chi,grad_chi      = interpolate_dof(xi_vec,node_phis,nodal_global_coords_current,nodal_global_coords_reference)
    N,grad_N_ref,detJ   = hex8.Hex8_get_shape_function_info(n,xi_vec,nodal_global_coords)
    PK2,SIGMA,M         = compute_stress(F,chi,grad_chi,params,state_variables)
    
    RBLM = compute_BLM_residual(N,grad_N_ref,detJ,PK2,SIGMA,M,RHO0)
    RMOM = compute_FMOM_residual(N,grad_N_ref,detJ,PK2,SIGMA,M,RHO0)
    
    return np.vstack([RBLM,RMOM])
    
def compute_BLM_residual_gpt(N,F,grad_N_ref,detJ,PK2,TRACTION,SIGMA,M,RHO0,ACCEL,BODYF):
    """Compute the residual of the balance of linear momentum"""
    RBLM = compute_resid_fint_fkin_bf(N,F,grad_N_ref,detJ,PK2,SIGMA,M,RHO0,ACCEL,BODYF)
    RBLM += compute_BLM_resid_tracation(N,F,grad_N_ref,detJ,TRACTION,SIGMA,M,RHO0,ACCEL,BODYF)
    return RBLM
    
def compute_resid_fint_fkin_bf(N,F,grad_N_ref,detJ,PK2,SIGMA,M,RHO0,ACCEL,BODYF):
    """Compute the residual resulting from the kinematic, body, and internal forces"""
    RBLM = np.zeros([3,])
    for j in range(3):
        RBLM[j] = N*RHO0*(BODYF[j]-ACCEL[j])*detJ
        for I in range(3):
            for J in range(3):
                RBLM[j] -= grad_N_ref[I]*PK2[I,J]*F[j,J]*detJ
    
    return RBLM
    
def compute_BLM_resid_traction(N,F,grad_N_ref,detJ,TRACTION,SIGMA,M,RHO0,ACCEL,BODYF):
    """Compute the residual from the applied traction"""
    ORDER = 2
    GP,W = hex8.get_face_gpw(ORDER)
    
    RES = np.zeros([3,])
    
    #TODO
    
    return RES
    
def compute_dBLMdU():
    """Compute the derivative of the balance of linear momentum with respect to the degree of freedom vector"""
    #TODO
    
    
def compute_FMOM_residual_gpt(N,F,grad_N_ref,detJ,PK2,SIGMA,M,RHO0,GYRATION,BODYCOUPLE):
    """Compute the residual of the balance of first moment of momentum"""
    RFMOM = compute_resid_Sint_Skin_bC(N,F,grad_N_ref,detJ,PK2,SIGMA,M,RHO0,GYRATION,BODYCOUPLE)
    RFMOM += compute_FMOM_resid_traction(N,F,grad_N_ref,detJ,COUPLE_TRACTION,M)
    return RFMOM
    
def compute_resid_Sint_Skin_bC(N,F,grad_N_ref,detJ,PK2,SIGMA,M,RHO0,GYRATION,BODYCOUPLE):
    """Compute the residual resulting from the kinematic, body, and internal forces"""
    RES = np.zeros([3,3])
    
    for i in range(3):
        for j in range(3):
            RES[i,j] = RHO0*(BODYCOUPLE[j,i]-GYRATION[j,i])*N*detJ
            
            for I in range(3):
                for J in range(3):
                    RES[i,j] += N*(F[i,I]*(PK2[I,J] - SIGMA[I,J])*F[j,J])*detJ
                    for K in range(3):
                        RES[i,j] -= grad_N_ref[K]*F[j,J]*chi[i,I]*M[K,J,I]
    return RES
    
def compute_FMOM_resid_traction(N,F,grad_N_ref,detJ,COUPLE_TRACTION,SIGMA,M,RHO0,ACCEL,BODYF):
    """Compute the residual from the applied traction"""
    ORDER = 2
    GP,W = hex8.get_face_gpw(ORDER)
    
    RES = np.zeros([3,])
    
    #TODO
    
    return RES
    
def compute_dFMOMdU():
    """Compute the derivative of the balance of first moment of momentum with respect to the degree of freedom vector"""
    #TODO
    
###### Stress/Strain Calculations ######
    
def compute_stress(F,chi,grad_chi,params,state_variables):
    """Compute the current value of the stress quantities. Returns the stress and required derivatives"""
    MODE = 1
    if MODE==1:
        return micro_LE.micromorphic_linear_elasticity(F,chi,grad_chi,params)
            
###### Compute Tangents ######
    
def compute_dpk2dU(dpk2dC,dpk2dPsi,dpk2dGamma,dCdU,dPsidU,dGammadU):
    """Compute the derivative of the second piola kirchhoff stress
    with respect to the degree of freedom vector"""
    
    dpk2dU = np.zeros([9*96])
    
    for n in range(8):
        dpk2dUn = compute_dpk2dUn(n,dpk2dC,dpk2dPsi,dpk2dGamma,dCdU,dPsidU,dGammadU)
        for I in range(3):
            for J in range(3):
                for K in range(12):
                    dpk2dU[T2V([I,J,K+n*12],[3,3,96])] = dpk2dUn[T2V([I,J,K],[3,3,12])]
                
    return dpk2dU
    
def compute_dpk2dUn(n,dpk2dC,dpk2dPsi,dpk2dGamma,dCdU,dPsidU,dGammadU):
    """Compute a submatrix of the derivative of the second piola kirchhoff stress
    with respect to the degree of freedom vector"""
    
    dpk2dUn = np.zeros([9*12,])
    
    for I in range(3):
        for J in range(3):
            for K in range(12):
                
                for L in range(3):
                    for M in range(3):
                        IJLM = T2V([I,J,L,M],[3,3,3,3])
                        LMK  = T2V([L,M,K+n*12],[3,3,96])
                        dpk2dU[IJK] += dpk2dC[IJLM]*dCdU[LMK] + dpk2dPsi[IJLM]*dPsidU[LMK]
                        
                for L in range(3):
                    for M in range(3):
                        for N in range(3):
                            IJLMN = T2V([I,J,L,M,N],[3,3,3,3,3])
                            LMNK  = T2V([L,M,N,K+n*12],[3,3,3,96])
                            dpk2dU[IJK] += dpk2dGamma[IJLMN]*dGammadU[LMNK]
    return dpk2dU
    
def compute_dsymmetric_stressdU(dSigmadC,dSigmadPsi,dSigmadGamma,dCdU,dPsidU,dGammadU):
    """Compute the derivative of the symmetric stress with respect 
    to the degree of freedom vector"""
    
    #The calculations for the tangent for the pk2 stress and the symmetric stress are identical
    #therefore
    
    return compute_dpk2dU(dSigmadC,dSigmadPsi,dSigmadGamma,dCdU,dPsidU,dGammadU)
    
def compute_dho_stressdU(dMdC,dMdPhi,dMdGamma,dCdU,dPsidU,dGammadU):
    """Compute the derivative of the symmetric stress with respect
    to the degree of freedom vector"""
    
    dMdU = np.zeros([3*3*3*96])
    
    for n in range(8):
        dMdUn = compute_dho_stressdUn(dMdC,dMdPsi,dMdGamma,dCdU,dPsidU,dGammadU)
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    for L in range(12):
                        dMdU[T2V([I,J,K,L+n*12],[3,3,3,96])] = dMdUn[T2V([I,J,K,L],[3,3,3,12])]
    return dMdU
      
def compute_dho_stressdUn(dMdC,dMdPsi,dMdGamma,dCdU,dPsidU,dGammadU):
    """Compute a submatrix of the derivative of the symmetric stress with respect
    to the degree of freedom vector"""
    
    dMdUn = np.zeros([3*3*3*12])
    for I in range(3):
        for J in range(3):
            for K in range(3):
                for L in range(12):
                    IJKL  = T2V([I,J,K,L],[3,3,3,12])
                    
                    for O in range(3):
                        for P in range(3):
                            IJKOP = T2V([I,J,K,O,P],[3,3,3,3,3])
                            OPL   = T2V([O,P,P+n*12],[3,3,96])
                            
                            dMdUn[IJKL] += dMdC[IJKOP]*dCdU[OPL] + dMdPsi[IJKOP]*dPsidU[OPL]
                            
                    for O in range(3):
                        for P in range(3):
                            for Q in range(3):
                                IJKOPQ = T2V([I,J,K,O,P,Q],[3,3,3,3,3,3])
                                OPQL   = T2V([O,P,Q,L+n*12],[3,3,3,96])
                                
                                dMdUn[IJKL] += dMdGamma[IJKOPQ]*dGammadU[OPQL]
    return dMdUn
    
class TestMicroElement(unittest.TestCase):

    f                    = None
    original_directory   = ""
    output_file_name     = r"micro_element_unittests.txt"
    output_file_location = r".\tests\unittests"
    currentResult        = None
    @classmethod
    def setUpClass(self):
        """Setup method"""
        output_file = os.path.join(self.output_file_location,self.output_file_name)
        if(os.path.isfile(output_file)):
            os.remove(output_file)
        self.f = open(output_file,"w+")
    @classmethod
    def tearDownClass(self):
        """Teardown method"""
        self.f.close()
        
    def setUp(self):
        pass
        
    def tearDown(self):
        ok = self.currentResult.wasSuccessful()
        tname = self.id().split(".")[-1]
        self.f.write(tname+"\t&\t"+str(ok)+"\n")
        
    def run(self, result=None):
        """Redefine run to keep track of results"""
        self.currentResult = result
        unittest.TestCase.run(self,result)

    def test_interpolate_dof(self):
        """Test the interpolate dof function"""
        #Get the analytic deformation gradient, the reference coordinates, and the current coordinates
        xi_vec = [0.,0.,0.]
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_loc = (0.,0.,0.)
                   
        Fanalytic,ccoords      = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vecs = self._get_chi_values(rcoords,[X,Y,Z])
        
        F,chi,grad_chi = interpolate_dof(xi_vec,phi_vecs,ccoords,rcoords)
        
        self.assertEqual(np.allclose(F,hex8.convert_M_to_V(hex8.reduce_tensor_to_matrix_form(Fanalytic),[3,3])),True)
        self.assertEqual(np.allclose(chi,hex8.convert_M_to_V(hex8.reduce_tensor_to_matrix_form(chia),[3,3])),True)
        self.assertEqual(np.allclose(grad_chi,hex8.reduce_tensor_to_vector_form(grad_chia)),True)
        
    def test_compute_F(self):
        """Test the computation of the deformation gradient"""
        
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_loc = (0.,0.,0.)
        
        #Get the analytic deformation gradient and the current coordinates
        Fanalytic,ccoords = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        
        F = compute_F(xi_loc,ccoords,rcoords)
            
        self.assertEqual(np.allclose(F,hex8.reduce_tensor_to_vector_form(Fanalytic)),True)
        
    def test_compute_chi(self):
        """Test compute_chi to take values of phi located at the nodes, interpolate, and assemble into chi"""
        
        #Define node coordinates
        rcoords = [[-2.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,.5,-1.],\
                   [-1.,-1., 1.],[1.,-7., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.zeros([3,])
        
        chia,grad_chi,phi_vecs = self._get_chi_values(rcoords,[X,Y,Z])
        chi = compute_chi(xi_vec,phi_vecs)
        
        self.assertEqual(np.allclose(chi,hex8.reduce_tensor_to_vector_form(chia)),True)
        
    def test_compute_grad_chi(self):
        """Test compute_grad_chi to take values of phi located at the nodes and compute their gradient at a location"""
        
        #Define node coordinates
        rcoords = [[-2.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,.5,-1.],\
                   [-1.,-1., 1.],[1.,-7., 1.],[1.,1., 1.],[-1.,1., 1.]]
                   
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.zeros([3,])
        
        chia,grad_chia,phi_vecs = self._get_chi_values(rcoords,[X,Y,Z])
        
        grad_chi = compute_grad_chi(xi_vec,phi_vecs,rcoords)
        
        self.assertEqual(np.allclose(grad_chi,hex8.reduce_tensor_to_vector_form(grad_chia)),True)
        
    def test_parse_dof_vector(self):
        """Test of the parsing of the DOF vector"""
        #Construct U vector
        DOF = [[["N"+str(n)+"U"+str(i+1) for i in range(3)],["N"+str(n)+"PHI"+str(i+1) for i in range(9)]] for n in range(8)]
        UA,PHIA = zip(*DOF)
        U = self._flatten_list(self._flatten_list(DOF))
        U,PHI = parse_dof_vector(U)
        
        resU   = sum([sum([Un[i]==UAn[i] for i in range(len(Un))]) for Un,UAn in zip(U,UA)])
        resPHI = sum([sum([PHIn[i]==PHIAn[i] for i in range(len(PHIn))]) for PHIn,PHIAn in zip(PHI,PHIA)])
        
        self.assertEqual(resU,3*8)
        self.assertEqual(resPHI,9*8)
        
    def test_get_deformation_measures(self):
        """Test get_deformation_measures to compute C, Psi, Gamma"""
        
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.,0.,0.])
        
                   
        Fanalytic,ccoords        = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute the fundamental values
        F,chi,grad_chi = interpolate_dof(xi_vec,phi_vectors,ccoords,rcoords)
        
        C,Psi,Gamma = get_deformation_measures(F,chi,grad_chi)
        
        C     = hex8.convert_V_to_T(C,      [3,3])
        Psi   = hex8.convert_V_to_T(Psi,    [3,3])
        Gamma = hex8.convert_V_to_T(Gamma,[3,3,3])
        
        F        = hex8.convert_V_to_T(F,         [3,3])
        chi      = hex8.convert_V_to_T(chi,       [3,3])
        grad_chi = hex8.convert_V_to_T(grad_chi,[3,3,3])
        
        CA     = np.einsum('iI,iJ',F,F)
        PsiA   = np.einsum('iI,iJ',F,chi)
        GammaA = np.einsum('iI,iJK',F,grad_chi)
        
        self.assertEqual(np.allclose(    C,     CA), True)
        self.assertEqual(np.allclose(  Psi,   PsiA), True)
        self.assertEqual(np.allclose(Gamma, GammaA), True)
        
    def test_compute_dFdU(self):
        """Test the computation of the matrix form of the derivative of the deformation gradient 
        with respect to the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def F_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            node_xs = [[u1+rc1,u2+rc2,u3+rc3] for (u1,u2,u3),(rc1,rc2,rc3) in zip(node_us,rcoords)]
            F = compute_F(xi_vec,node_xs,rcoords)
            F = hex8.convert_V_to_M(F,[3,3])
            return np.reshape(F,[9,])
            
        #Compute the gradients
        dFdUn = fd.numeric_gradient(F_parser,U,1e-6)
        dFdU  = compute_dFdU(xi_vec,rcoords)
        
        #tmp = 4
        #print "Numeric"
        #print dFdUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print hex8.convert_V_to_M(dFdU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dFdUn.T[:,(12*tmp):12*(tmp+1)] - hex8.convert_V_to_M(dFdU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(hex8.convert_M_to_V(dFdUn.T,[3,3,96]),dFdU),True)
        
    def test_compute_dCdU(self):
        """Test the computation of the matrix form of the derivative of the 
        right Cauchy-Green deformation tensor with respect to the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def C_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            node_xs = [[u1+rc1,u2+rc2,u3+rc3] for (u1,u2,u3),(rc1,rc2,rc3) in zip(node_us,rcoords)]
            F = compute_F(xi_vec,node_xs,rcoords)
            C = hex8.matrix_Tdot_V(F,F)
            C = hex8.convert_V_to_M(C,[3,3])
            return np.reshape(C,[9,])
        
        #Compute the gradients
        dCdUn = fd.numeric_gradient(C_parser,U,1e-6)
        F     = compute_F(xi_vec,ccoords,rcoords)
        dFdU  = compute_dFdU(xi_vec,rcoords)
        dCdU  = compute_dCdU(F,dFdU)
        
        #tmp = 0
        #print "Numeric"
        #print dCdUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print hex8.convert_V_to_M(dCdU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dCdUn.T[:,(12*tmp):12*(tmp+1)] - hex8.convert_V_to_M(dCdU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dCdUn.T,hex8.convert_V_to_M(dCdU,[3,3,96])),True)
        
    def test_compute_dCinvdU(self):
        """Test the computation of the matrix form of the derivative of the 
        inverse of the right Cauchy-Green deformation tensor with respect to 
        the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def Cinv_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            node_xs = [[u1+rc1,u2+rc2,u3+rc3] for (u1,u2,u3),(rc1,rc2,rc3) in zip(node_us,rcoords)]
            F = compute_F(xi_vec,node_xs,rcoords)
            Cinv = hex8.invert_3x3_matrix_V(hex8.matrix_Tdot_V(F,F))[1]
            Cinv = hex8.convert_V_to_M(Cinv,[3,3])
            return np.reshape(Cinv,[9,])
            
        #Compute the gradients
        dCinvdUn = fd.numeric_gradient(Cinv_parser,U,1e-6)
        F        = compute_F(xi_vec,ccoords,rcoords)
        Cinv     = hex8.invert_3x3_matrix_V(hex8.matrix_Tdot_V(F,F))[1]
        dFdU     = compute_dFdU(xi_vec,rcoords)
        dCdU     = compute_dCdU(F,dFdU)
        dCinvdU  = compute_dCinvdU(Cinv,dCdU)
        
        #tmp = 0
        #print "Numeric"
        #print dCinvdUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print hex8.convert_V_to_M(dCinvdU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dCinvdUn.T[:,(12*tmp):12*(tmp+1)] - hex8.convert_V_to_M(dCinvdU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dCinvdUn.T,hex8.convert_V_to_M(dCinvdU,[3,3,96])),True)
        
    def test_compute_dchidU(self):
        """Test the computation of the matrix form of the derivative of the micro-displacement 
        with respect to the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def chi_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            chi = compute_chi(xi_vec,node_phis)
            chi = hex8.convert_V_to_M(chi,[3,3])
            return np.reshape(chi,[9,])
            
        #Compute the gradients
        dchidUn = fd.numeric_gradient(chi_parser,U,1e-6)
        dchidU  = compute_dchidU(xi_vec)
        
        #tmp = 0
        #print "Numeric"
        #print dchidUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print hex8.convert_V_to_M(dchidU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dchidUn.T[:,(12*tmp):12*(tmp+1)] - hex8.convert_V_to_M(dchidU,[3,3,96])[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dchidUn.T,hex8.convert_V_to_M(dchidU,[3,3,96])),True)
        
    def test_compute_dPsidU(self):
        """Test the computation of the matrix form of the derivative of the
        deformation measure Psi with respect to the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def Phi_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            node_xs = [[u1+rc1,u2+rc2,u3+rc3] for (u1,u2,u3),(rc1,rc2,rc3) in zip(node_us,rcoords)]
            F       = compute_F(xi_vec,node_xs,rcoords)
            chi     = compute_chi(xi_vec,node_phis)
            Phi     = hex8.matrix_Tdot_V(F,chi)
            Phi     = hex8.convert_V_to_M(Phi,[3,3])
            return np.reshape(Phi,[9,])
            
        #Compute required measures
        F      = compute_F(xi_vec,ccoords,rcoords)
        dFdU   = compute_dFdU(xi_vec,rcoords)
        chi    = compute_chi(xi_vec,phi_vectors)
        dchidU = compute_dchidU(xi_vec)
            
        #Compute the gradients
        dPsidUn = fd.numeric_gradient(Phi_parser,U,1e-6)
        dPsidU  = compute_dPsidU(F,chi,dFdU,dchidU)
        
        self.assertEqual(np.allclose(dPsidUn.T,hex8.convert_V_to_M(dPsidU,[3,3,96])),True)
        
    def test_compute_dgrad_chidU(self):
        """Test the computation of the matrix form of the derivative of the 
        gradient of the micro-displacement with respect to the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def grad_chi_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            grad_chi = compute_grad_chi(xi_vec,node_phis,rcoords)
            grad_chi = hex8.convert_V_to_M(grad_chi,[3,3,3])
            grad_chi = np.vstack([grad_chi[:,0],grad_chi[:,1],grad_chi[:,2]])
            return np.reshape(grad_chi,[27,])
            
        #Compute the gradients
        dgrad_chidUn = fd.numeric_gradient(grad_chi_parser,U,1e-6)
        dgrad_chidU  = compute_dgrad_chidU(xi_vec,rcoords)
        
        #tmp = 0
        #print "Numeric"
        #print dgrad_chidUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print hex8.convert_V_to_M(dgrad_chidU,[3,3,3,96])[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dgrad_chidUn.T[:,(12*tmp):12*(tmp+1)] - hex8.convert_V_to_M(dgrad_chidU,[3,3,3,96])[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dgrad_chidUn.T,hex8.convert_V_to_M(dgrad_chidU,[3,3,3,96])),True)
        
    def test_compute_dGammadU(self):
        """Test the computation of the matrix form of the derivative of the 
        degormation measure Gamma with respect to the degree of freedom vector"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Define a helper function for the gradient calculation
        def Gamma_parser(Uin):
            """Function which parses U for the computation of F"""
            node_us,node_phis = parse_dof_vector(Uin)
            node_xs  = [[u1+rc1,u2+rc2,u3+rc3] for (u1,u2,u3),(rc1,rc2,rc3) in zip(node_us,rcoords)]
            F        = compute_F(xi_vec,node_xs,rcoords)
            grad_chi = compute_grad_chi(xi_vec,node_phis,rcoords)
            
            Gamma = np.zeros([27,])
            for I in range(3):
                for J in range(3):
                    for K in range(3):
                        index = T2V([I,J,K],[3,3,3]) #Identify the index of the Gamma vector
                        for i in range(3):
                            Findx = T2V([i,I],[3,3])      #Identify the index of the F vector
                            Gindx = T2V([i,J,K],[3,3,3]) #Identify the index of the grad_chi vector
                            Gamma[index] += F[Findx]*grad_chi[Gindx]
            
            Gamma = hex8.convert_V_to_M(Gamma,[3,3,3])
            Gamma = np.vstack([Gamma[:,0],Gamma[:,1],Gamma[:,2]])
            return np.reshape(Gamma,[27,])
        
        #Compute required measures
        F           = compute_F(xi_vec,ccoords,rcoords)
        dFdU        = compute_dFdU(xi_vec,rcoords)
        grad_chi    = compute_grad_chi(xi_vec,phi_vectors,rcoords)
        dgrad_chidU = compute_dgrad_chidU(xi_vec,rcoords)
        
        #Compute the gradients
        dGammadUn = fd.numeric_gradient(Gamma_parser,U,1e-6)
        dGammadU  = compute_dGammadU(F,grad_chi,dFdU,dgrad_chidU)
        
        #tmp = 0
        #print "Numeric"
        #print dGammadUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dGammadU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dGammadUn.T[:,(12*tmp):12*(tmp+1)] - dGammadU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dGammadUn.T,hex8.convert_V_to_M(dGammadU,[3,3,3,96])),True)
        
    def test_compute_fundamental_derivatives(self):
        """Test the function compute_fundamental_derivatives"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        #Compute the derivatives separately
        dFdUt        = compute_dFdU(xi_vec,rcoords)
        dchidUt      = compute_dchidU(xi_vec)
        dgrad_chidUt = compute_dgrad_chidU(xi_vec,rcoords)
        
        #Compute the derivatives together
        dFdU,dchidU,dgrad_chidU = compute_fundamental_derivatives(xi_vec,rcoords)
        
        self.assertEqual(np.allclose(dFdUt,dFdU),True)
        self.assertEqual(np.allclose(dchidUt,dchidU),True)
        self.assertEqual(np.allclose(dgrad_chidUt,dgrad_chidU),True)
        
    def test_compute_DM_derivatives(self):
        """Test the computation of the deformation measure derivatives"""
        #Define the node coordinates
        rcoords = [[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],\
                   [-1.,-1., 1.],[1.,-1., 1.],[1.,1., 1.],[-1.,1., 1.]]
        #Identify a point
        Xs,Ys,Zs = zip(*rcoords)
        X=sum(Xs)/len(Xs)
        Y=sum(Ys)/len(Ys)
        Z=sum(Zs)/len(Zs)
        xi_vec = np.array([0.1,-0.27,0.3])
        
        #Get quantities of interest
        Fanalytic,ccoords          = self._get_deformation_gradient_values(rcoords,[X,Y,Z])
        chia,grad_chia,phi_vectors = self._get_chi_values(rcoords,[X,Y,Z])
        
        #Compute u
        u_vecs = [[cc1-rc1,cc2-rc2,cc3-rc3] for (cc1,cc2,cc3),(rc1,rc2,rc3) in zip(ccoords,rcoords)]
        #Create the dof vector
        U = np.concatenate([np.concatenate((u_vec,phi_vec)) for u_vec,phi_vec in zip(u_vecs,phi_vectors)])
        
        F,chi,grad_chi          = interpolate_dof(xi_vec,phi_vectors,ccoords,rcoords)
        dFdU,dchidU,dgrad_chidU = compute_fundamental_derivatives(xi_vec,rcoords)
        dCdU,dPsidU,dGammadU = compute_DM_derivatives(F,chi,grad_chi,dFdU,dchidU,dgrad_chidU)
        
        dCdUt     = compute_dCdU(F,dFdU)
        dPsidUt   = compute_dPsidU(F,chi,dFdU,dchidU)
        dGammadUt = compute_dGammadU(F,grad_chi,dFdU,dgrad_chidU)
        
        self.assertEqual(np.allclose(dCdUt,        dCdU),True)
        self.assertEqual(np.allclose(dPsidUt,    dPsidU),True)
        self.assertEqual(np.allclose(dGammadUt,dGammadU),True)
        
    def _get_deformation_gradient_values(self,rcoords,X_vec):
        """Get the values required to compute the deformation gradient for testing"""
        #Define the deformation gradient
        x    = lambda X,Y,Z: 1.300*X-0.375*Y+1.2*Z
        dxdX = lambda X,Y,Z:  1.300
        dxdY = lambda X,Y,Z: -0.375
        dxdZ = lambda X,Y,Z: 1.2
        
        y    = lambda X,Y,Z: 0.75*X + 0.650*Y - .31*Z
        dydX = lambda X,Y,Z: 0.75
        dydY = lambda X,Y,Z: 0.650
        dydZ = lambda X,Y,Z: -0.31
        
        z    = lambda X,Y,Z: -2.3*X + 1.4*Y + .44*Z
        dzdX = lambda X,Y,Z: -2.3
        dzdY = lambda X,Y,Z: 1.4
        dzdZ = lambda X,Y,Z: 0.44
        
        px,py,pz = X_vec
        
        #Construct the analytic deformation gradient
        
        Fanalytic = np.reshape(np.array([dxdX(px,py,pz),dxdY(px,py,pz),dxdZ(px,py,pz),\
                                         dydX(px,py,pz),dydY(px,py,pz),dydZ(px,py,pz),\
                                         dzdX(px,py,pz),dzdY(px,py,pz),dzdZ(px,py,pz)]),[3,3])
                   
        ccoords = [[x(c1,c2,c3),y(c1,c2,c3),z(c1,c2,c3)] for c1,c2,c3 in rcoords]
        return Fanalytic,ccoords
        
    def _get_chi_values(self,rcoords,X_vec):
        """Get the values required to compute chi and grad_chi for testing"""
        
        #Define chi field
        chi11 = lambda X,Y,Z:  2.30*X-5.10*Y+0.93*Z+.8
        chi22 = lambda X,Y,Z:  0.10*X-0.71*Y-0.13*Z-.2
        chi33 = lambda X,Y,Z:  0.30*X+2.31*Y+7.30*Z
        chi23 = lambda X,Y,Z:  0.76*X-3.10*Y+2.30*Z
        chi13 = lambda X,Y,Z: -0.21*X+8.70*Y+0.30*Z
        chi12 = lambda X,Y,Z:  1.10*X-1.21*Y+0.63*Z
        chi32 = lambda X,Y,Z: -8.02*X-5.21*Y+0.23*Z
        chi31 = lambda X,Y,Z:  4.20*X+0.29*Y+3.30*Z
        chi21 = lambda X,Y,Z:  9.40*X-1.10*Y+0.30*Z

        #Compute the gradient values        
        dchi11dX = lambda X,Y,Z:  2.3
        dchi22dX = lambda X,Y,Z:  0.10
        dchi33dX = lambda X,Y,Z:  0.30
        dchi23dX = lambda X,Y,Z:  0.76
        dchi13dX = lambda X,Y,Z: -0.21
        dchi12dX = lambda X,Y,Z:  1.10
        dchi32dX = lambda X,Y,Z: -8.02
        dchi31dX = lambda X,Y,Z:  4.20
        dchi21dX = lambda X,Y,Z:  9.40
        
        dchi11dY = lambda X,Y,Z: -5.10
        dchi22dY = lambda X,Y,Z: -0.71
        dchi33dY = lambda X,Y,Z:  2.31
        dchi23dY = lambda X,Y,Z: -3.10
        dchi13dY = lambda X,Y,Z:  8.70
        dchi12dY = lambda X,Y,Z: -1.21
        dchi32dY = lambda X,Y,Z: -5.21
        dchi31dY = lambda X,Y,Z:  0.29
        dchi21dY = lambda X,Y,Z: -1.10
        
        dchi11dZ = lambda X,Y,Z:  0.93
        dchi22dZ = lambda X,Y,Z: -0.13
        dchi33dZ = lambda X,Y,Z:  7.30
        dchi23dZ = lambda X,Y,Z:  2.30
        dchi13dZ = lambda X,Y,Z:  0.30
        dchi12dZ = lambda X,Y,Z:  0.63
        dchi32dZ = lambda X,Y,Z:  0.23
        dchi31dZ = lambda X,Y,Z:  3.30
        dchi21dZ = lambda X,Y,Z:  0.30
        
        #Compute chi at the nodes
        chi = np.zeros([3,3])
        chi[0,0] = chi11(*X_vec)
        chi[1,1] = chi22(*X_vec)
        chi[2,2] = chi33(*X_vec)
        chi[1,2] = chi23(*X_vec)
        chi[0,2] = chi13(*X_vec)
        chi[0,1] = chi12(*X_vec)
        chi[2,1] = chi32(*X_vec)
        chi[2,0] = chi31(*X_vec)
        chi[1,0] = chi21(*X_vec)
        
        #print chi

        #Form the gradient
        grad_chi = np.zeros([3,3,3])
        
        grad_chi[0,0,0] = dchi11dX(*X_vec)
        grad_chi[1,1,0] = dchi22dX(*X_vec)
        grad_chi[2,2,0] = dchi33dX(*X_vec)
        grad_chi[1,2,0] = dchi23dX(*X_vec)
        grad_chi[0,2,0] = dchi13dX(*X_vec)
        grad_chi[0,1,0] = dchi12dX(*X_vec)
        grad_chi[2,1,0] = dchi32dX(*X_vec)
        grad_chi[2,0,0] = dchi31dX(*X_vec)
        grad_chi[1,0,0] = dchi21dX(*X_vec)
        
        grad_chi[0,0,1] = dchi11dY(*X_vec)
        grad_chi[1,1,1] = dchi22dY(*X_vec)
        grad_chi[2,2,1] = dchi33dY(*X_vec)
        grad_chi[1,2,1] = dchi23dY(*X_vec)
        grad_chi[0,2,1] = dchi13dY(*X_vec)
        grad_chi[0,1,1] = dchi12dY(*X_vec)
        grad_chi[2,1,1] = dchi32dY(*X_vec)
        grad_chi[2,0,1] = dchi31dY(*X_vec)
        grad_chi[1,0,1] = dchi21dY(*X_vec)
        
        grad_chi[0,0,2] = dchi11dZ(*X_vec)
        grad_chi[1,1,2] = dchi22dZ(*X_vec)
        grad_chi[2,2,2] = dchi33dZ(*X_vec)
        grad_chi[1,2,2] = dchi23dZ(*X_vec)
        grad_chi[0,2,2] = dchi13dZ(*X_vec)
        grad_chi[0,1,2] = dchi12dZ(*X_vec)
        grad_chi[2,1,2] = dchi32dZ(*X_vec)
        grad_chi[2,0,2] = dchi31dZ(*X_vec)
        grad_chi[1,0,2] = dchi21dZ(*X_vec)
        
        #print grad_chi
        
        #Form phi vectors
        
        phi_vectors = [[chi11(c1,c2,c3)-1,chi22(c1,c2,c3)-1,chi33(c1,c2,c3)-1,\
                        chi23(c1,c2,c3),  chi13(c1,c2,c3),  chi12(c1,c2,c3),\
                        chi32(c1,c2,c3),  chi31(c1,c2,c3),  chi21(c1,c2,c3)] for c1,c2,c3 in rcoords]
        
        return chi,grad_chi,phi_vectors
        
    def _flatten_list(self,l):
        """Flatten list of lists"""
        return [item for sublist in l for item in sublist]
        
    def SOTtensor_to_vector(self,T):
        """Convert a second order tensor to vector form"""
        V = np.zeros([9,])
        V[0] = T[0,0]
        V[1] = T[1,1]
        V[2] = T[2,2]
        V[3] = T[1,2]
        V[4] = T[0,2]
        V[5] = T[0,1]
        V[6] = T[2,1]
        V[7] = T[2,0]
        V[8] = T[1,0]
        return V
        
    def TOTtensor_to_vector(self,T):
        """Convert a third order tensor to vector form"""
        V = np.zeros([27,])
        for i in range(3):
            V[0+9*i] = T[0,0,i]
            V[1+9*i] = T[1,1,i]
            V[2+9*i] = T[2,2,i]
            V[3+9*i] = T[1,2,i]
            V[4+9*i] = T[0,2,i]
            V[5+9*i] = T[0,1,i]
            V[6+9*i] = T[2,1,i]
            V[7+9*i] = T[2,0,i]
            V[8+9*i] = T[1,0,i]
        return V
            

if __name__ == '__main__':
    unittest.main()
    