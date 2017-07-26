import numpy as np
import hex8
import unittest
import finite_difference as fd
import micromorphic_linear_elasticity as micro_LE
from hex8 import T_to_V_mapping as T2V
from hex8 import V_to_T_mapping as V2T

"""Definition of a 3D Micromorphic Element"""

###### Main Subroutine ######
# Intended to mimic an      #
# Abaqus UEL to allow for   #
# easy transfer between     #
# codes                     #
#############################

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
    return hex8.matrix_dot(dxdxi,dxidX)
    
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
    chi = np.empty([3,3])
    #Populate array
    chi[0,0] = 1.+phi_vector[0]
    chi[1,1] = 1.+phi_vector[1]
    chi[2,2] = 1.+phi_vector[2]
    chi[1,2] =    phi_vector[3]
    chi[0,2] =    phi_vector[4]
    chi[0,1] =    phi_vector[5]
    chi[2,1] =    phi_vector[6]
    chi[2,0] =    phi_vector[7]
    chi[1,0] =    phi_vector[8]
    return chi
    
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
                    grad_chi[i,j,k] += chin[i][j]*dNdxs[n][k]
    return grad_chi
    
###### Compute the derivatives of the fundamental quantities ######
    
def compute_fundamental_derivatives(xi_vec,nodal_global_coords_reference): #Test function written
    """Compute the derivative of the fundamental expressions with respect to the degree of freedom vector"""
    dFdU        = compute_dFdU(xi_vec,nodal_global_coords_reference)
    dchidU      = compute_dchidU(xi_vec)
    dgrad_chidU = compute_dgrad_chidU(xi_vec,nodal_global_coords_reference)
    return dFdU,dchidU,dgrad_chidU
    
def compute_dFdU(xi_vec,nodal_global_coords_reference): #Test function written
    """Compute the derivative of the deformation gradient with respect to the degree of freedom vector"""
    dFdU = np.zeros([9,96])
    for i in range(8):
        dFdU[:,(12*i):(12*(i+1))] = compute_dFndUn(i,xi_vec,nodal_global_coords_reference)
    return dFdU
    
def compute_dFndUn(n,xi_vec,nodal_global_coords_reference): #Test function written (part of test_compute_dFdU)
    """Compute the nth matrix in dFdU"""
    dFdUn   = np.zeros([9,12])
    dNdX = hex8.Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords_reference)[0]
    
    #Set the indices
    I1 = 0
    I2 = 1
    I3 = 2
    I4 = 3
    I5 = 4
    I6 = 5
    I7 = 6
    I8 = 7
    I9 = 8
    
    #Assemble the matrix
    dFdUn[I1,I1]  = dNdX[I1]
    dFdUn[I2,I2]  = dNdX[I2]
    dFdUn[I3,I3]  = dNdX[I3]
    dFdUn[I4,I2]  = dNdX[I3]
    dFdUn[I5,I1]  = dNdX[I3]
    dFdUn[I6,I1]  = dNdX[I2]
    dFdUn[I7,I3]  = dNdX[I2]
    dFdUn[I8,I3]  = dNdX[I1]
    dFdUn[I9,I2]  = dNdX[I1]
    return dFdUn
    
def compute_dchidU(xi_vec): #Test function written
    """Compute the derivative of the micro-displacement tensor with respect to the degree of freedom vector"""
    X = np.zeros([9,96])
    for i in range(8):
        X[:,(12*i):(12*(i+1))] = compute_dchindU(i,xi_vec)
    return X
    
def compute_dchindU(n,xi_vec): #Test function written (part of compute_dchidU)
    """Compute the nth matrix in dchidU"""
    return np.hstack((np.zeros([9,3]),hex8.Hex8_shape_function(n,xi_vec)*np.eye(9)))
    
def compute_dgrad_chidU(xi_vec,nodal_global_coords_reference): #Test function written
    """Compute the derivative of the gradient of the micro-displacement tensor with respect to the degree of freedom vector"""
    G = np.zeros([27,96])
    for n in range(8):
        G[:,(12*n):(12*(n+1))] = compute_dgrad_chindUn(n,xi_vec,nodal_global_coords_reference)
        
    return G
    
def compute_dgrad_chindUn(n,xi_vec,nodal_global_coords_reference): #Test function written (part of compute_dgrad_chidU)
    """Compute the nth matrix in dgrad_chindU"""
    #Initialize matrix
    Gn = np.zeros([27,12])
    
    #Set constants
    one = 1.
    
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Form the identity matrix
    I = np.zeros([9,9])
    
    #Compute the shape function gradient
    dNdX = hex8.Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords_reference)[0]
    
    for i in range(9):
        Gn[i,    i+3] = dNdX[I1]
        Gn[i+9,  i+3] = dNdX[I2]
        Gn[i+18, i+3] = dNdX[I3]
        
    return Gn
    
###### Compute Deformation Measures ######
                    
def get_deformation_measures(F,chi,grad_chi): #Test function written
    """Compute and return the deformation measures"""
    C   = hex8.matrix_Tdot(F,F)
    Phi = hex8.matrix_Tdot(F,chi)
    Gamma = hex8.matrix_Tdot_TOT(F,grad_chi)
    return C,Phi,Gamma
    
###### Compute Derivatives of Deformation Measures ######
def compute_DM_derivatives(F,chi,grad_chi,dFdU,dchidU,dgrad_chidU):
    """Compute the derivatives of the deformation measures"""
    dCdU     = compute_dCdU(F,dFdU)
    dPhidU   = compute_dPhidU(F,chi,dFdU,dchidU)
    dGammadU = compute_dGammadU(F,grad_chi,dFdU,dgrad_chidU)
    return dCdU,dPhidU,dGammadU
    
def compute_dCdU(F,dFdU): #Test function written
    """Compute the derivative of the right Cauchy-Green Deformation Tensor"""
    dCdU = np.zeros([9,96])
    
    for n in range(8):
        dCdUn = compute_dCdUn(n,F,dFdU)
        for i in range(9):
            for j in range(3):
                dCdU[i,j+(12*n)] = dCdUn[i,j]
    return dCdU

def compute_dCdUn(n,F,dFdU): #Test function written (part of test_compute_dCdU)
    """Compute a submatrix of the derivative of the right Cauchy-Green Deformation Tensor
    n goes from 0 to 7
    """
    dCdUn = np.zeros([9,12])
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Set the constants
    two = 2.
    
    #Extract the required values
    F11 = F[I1,I1]
    F22 = F[I2,I2]
    F33 = F[I3,I3]
    F23 = F[I2,I3]
    F13 = F[I1,I3]
    F12 = F[I1,I2]
    F32 = F[I3,I2]
    F31 = F[I3,I1]
    F21 = F[I2,I1]
    
    dFdU11 = dFdU[I1,  12*n]
    dFdU22 = dFdU[I2,12*n+1]
    dFdU33 = dFdU[I3,12*n+2]
    dFdU23 = dFdU[I4,12*n+1]
    dFdU13 = dFdU[I5,  12*n]
    dFdU12 = dFdU[I6,  12*n]
    dFdU32 = dFdU[I7,12*n+2]
    dFdU31 = dFdU[I8,12*n+2]
    dFdU21 = dFdU[I9,12*n+1]
    
    #Assemble the submatrix
    #Column 1
    dCdUn[I1,I1] = two*F11*dFdU11
    dCdUn[I2,I1] = two*F12*dFdU12
    dCdUn[I3,I1] = two*F13*dFdU13
    dCdUn[I4,I1] = F12*dFdU13 + F13*dFdU12
    dCdUn[I5,I1] = F11*dFdU13 + F13*dFdU11
    dCdUn[I6,I1] = F11*dFdU12 + F12*dFdU11
    dCdUn[I7,I1] = F12*dFdU13 + F13*dFdU12
    dCdUn[I8,I1] = F11*dFdU13 + F13*dFdU11
    dCdUn[I9,I1] = F11*dFdU12 + F12*dFdU11
    #Column 2
    dCdUn[I1,I2] = two*F21*dFdU21
    dCdUn[I2,I2] = two*F22*dFdU22
    dCdUn[I3,I2] = two*F23*dFdU23
    dCdUn[I4,I2] = F22*dFdU23 + F23*dFdU22
    dCdUn[I5,I2] = F21*dFdU23 + F23*dFdU21
    dCdUn[I6,I2] = F21*dFdU22 + F22*dFdU21
    dCdUn[I7,I2] = F22*dFdU23 + F23*dFdU22
    dCdUn[I8,I2] = F21*dFdU23 + F23*dFdU21
    dCdUn[I9,I2] = F21*dFdU22 + F22*dFdU21
    #Column 3
    dCdUn[I1,I3] = two*F31*dFdU31
    dCdUn[I2,I3] = two*F32*dFdU32
    dCdUn[I3,I3] = two*F33*dFdU33
    dCdUn[I4,I3] = F32*dFdU33 + F33*dFdU32
    dCdUn[I5,I3] = F31*dFdU33 + F33*dFdU31
    dCdUn[I6,I3] = F31*dFdU32 + F32*dFdU31
    dCdUn[I7,I3] = F32*dFdU33 + F33*dFdU32
    dCdUn[I8,I3] = F31*dFdU33 + F33*dFdU31
    dCdUn[I9,I3] = F31*dFdU32 + F32*dFdU31
    return dCdUn
    
def compute_dCinvdU(Cinv,dCdU): #Test function written
    """Compute the derivative of the inverse of the 
    right Cauchy-Green deformation tensor with respect 
    to the degree of freedom vector"""
    dCinvdU = np.zeros([9,96])
    
    for n in range(8):
        dCinvdUn = compute_dCinvdUn(n,Cinv,dCdU)
        for i in range(9):
            for j in range(3):
                dCinvdU[i,j+(12*n)] = dCinvdUn[i,j]
    return dCinvdU
    
def compute_dCinvdUn(n,Cinv,dCdU): #Test function written (part of compute_dCinvdU)
    """Compute a submatrix of the derivative of the inverse 
    of the right Cauchy-Green deformation tensor with respect 
    to the degree of freedom vector. n goes from 0 to 7"""
    #Initialize the output matrix
    dCinvdUn = np.zeros([9,12])
    
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Extract the required values
    Cinv11 = Cinv[I1,I1]
    Cinv22 = Cinv[I2,I2]
    Cinv33 = Cinv[I3,I3]
    Cinv23 = Cinv[I2,I3]
    Cinv13 = Cinv[I1,I3]
    Cinv12 = Cinv[I1,I2]
    Cinv32 = Cinv[I3,I2]
    Cinv31 = Cinv[I3,I1]
    Cinv21 = Cinv[I2,I1]
    
    dCdU111 = dCdU[I1,I1+12*n]
    dCdU221 = dCdU[I2,I1+12*n]
    dCdU331 = dCdU[I3,I1+12*n]
    dCdU231 = dCdU[I4,I1+12*n]
    dCdU131 = dCdU[I5,I1+12*n]
    dCdU121 = dCdU[I6,I1+12*n]
    dCdU321 = dCdU[I7,I1+12*n]
    dCdU311 = dCdU[I8,I1+12*n]
    dCdU211 = dCdU[I9,I1+12*n]
    
    dCdU112 = dCdU[I1,I2+12*n]
    dCdU222 = dCdU[I2,I2+12*n]
    dCdU332 = dCdU[I3,I2+12*n]
    dCdU232 = dCdU[I4,I2+12*n]
    dCdU132 = dCdU[I5,I2+12*n]
    dCdU122 = dCdU[I6,I2+12*n]
    dCdU322 = dCdU[I7,I2+12*n]
    dCdU312 = dCdU[I8,I2+12*n]
    dCdU212 = dCdU[I9,I2+12*n]
    
    dCdU113 = dCdU[I1,I3+12*n]
    dCdU223 = dCdU[I2,I3+12*n]
    dCdU333 = dCdU[I3,I3+12*n]
    dCdU233 = dCdU[I4,I3+12*n]
    dCdU133 = dCdU[I5,I3+12*n]
    dCdU123 = dCdU[I6,I3+12*n]
    dCdU323 = dCdU[I7,I3+12*n]
    dCdU313 = dCdU[I8,I3+12*n]
    dCdU213 = dCdU[I9,I3+12*n]
    
    #Compute and store the derivative
    
    #Column 1

    dCinvdUn[I1,I1] = -Cinv11**2*dCdU111     - Cinv11*Cinv12*dCdU211 - Cinv11*Cinv13*dCdU311 - Cinv11*Cinv21*dCdU121 - Cinv11*Cinv31*dCdU131 - Cinv12*Cinv21*dCdU221 - Cinv12*Cinv31*dCdU231 - Cinv13*Cinv21*dCdU321 - Cinv13*Cinv31*dCdU331
    dCinvdUn[I2,I1] = -Cinv12*Cinv21*dCdU111 - Cinv12*Cinv22*dCdU211 - Cinv12*Cinv23*dCdU311 - Cinv21*Cinv22*dCdU121 - Cinv21*Cinv32*dCdU131 - Cinv22**2*dCdU221 - Cinv22*Cinv23*dCdU321 - Cinv22*Cinv32*dCdU231 - Cinv23*Cinv32*dCdU331
    dCinvdUn[I3,I1] = -Cinv13*Cinv31*dCdU111 - Cinv13*Cinv32*dCdU211 - Cinv13*Cinv33*dCdU311 - Cinv23*Cinv31*dCdU121 - Cinv23*Cinv32*dCdU221 - Cinv23*Cinv33*dCdU321 - Cinv31*Cinv33*dCdU131 - Cinv32*Cinv33*dCdU231 - Cinv33**2*dCdU331
    dCinvdUn[I4,I1] = -Cinv13*Cinv21*dCdU111 - Cinv13*Cinv22*dCdU211 - Cinv13*Cinv23*dCdU311 - Cinv21*Cinv23*dCdU121 - Cinv21*Cinv33*dCdU131 - Cinv22*Cinv23*dCdU221 - Cinv22*Cinv33*dCdU231 - Cinv23**2*dCdU321 - Cinv23*Cinv33*dCdU331
    dCinvdUn[I5,I1] = -Cinv11*Cinv13*dCdU111 - Cinv11*Cinv23*dCdU121 - Cinv11*Cinv33*dCdU131 - Cinv12*Cinv13*dCdU211 - Cinv12*Cinv23*dCdU221 - Cinv12*Cinv33*dCdU231 - Cinv13**2*dCdU311 - Cinv13*Cinv23*dCdU321 - Cinv13*Cinv33*dCdU331
    dCinvdUn[I6,I1] = -Cinv11*Cinv12*dCdU111 - Cinv11*Cinv22*dCdU121 - Cinv11*Cinv32*dCdU131 - Cinv12**2*dCdU211 - Cinv12*Cinv13*dCdU311 - Cinv12*Cinv22*dCdU221 - Cinv12*Cinv32*dCdU231 - Cinv13*Cinv22*dCdU321 - Cinv13*Cinv32*dCdU331
    dCinvdUn[I7,I1] = -Cinv12*Cinv31*dCdU111 - Cinv12*Cinv32*dCdU211 - Cinv12*Cinv33*dCdU311 - Cinv22*Cinv31*dCdU121 - Cinv22*Cinv32*dCdU221 - Cinv22*Cinv33*dCdU321 - Cinv31*Cinv32*dCdU131 - Cinv32**2*dCdU231 - Cinv32*Cinv33*dCdU331
    dCinvdUn[I8,I1] = -Cinv11*Cinv31*dCdU111 - Cinv11*Cinv32*dCdU211 - Cinv11*Cinv33*dCdU311 - Cinv21*Cinv31*dCdU121 - Cinv21*Cinv32*dCdU221 - Cinv21*Cinv33*dCdU321 - Cinv31**2*dCdU131 - Cinv31*Cinv32*dCdU231 - Cinv31*Cinv33*dCdU331
    dCinvdUn[I9,I1] = -Cinv11*Cinv21*dCdU111 - Cinv11*Cinv22*dCdU211 - Cinv11*Cinv23*dCdU311 - Cinv21**2*dCdU121 - Cinv21*Cinv22*dCdU221 - Cinv21*Cinv23*dCdU321 - Cinv21*Cinv31*dCdU131 - Cinv22*Cinv31*dCdU231 - Cinv23*Cinv31*dCdU331

    #Column 2

    dCinvdUn[I1,I2] = -Cinv11**2*dCdU112     - Cinv11*Cinv12*dCdU212 - Cinv11*Cinv13*dCdU312 - Cinv11*Cinv21*dCdU122 - Cinv11*Cinv31*dCdU132 - Cinv12*Cinv21*dCdU222 - Cinv12*Cinv31*dCdU232 - Cinv13*Cinv21*dCdU322 - Cinv13*Cinv31*dCdU332
    dCinvdUn[I2,I2] = -Cinv12*Cinv21*dCdU112 - Cinv12*Cinv22*dCdU212 - Cinv12*Cinv23*dCdU312 - Cinv21*Cinv22*dCdU122 - Cinv21*Cinv32*dCdU132 - Cinv22**2*dCdU222 - Cinv22*Cinv23*dCdU322 - Cinv22*Cinv32*dCdU232 - Cinv23*Cinv32*dCdU332
    dCinvdUn[I3,I2] = -Cinv13*Cinv31*dCdU112 - Cinv13*Cinv32*dCdU212 - Cinv13*Cinv33*dCdU312 - Cinv23*Cinv31*dCdU122 - Cinv23*Cinv32*dCdU222 - Cinv23*Cinv33*dCdU322 - Cinv31*Cinv33*dCdU132 - Cinv32*Cinv33*dCdU232 - Cinv33**2*dCdU332
    dCinvdUn[I4,I2] = -Cinv13*Cinv21*dCdU112 - Cinv13*Cinv22*dCdU212 - Cinv13*Cinv23*dCdU312 - Cinv21*Cinv23*dCdU122 - Cinv21*Cinv33*dCdU132 - Cinv22*Cinv23*dCdU222 - Cinv22*Cinv33*dCdU232 - Cinv23**2*dCdU322 - Cinv23*Cinv33*dCdU332
    dCinvdUn[I5,I2] = -Cinv11*Cinv13*dCdU112 - Cinv11*Cinv23*dCdU122 - Cinv11*Cinv33*dCdU132 - Cinv12*Cinv13*dCdU212 - Cinv12*Cinv23*dCdU222 - Cinv12*Cinv33*dCdU232 - Cinv13**2*dCdU312 - Cinv13*Cinv23*dCdU322 - Cinv13*Cinv33*dCdU332
    dCinvdUn[I6,I2] = -Cinv11*Cinv12*dCdU112 - Cinv11*Cinv22*dCdU122 - Cinv11*Cinv32*dCdU132 - Cinv12**2*dCdU212 - Cinv12*Cinv13*dCdU312 - Cinv12*Cinv22*dCdU222 - Cinv12*Cinv32*dCdU232 - Cinv13*Cinv22*dCdU322 - Cinv13*Cinv32*dCdU332
    dCinvdUn[I7,I2] = -Cinv12*Cinv31*dCdU112 - Cinv12*Cinv32*dCdU212 - Cinv12*Cinv33*dCdU312 - Cinv22*Cinv31*dCdU122 - Cinv22*Cinv32*dCdU222 - Cinv22*Cinv33*dCdU322 - Cinv31*Cinv32*dCdU132 - Cinv32**2*dCdU232 - Cinv32*Cinv33*dCdU332
    dCinvdUn[I8,I2] = -Cinv11*Cinv31*dCdU112 - Cinv11*Cinv32*dCdU212 - Cinv11*Cinv33*dCdU312 - Cinv21*Cinv31*dCdU122 - Cinv21*Cinv32*dCdU222 - Cinv21*Cinv33*dCdU322 - Cinv31**2*dCdU132 - Cinv31*Cinv32*dCdU232 - Cinv31*Cinv33*dCdU332
    dCinvdUn[I9,I2] = -Cinv11*Cinv21*dCdU112 - Cinv11*Cinv22*dCdU212 - Cinv11*Cinv23*dCdU312 - Cinv21**2*dCdU122 - Cinv21*Cinv22*dCdU222 - Cinv21*Cinv23*dCdU322 - Cinv21*Cinv31*dCdU132 - Cinv22*Cinv31*dCdU232 - Cinv23*Cinv31*dCdU332

    #Column 3

    dCinvdUn[I1,I3] = -Cinv11**2*dCdU113     - Cinv11*Cinv12*dCdU213 - Cinv11*Cinv13*dCdU313 - Cinv11*Cinv21*dCdU123 - Cinv11*Cinv31*dCdU133 - Cinv12*Cinv21*dCdU223 - Cinv12*Cinv31*dCdU233 - Cinv13*Cinv21*dCdU323 - Cinv13*Cinv31*dCdU333
    dCinvdUn[I2,I3] = -Cinv12*Cinv21*dCdU113 - Cinv12*Cinv22*dCdU213 - Cinv12*Cinv23*dCdU313 - Cinv21*Cinv22*dCdU123 - Cinv21*Cinv32*dCdU133 - Cinv22**2*dCdU223     - Cinv22*Cinv23*dCdU323 - Cinv22*Cinv32*dCdU233 - Cinv23*Cinv32*dCdU333
    dCinvdUn[I3,I3] = -Cinv13*Cinv31*dCdU113 - Cinv13*Cinv32*dCdU213 - Cinv13*Cinv33*dCdU313 - Cinv23*Cinv31*dCdU123 - Cinv23*Cinv32*dCdU223 - Cinv23*Cinv33*dCdU323 - Cinv31*Cinv33*dCdU133 - Cinv32*Cinv33*dCdU233 - Cinv33**2*dCdU333
    dCinvdUn[I4,I3] = -Cinv13*Cinv21*dCdU113 - Cinv13*Cinv22*dCdU213 - Cinv13*Cinv23*dCdU313 - Cinv21*Cinv23*dCdU123 - Cinv21*Cinv33*dCdU133 - Cinv22*Cinv23*dCdU223 - Cinv22*Cinv33*dCdU233 - Cinv23**2*dCdU323     - Cinv23*Cinv33*dCdU333
    dCinvdUn[I5,I3] = -Cinv11*Cinv13*dCdU113 - Cinv11*Cinv23*dCdU123 - Cinv11*Cinv33*dCdU133 - Cinv12*Cinv13*dCdU213 - Cinv12*Cinv23*dCdU223 - Cinv12*Cinv33*dCdU233 - Cinv13**2*dCdU313     - Cinv13*Cinv23*dCdU323 - Cinv13*Cinv33*dCdU333
    dCinvdUn[I6,I3] = -Cinv11*Cinv12*dCdU113 - Cinv11*Cinv22*dCdU123 - Cinv11*Cinv32*dCdU133 - Cinv12**2*dCdU213     - Cinv12*Cinv13*dCdU313 - Cinv12*Cinv22*dCdU223 - Cinv12*Cinv32*dCdU233 - Cinv13*Cinv22*dCdU323 - Cinv13*Cinv32*dCdU333
    dCinvdUn[I7,I3] = -Cinv12*Cinv31*dCdU113 - Cinv12*Cinv32*dCdU213 - Cinv12*Cinv33*dCdU313 - Cinv22*Cinv31*dCdU123 - Cinv22*Cinv32*dCdU223 - Cinv22*Cinv33*dCdU323 - Cinv31*Cinv32*dCdU133 - Cinv32**2*dCdU233     - Cinv32*Cinv33*dCdU333
    dCinvdUn[I8,I3] = -Cinv11*Cinv31*dCdU113 - Cinv11*Cinv32*dCdU213 - Cinv11*Cinv33*dCdU313 - Cinv21*Cinv31*dCdU123 - Cinv21*Cinv32*dCdU223 - Cinv21*Cinv33*dCdU323 - Cinv31**2*dCdU133     - Cinv31*Cinv32*dCdU233 - Cinv31*Cinv33*dCdU333
    dCinvdUn[I9,I3] = -Cinv11*Cinv21*dCdU113 - Cinv11*Cinv22*dCdU213 - Cinv11*Cinv23*dCdU313 - Cinv21**2*dCdU123     - Cinv21*Cinv22*dCdU223 - Cinv21*Cinv23*dCdU323 - Cinv21*Cinv31*dCdU133 - Cinv22*Cinv31*dCdU233 - Cinv23*Cinv31*dCdU333
    
    return dCinvdUn
    
def compute_dPhidU(F,chi,dFdU,dchidU): #Test function written
    """Compute the derivative of the deformation measure
    Phi with respect to the deformation gradient"""
    dPhidU = np.zeros([9,96])
    
    for n in range(8):
        dPhidUn = compute_dPhidUn(n,F,chi,dFdU,dchidU)
        for i in range(9):
            for j in range(12):
                dPhidU[i,j+(12*n)] = dPhidUn[i,j]
    return dPhidU
    
def compute_dPhidUn(n,F,chi,dFdU,dchidU): #Test function written (part of test_compute_dPhidU)
    """Compute a submatrix of the derivative of the deformation measure Phi with 
    respect to the dof vector n goes from 0 to 7
    """
    dPhidUn = np.zeros([9,12])
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Extract the required values
    F11 = F[I1,I1]
    F22 = F[I2,I2]
    F33 = F[I3,I3]
    F23 = F[I2,I3]
    F13 = F[I1,I3]
    F12 = F[I1,I2]
    F32 = F[I3,I2]
    F31 = F[I3,I1]
    F21 = F[I2,I1]
    
    chi11 = chi[I1,I1]
    chi22 = chi[I2,I2]
    chi33 = chi[I3,I3]
    chi23 = chi[I2,I3]
    chi13 = chi[I1,I3]
    chi12 = chi[I1,I2]
    chi32 = chi[I3,I2]
    chi31 = chi[I3,I1]
    chi21 = chi[I2,I1]
    
    dFdU11 = dFdU[I1,  12*n]
    dFdU22 = dFdU[I2,12*n+1]
    dFdU33 = dFdU[I3,12*n+2]
    dFdU23 = dFdU[I4,12*n+1]
    dFdU13 = dFdU[I5,  12*n]
    dFdU12 = dFdU[I6,  12*n]
    dFdU32 = dFdU[I7,12*n+2]
    dFdU31 = dFdU[I8,12*n+2]
    dFdU21 = dFdU[I9,12*n+1]
    
    dchidU11 = dchidU[I1, I4+12*n]
    dchidU22 = dchidU[I2, I5+12*n]
    dchidU33 = dchidU[I3, I6+12*n]
    dchidU23 = dchidU[I4, I7+12*n]
    dchidU13 = dchidU[I5, I8+12*n]
    dchidU12 = dchidU[I6, I9+12*n]
    dchidU32 = dchidU[I7,I10+12*n]
    dchidU31 = dchidU[I8,I11+12*n]
    dchidU21 = dchidU[I9,I12+12*n]
    
    #Assemble the submatrix
    #Column 1
    dPhidUn[I1,I1]  = chi11*dFdU11
    dPhidUn[I2,I1]  = chi12*dFdU12
    dPhidUn[I3,I1]  = chi13*dFdU13
    dPhidUn[I4,I1]  = chi13*dFdU12
    dPhidUn[I5,I1]  = chi13*dFdU11
    dPhidUn[I6,I1]  = chi12*dFdU11
    dPhidUn[I7,I1]  = chi12*dFdU13
    dPhidUn[I8,I1]  = chi11*dFdU13
    dPhidUn[I9,I1]  = chi11*dFdU12
    #Column 2
    dPhidUn[I1,I2]  = chi21*dFdU21 
    dPhidUn[I2,I2]  = chi22*dFdU22
    dPhidUn[I3,I2]  = chi23*dFdU23
    dPhidUn[I4,I2]  = chi23*dFdU22
    dPhidUn[I5,I2]  = chi23*dFdU21
    dPhidUn[I6,I2]  = chi22*dFdU21
    dPhidUn[I7,I2]  = chi22*dFdU23
    dPhidUn[I8,I2]  = chi21*dFdU23
    dPhidUn[I9,I2]  = chi21*dFdU22
    #Column 3
    dPhidUn[I1,I3]  = chi31*dFdU31
    dPhidUn[I2,I3]  = chi32*dFdU32
    dPhidUn[I3,I3]  = chi33*dFdU33
    dPhidUn[I4,I3]  = chi33*dFdU32
    dPhidUn[I5,I3]  = chi33*dFdU31
    dPhidUn[I6,I3]  = chi32*dFdU31
    dPhidUn[I7,I3]  = chi32*dFdU33
    dPhidUn[I8,I3]  = chi31*dFdU33
    dPhidUn[I9,I3]  = chi31*dFdU32
    #Column  4
    dPhidUn[I1,I4]  = F11*dchidU11
    dPhidUn[I8,I4]  = F13*dchidU11
    dPhidUn[I9,I4]  = F12*dchidU11
    #Column  5
    dPhidUn[I2,I5]  = F22*dchidU22
    dPhidUn[I6,I5]  = F21*dchidU22
    dPhidUn[I7,I5]  = F23*dchidU22
    #Column  6
    dPhidUn[I3,I6]  = F33*dchidU33
    dPhidUn[I4,I6]  = F32*dchidU33
    dPhidUn[I5,I6]  = F31*dchidU33
    #Column  7
    dPhidUn[I3,I7]  = F23*dchidU23
    dPhidUn[I4,I7]  = F22*dchidU23
    dPhidUn[I5,I7]  = F21*dchidU23
    #Column  8
    dPhidUn[I3,I8]  = F13*dchidU13
    dPhidUn[I4,I8]  = F12*dchidU13
    dPhidUn[I5,I8]  = F11*dchidU13
    #Column  9
    dPhidUn[I2,I9]  = F12*dchidU12
    dPhidUn[I6,I9]  = F11*dchidU12
    dPhidUn[I7,I9]  = F13*dchidU12
    #Column 10
    dPhidUn[I2,I10] = F32*dchidU32
    dPhidUn[I6,I10] = F31*dchidU32
    dPhidUn[I7,I10] = F33*dchidU32
    #Column 11
    dPhidUn[I1,I11] = F31*dchidU31
    dPhidUn[I8,I11] = F33*dchidU31
    dPhidUn[I9,I11] = F32*dchidU31
    #Column 12
    dPhidUn[I1,I12] = F21*dchidU21
    dPhidUn[I8,I12] = F23*dchidU21
    dPhidUn[I9,I12] = F22*dchidU21
    
    return dPhidUn
    
def compute_dGammadU(F,grad_chi,dFdU,dgrad_chidU): #Test function written
    """Compute the derivative of the deformation measure
    Gamma with respect to the deformation gradient"""
    dGammadU = np.zeros([27,96])
    
    for n in range(8):
        dGammadUn = compute_dGammadUn(n,F,grad_chi,dFdU,dgrad_chidU)
        for i in range(27):
            for j in range(12):
                dGammadU[   i,j+(12*n)] = dGammadUn[   i,j]
    return dGammadU
    
def compute_dGammadUn(n,F,grad_chi,dFdU,dgrad_chidU): #Test function written (contained in test_compute_dGammadU)
    """Compute a submatrix of the derivative of the
    deformation measure Gamma with respect to the
    degree of freedom vector. n goes from 0 to 7"""
    #Initialize the array
    dGammadUn = np.zeros([27,12])
    
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Extract the required values
    
    #Extract F
    F11 = F[I1,I1]
    F22 = F[I2,I2]
    F33 = F[I3,I3]
    F23 = F[I2,I3]
    F13 = F[I1,I3]
    F12 = F[I1,I2]
    F32 = F[I3,I2]
    F31 = F[I3,I1]
    F21 = F[I2,I1]
    
    #Extract the gradient of chi
    grad_chi111 = grad_chi[I1,I1,I1]
    grad_chi221 = grad_chi[I2,I2,I1]
    grad_chi331 = grad_chi[I3,I3,I1]
    grad_chi231 = grad_chi[I2,I3,I1]
    grad_chi131 = grad_chi[I1,I3,I1]
    grad_chi121 = grad_chi[I1,I2,I1]
    grad_chi321 = grad_chi[I3,I2,I1]
    grad_chi311 = grad_chi[I3,I1,I1]
    grad_chi211 = grad_chi[I2,I1,I1]
    
    grad_chi112 = grad_chi[I1,I1,I2]
    grad_chi222 = grad_chi[I2,I2,I2]
    grad_chi332 = grad_chi[I3,I3,I2]
    grad_chi232 = grad_chi[I2,I3,I2]
    grad_chi132 = grad_chi[I1,I3,I2]
    grad_chi122 = grad_chi[I1,I2,I2]
    grad_chi322 = grad_chi[I3,I2,I2]
    grad_chi312 = grad_chi[I3,I1,I2]
    grad_chi212 = grad_chi[I2,I1,I2]
    
    grad_chi113 = grad_chi[I1,I1,I3]
    grad_chi223 = grad_chi[I2,I2,I3]
    grad_chi333 = grad_chi[I3,I3,I3]
    grad_chi233 = grad_chi[I2,I3,I3]
    grad_chi133 = grad_chi[I1,I3,I3]
    grad_chi123 = grad_chi[I1,I2,I3]
    grad_chi323 = grad_chi[I3,I2,I3]
    grad_chi313 = grad_chi[I3,I1,I3]
    grad_chi213 = grad_chi[I2,I1,I3]
    
    #Extract dFdU
    dFdU11 = dFdU[I1,  12*n]
    dFdU22 = dFdU[I2,12*n+1]
    dFdU33 = dFdU[I3,12*n+2]
    dFdU23 = dFdU[I4,12*n+1]
    dFdU13 = dFdU[I5,  12*n]
    dFdU12 = dFdU[I6,  12*n]
    dFdU32 = dFdU[I7,12*n+2]
    dFdU31 = dFdU[I8,12*n+2]
    dFdU21 = dFdU[I9,12*n+1]
    
    #Extract dgrad_chidU
    dgrad_chidU111 = dgrad_chidU[   I1, I4+12*n]
    dgrad_chidU221 = dgrad_chidU[   I2, I5+12*n]
    dgrad_chidU331 = dgrad_chidU[   I3, I6+12*n]
    dgrad_chidU231 = dgrad_chidU[   I4, I7+12*n]
    dgrad_chidU131 = dgrad_chidU[   I5, I8+12*n]
    dgrad_chidU121 = dgrad_chidU[   I6, I9+12*n]
    dgrad_chidU321 = dgrad_chidU[   I7,I10+12*n]
    dgrad_chidU311 = dgrad_chidU[   I8,I11+12*n]
    dgrad_chidU211 = dgrad_chidU[   I9,I12+12*n]
    
    dgrad_chidU112 = dgrad_chidU[ I1+9, I4+12*n]
    dgrad_chidU222 = dgrad_chidU[ I2+9, I5+12*n]
    dgrad_chidU332 = dgrad_chidU[ I3+9, I6+12*n]
    dgrad_chidU232 = dgrad_chidU[ I4+9, I7+12*n]
    dgrad_chidU132 = dgrad_chidU[ I5+9, I8+12*n]
    dgrad_chidU122 = dgrad_chidU[ I6+9, I9+12*n]
    dgrad_chidU322 = dgrad_chidU[ I7+9,I10+12*n]
    dgrad_chidU312 = dgrad_chidU[ I8+9,I11+12*n]
    dgrad_chidU212 = dgrad_chidU[ I9+9,I12+12*n]
    
    dgrad_chidU113 = dgrad_chidU[I1+18, I4+12*n]
    dgrad_chidU223 = dgrad_chidU[I2+18, I5+12*n]
    dgrad_chidU333 = dgrad_chidU[I3+18, I6+12*n]
    dgrad_chidU233 = dgrad_chidU[I4+18, I7+12*n]
    dgrad_chidU133 = dgrad_chidU[I5+18, I8+12*n]
    dgrad_chidU123 = dgrad_chidU[I6+18, I9+12*n]
    dgrad_chidU323 = dgrad_chidU[I7+18,I10+12*n]
    dgrad_chidU313 = dgrad_chidU[I8+18,I11+12*n]
    dgrad_chidU213 = dgrad_chidU[I9+18,I12+12*n]
    
    #Assemble the submatrix
    #Column 1

    #print I1+9
    #print I1+18
    #print dFdU11*grad_chi112
    
    dGammadUn[   I1,I1] = dFdU11*grad_chi111
    dGammadUn[   I2,I1] = dFdU12*grad_chi121
    dGammadUn[   I3,I1] = dFdU13*grad_chi131
    dGammadUn[   I4,I1] = dFdU12*grad_chi131
    dGammadUn[   I5,I1] = dFdU11*grad_chi131
    dGammadUn[   I6,I1] = dFdU11*grad_chi121
    dGammadUn[   I7,I1] = dFdU13*grad_chi121
    dGammadUn[   I8,I1] = dFdU13*grad_chi111
    dGammadUn[   I9,I1] = dFdU12*grad_chi111
    dGammadUn[ I1+9,I1] = dFdU11*grad_chi112
    dGammadUn[ I2+9,I1] = dFdU12*grad_chi122
    dGammadUn[ I3+9,I1] = dFdU13*grad_chi132
    dGammadUn[ I4+9,I1] = dFdU12*grad_chi132
    dGammadUn[ I5+9,I1] = dFdU11*grad_chi132
    dGammadUn[ I6+9,I1] = dFdU11*grad_chi122
    dGammadUn[ I7+9,I1] = dFdU13*grad_chi122
    dGammadUn[ I8+9,I1] = dFdU13*grad_chi112
    dGammadUn[ I9+9,I1] = dFdU12*grad_chi112
    dGammadUn[I1+18,I1] = dFdU11*grad_chi113
    dGammadUn[I2+18,I1] = dFdU12*grad_chi123
    dGammadUn[I3+18,I1] = dFdU13*grad_chi133
    dGammadUn[I4+18,I1] = dFdU12*grad_chi133
    dGammadUn[I5+18,I1] = dFdU11*grad_chi133
    dGammadUn[I6+18,I1] = dFdU11*grad_chi123
    dGammadUn[I7+18,I1] = dFdU13*grad_chi123
    dGammadUn[I8+18,I1] = dFdU13*grad_chi113
    dGammadUn[I9+18,I1] = dFdU12*grad_chi113

    #Column 2

    dGammadUn[   I1,I2] = dFdU21*grad_chi211
    dGammadUn[   I2,I2] = dFdU22*grad_chi221
    dGammadUn[   I3,I2] = dFdU23*grad_chi231
    dGammadUn[   I4,I2] = dFdU22*grad_chi231
    dGammadUn[   I5,I2] = dFdU21*grad_chi231
    dGammadUn[   I6,I2] = dFdU21*grad_chi221
    dGammadUn[   I7,I2] = dFdU23*grad_chi221
    dGammadUn[   I8,I2] = dFdU23*grad_chi211
    dGammadUn[   I9,I2] = dFdU22*grad_chi211
    dGammadUn[ I1+9,I2] = dFdU21*grad_chi212
    dGammadUn[ I2+9,I2] = dFdU22*grad_chi222
    dGammadUn[ I3+9,I2] = dFdU23*grad_chi232
    dGammadUn[ I4+9,I2] = dFdU22*grad_chi232
    dGammadUn[ I5+9,I2] = dFdU21*grad_chi232
    dGammadUn[ I6+9,I2] = dFdU21*grad_chi222
    dGammadUn[ I7+9,I2] = dFdU23*grad_chi222
    dGammadUn[ I8+9,I2] = dFdU23*grad_chi212
    dGammadUn[ I9+9,I2] = dFdU22*grad_chi212
    dGammadUn[I1+18,I2] = dFdU21*grad_chi213
    dGammadUn[I2+18,I2] = dFdU22*grad_chi223
    dGammadUn[I3+18,I2] = dFdU23*grad_chi233
    dGammadUn[I4+18,I2] = dFdU22*grad_chi233
    dGammadUn[I5+18,I2] = dFdU21*grad_chi233
    dGammadUn[I6+18,I2] = dFdU21*grad_chi223
    dGammadUn[I7+18,I2] = dFdU23*grad_chi223
    dGammadUn[I8+18,I2] = dFdU23*grad_chi213
    dGammadUn[I9+18,I2] = dFdU22*grad_chi213

    #Column 3

    dGammadUn[   I1,I3] = dFdU31*grad_chi311
    dGammadUn[   I2,I3] = dFdU32*grad_chi321
    dGammadUn[   I3,I3] = dFdU33*grad_chi331
    dGammadUn[   I4,I3] = dFdU32*grad_chi331
    dGammadUn[   I5,I3] = dFdU31*grad_chi331
    dGammadUn[   I6,I3] = dFdU31*grad_chi321
    dGammadUn[   I7,I3] = dFdU33*grad_chi321
    dGammadUn[   I8,I3] = dFdU33*grad_chi311
    dGammadUn[   I9,I3] = dFdU32*grad_chi311
    dGammadUn[ I1+9,I3] = dFdU31*grad_chi312
    dGammadUn[ I2+9,I3] = dFdU32*grad_chi322
    dGammadUn[ I3+9,I3] = dFdU33*grad_chi332
    dGammadUn[ I4+9,I3] = dFdU32*grad_chi332
    dGammadUn[ I5+9,I3] = dFdU31*grad_chi332
    dGammadUn[ I6+9,I3] = dFdU31*grad_chi322
    dGammadUn[ I7+9,I3] = dFdU33*grad_chi322
    dGammadUn[ I8+9,I3] = dFdU33*grad_chi312
    dGammadUn[ I9+9,I3] = dFdU32*grad_chi312
    dGammadUn[I1+18,I3] = dFdU31*grad_chi313
    dGammadUn[I2+18,I3] = dFdU32*grad_chi323
    dGammadUn[I3+18,I3] = dFdU33*grad_chi333
    dGammadUn[I4+18,I3] = dFdU32*grad_chi333
    dGammadUn[I5+18,I3] = dFdU31*grad_chi333
    dGammadUn[I6+18,I3] = dFdU31*grad_chi323
    dGammadUn[I7+18,I3] = dFdU33*grad_chi323
    dGammadUn[I8+18,I3] = dFdU33*grad_chi313
    dGammadUn[I9+18,I3] = dFdU32*grad_chi313

    #Column 4

    dGammadUn[   I1,I4] = F11*dgrad_chidU111
    dGammadUn[   I8,I4] = F13*dgrad_chidU111
    dGammadUn[   I9,I4] = F12*dgrad_chidU111
    dGammadUn[ I1+9,I4] = F11*dgrad_chidU112
    dGammadUn[ I8+9,I4] = F13*dgrad_chidU112
    dGammadUn[ I9+9,I4] = F12*dgrad_chidU112
    dGammadUn[I1+18,I4] = F11*dgrad_chidU113
    dGammadUn[I8+18,I4] = F13*dgrad_chidU113
    dGammadUn[I9+18,I4] = F12*dgrad_chidU113

    #Column 5

    dGammadUn[   I2,I5] = F22*dgrad_chidU221
    dGammadUn[   I6,I5] = F21*dgrad_chidU221
    dGammadUn[   I7,I5] = F23*dgrad_chidU221
    dGammadUn[ I2+9,I5] = F22*dgrad_chidU222
    dGammadUn[ I6+9,I5] = F21*dgrad_chidU222
    dGammadUn[ I7+9,I5] = F23*dgrad_chidU222
    dGammadUn[I2+18,I5] = F22*dgrad_chidU223
    dGammadUn[I6+18,I5] = F21*dgrad_chidU223
    dGammadUn[I7+18,I5] = F23*dgrad_chidU223

    #Column 6

    dGammadUn[   I3,I6] = F33*dgrad_chidU331
    dGammadUn[   I4,I6] = F32*dgrad_chidU331
    dGammadUn[   I5,I6] = F31*dgrad_chidU331
    dGammadUn[ I3+9,I6] = F33*dgrad_chidU332
    dGammadUn[ I4+9,I6] = F32*dgrad_chidU332
    dGammadUn[ I5+9,I6] = F31*dgrad_chidU332
    dGammadUn[I3+18,I6] = F33*dgrad_chidU333
    dGammadUn[I4+18,I6] = F32*dgrad_chidU333
    dGammadUn[I5+18,I6] = F31*dgrad_chidU333

    #Column 7

    dGammadUn[   I3,I7] = F23*dgrad_chidU231
    dGammadUn[   I4,I7] = F22*dgrad_chidU231
    dGammadUn[   I5,I7] = F21*dgrad_chidU231
    dGammadUn[ I3+9,I7] = F23*dgrad_chidU232
    dGammadUn[ I4+9,I7] = F22*dgrad_chidU232
    dGammadUn[ I5+9,I7] = F21*dgrad_chidU232
    dGammadUn[I3+18,I7] = F23*dgrad_chidU233
    dGammadUn[I4+18,I7] = F22*dgrad_chidU233
    dGammadUn[I5+18,I7] = F21*dgrad_chidU233

    #Column 8

    dGammadUn[   I3,I8] = F13*dgrad_chidU131
    dGammadUn[   I4,I8] = F12*dgrad_chidU131
    dGammadUn[   I5,I8] = F11*dgrad_chidU131
    dGammadUn[ I3+9,I8] = F13*dgrad_chidU132
    dGammadUn[ I4+9,I8] = F12*dgrad_chidU132
    dGammadUn[ I5+9,I8] = F11*dgrad_chidU132
    dGammadUn[I3+18,I8] = F13*dgrad_chidU133
    dGammadUn[I4+18,I8] = F12*dgrad_chidU133
    dGammadUn[I5+18,I8] = F11*dgrad_chidU133

    #Column 9

    dGammadUn[   I2,I9] = F12*dgrad_chidU121
    dGammadUn[   I6,I9] = F11*dgrad_chidU121
    dGammadUn[   I7,I9] = F13*dgrad_chidU121
    dGammadUn[ I2+9,I9] = F12*dgrad_chidU122
    dGammadUn[ I6+9,I9] = F11*dgrad_chidU122
    dGammadUn[ I7+9,I9] = F13*dgrad_chidU122
    dGammadUn[I2+18,I9] = F12*dgrad_chidU123
    dGammadUn[I6+18,I9] = F11*dgrad_chidU123
    dGammadUn[I7+18,I9] = F13*dgrad_chidU123

    #Column 10

    dGammadUn[   I2,I10] = F32*dgrad_chidU321
    dGammadUn[   I6,I10] = F31*dgrad_chidU321
    dGammadUn[   I7,I10] = F33*dgrad_chidU321
    dGammadUn[ I2+9,I10] = F32*dgrad_chidU322
    dGammadUn[ I6+9,I10] = F31*dgrad_chidU322
    dGammadUn[ I7+9,I10] = F33*dgrad_chidU322
    dGammadUn[I2+18,I10] = F32*dgrad_chidU323
    dGammadUn[I6+18,I10] = F31*dgrad_chidU323
    dGammadUn[I7+18,I10] = F33*dgrad_chidU323

    #Column 11

    dGammadUn[   I1,I11] = F31*dgrad_chidU311
    dGammadUn[   I8,I11] = F33*dgrad_chidU311
    dGammadUn[   I9,I11] = F32*dgrad_chidU311
    dGammadUn[ I1+9,I11] = F31*dgrad_chidU312
    dGammadUn[ I8+9,I11] = F33*dgrad_chidU312
    dGammadUn[ I9+9,I11] = F32*dgrad_chidU312
    dGammadUn[I1+18,I11] = F31*dgrad_chidU313
    dGammadUn[I8+18,I11] = F33*dgrad_chidU313
    dGammadUn[I9+18,I11] = F32*dgrad_chidU313

    #Column 12

    dGammadUn[   I1,I12] = F21*dgrad_chidU211
    dGammadUn[   I8,I12] = F23*dgrad_chidU211
    dGammadUn[   I9,I12] = F22*dgrad_chidU211
    dGammadUn[ I1+9,I12] = F21*dgrad_chidU212
    dGammadUn[ I8+9,I12] = F23*dgrad_chidU212
    dGammadUn[ I9+9,I12] = F22*dgrad_chidU212
    dGammadUn[I1+18,I12] = F21*dgrad_chidU213
    dGammadUn[I8+18,I12] = F23*dgrad_chidU213
    dGammadUn[I9+18,I12] = F22*dgrad_chidU213
    
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
    
def compute_dpk2dU(dSdC,dSdPhi,dSdGamma,dCdU,dPhidU,dGammadU):
    """Compute the derivative of the second piola kirchhoff stress
    with respect to the degree of freedom vector"""
    
    #Initialize
    dpk2dU = np.zeros([9,12*8])
    
    for n in range(8):
        dpk2dUn = compute_dpk2dUn(dSdC[:,(n*12):((n+1)*12)],dSdPhi[:,(n*12):((n+1)*12)],dSdGamma[:,(n*12):((n+1)*12)],\
                                  dCdU[:,(n*12):((n+1)*12)],dPhidU[:,(n*12):((n+1)*12)],dGammadU[:,(n*12):((n+1)*12)])
        for i in range(9):
            for j in range(12):
                dpk2dU[i,j+n*12] = dpk2dUn[i,j]
                
    return dpk2dU
    
def compute_dpk2dUn(dSdC,dSdPhi,dSdGamma,dCdU,dPhidU,dGammadU):
    """Compute a submatrix of the derivative of the second piola kirchhoff stress
    with respect to the degree of freedom vector"""
    #Initialize
    dpk2dU = np.zeros([9,12])
    
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Extract values
    
    #Extract components of dCdU
    dCdU111 = dCdU[I1,I1]
    dCdU221 = dCdU[I2,I1]
    dCdU331 = dCdU[I3,I1]
    dCdU231 = dCdU[I4,I1]
    dCdU131 = dCdU[I5,I1]
    dCdU121 = dCdU[I6,I1]
    dCdU321 = dCdU[I7,I1]
    dCdU311 = dCdU[I8,I1]
    dCdU211 = dCdU[I9,I1]
    dCdU112 = dCdU[I1,I2]
    dCdU222 = dCdU[I2,I2]
    dCdU332 = dCdU[I3,I2]
    dCdU232 = dCdU[I4,I2]
    dCdU132 = dCdU[I5,I2]
    dCdU122 = dCdU[I6,I2]
    dCdU322 = dCdU[I7,I2]
    dCdU312 = dCdU[I8,I2]
    dCdU212 = dCdU[I9,I2]
    dCdU113 = dCdU[I1,I3]
    dCdU223 = dCdU[I2,I3]
    dCdU333 = dCdU[I3,I3]
    dCdU233 = dCdU[I4,I3]
    dCdU133 = dCdU[I5,I3]
    dCdU123 = dCdU[I6,I3]
    dCdU323 = dCdU[I7,I3]
    dCdU313 = dCdU[I8,I3]
    dCdU213 = dCdU[I9,I3]

    #Extract components of dPsidU
    dPhidU111 = dPsidU[I1,I1]
    dPhidU221 = dPsidU[I2,I1]
    dPhidU331 = dPsidU[I3,I1]
    dPhidU231 = dPsidU[I4,I1]
    dPhidU131 = dPsidU[I5,I1]
    dPhidU121 = dPsidU[I6,I1]
    dPhidU321 = dPsidU[I7,I1]
    dPhidU311 = dPsidU[I8,I1]
    dPhidU211 = dPsidU[I9,I1]
    dPhidU112 = dPsidU[I1,I2]
    dPhidU222 = dPsidU[I2,I2]
    dPhidU332 = dPsidU[I3,I2]
    dPhidU232 = dPsidU[I4,I2]
    dPhidU132 = dPsidU[I5,I2]
    dPhidU122 = dPsidU[I6,I2]
    dPhidU322 = dPsidU[I7,I2]
    dPhidU312 = dPsidU[I8,I2]
    dPhidU212 = dPsidU[I9,I2]
    dPhidU113 = dPsidU[I1,I3]
    dPhidU223 = dPsidU[I2,I3]
    dPhidU333 = dPsidU[I3,I3]
    dPhidU233 = dPsidU[I4,I3]
    dPhidU133 = dPsidU[I5,I3]
    dPhidU123 = dPsidU[I6,I3]
    dPhidU323 = dPsidU[I7,I3]
    dPhidU313 = dPsidU[I8,I3]
    dPhidU213 = dPsidU[I9,I3]
    dPhidU114 = dPsidU[I1,I4]
    dPhidU314 = dPsidU[I8,I4]
    dPhidU214 = dPsidU[I9,I4]
    dPhidU225 = dPsidU[I2,I5]
    dPhidU125 = dPsidU[I6,I5]
    dPhidU325 = dPsidU[I7,I5]
    dPhidU336 = dPsidU[I3,I6]
    dPhidU236 = dPsidU[I4,I6]
    dPhidU136 = dPsidU[I5,I6]
    dPhidU337 = dPsidU[I3,I7]
    dPhidU237 = dPsidU[I4,I7]
    dPhidU137 = dPsidU[I5,I7]
    dPhidU338 = dPsidU[I3,I8]
    dPhidU238 = dPsidU[I4,I8]
    dPhidU138 = dPsidU[I5,I8]
    dPhidU229 = dPsidU[I2,I9]
    dPhidU129 = dPsidU[I6,I9]
    dPhidU329 = dPsidU[I7,I9]
    dPhidU2210 = dPsidU[I2,I10]
    dPhidU1210 = dPsidU[I6,I10]
    dPhidU3210 = dPsidU[I7,I10]
    dPhidU1111 = dPsidU[I1,I11]
    dPhidU3111 = dPsidU[I8,I11]
    dPhidU2111 = dPsidU[I9,I11]
    dPhidU1112 = dPsidU[I1,I12]
    dPhidU3112 = dPsidU[I8,I12]
    dPhidU2112 = dPsidU[I9,I12]

    #Extract components of dGammadU
    dGammadU1111 = dGammadU[I1,I1]
    dGammadU2211 = dGammadU[I2,I1]
    dGammadU3311 = dGammadU[I3,I1]
    dGammadU2311 = dGammadU[I4,I1]
    dGammadU1311 = dGammadU[I5,I1]
    dGammadU1211 = dGammadU[I6,I1]
    dGammadU3211 = dGammadU[I7,I1]
    dGammadU3111 = dGammadU[I8,I1]
    dGammadU2111 = dGammadU[I9,I1]
    dGammadU1121 = dGammadU[I1+9,I1]
    dGammadU2221 = dGammadU[I2+9,I1]
    dGammadU3321 = dGammadU[I3+9,I1]
    dGammadU2321 = dGammadU[I4+9,I1]
    dGammadU1321 = dGammadU[I5+9,I1]
    dGammadU1221 = dGammadU[I6+9,I1]
    dGammadU3221 = dGammadU[I7+9,I1]
    dGammadU3121 = dGammadU[I8+9,I1]
    dGammadU2121 = dGammadU[I9+9,I1]
    dGammadU1131 = dGammadU[I1+18,I1]
    dGammadU2231 = dGammadU[I2+18,I1]
    dGammadU3331 = dGammadU[I3+18,I1]
    dGammadU2331 = dGammadU[I4+18,I1]
    dGammadU1331 = dGammadU[I5+18,I1]
    dGammadU1231 = dGammadU[I6+18,I1]
    dGammadU3231 = dGammadU[I7+18,I1]
    dGammadU3131 = dGammadU[I8+18,I1]
    dGammadU2131 = dGammadU[I9+18,I1]
    dGammadU1112 = dGammadU[I1,I2]
    dGammadU2212 = dGammadU[I2,I2]
    dGammadU3312 = dGammadU[I3,I2]
    dGammadU2312 = dGammadU[I4,I2]
    dGammadU1312 = dGammadU[I5,I2]
    dGammadU1212 = dGammadU[I6,I2]
    dGammadU3212 = dGammadU[I7,I2]
    dGammadU3112 = dGammadU[I8,I2]
    dGammadU2112 = dGammadU[I9,I2]
    dGammadU1122 = dGammadU[I1+9,I2]
    dGammadU2222 = dGammadU[I2+9,I2]
    dGammadU3322 = dGammadU[I3+9,I2]
    dGammadU2322 = dGammadU[I4+9,I2]
    dGammadU1322 = dGammadU[I5+9,I2]
    dGammadU1222 = dGammadU[I6+9,I2]
    dGammadU3222 = dGammadU[I7+9,I2]
    dGammadU3122 = dGammadU[I8+9,I2]
    dGammadU2122 = dGammadU[I9+9,I2]
    dGammadU1132 = dGammadU[I1+18,I2]
    dGammadU2232 = dGammadU[I2+18,I2]
    dGammadU3332 = dGammadU[I3+18,I2]
    dGammadU2332 = dGammadU[I4+18,I2]
    dGammadU1332 = dGammadU[I5+18,I2]
    dGammadU1232 = dGammadU[I6+18,I2]
    dGammadU3232 = dGammadU[I7+18,I2]
    dGammadU3132 = dGammadU[I8+18,I2]
    dGammadU2132 = dGammadU[I9+18,I2]
    dGammadU1113 = dGammadU[I1,I3]
    dGammadU2213 = dGammadU[I2,I3]
    dGammadU3313 = dGammadU[I3,I3]
    dGammadU2313 = dGammadU[I4,I3]
    dGammadU1313 = dGammadU[I5,I3]
    dGammadU1213 = dGammadU[I6,I3]
    dGammadU3213 = dGammadU[I7,I3]
    dGammadU3113 = dGammadU[I8,I3]
    dGammadU2113 = dGammadU[I9,I3]
    dGammadU1123 = dGammadU[I1+9,I3]
    dGammadU2223 = dGammadU[I2+9,I3]
    dGammadU3323 = dGammadU[I3+9,I3]
    dGammadU2323 = dGammadU[I4+9,I3]
    dGammadU1323 = dGammadU[I5+9,I3]
    dGammadU1223 = dGammadU[I6+9,I3]
    dGammadU3223 = dGammadU[I7+9,I3]
    dGammadU3123 = dGammadU[I8+9,I3]
    dGammadU2123 = dGammadU[I9+9,I3]
    dGammadU1133 = dGammadU[I1+18,I3]
    dGammadU2233 = dGammadU[I2+18,I3]
    dGammadU3333 = dGammadU[I3+18,I3]
    dGammadU2333 = dGammadU[I4+18,I3]
    dGammadU1333 = dGammadU[I5+18,I3]
    dGammadU1233 = dGammadU[I6+18,I3]
    dGammadU3233 = dGammadU[I7+18,I3]
    dGammadU3133 = dGammadU[I8+18,I3]
    dGammadU2133 = dGammadU[I9+18,I3]
    dGammadU1114 = dGammadU[I1,I4]
    dGammadU3114 = dGammadU[I8,I4]
    dGammadU2114 = dGammadU[I9,I4]
    dGammadU1124 = dGammadU[I1+9,I4]
    dGammadU3124 = dGammadU[I8+9,I4]
    dGammadU2124 = dGammadU[I9+9,I4]
    dGammadU1134 = dGammadU[I1+18,I4]
    dGammadU3134 = dGammadU[I8+18,I4]
    dGammadU2134 = dGammadU[I9+18,I4]
    dGammadU2215 = dGammadU[I2,I5]
    dGammadU1215 = dGammadU[I6,I5]
    dGammadU3215 = dGammadU[I7,I5]
    dGammadU2225 = dGammadU[I2+9,I5]
    dGammadU1225 = dGammadU[I6+9,I5]
    dGammadU3225 = dGammadU[I7+9,I5]
    dGammadU2235 = dGammadU[I2+18,I5]
    dGammadU1235 = dGammadU[I6+18,I5]
    dGammadU3235 = dGammadU[I7+18,I5]
    dGammadU3316 = dGammadU[I3,I6]
    dGammadU2316 = dGammadU[I4,I6]
    dGammadU1316 = dGammadU[I5,I6]
    dGammadU3326 = dGammadU[I3+9,I6]
    dGammadU2326 = dGammadU[I4+9,I6]
    dGammadU1326 = dGammadU[I5+9,I6]
    dGammadU3336 = dGammadU[I3+18,I6]
    dGammadU2336 = dGammadU[I4+18,I6]
    dGammadU1336 = dGammadU[I5+18,I6]
    dGammadU3317 = dGammadU[I3,I7]
    dGammadU2317 = dGammadU[I4,I7]
    dGammadU1317 = dGammadU[I5,I7]
    dGammadU3327 = dGammadU[I3+9,I7]
    dGammadU2327 = dGammadU[I4+9,I7]
    dGammadU1327 = dGammadU[I5+9,I7]
    dGammadU3337 = dGammadU[I3+18,I7]
    dGammadU2337 = dGammadU[I4+18,I7]
    dGammadU1337 = dGammadU[I5+18,I7]
    dGammadU3318 = dGammadU[I3,I8]
    dGammadU2318 = dGammadU[I4,I8]
    dGammadU1318 = dGammadU[I5,I8]
    dGammadU3328 = dGammadU[I3+9,I8]
    dGammadU2328 = dGammadU[I4+9,I8]
    dGammadU1328 = dGammadU[I5+9,I8]
    dGammadU3338 = dGammadU[I3+18,I8]
    dGammadU2338 = dGammadU[I4+18,I8]
    dGammadU1338 = dGammadU[I5+18,I8]
    dGammadU2219 = dGammadU[I2,I9]
    dGammadU1219 = dGammadU[I6,I9]
    dGammadU3219 = dGammadU[I7,I9]
    dGammadU2229 = dGammadU[I2+9,I9]
    dGammadU1229 = dGammadU[I6+9,I9]
    dGammadU3229 = dGammadU[I7+9,I9]
    dGammadU2239 = dGammadU[I2+18,I9]
    dGammadU1239 = dGammadU[I6+18,I9]
    dGammadU3239 = dGammadU[I7+18,I9]
    dGammadU22110 = dGammadU[I2,I10]
    dGammadU12110 = dGammadU[I6,I10]
    dGammadU32110 = dGammadU[I7,I10]
    dGammadU22210 = dGammadU[I2+9,I10]
    dGammadU12210 = dGammadU[I6+9,I10]
    dGammadU32210 = dGammadU[I7+9,I10]
    dGammadU22310 = dGammadU[I2+18,I10]
    dGammadU12310 = dGammadU[I6+18,I10]
    dGammadU32310 = dGammadU[I7+18,I10]
    dGammadU11111 = dGammadU[I1,I11]
    dGammadU31111 = dGammadU[I8,I11]
    dGammadU21111 = dGammadU[I9,I11]
    dGammadU11211 = dGammadU[I1+9,I11]
    dGammadU31211 = dGammadU[I8+9,I11]
    dGammadU21211 = dGammadU[I9+9,I11]
    dGammadU11311 = dGammadU[I1+18,I11]
    dGammadU31311 = dGammadU[I8+18,I11]
    dGammadU21311 = dGammadU[I9+18,I11]
    dGammadU11112 = dGammadU[I1,I12]
    dGammadU31112 = dGammadU[I8,I12]
    dGammadU21112 = dGammadU[I9,I12]
    dGammadU11212 = dGammadU[I1+9,I12]
    dGammadU31212 = dGammadU[I8+9,I12]
    dGammadU21212 = dGammadU[I9+9,I12]
    dGammadU11312 = dGammadU[I1+18,I12]
    dGammadU31312 = dGammadU[I8+18,I12]
    dGammadU21312 = dGammadU[I9+18,I12]

    #Extract components of dpk2dC
    dpk2dC1111 = dpk2dC[I1,I1]
    dpk2dC2211 = dpk2dC[I2,I1]
    dpk2dC3311 = dpk2dC[I3,I1]
    dpk2dC2311 = dpk2dC[I4,I1]
    dpk2dC1311 = dpk2dC[I5,I1]
    dpk2dC1211 = dpk2dC[I6,I1]
    dpk2dC3211 = dpk2dC[I7,I1]
    dpk2dC3111 = dpk2dC[I8,I1]
    dpk2dC2111 = dpk2dC[I9,I1]
    dpk2dC1121 = dpk2dC[I1+9,I1]
    dpk2dC2221 = dpk2dC[I2+9,I1]
    dpk2dC3321 = dpk2dC[I3+9,I1]
    dpk2dC2321 = dpk2dC[I4+9,I1]
    dpk2dC1321 = dpk2dC[I5+9,I1]
    dpk2dC1221 = dpk2dC[I6+9,I1]
    dpk2dC3221 = dpk2dC[I7+9,I1]
    dpk2dC3121 = dpk2dC[I8+9,I1]
    dpk2dC2121 = dpk2dC[I9+9,I1]
    dpk2dC1131 = dpk2dC[I1+18,I1]
    dpk2dC2231 = dpk2dC[I2+18,I1]
    dpk2dC3331 = dpk2dC[I3+18,I1]
    dpk2dC2331 = dpk2dC[I4+18,I1]
    dpk2dC1331 = dpk2dC[I5+18,I1]
    dpk2dC1231 = dpk2dC[I6+18,I1]
    dpk2dC3231 = dpk2dC[I7+18,I1]
    dpk2dC3131 = dpk2dC[I8+18,I1]
    dpk2dC2131 = dpk2dC[I9+18,I1]
    dpk2dC1112 = dpk2dC[I1,I2]
    dpk2dC2212 = dpk2dC[I2,I2]
    dpk2dC3312 = dpk2dC[I3,I2]
    dpk2dC2312 = dpk2dC[I4,I2]
    dpk2dC1312 = dpk2dC[I5,I2]
    dpk2dC1212 = dpk2dC[I6,I2]
    dpk2dC3212 = dpk2dC[I7,I2]
    dpk2dC3112 = dpk2dC[I8,I2]
    dpk2dC2112 = dpk2dC[I9,I2]
    dpk2dC1122 = dpk2dC[I1+9,I2]
    dpk2dC2222 = dpk2dC[I2+9,I2]
    dpk2dC3322 = dpk2dC[I3+9,I2]
    dpk2dC2322 = dpk2dC[I4+9,I2]
    dpk2dC1322 = dpk2dC[I5+9,I2]
    dpk2dC1222 = dpk2dC[I6+9,I2]
    dpk2dC3222 = dpk2dC[I7+9,I2]
    dpk2dC3122 = dpk2dC[I8+9,I2]
    dpk2dC2122 = dpk2dC[I9+9,I2]
    dpk2dC1132 = dpk2dC[I1+18,I2]
    dpk2dC2232 = dpk2dC[I2+18,I2]
    dpk2dC3332 = dpk2dC[I3+18,I2]
    dpk2dC2332 = dpk2dC[I4+18,I2]
    dpk2dC1332 = dpk2dC[I5+18,I2]
    dpk2dC1232 = dpk2dC[I6+18,I2]
    dpk2dC3232 = dpk2dC[I7+18,I2]
    dpk2dC3132 = dpk2dC[I8+18,I2]
    dpk2dC2132 = dpk2dC[I9+18,I2]
    dpk2dC1113 = dpk2dC[I1,I3]
    dpk2dC2213 = dpk2dC[I2,I3]
    dpk2dC3313 = dpk2dC[I3,I3]
    dpk2dC2313 = dpk2dC[I4,I3]
    dpk2dC1313 = dpk2dC[I5,I3]
    dpk2dC1213 = dpk2dC[I6,I3]
    dpk2dC3213 = dpk2dC[I7,I3]
    dpk2dC3113 = dpk2dC[I8,I3]
    dpk2dC2113 = dpk2dC[I9,I3]
    dpk2dC1123 = dpk2dC[I1+9,I3]
    dpk2dC2223 = dpk2dC[I2+9,I3]
    dpk2dC3323 = dpk2dC[I3+9,I3]
    dpk2dC2323 = dpk2dC[I4+9,I3]
    dpk2dC1323 = dpk2dC[I5+9,I3]
    dpk2dC1223 = dpk2dC[I6+9,I3]
    dpk2dC3223 = dpk2dC[I7+9,I3]
    dpk2dC3123 = dpk2dC[I8+9,I3]
    dpk2dC2123 = dpk2dC[I9+9,I3]
    dpk2dC1133 = dpk2dC[I1+18,I3]
    dpk2dC2233 = dpk2dC[I2+18,I3]
    dpk2dC3333 = dpk2dC[I3+18,I3]
    dpk2dC2333 = dpk2dC[I4+18,I3]
    dpk2dC1333 = dpk2dC[I5+18,I3]
    dpk2dC1233 = dpk2dC[I6+18,I3]
    dpk2dC3233 = dpk2dC[I7+18,I3]
    dpk2dC3133 = dpk2dC[I8+18,I3]
    dpk2dC2133 = dpk2dC[I9+18,I3]

    #Extract components of dpk2dPsi
    dpk2dPsi1111 = dpk2dPsi[I1,I1]
    dpk2dPsi2211 = dpk2dPsi[I2,I1]
    dpk2dPsi3311 = dpk2dPsi[I3,I1]
    dpk2dPsi2311 = dpk2dPsi[I4,I1]
    dpk2dPsi1311 = dpk2dPsi[I5,I1]
    dpk2dPsi1211 = dpk2dPsi[I6,I1]
    dpk2dPsi3211 = dpk2dPsi[I7,I1]
    dpk2dPsi3111 = dpk2dPsi[I8,I1]
    dpk2dPsi2111 = dpk2dPsi[I9,I1]
    dpk2dPsi1121 = dpk2dPsi[I1+9,I1]
    dpk2dPsi2221 = dpk2dPsi[I2+9,I1]
    dpk2dPsi3321 = dpk2dPsi[I3+9,I1]
    dpk2dPsi2321 = dpk2dPsi[I4+9,I1]
    dpk2dPsi1321 = dpk2dPsi[I5+9,I1]
    dpk2dPsi1221 = dpk2dPsi[I6+9,I1]
    dpk2dPsi3221 = dpk2dPsi[I7+9,I1]
    dpk2dPsi3121 = dpk2dPsi[I8+9,I1]
    dpk2dPsi2121 = dpk2dPsi[I9+9,I1]
    dpk2dPsi1131 = dpk2dPsi[I1+18,I1]
    dpk2dPsi2231 = dpk2dPsi[I2+18,I1]
    dpk2dPsi3331 = dpk2dPsi[I3+18,I1]
    dpk2dPsi2331 = dpk2dPsi[I4+18,I1]
    dpk2dPsi1331 = dpk2dPsi[I5+18,I1]
    dpk2dPsi1231 = dpk2dPsi[I6+18,I1]
    dpk2dPsi3231 = dpk2dPsi[I7+18,I1]
    dpk2dPsi3131 = dpk2dPsi[I8+18,I1]
    dpk2dPsi2131 = dpk2dPsi[I9+18,I1]
    dpk2dPsi1112 = dpk2dPsi[I1,I2]
    dpk2dPsi2212 = dpk2dPsi[I2,I2]
    dpk2dPsi3312 = dpk2dPsi[I3,I2]
    dpk2dPsi2312 = dpk2dPsi[I4,I2]
    dpk2dPsi1312 = dpk2dPsi[I5,I2]
    dpk2dPsi1212 = dpk2dPsi[I6,I2]
    dpk2dPsi3212 = dpk2dPsi[I7,I2]
    dpk2dPsi3112 = dpk2dPsi[I8,I2]
    dpk2dPsi2112 = dpk2dPsi[I9,I2]
    dpk2dPsi1122 = dpk2dPsi[I1+9,I2]
    dpk2dPsi2222 = dpk2dPsi[I2+9,I2]
    dpk2dPsi3322 = dpk2dPsi[I3+9,I2]
    dpk2dPsi2322 = dpk2dPsi[I4+9,I2]
    dpk2dPsi1322 = dpk2dPsi[I5+9,I2]
    dpk2dPsi1222 = dpk2dPsi[I6+9,I2]
    dpk2dPsi3222 = dpk2dPsi[I7+9,I2]
    dpk2dPsi3122 = dpk2dPsi[I8+9,I2]
    dpk2dPsi2122 = dpk2dPsi[I9+9,I2]
    dpk2dPsi1132 = dpk2dPsi[I1+18,I2]
    dpk2dPsi2232 = dpk2dPsi[I2+18,I2]
    dpk2dPsi3332 = dpk2dPsi[I3+18,I2]
    dpk2dPsi2332 = dpk2dPsi[I4+18,I2]
    dpk2dPsi1332 = dpk2dPsi[I5+18,I2]
    dpk2dPsi1232 = dpk2dPsi[I6+18,I2]
    dpk2dPsi3232 = dpk2dPsi[I7+18,I2]
    dpk2dPsi3132 = dpk2dPsi[I8+18,I2]
    dpk2dPsi2132 = dpk2dPsi[I9+18,I2]
    dpk2dPsi1113 = dpk2dPsi[I1,I3]
    dpk2dPsi2213 = dpk2dPsi[I2,I3]
    dpk2dPsi3313 = dpk2dPsi[I3,I3]
    dpk2dPsi2313 = dpk2dPsi[I4,I3]
    dpk2dPsi1313 = dpk2dPsi[I5,I3]
    dpk2dPsi1213 = dpk2dPsi[I6,I3]
    dpk2dPsi3213 = dpk2dPsi[I7,I3]
    dpk2dPsi3113 = dpk2dPsi[I8,I3]
    dpk2dPsi2113 = dpk2dPsi[I9,I3]
    dpk2dPsi1123 = dpk2dPsi[I1+9,I3]
    dpk2dPsi2223 = dpk2dPsi[I2+9,I3]
    dpk2dPsi3323 = dpk2dPsi[I3+9,I3]
    dpk2dPsi2323 = dpk2dPsi[I4+9,I3]
    dpk2dPsi1323 = dpk2dPsi[I5+9,I3]
    dpk2dPsi1223 = dpk2dPsi[I6+9,I3]
    dpk2dPsi3223 = dpk2dPsi[I7+9,I3]
    dpk2dPsi3123 = dpk2dPsi[I8+9,I3]
    dpk2dPsi2123 = dpk2dPsi[I9+9,I3]
    dpk2dPsi1133 = dpk2dPsi[I1+18,I3]
    dpk2dPsi2233 = dpk2dPsi[I2+18,I3]
    dpk2dPsi3333 = dpk2dPsi[I3+18,I3]
    dpk2dPsi2333 = dpk2dPsi[I4+18,I3]
    dpk2dPsi1333 = dpk2dPsi[I5+18,I3]
    dpk2dPsi1233 = dpk2dPsi[I6+18,I3]
    dpk2dPsi3233 = dpk2dPsi[I7+18,I3]
    dpk2dPsi3133 = dpk2dPsi[I8+18,I3]
    dpk2dPsi2133 = dpk2dPsi[I9+18,I3]

    #Extract components of dpk2dGamma
    dpk2dGamma11111 = dpk2dGamma[I1,I1]
    dpk2dGamma22111 = dpk2dGamma[I2,I1]
    dpk2dGamma33111 = dpk2dGamma[I3,I1]
    dpk2dGamma23111 = dpk2dGamma[I4,I1]
    dpk2dGamma13111 = dpk2dGamma[I5,I1]
    dpk2dGamma12111 = dpk2dGamma[I6,I1]
    dpk2dGamma32111 = dpk2dGamma[I7,I1]
    dpk2dGamma31111 = dpk2dGamma[I8,I1]
    dpk2dGamma21111 = dpk2dGamma[I9,I1]
    dpk2dGamma11121 = dpk2dGamma[I1+9,I1]
    dpk2dGamma22121 = dpk2dGamma[I2+9,I1]
    dpk2dGamma33121 = dpk2dGamma[I3+9,I1]
    dpk2dGamma23121 = dpk2dGamma[I4+9,I1]
    dpk2dGamma13121 = dpk2dGamma[I5+9,I1]
    dpk2dGamma12121 = dpk2dGamma[I6+9,I1]
    dpk2dGamma32121 = dpk2dGamma[I7+9,I1]
    dpk2dGamma31121 = dpk2dGamma[I8+9,I1]
    dpk2dGamma21121 = dpk2dGamma[I9+9,I1]
    dpk2dGamma11131 = dpk2dGamma[I1+18,I1]
    dpk2dGamma22131 = dpk2dGamma[I2+18,I1]
    dpk2dGamma33131 = dpk2dGamma[I3+18,I1]
    dpk2dGamma23131 = dpk2dGamma[I4+18,I1]
    dpk2dGamma13131 = dpk2dGamma[I5+18,I1]
    dpk2dGamma12131 = dpk2dGamma[I6+18,I1]
    dpk2dGamma32131 = dpk2dGamma[I7+18,I1]
    dpk2dGamma31131 = dpk2dGamma[I8+18,I1]
    dpk2dGamma21131 = dpk2dGamma[I9+18,I1]
    dpk2dGamma11211 = dpk2dGamma[I1+27,I1]
    dpk2dGamma22211 = dpk2dGamma[I2+27,I1]
    dpk2dGamma33211 = dpk2dGamma[I3+27,I1]
    dpk2dGamma23211 = dpk2dGamma[I4+27,I1]
    dpk2dGamma13211 = dpk2dGamma[I5+27,I1]
    dpk2dGamma12211 = dpk2dGamma[I6+27,I1]
    dpk2dGamma32211 = dpk2dGamma[I7+27,I1]
    dpk2dGamma31211 = dpk2dGamma[I8+27,I1]
    dpk2dGamma21211 = dpk2dGamma[I9+27,I1]
    dpk2dGamma11221 = dpk2dGamma[I1+36,I1]
    dpk2dGamma22221 = dpk2dGamma[I2+36,I1]
    dpk2dGamma33221 = dpk2dGamma[I3+36,I1]
    dpk2dGamma23221 = dpk2dGamma[I4+36,I1]
    dpk2dGamma13221 = dpk2dGamma[I5+36,I1]
    dpk2dGamma12221 = dpk2dGamma[I6+36,I1]
    dpk2dGamma32221 = dpk2dGamma[I7+36,I1]
    dpk2dGamma31221 = dpk2dGamma[I8+36,I1]
    dpk2dGamma21221 = dpk2dGamma[I9+36,I1]
    dpk2dGamma11231 = dpk2dGamma[I1+45,I1]
    dpk2dGamma22231 = dpk2dGamma[I2+45,I1]
    dpk2dGamma33231 = dpk2dGamma[I3+45,I1]
    dpk2dGamma23231 = dpk2dGamma[I4+45,I1]
    dpk2dGamma13231 = dpk2dGamma[I5+45,I1]
    dpk2dGamma12231 = dpk2dGamma[I6+45,I1]
    dpk2dGamma32231 = dpk2dGamma[I7+45,I1]
    dpk2dGamma31231 = dpk2dGamma[I8+45,I1]
    dpk2dGamma21231 = dpk2dGamma[I9+45,I1]
    dpk2dGamma11311 = dpk2dGamma[I1+54,I1]
    dpk2dGamma22311 = dpk2dGamma[I2+54,I1]
    dpk2dGamma33311 = dpk2dGamma[I3+54,I1]
    dpk2dGamma23311 = dpk2dGamma[I4+54,I1]
    dpk2dGamma13311 = dpk2dGamma[I5+54,I1]
    dpk2dGamma12311 = dpk2dGamma[I6+54,I1]
    dpk2dGamma32311 = dpk2dGamma[I7+54,I1]
    dpk2dGamma31311 = dpk2dGamma[I8+54,I1]
    dpk2dGamma21311 = dpk2dGamma[I9+54,I1]
    dpk2dGamma11321 = dpk2dGamma[I1+63,I1]
    dpk2dGamma22321 = dpk2dGamma[I2+63,I1]
    dpk2dGamma33321 = dpk2dGamma[I3+63,I1]
    dpk2dGamma23321 = dpk2dGamma[I4+63,I1]
    dpk2dGamma13321 = dpk2dGamma[I5+63,I1]
    dpk2dGamma12321 = dpk2dGamma[I6+63,I1]
    dpk2dGamma32321 = dpk2dGamma[I7+63,I1]
    dpk2dGamma31321 = dpk2dGamma[I8+63,I1]
    dpk2dGamma21321 = dpk2dGamma[I9+63,I1]
    dpk2dGamma11331 = dpk2dGamma[I1+72,I1]
    dpk2dGamma22331 = dpk2dGamma[I2+72,I1]
    dpk2dGamma33331 = dpk2dGamma[I3+72,I1]
    dpk2dGamma23331 = dpk2dGamma[I4+72,I1]
    dpk2dGamma13331 = dpk2dGamma[I5+72,I1]
    dpk2dGamma12331 = dpk2dGamma[I6+72,I1]
    dpk2dGamma32331 = dpk2dGamma[I7+72,I1]
    dpk2dGamma31331 = dpk2dGamma[I8+72,I1]
    dpk2dGamma21331 = dpk2dGamma[I9+72,I1]
    dpk2dGamma11112 = dpk2dGamma[I1,I2]
    dpk2dGamma22112 = dpk2dGamma[I2,I2]
    dpk2dGamma33112 = dpk2dGamma[I3,I2]
    dpk2dGamma23112 = dpk2dGamma[I4,I2]
    dpk2dGamma13112 = dpk2dGamma[I5,I2]
    dpk2dGamma12112 = dpk2dGamma[I6,I2]
    dpk2dGamma32112 = dpk2dGamma[I7,I2]
    dpk2dGamma31112 = dpk2dGamma[I8,I2]
    dpk2dGamma21112 = dpk2dGamma[I9,I2]
    dpk2dGamma11122 = dpk2dGamma[I1+9,I2]
    dpk2dGamma22122 = dpk2dGamma[I2+9,I2]
    dpk2dGamma33122 = dpk2dGamma[I3+9,I2]
    dpk2dGamma23122 = dpk2dGamma[I4+9,I2]
    dpk2dGamma13122 = dpk2dGamma[I5+9,I2]
    dpk2dGamma12122 = dpk2dGamma[I6+9,I2]
    dpk2dGamma32122 = dpk2dGamma[I7+9,I2]
    dpk2dGamma31122 = dpk2dGamma[I8+9,I2]
    dpk2dGamma21122 = dpk2dGamma[I9+9,I2]
    dpk2dGamma11132 = dpk2dGamma[I1+18,I2]
    dpk2dGamma22132 = dpk2dGamma[I2+18,I2]
    dpk2dGamma33132 = dpk2dGamma[I3+18,I2]
    dpk2dGamma23132 = dpk2dGamma[I4+18,I2]
    dpk2dGamma13132 = dpk2dGamma[I5+18,I2]
    dpk2dGamma12132 = dpk2dGamma[I6+18,I2]
    dpk2dGamma32132 = dpk2dGamma[I7+18,I2]
    dpk2dGamma31132 = dpk2dGamma[I8+18,I2]
    dpk2dGamma21132 = dpk2dGamma[I9+18,I2]
    dpk2dGamma11212 = dpk2dGamma[I1+27,I2]
    dpk2dGamma22212 = dpk2dGamma[I2+27,I2]
    dpk2dGamma33212 = dpk2dGamma[I3+27,I2]
    dpk2dGamma23212 = dpk2dGamma[I4+27,I2]
    dpk2dGamma13212 = dpk2dGamma[I5+27,I2]
    dpk2dGamma12212 = dpk2dGamma[I6+27,I2]
    dpk2dGamma32212 = dpk2dGamma[I7+27,I2]
    dpk2dGamma31212 = dpk2dGamma[I8+27,I2]
    dpk2dGamma21212 = dpk2dGamma[I9+27,I2]
    dpk2dGamma11222 = dpk2dGamma[I1+36,I2]
    dpk2dGamma22222 = dpk2dGamma[I2+36,I2]
    dpk2dGamma33222 = dpk2dGamma[I3+36,I2]
    dpk2dGamma23222 = dpk2dGamma[I4+36,I2]
    dpk2dGamma13222 = dpk2dGamma[I5+36,I2]
    dpk2dGamma12222 = dpk2dGamma[I6+36,I2]
    dpk2dGamma32222 = dpk2dGamma[I7+36,I2]
    dpk2dGamma31222 = dpk2dGamma[I8+36,I2]
    dpk2dGamma21222 = dpk2dGamma[I9+36,I2]
    dpk2dGamma11232 = dpk2dGamma[I1+45,I2]
    dpk2dGamma22232 = dpk2dGamma[I2+45,I2]
    dpk2dGamma33232 = dpk2dGamma[I3+45,I2]
    dpk2dGamma23232 = dpk2dGamma[I4+45,I2]
    dpk2dGamma13232 = dpk2dGamma[I5+45,I2]
    dpk2dGamma12232 = dpk2dGamma[I6+45,I2]
    dpk2dGamma32232 = dpk2dGamma[I7+45,I2]
    dpk2dGamma31232 = dpk2dGamma[I8+45,I2]
    dpk2dGamma21232 = dpk2dGamma[I9+45,I2]
    dpk2dGamma11312 = dpk2dGamma[I1+54,I2]
    dpk2dGamma22312 = dpk2dGamma[I2+54,I2]
    dpk2dGamma33312 = dpk2dGamma[I3+54,I2]
    dpk2dGamma23312 = dpk2dGamma[I4+54,I2]
    dpk2dGamma13312 = dpk2dGamma[I5+54,I2]
    dpk2dGamma12312 = dpk2dGamma[I6+54,I2]
    dpk2dGamma32312 = dpk2dGamma[I7+54,I2]
    dpk2dGamma31312 = dpk2dGamma[I8+54,I2]
    dpk2dGamma21312 = dpk2dGamma[I9+54,I2]
    dpk2dGamma11322 = dpk2dGamma[I1+63,I2]
    dpk2dGamma22322 = dpk2dGamma[I2+63,I2]
    dpk2dGamma33322 = dpk2dGamma[I3+63,I2]
    dpk2dGamma23322 = dpk2dGamma[I4+63,I2]
    dpk2dGamma13322 = dpk2dGamma[I5+63,I2]
    dpk2dGamma12322 = dpk2dGamma[I6+63,I2]
    dpk2dGamma32322 = dpk2dGamma[I7+63,I2]
    dpk2dGamma31322 = dpk2dGamma[I8+63,I2]
    dpk2dGamma21322 = dpk2dGamma[I9+63,I2]
    dpk2dGamma11332 = dpk2dGamma[I1+72,I2]
    dpk2dGamma22332 = dpk2dGamma[I2+72,I2]
    dpk2dGamma33332 = dpk2dGamma[I3+72,I2]
    dpk2dGamma23332 = dpk2dGamma[I4+72,I2]
    dpk2dGamma13332 = dpk2dGamma[I5+72,I2]
    dpk2dGamma12332 = dpk2dGamma[I6+72,I2]
    dpk2dGamma32332 = dpk2dGamma[I7+72,I2]
    dpk2dGamma31332 = dpk2dGamma[I8+72,I2]
    dpk2dGamma21332 = dpk2dGamma[I9+72,I2]
    dpk2dGamma11113 = dpk2dGamma[I1,I3]
    dpk2dGamma22113 = dpk2dGamma[I2,I3]
    dpk2dGamma33113 = dpk2dGamma[I3,I3]
    dpk2dGamma23113 = dpk2dGamma[I4,I3]
    dpk2dGamma13113 = dpk2dGamma[I5,I3]
    dpk2dGamma12113 = dpk2dGamma[I6,I3]
    dpk2dGamma32113 = dpk2dGamma[I7,I3]
    dpk2dGamma31113 = dpk2dGamma[I8,I3]
    dpk2dGamma21113 = dpk2dGamma[I9,I3]
    dpk2dGamma11123 = dpk2dGamma[I1+9,I3]
    dpk2dGamma22123 = dpk2dGamma[I2+9,I3]
    dpk2dGamma33123 = dpk2dGamma[I3+9,I3]
    dpk2dGamma23123 = dpk2dGamma[I4+9,I3]
    dpk2dGamma13123 = dpk2dGamma[I5+9,I3]
    dpk2dGamma12123 = dpk2dGamma[I6+9,I3]
    dpk2dGamma32123 = dpk2dGamma[I7+9,I3]
    dpk2dGamma31123 = dpk2dGamma[I8+9,I3]
    dpk2dGamma21123 = dpk2dGamma[I9+9,I3]
    dpk2dGamma11133 = dpk2dGamma[I1+18,I3]
    dpk2dGamma22133 = dpk2dGamma[I2+18,I3]
    dpk2dGamma33133 = dpk2dGamma[I3+18,I3]
    dpk2dGamma23133 = dpk2dGamma[I4+18,I3]
    dpk2dGamma13133 = dpk2dGamma[I5+18,I3]
    dpk2dGamma12133 = dpk2dGamma[I6+18,I3]
    dpk2dGamma32133 = dpk2dGamma[I7+18,I3]
    dpk2dGamma31133 = dpk2dGamma[I8+18,I3]
    dpk2dGamma21133 = dpk2dGamma[I9+18,I3]
    dpk2dGamma11213 = dpk2dGamma[I1+27,I3]
    dpk2dGamma22213 = dpk2dGamma[I2+27,I3]
    dpk2dGamma33213 = dpk2dGamma[I3+27,I3]
    dpk2dGamma23213 = dpk2dGamma[I4+27,I3]
    dpk2dGamma13213 = dpk2dGamma[I5+27,I3]
    dpk2dGamma12213 = dpk2dGamma[I6+27,I3]
    dpk2dGamma32213 = dpk2dGamma[I7+27,I3]
    dpk2dGamma31213 = dpk2dGamma[I8+27,I3]
    dpk2dGamma21213 = dpk2dGamma[I9+27,I3]
    dpk2dGamma11223 = dpk2dGamma[I1+36,I3]
    dpk2dGamma22223 = dpk2dGamma[I2+36,I3]
    dpk2dGamma33223 = dpk2dGamma[I3+36,I3]
    dpk2dGamma23223 = dpk2dGamma[I4+36,I3]
    dpk2dGamma13223 = dpk2dGamma[I5+36,I3]
    dpk2dGamma12223 = dpk2dGamma[I6+36,I3]
    dpk2dGamma32223 = dpk2dGamma[I7+36,I3]
    dpk2dGamma31223 = dpk2dGamma[I8+36,I3]
    dpk2dGamma21223 = dpk2dGamma[I9+36,I3]
    dpk2dGamma11233 = dpk2dGamma[I1+45,I3]
    dpk2dGamma22233 = dpk2dGamma[I2+45,I3]
    dpk2dGamma33233 = dpk2dGamma[I3+45,I3]
    dpk2dGamma23233 = dpk2dGamma[I4+45,I3]
    dpk2dGamma13233 = dpk2dGamma[I5+45,I3]
    dpk2dGamma12233 = dpk2dGamma[I6+45,I3]
    dpk2dGamma32233 = dpk2dGamma[I7+45,I3]
    dpk2dGamma31233 = dpk2dGamma[I8+45,I3]
    dpk2dGamma21233 = dpk2dGamma[I9+45,I3]
    dpk2dGamma11313 = dpk2dGamma[I1+54,I3]
    dpk2dGamma22313 = dpk2dGamma[I2+54,I3]
    dpk2dGamma33313 = dpk2dGamma[I3+54,I3]
    dpk2dGamma23313 = dpk2dGamma[I4+54,I3]
    dpk2dGamma13313 = dpk2dGamma[I5+54,I3]
    dpk2dGamma12313 = dpk2dGamma[I6+54,I3]
    dpk2dGamma32313 = dpk2dGamma[I7+54,I3]
    dpk2dGamma31313 = dpk2dGamma[I8+54,I3]
    dpk2dGamma21313 = dpk2dGamma[I9+54,I3]
    dpk2dGamma11323 = dpk2dGamma[I1+63,I3]
    dpk2dGamma22323 = dpk2dGamma[I2+63,I3]
    dpk2dGamma33323 = dpk2dGamma[I3+63,I3]
    dpk2dGamma23323 = dpk2dGamma[I4+63,I3]
    dpk2dGamma13323 = dpk2dGamma[I5+63,I3]
    dpk2dGamma12323 = dpk2dGamma[I6+63,I3]
    dpk2dGamma32323 = dpk2dGamma[I7+63,I3]
    dpk2dGamma31323 = dpk2dGamma[I8+63,I3]
    dpk2dGamma21323 = dpk2dGamma[I9+63,I3]
    dpk2dGamma11333 = dpk2dGamma[I1+72,I3]
    dpk2dGamma22333 = dpk2dGamma[I2+72,I3]
    dpk2dGamma33333 = dpk2dGamma[I3+72,I3]
    dpk2dGamma23333 = dpk2dGamma[I4+72,I3]
    dpk2dGamma13333 = dpk2dGamma[I5+72,I3]
    dpk2dGamma12333 = dpk2dGamma[I6+72,I3]
    dpk2dGamma32333 = dpk2dGamma[I7+72,I3]
    dpk2dGamma31333 = dpk2dGamma[I8+72,I3]
    dpk2dGamma21333 = dpk2dGamma[I9+72,I3]

    #Compute Tangent
    
    #Column 1
    dpk2dU[I1,I1] = dCdU111*dpk2dC1111 + dCdU121*dpk2dC1112 + dCdU131*dpk2dC1113 + dCdU211*dpk2dC1121 + dCdU221*dpk2dC1122 + dCdU231*dpk2dC1123 + dCdU311*dpk2dC1131 + dCdU321*dpk2dC1132 + dCdU331*dpk2dC1133 + dGammadU1111*dpk2dGamma11111 + dGammadU1121*dpk2dGamma11112 + dGammadU1131*dpk2dGamma11113 + dGammadU1211*dpk2dGamma11121 + dGammadU1221*dpk2dGamma11122 + dGammadU1231*dpk2dGamma11123 + dGammadU1311*dpk2dGamma11131 + dGammadU1321*dpk2dGamma11132 + dGammadU1331*dpk2dGamma11133 + dGammadU2111*dpk2dGamma11211 + dGammadU2121*dpk2dGamma11212 + dGammadU2131*dpk2dGamma11213 + dGammadU2211*dpk2dGamma11221 + dGammadU2221*dpk2dGamma11222 + dGammadU2231*dpk2dGamma11223 + dGammadU2311*dpk2dGamma11231 + dGammadU2321*dpk2dGamma11232 + dGammadU2331*dpk2dGamma11233 + dGammadU3111*dpk2dGamma11311 + dGammadU3121*dpk2dGamma11312 + dGammadU3131*dpk2dGamma11313 + dGammadU3211*dpk2dGamma11321 + dGammadU3221*dpk2dGamma11322 + dGammadU3231*dpk2dGamma11323 + dGammadU3311*dpk2dGamma11331 + dGammadU3321*dpk2dGamma11332 + dGammadU3331*dpk2dGamma11333 + dPhidU111*dpk2dPsi1111 + dPhidU121*dpk2dPsi1112 + dPhidU131*dpk2dPsi1113 + dPhidU211*dpk2dPsi1121 + dPhidU221*dpk2dPsi1122 + dPhidU231*dpk2dPsi1123 + dPhidU311*dpk2dPsi1131 + dPhidU321*dpk2dPsi1132 + dPhidU331*dpk2dPsi1133
    dpk2dU[I2,I1] = dCdU111*dpk2dC2211 + dCdU121*dpk2dC2212 + dCdU131*dpk2dC2213 + dCdU211*dpk2dC2221 + dCdU221*dpk2dC2222 + dCdU231*dpk2dC2223 + dCdU311*dpk2dC2231 + dCdU321*dpk2dC2232 + dCdU331*dpk2dC2233 + dGammadU1111*dpk2dGamma22111 + dGammadU1121*dpk2dGamma22112 + dGammadU1131*dpk2dGamma22113 + dGammadU1211*dpk2dGamma22121 + dGammadU1221*dpk2dGamma22122 + dGammadU1231*dpk2dGamma22123 + dGammadU1311*dpk2dGamma22131 + dGammadU1321*dpk2dGamma22132 + dGammadU1331*dpk2dGamma22133 + dGammadU2111*dpk2dGamma22211 + dGammadU2121*dpk2dGamma22212 + dGammadU2131*dpk2dGamma22213 + dGammadU2211*dpk2dGamma22221 + dGammadU2221*dpk2dGamma22222 + dGammadU2231*dpk2dGamma22223 + dGammadU2311*dpk2dGamma22231 + dGammadU2321*dpk2dGamma22232 + dGammadU2331*dpk2dGamma22233 + dGammadU3111*dpk2dGamma22311 + dGammadU3121*dpk2dGamma22312 + dGammadU3131*dpk2dGamma22313 + dGammadU3211*dpk2dGamma22321 + dGammadU3221*dpk2dGamma22322 + dGammadU3231*dpk2dGamma22323 + dGammadU3311*dpk2dGamma22331 + dGammadU3321*dpk2dGamma22332 + dGammadU3331*dpk2dGamma22333 + dPhidU111*dpk2dPsi2211 + dPhidU121*dpk2dPsi2212 + dPhidU131*dpk2dPsi2213 + dPhidU211*dpk2dPsi2221 + dPhidU221*dpk2dPsi2222 + dPhidU231*dpk2dPsi2223 + dPhidU311*dpk2dPsi2231 + dPhidU321*dpk2dPsi2232 + dPhidU331*dpk2dPsi2233
    dpk2dU[I3,I1] = dCdU111*dpk2dC3311 + dCdU121*dpk2dC3312 + dCdU131*dpk2dC3313 + dCdU211*dpk2dC3321 + dCdU221*dpk2dC3322 + dCdU231*dpk2dC3323 + dCdU311*dpk2dC3331 + dCdU321*dpk2dC3332 + dCdU331*dpk2dC3333 + dGammadU1111*dpk2dGamma33111 + dGammadU1121*dpk2dGamma33112 + dGammadU1131*dpk2dGamma33113 + dGammadU1211*dpk2dGamma33121 + dGammadU1221*dpk2dGamma33122 + dGammadU1231*dpk2dGamma33123 + dGammadU1311*dpk2dGamma33131 + dGammadU1321*dpk2dGamma33132 + dGammadU1331*dpk2dGamma33133 + dGammadU2111*dpk2dGamma33211 + dGammadU2121*dpk2dGamma33212 + dGammadU2131*dpk2dGamma33213 + dGammadU2211*dpk2dGamma33221 + dGammadU2221*dpk2dGamma33222 + dGammadU2231*dpk2dGamma33223 + dGammadU2311*dpk2dGamma33231 + dGammadU2321*dpk2dGamma33232 + dGammadU2331*dpk2dGamma33233 + dGammadU3111*dpk2dGamma33311 + dGammadU3121*dpk2dGamma33312 + dGammadU3131*dpk2dGamma33313 + dGammadU3211*dpk2dGamma33321 + dGammadU3221*dpk2dGamma33322 + dGammadU3231*dpk2dGamma33323 + dGammadU3311*dpk2dGamma33331 + dGammadU3321*dpk2dGamma33332 + dGammadU3331*dpk2dGamma33333 + dPhidU111*dpk2dPsi3311 + dPhidU121*dpk2dPsi3312 + dPhidU131*dpk2dPsi3313 + dPhidU211*dpk2dPsi3321 + dPhidU221*dpk2dPsi3322 + dPhidU231*dpk2dPsi3323 + dPhidU311*dpk2dPsi3331 + dPhidU321*dpk2dPsi3332 + dPhidU331*dpk2dPsi3333
    dpk2dU[I4,I1] = dCdU111*dpk2dC2311 + dCdU121*dpk2dC2312 + dCdU131*dpk2dC2313 + dCdU211*dpk2dC2321 + dCdU221*dpk2dC2322 + dCdU231*dpk2dC2323 + dCdU311*dpk2dC2331 + dCdU321*dpk2dC2332 + dCdU331*dpk2dC2333 + dGammadU1111*dpk2dGamma23111 + dGammadU1121*dpk2dGamma23112 + dGammadU1131*dpk2dGamma23113 + dGammadU1211*dpk2dGamma23121 + dGammadU1221*dpk2dGamma23122 + dGammadU1231*dpk2dGamma23123 + dGammadU1311*dpk2dGamma23131 + dGammadU1321*dpk2dGamma23132 + dGammadU1331*dpk2dGamma23133 + dGammadU2111*dpk2dGamma23211 + dGammadU2121*dpk2dGamma23212 + dGammadU2131*dpk2dGamma23213 + dGammadU2211*dpk2dGamma23221 + dGammadU2221*dpk2dGamma23222 + dGammadU2231*dpk2dGamma23223 + dGammadU2311*dpk2dGamma23231 + dGammadU2321*dpk2dGamma23232 + dGammadU2331*dpk2dGamma23233 + dGammadU3111*dpk2dGamma23311 + dGammadU3121*dpk2dGamma23312 + dGammadU3131*dpk2dGamma23313 + dGammadU3211*dpk2dGamma23321 + dGammadU3221*dpk2dGamma23322 + dGammadU3231*dpk2dGamma23323 + dGammadU3311*dpk2dGamma23331 + dGammadU3321*dpk2dGamma23332 + dGammadU3331*dpk2dGamma23333 + dPhidU111*dpk2dPsi2311 + dPhidU121*dpk2dPsi2312 + dPhidU131*dpk2dPsi2313 + dPhidU211*dpk2dPsi2321 + dPhidU221*dpk2dPsi2322 + dPhidU231*dpk2dPsi2323 + dPhidU311*dpk2dPsi2331 + dPhidU321*dpk2dPsi2332 + dPhidU331*dpk2dPsi2333
    dpk2dU[I5,I1] = dCdU111*dpk2dC1311 + dCdU121*dpk2dC1312 + dCdU131*dpk2dC1313 + dCdU211*dpk2dC1321 + dCdU221*dpk2dC1322 + dCdU231*dpk2dC1323 + dCdU311*dpk2dC1331 + dCdU321*dpk2dC1332 + dCdU331*dpk2dC1333 + dGammadU1111*dpk2dGamma13111 + dGammadU1121*dpk2dGamma13112 + dGammadU1131*dpk2dGamma13113 + dGammadU1211*dpk2dGamma13121 + dGammadU1221*dpk2dGamma13122 + dGammadU1231*dpk2dGamma13123 + dGammadU1311*dpk2dGamma13131 + dGammadU1321*dpk2dGamma13132 + dGammadU1331*dpk2dGamma13133 + dGammadU2111*dpk2dGamma13211 + dGammadU2121*dpk2dGamma13212 + dGammadU2131*dpk2dGamma13213 + dGammadU2211*dpk2dGamma13221 + dGammadU2221*dpk2dGamma13222 + dGammadU2231*dpk2dGamma13223 + dGammadU2311*dpk2dGamma13231 + dGammadU2321*dpk2dGamma13232 + dGammadU2331*dpk2dGamma13233 + dGammadU3111*dpk2dGamma13311 + dGammadU3121*dpk2dGamma13312 + dGammadU3131*dpk2dGamma13313 + dGammadU3211*dpk2dGamma13321 + dGammadU3221*dpk2dGamma13322 + dGammadU3231*dpk2dGamma13323 + dGammadU3311*dpk2dGamma13331 + dGammadU3321*dpk2dGamma13332 + dGammadU3331*dpk2dGamma13333 + dPhidU111*dpk2dPsi1311 + dPhidU121*dpk2dPsi1312 + dPhidU131*dpk2dPsi1313 + dPhidU211*dpk2dPsi1321 + dPhidU221*dpk2dPsi1322 + dPhidU231*dpk2dPsi1323 + dPhidU311*dpk2dPsi1331 + dPhidU321*dpk2dPsi1332 + dPhidU331*dpk2dPsi1333
    dpk2dU[I6,I1] = dCdU111*dpk2dC1211 + dCdU121*dpk2dC1212 + dCdU131*dpk2dC1213 + dCdU211*dpk2dC1221 + dCdU221*dpk2dC1222 + dCdU231*dpk2dC1223 + dCdU311*dpk2dC1231 + dCdU321*dpk2dC1232 + dCdU331*dpk2dC1233 + dGammadU1111*dpk2dGamma12111 + dGammadU1121*dpk2dGamma12112 + dGammadU1131*dpk2dGamma12113 + dGammadU1211*dpk2dGamma12121 + dGammadU1221*dpk2dGamma12122 + dGammadU1231*dpk2dGamma12123 + dGammadU1311*dpk2dGamma12131 + dGammadU1321*dpk2dGamma12132 + dGammadU1331*dpk2dGamma12133 + dGammadU2111*dpk2dGamma12211 + dGammadU2121*dpk2dGamma12212 + dGammadU2131*dpk2dGamma12213 + dGammadU2211*dpk2dGamma12221 + dGammadU2221*dpk2dGamma12222 + dGammadU2231*dpk2dGamma12223 + dGammadU2311*dpk2dGamma12231 + dGammadU2321*dpk2dGamma12232 + dGammadU2331*dpk2dGamma12233 + dGammadU3111*dpk2dGamma12311 + dGammadU3121*dpk2dGamma12312 + dGammadU3131*dpk2dGamma12313 + dGammadU3211*dpk2dGamma12321 + dGammadU3221*dpk2dGamma12322 + dGammadU3231*dpk2dGamma12323 + dGammadU3311*dpk2dGamma12331 + dGammadU3321*dpk2dGamma12332 + dGammadU3331*dpk2dGamma12333 + dPhidU111*dpk2dPsi1211 + dPhidU121*dpk2dPsi1212 + dPhidU131*dpk2dPsi1213 + dPhidU211*dpk2dPsi1221 + dPhidU221*dpk2dPsi1222 + dPhidU231*dpk2dPsi1223 + dPhidU311*dpk2dPsi1231 + dPhidU321*dpk2dPsi1232 + dPhidU331*dpk2dPsi1233
    dpk2dU[I7,I1] = dCdU111*dpk2dC3211 + dCdU121*dpk2dC3212 + dCdU131*dpk2dC3213 + dCdU211*dpk2dC3221 + dCdU221*dpk2dC3222 + dCdU231*dpk2dC3223 + dCdU311*dpk2dC3231 + dCdU321*dpk2dC3232 + dCdU331*dpk2dC3233 + dGammadU1111*dpk2dGamma32111 + dGammadU1121*dpk2dGamma32112 + dGammadU1131*dpk2dGamma32113 + dGammadU1211*dpk2dGamma32121 + dGammadU1221*dpk2dGamma32122 + dGammadU1231*dpk2dGamma32123 + dGammadU1311*dpk2dGamma32131 + dGammadU1321*dpk2dGamma32132 + dGammadU1331*dpk2dGamma32133 + dGammadU2111*dpk2dGamma32211 + dGammadU2121*dpk2dGamma32212 + dGammadU2131*dpk2dGamma32213 + dGammadU2211*dpk2dGamma32221 + dGammadU2221*dpk2dGamma32222 + dGammadU2231*dpk2dGamma32223 + dGammadU2311*dpk2dGamma32231 + dGammadU2321*dpk2dGamma32232 + dGammadU2331*dpk2dGamma32233 + dGammadU3111*dpk2dGamma32311 + dGammadU3121*dpk2dGamma32312 + dGammadU3131*dpk2dGamma32313 + dGammadU3211*dpk2dGamma32321 + dGammadU3221*dpk2dGamma32322 + dGammadU3231*dpk2dGamma32323 + dGammadU3311*dpk2dGamma32331 + dGammadU3321*dpk2dGamma32332 + dGammadU3331*dpk2dGamma32333 + dPhidU111*dpk2dPsi3211 + dPhidU121*dpk2dPsi3212 + dPhidU131*dpk2dPsi3213 + dPhidU211*dpk2dPsi3221 + dPhidU221*dpk2dPsi3222 + dPhidU231*dpk2dPsi3223 + dPhidU311*dpk2dPsi3231 + dPhidU321*dpk2dPsi3232 + dPhidU331*dpk2dPsi3233
    dpk2dU[I8,I1] = dCdU111*dpk2dC3111 + dCdU121*dpk2dC3112 + dCdU131*dpk2dC3113 + dCdU211*dpk2dC3121 + dCdU221*dpk2dC3122 + dCdU231*dpk2dC3123 + dCdU311*dpk2dC3131 + dCdU321*dpk2dC3132 + dCdU331*dpk2dC3133 + dGammadU1111*dpk2dGamma31111 + dGammadU1121*dpk2dGamma31112 + dGammadU1131*dpk2dGamma31113 + dGammadU1211*dpk2dGamma31121 + dGammadU1221*dpk2dGamma31122 + dGammadU1231*dpk2dGamma31123 + dGammadU1311*dpk2dGamma31131 + dGammadU1321*dpk2dGamma31132 + dGammadU1331*dpk2dGamma31133 + dGammadU2111*dpk2dGamma31211 + dGammadU2121*dpk2dGamma31212 + dGammadU2131*dpk2dGamma31213 + dGammadU2211*dpk2dGamma31221 + dGammadU2221*dpk2dGamma31222 + dGammadU2231*dpk2dGamma31223 + dGammadU2311*dpk2dGamma31231 + dGammadU2321*dpk2dGamma31232 + dGammadU2331*dpk2dGamma31233 + dGammadU3111*dpk2dGamma31311 + dGammadU3121*dpk2dGamma31312 + dGammadU3131*dpk2dGamma31313 + dGammadU3211*dpk2dGamma31321 + dGammadU3221*dpk2dGamma31322 + dGammadU3231*dpk2dGamma31323 + dGammadU3311*dpk2dGamma31331 + dGammadU3321*dpk2dGamma31332 + dGammadU3331*dpk2dGamma31333 + dPhidU111*dpk2dPsi3111 + dPhidU121*dpk2dPsi3112 + dPhidU131*dpk2dPsi3113 + dPhidU211*dpk2dPsi3121 + dPhidU221*dpk2dPsi3122 + dPhidU231*dpk2dPsi3123 + dPhidU311*dpk2dPsi3131 + dPhidU321*dpk2dPsi3132 + dPhidU331*dpk2dPsi3133
    dpk2dU[I9,I1] = dCdU111*dpk2dC2111 + dCdU121*dpk2dC2112 + dCdU131*dpk2dC2113 + dCdU211*dpk2dC2121 + dCdU221*dpk2dC2122 + dCdU231*dpk2dC2123 + dCdU311*dpk2dC2131 + dCdU321*dpk2dC2132 + dCdU331*dpk2dC2133 + dGammadU1111*dpk2dGamma21111 + dGammadU1121*dpk2dGamma21112 + dGammadU1131*dpk2dGamma21113 + dGammadU1211*dpk2dGamma21121 + dGammadU1221*dpk2dGamma21122 + dGammadU1231*dpk2dGamma21123 + dGammadU1311*dpk2dGamma21131 + dGammadU1321*dpk2dGamma21132 + dGammadU1331*dpk2dGamma21133 + dGammadU2111*dpk2dGamma21211 + dGammadU2121*dpk2dGamma21212 + dGammadU2131*dpk2dGamma21213 + dGammadU2211*dpk2dGamma21221 + dGammadU2221*dpk2dGamma21222 + dGammadU2231*dpk2dGamma21223 + dGammadU2311*dpk2dGamma21231 + dGammadU2321*dpk2dGamma21232 + dGammadU2331*dpk2dGamma21233 + dGammadU3111*dpk2dGamma21311 + dGammadU3121*dpk2dGamma21312 + dGammadU3131*dpk2dGamma21313 + dGammadU3211*dpk2dGamma21321 + dGammadU3221*dpk2dGamma21322 + dGammadU3231*dpk2dGamma21323 + dGammadU3311*dpk2dGamma21331 + dGammadU3321*dpk2dGamma21332 + dGammadU3331*dpk2dGamma21333 + dPhidU111*dpk2dPsi2111 + dPhidU121*dpk2dPsi2112 + dPhidU131*dpk2dPsi2113 + dPhidU211*dpk2dPsi2121 + dPhidU221*dpk2dPsi2122 + dPhidU231*dpk2dPsi2123 + dPhidU311*dpk2dPsi2131 + dPhidU321*dpk2dPsi2132 + dPhidU331*dpk2dPsi2133

    #Column 2
    dpk2dU[I1,I2] = dCdU112*dpk2dC1111 + dCdU122*dpk2dC1112 + dCdU132*dpk2dC1113 + dCdU212*dpk2dC1121 + dCdU222*dpk2dC1122 + dCdU232*dpk2dC1123 + dCdU312*dpk2dC1131 + dCdU322*dpk2dC1132 + dCdU332*dpk2dC1133 + dGammadU1112*dpk2dGamma11111 + dGammadU1122*dpk2dGamma11112 + dGammadU1132*dpk2dGamma11113 + dGammadU1212*dpk2dGamma11121 + dGammadU1222*dpk2dGamma11122 + dGammadU1232*dpk2dGamma11123 + dGammadU1312*dpk2dGamma11131 + dGammadU1322*dpk2dGamma11132 + dGammadU1332*dpk2dGamma11133 + dGammadU2112*dpk2dGamma11211 + dGammadU2122*dpk2dGamma11212 + dGammadU2132*dpk2dGamma11213 + dGammadU2212*dpk2dGamma11221 + dGammadU2222*dpk2dGamma11222 + dGammadU2232*dpk2dGamma11223 + dGammadU2312*dpk2dGamma11231 + dGammadU2322*dpk2dGamma11232 + dGammadU2332*dpk2dGamma11233 + dGammadU3112*dpk2dGamma11311 + dGammadU3122*dpk2dGamma11312 + dGammadU3132*dpk2dGamma11313 + dGammadU3212*dpk2dGamma11321 + dGammadU3222*dpk2dGamma11322 + dGammadU3232*dpk2dGamma11323 + dGammadU3312*dpk2dGamma11331 + dGammadU3322*dpk2dGamma11332 + dGammadU3332*dpk2dGamma11333 + dPhidU112*dpk2dPsi1111 + dPhidU122*dpk2dPsi1112 + dPhidU132*dpk2dPsi1113 + dPhidU212*dpk2dPsi1121 + dPhidU222*dpk2dPsi1122 + dPhidU232*dpk2dPsi1123 + dPhidU312*dpk2dPsi1131 + dPhidU322*dpk2dPsi1132 + dPhidU332*dpk2dPsi1133
    dpk2dU[I2,I2] = dCdU112*dpk2dC2211 + dCdU122*dpk2dC2212 + dCdU132*dpk2dC2213 + dCdU212*dpk2dC2221 + dCdU222*dpk2dC2222 + dCdU232*dpk2dC2223 + dCdU312*dpk2dC2231 + dCdU322*dpk2dC2232 + dCdU332*dpk2dC2233 + dGammadU1112*dpk2dGamma22111 + dGammadU1122*dpk2dGamma22112 + dGammadU1132*dpk2dGamma22113 + dGammadU1212*dpk2dGamma22121 + dGammadU1222*dpk2dGamma22122 + dGammadU1232*dpk2dGamma22123 + dGammadU1312*dpk2dGamma22131 + dGammadU1322*dpk2dGamma22132 + dGammadU1332*dpk2dGamma22133 + dGammadU2112*dpk2dGamma22211 + dGammadU2122*dpk2dGamma22212 + dGammadU2132*dpk2dGamma22213 + dGammadU2212*dpk2dGamma22221 + dGammadU2222*dpk2dGamma22222 + dGammadU2232*dpk2dGamma22223 + dGammadU2312*dpk2dGamma22231 + dGammadU2322*dpk2dGamma22232 + dGammadU2332*dpk2dGamma22233 + dGammadU3112*dpk2dGamma22311 + dGammadU3122*dpk2dGamma22312 + dGammadU3132*dpk2dGamma22313 + dGammadU3212*dpk2dGamma22321 + dGammadU3222*dpk2dGamma22322 + dGammadU3232*dpk2dGamma22323 + dGammadU3312*dpk2dGamma22331 + dGammadU3322*dpk2dGamma22332 + dGammadU3332*dpk2dGamma22333 + dPhidU112*dpk2dPsi2211 + dPhidU122*dpk2dPsi2212 + dPhidU132*dpk2dPsi2213 + dPhidU212*dpk2dPsi2221 + dPhidU222*dpk2dPsi2222 + dPhidU232*dpk2dPsi2223 + dPhidU312*dpk2dPsi2231 + dPhidU322*dpk2dPsi2232 + dPhidU332*dpk2dPsi2233
    dpk2dU[I3,I2] = dCdU112*dpk2dC3311 + dCdU122*dpk2dC3312 + dCdU132*dpk2dC3313 + dCdU212*dpk2dC3321 + dCdU222*dpk2dC3322 + dCdU232*dpk2dC3323 + dCdU312*dpk2dC3331 + dCdU322*dpk2dC3332 + dCdU332*dpk2dC3333 + dGammadU1112*dpk2dGamma33111 + dGammadU1122*dpk2dGamma33112 + dGammadU1132*dpk2dGamma33113 + dGammadU1212*dpk2dGamma33121 + dGammadU1222*dpk2dGamma33122 + dGammadU1232*dpk2dGamma33123 + dGammadU1312*dpk2dGamma33131 + dGammadU1322*dpk2dGamma33132 + dGammadU1332*dpk2dGamma33133 + dGammadU2112*dpk2dGamma33211 + dGammadU2122*dpk2dGamma33212 + dGammadU2132*dpk2dGamma33213 + dGammadU2212*dpk2dGamma33221 + dGammadU2222*dpk2dGamma33222 + dGammadU2232*dpk2dGamma33223 + dGammadU2312*dpk2dGamma33231 + dGammadU2322*dpk2dGamma33232 + dGammadU2332*dpk2dGamma33233 + dGammadU3112*dpk2dGamma33311 + dGammadU3122*dpk2dGamma33312 + dGammadU3132*dpk2dGamma33313 + dGammadU3212*dpk2dGamma33321 + dGammadU3222*dpk2dGamma33322 + dGammadU3232*dpk2dGamma33323 + dGammadU3312*dpk2dGamma33331 + dGammadU3322*dpk2dGamma33332 + dGammadU3332*dpk2dGamma33333 + dPhidU112*dpk2dPsi3311 + dPhidU122*dpk2dPsi3312 + dPhidU132*dpk2dPsi3313 + dPhidU212*dpk2dPsi3321 + dPhidU222*dpk2dPsi3322 + dPhidU232*dpk2dPsi3323 + dPhidU312*dpk2dPsi3331 + dPhidU322*dpk2dPsi3332 + dPhidU332*dpk2dPsi3333
    dpk2dU[I4,I2] = dCdU112*dpk2dC2311 + dCdU122*dpk2dC2312 + dCdU132*dpk2dC2313 + dCdU212*dpk2dC2321 + dCdU222*dpk2dC2322 + dCdU232*dpk2dC2323 + dCdU312*dpk2dC2331 + dCdU322*dpk2dC2332 + dCdU332*dpk2dC2333 + dGammadU1112*dpk2dGamma23111 + dGammadU1122*dpk2dGamma23112 + dGammadU1132*dpk2dGamma23113 + dGammadU1212*dpk2dGamma23121 + dGammadU1222*dpk2dGamma23122 + dGammadU1232*dpk2dGamma23123 + dGammadU1312*dpk2dGamma23131 + dGammadU1322*dpk2dGamma23132 + dGammadU1332*dpk2dGamma23133 + dGammadU2112*dpk2dGamma23211 + dGammadU2122*dpk2dGamma23212 + dGammadU2132*dpk2dGamma23213 + dGammadU2212*dpk2dGamma23221 + dGammadU2222*dpk2dGamma23222 + dGammadU2232*dpk2dGamma23223 + dGammadU2312*dpk2dGamma23231 + dGammadU2322*dpk2dGamma23232 + dGammadU2332*dpk2dGamma23233 + dGammadU3112*dpk2dGamma23311 + dGammadU3122*dpk2dGamma23312 + dGammadU3132*dpk2dGamma23313 + dGammadU3212*dpk2dGamma23321 + dGammadU3222*dpk2dGamma23322 + dGammadU3232*dpk2dGamma23323 + dGammadU3312*dpk2dGamma23331 + dGammadU3322*dpk2dGamma23332 + dGammadU3332*dpk2dGamma23333 + dPhidU112*dpk2dPsi2311 + dPhidU122*dpk2dPsi2312 + dPhidU132*dpk2dPsi2313 + dPhidU212*dpk2dPsi2321 + dPhidU222*dpk2dPsi2322 + dPhidU232*dpk2dPsi2323 + dPhidU312*dpk2dPsi2331 + dPhidU322*dpk2dPsi2332 + dPhidU332*dpk2dPsi2333
    dpk2dU[I5,I2] = dCdU112*dpk2dC1311 + dCdU122*dpk2dC1312 + dCdU132*dpk2dC1313 + dCdU212*dpk2dC1321 + dCdU222*dpk2dC1322 + dCdU232*dpk2dC1323 + dCdU312*dpk2dC1331 + dCdU322*dpk2dC1332 + dCdU332*dpk2dC1333 + dGammadU1112*dpk2dGamma13111 + dGammadU1122*dpk2dGamma13112 + dGammadU1132*dpk2dGamma13113 + dGammadU1212*dpk2dGamma13121 + dGammadU1222*dpk2dGamma13122 + dGammadU1232*dpk2dGamma13123 + dGammadU1312*dpk2dGamma13131 + dGammadU1322*dpk2dGamma13132 + dGammadU1332*dpk2dGamma13133 + dGammadU2112*dpk2dGamma13211 + dGammadU2122*dpk2dGamma13212 + dGammadU2132*dpk2dGamma13213 + dGammadU2212*dpk2dGamma13221 + dGammadU2222*dpk2dGamma13222 + dGammadU2232*dpk2dGamma13223 + dGammadU2312*dpk2dGamma13231 + dGammadU2322*dpk2dGamma13232 + dGammadU2332*dpk2dGamma13233 + dGammadU3112*dpk2dGamma13311 + dGammadU3122*dpk2dGamma13312 + dGammadU3132*dpk2dGamma13313 + dGammadU3212*dpk2dGamma13321 + dGammadU3222*dpk2dGamma13322 + dGammadU3232*dpk2dGamma13323 + dGammadU3312*dpk2dGamma13331 + dGammadU3322*dpk2dGamma13332 + dGammadU3332*dpk2dGamma13333 + dPhidU112*dpk2dPsi1311 + dPhidU122*dpk2dPsi1312 + dPhidU132*dpk2dPsi1313 + dPhidU212*dpk2dPsi1321 + dPhidU222*dpk2dPsi1322 + dPhidU232*dpk2dPsi1323 + dPhidU312*dpk2dPsi1331 + dPhidU322*dpk2dPsi1332 + dPhidU332*dpk2dPsi1333
    dpk2dU[I6,I2] = dCdU112*dpk2dC1211 + dCdU122*dpk2dC1212 + dCdU132*dpk2dC1213 + dCdU212*dpk2dC1221 + dCdU222*dpk2dC1222 + dCdU232*dpk2dC1223 + dCdU312*dpk2dC1231 + dCdU322*dpk2dC1232 + dCdU332*dpk2dC1233 + dGammadU1112*dpk2dGamma12111 + dGammadU1122*dpk2dGamma12112 + dGammadU1132*dpk2dGamma12113 + dGammadU1212*dpk2dGamma12121 + dGammadU1222*dpk2dGamma12122 + dGammadU1232*dpk2dGamma12123 + dGammadU1312*dpk2dGamma12131 + dGammadU1322*dpk2dGamma12132 + dGammadU1332*dpk2dGamma12133 + dGammadU2112*dpk2dGamma12211 + dGammadU2122*dpk2dGamma12212 + dGammadU2132*dpk2dGamma12213 + dGammadU2212*dpk2dGamma12221 + dGammadU2222*dpk2dGamma12222 + dGammadU2232*dpk2dGamma12223 + dGammadU2312*dpk2dGamma12231 + dGammadU2322*dpk2dGamma12232 + dGammadU2332*dpk2dGamma12233 + dGammadU3112*dpk2dGamma12311 + dGammadU3122*dpk2dGamma12312 + dGammadU3132*dpk2dGamma12313 + dGammadU3212*dpk2dGamma12321 + dGammadU3222*dpk2dGamma12322 + dGammadU3232*dpk2dGamma12323 + dGammadU3312*dpk2dGamma12331 + dGammadU3322*dpk2dGamma12332 + dGammadU3332*dpk2dGamma12333 + dPhidU112*dpk2dPsi1211 + dPhidU122*dpk2dPsi1212 + dPhidU132*dpk2dPsi1213 + dPhidU212*dpk2dPsi1221 + dPhidU222*dpk2dPsi1222 + dPhidU232*dpk2dPsi1223 + dPhidU312*dpk2dPsi1231 + dPhidU322*dpk2dPsi1232 + dPhidU332*dpk2dPsi1233
    dpk2dU[I7,I2] = dCdU112*dpk2dC3211 + dCdU122*dpk2dC3212 + dCdU132*dpk2dC3213 + dCdU212*dpk2dC3221 + dCdU222*dpk2dC3222 + dCdU232*dpk2dC3223 + dCdU312*dpk2dC3231 + dCdU322*dpk2dC3232 + dCdU332*dpk2dC3233 + dGammadU1112*dpk2dGamma32111 + dGammadU1122*dpk2dGamma32112 + dGammadU1132*dpk2dGamma32113 + dGammadU1212*dpk2dGamma32121 + dGammadU1222*dpk2dGamma32122 + dGammadU1232*dpk2dGamma32123 + dGammadU1312*dpk2dGamma32131 + dGammadU1322*dpk2dGamma32132 + dGammadU1332*dpk2dGamma32133 + dGammadU2112*dpk2dGamma32211 + dGammadU2122*dpk2dGamma32212 + dGammadU2132*dpk2dGamma32213 + dGammadU2212*dpk2dGamma32221 + dGammadU2222*dpk2dGamma32222 + dGammadU2232*dpk2dGamma32223 + dGammadU2312*dpk2dGamma32231 + dGammadU2322*dpk2dGamma32232 + dGammadU2332*dpk2dGamma32233 + dGammadU3112*dpk2dGamma32311 + dGammadU3122*dpk2dGamma32312 + dGammadU3132*dpk2dGamma32313 + dGammadU3212*dpk2dGamma32321 + dGammadU3222*dpk2dGamma32322 + dGammadU3232*dpk2dGamma32323 + dGammadU3312*dpk2dGamma32331 + dGammadU3322*dpk2dGamma32332 + dGammadU3332*dpk2dGamma32333 + dPhidU112*dpk2dPsi3211 + dPhidU122*dpk2dPsi3212 + dPhidU132*dpk2dPsi3213 + dPhidU212*dpk2dPsi3221 + dPhidU222*dpk2dPsi3222 + dPhidU232*dpk2dPsi3223 + dPhidU312*dpk2dPsi3231 + dPhidU322*dpk2dPsi3232 + dPhidU332*dpk2dPsi3233
    dpk2dU[I8,I2] = dCdU112*dpk2dC3111 + dCdU122*dpk2dC3112 + dCdU132*dpk2dC3113 + dCdU212*dpk2dC3121 + dCdU222*dpk2dC3122 + dCdU232*dpk2dC3123 + dCdU312*dpk2dC3131 + dCdU322*dpk2dC3132 + dCdU332*dpk2dC3133 + dGammadU1112*dpk2dGamma31111 + dGammadU1122*dpk2dGamma31112 + dGammadU1132*dpk2dGamma31113 + dGammadU1212*dpk2dGamma31121 + dGammadU1222*dpk2dGamma31122 + dGammadU1232*dpk2dGamma31123 + dGammadU1312*dpk2dGamma31131 + dGammadU1322*dpk2dGamma31132 + dGammadU1332*dpk2dGamma31133 + dGammadU2112*dpk2dGamma31211 + dGammadU2122*dpk2dGamma31212 + dGammadU2132*dpk2dGamma31213 + dGammadU2212*dpk2dGamma31221 + dGammadU2222*dpk2dGamma31222 + dGammadU2232*dpk2dGamma31223 + dGammadU2312*dpk2dGamma31231 + dGammadU2322*dpk2dGamma31232 + dGammadU2332*dpk2dGamma31233 + dGammadU3112*dpk2dGamma31311 + dGammadU3122*dpk2dGamma31312 + dGammadU3132*dpk2dGamma31313 + dGammadU3212*dpk2dGamma31321 + dGammadU3222*dpk2dGamma31322 + dGammadU3232*dpk2dGamma31323 + dGammadU3312*dpk2dGamma31331 + dGammadU3322*dpk2dGamma31332 + dGammadU3332*dpk2dGamma31333 + dPhidU112*dpk2dPsi3111 + dPhidU122*dpk2dPsi3112 + dPhidU132*dpk2dPsi3113 + dPhidU212*dpk2dPsi3121 + dPhidU222*dpk2dPsi3122 + dPhidU232*dpk2dPsi3123 + dPhidU312*dpk2dPsi3131 + dPhidU322*dpk2dPsi3132 + dPhidU332*dpk2dPsi3133
    dpk2dU[I9,I2] = dCdU112*dpk2dC2111 + dCdU122*dpk2dC2112 + dCdU132*dpk2dC2113 + dCdU212*dpk2dC2121 + dCdU222*dpk2dC2122 + dCdU232*dpk2dC2123 + dCdU312*dpk2dC2131 + dCdU322*dpk2dC2132 + dCdU332*dpk2dC2133 + dGammadU1112*dpk2dGamma21111 + dGammadU1122*dpk2dGamma21112 + dGammadU1132*dpk2dGamma21113 + dGammadU1212*dpk2dGamma21121 + dGammadU1222*dpk2dGamma21122 + dGammadU1232*dpk2dGamma21123 + dGammadU1312*dpk2dGamma21131 + dGammadU1322*dpk2dGamma21132 + dGammadU1332*dpk2dGamma21133 + dGammadU2112*dpk2dGamma21211 + dGammadU2122*dpk2dGamma21212 + dGammadU2132*dpk2dGamma21213 + dGammadU2212*dpk2dGamma21221 + dGammadU2222*dpk2dGamma21222 + dGammadU2232*dpk2dGamma21223 + dGammadU2312*dpk2dGamma21231 + dGammadU2322*dpk2dGamma21232 + dGammadU2332*dpk2dGamma21233 + dGammadU3112*dpk2dGamma21311 + dGammadU3122*dpk2dGamma21312 + dGammadU3132*dpk2dGamma21313 + dGammadU3212*dpk2dGamma21321 + dGammadU3222*dpk2dGamma21322 + dGammadU3232*dpk2dGamma21323 + dGammadU3312*dpk2dGamma21331 + dGammadU3322*dpk2dGamma21332 + dGammadU3332*dpk2dGamma21333 + dPhidU112*dpk2dPsi2111 + dPhidU122*dpk2dPsi2112 + dPhidU132*dpk2dPsi2113 + dPhidU212*dpk2dPsi2121 + dPhidU222*dpk2dPsi2122 + dPhidU232*dpk2dPsi2123 + dPhidU312*dpk2dPsi2131 + dPhidU322*dpk2dPsi2132 + dPhidU332*dpk2dPsi2133

    #Column 3
    dpk2dU[I1,I3] = dCdU113*dpk2dC1111 + dCdU123*dpk2dC1112 + dCdU133*dpk2dC1113 + dCdU213*dpk2dC1121 + dCdU223*dpk2dC1122 + dCdU233*dpk2dC1123 + dCdU313*dpk2dC1131 + dCdU323*dpk2dC1132 + dCdU333*dpk2dC1133 + dGammadU1113*dpk2dGamma11111 + dGammadU1123*dpk2dGamma11112 + dGammadU1133*dpk2dGamma11113 + dGammadU1213*dpk2dGamma11121 + dGammadU1223*dpk2dGamma11122 + dGammadU1233*dpk2dGamma11123 + dGammadU1313*dpk2dGamma11131 + dGammadU1323*dpk2dGamma11132 + dGammadU1333*dpk2dGamma11133 + dGammadU2113*dpk2dGamma11211 + dGammadU2123*dpk2dGamma11212 + dGammadU2133*dpk2dGamma11213 + dGammadU2213*dpk2dGamma11221 + dGammadU2223*dpk2dGamma11222 + dGammadU2233*dpk2dGamma11223 + dGammadU2313*dpk2dGamma11231 + dGammadU2323*dpk2dGamma11232 + dGammadU2333*dpk2dGamma11233 + dGammadU3113*dpk2dGamma11311 + dGammadU3123*dpk2dGamma11312 + dGammadU3133*dpk2dGamma11313 + dGammadU3213*dpk2dGamma11321 + dGammadU3223*dpk2dGamma11322 + dGammadU3233*dpk2dGamma11323 + dGammadU3313*dpk2dGamma11331 + dGammadU3323*dpk2dGamma11332 + dGammadU3333*dpk2dGamma11333 + dPhidU113*dpk2dPsi1111 + dPhidU123*dpk2dPsi1112 + dPhidU133*dpk2dPsi1113 + dPhidU213*dpk2dPsi1121 + dPhidU223*dpk2dPsi1122 + dPhidU233*dpk2dPsi1123 + dPhidU313*dpk2dPsi1131 + dPhidU323*dpk2dPsi1132 + dPhidU333*dpk2dPsi1133
    dpk2dU[I2,I3] = dCdU113*dpk2dC2211 + dCdU123*dpk2dC2212 + dCdU133*dpk2dC2213 + dCdU213*dpk2dC2221 + dCdU223*dpk2dC2222 + dCdU233*dpk2dC2223 + dCdU313*dpk2dC2231 + dCdU323*dpk2dC2232 + dCdU333*dpk2dC2233 + dGammadU1113*dpk2dGamma22111 + dGammadU1123*dpk2dGamma22112 + dGammadU1133*dpk2dGamma22113 + dGammadU1213*dpk2dGamma22121 + dGammadU1223*dpk2dGamma22122 + dGammadU1233*dpk2dGamma22123 + dGammadU1313*dpk2dGamma22131 + dGammadU1323*dpk2dGamma22132 + dGammadU1333*dpk2dGamma22133 + dGammadU2113*dpk2dGamma22211 + dGammadU2123*dpk2dGamma22212 + dGammadU2133*dpk2dGamma22213 + dGammadU2213*dpk2dGamma22221 + dGammadU2223*dpk2dGamma22222 + dGammadU2233*dpk2dGamma22223 + dGammadU2313*dpk2dGamma22231 + dGammadU2323*dpk2dGamma22232 + dGammadU2333*dpk2dGamma22233 + dGammadU3113*dpk2dGamma22311 + dGammadU3123*dpk2dGamma22312 + dGammadU3133*dpk2dGamma22313 + dGammadU3213*dpk2dGamma22321 + dGammadU3223*dpk2dGamma22322 + dGammadU3233*dpk2dGamma22323 + dGammadU3313*dpk2dGamma22331 + dGammadU3323*dpk2dGamma22332 + dGammadU3333*dpk2dGamma22333 + dPhidU113*dpk2dPsi2211 + dPhidU123*dpk2dPsi2212 + dPhidU133*dpk2dPsi2213 + dPhidU213*dpk2dPsi2221 + dPhidU223*dpk2dPsi2222 + dPhidU233*dpk2dPsi2223 + dPhidU313*dpk2dPsi2231 + dPhidU323*dpk2dPsi2232 + dPhidU333*dpk2dPsi2233
    dpk2dU[I3,I3] = dCdU113*dpk2dC3311 + dCdU123*dpk2dC3312 + dCdU133*dpk2dC3313 + dCdU213*dpk2dC3321 + dCdU223*dpk2dC3322 + dCdU233*dpk2dC3323 + dCdU313*dpk2dC3331 + dCdU323*dpk2dC3332 + dCdU333*dpk2dC3333 + dGammadU1113*dpk2dGamma33111 + dGammadU1123*dpk2dGamma33112 + dGammadU1133*dpk2dGamma33113 + dGammadU1213*dpk2dGamma33121 + dGammadU1223*dpk2dGamma33122 + dGammadU1233*dpk2dGamma33123 + dGammadU1313*dpk2dGamma33131 + dGammadU1323*dpk2dGamma33132 + dGammadU1333*dpk2dGamma33133 + dGammadU2113*dpk2dGamma33211 + dGammadU2123*dpk2dGamma33212 + dGammadU2133*dpk2dGamma33213 + dGammadU2213*dpk2dGamma33221 + dGammadU2223*dpk2dGamma33222 + dGammadU2233*dpk2dGamma33223 + dGammadU2313*dpk2dGamma33231 + dGammadU2323*dpk2dGamma33232 + dGammadU2333*dpk2dGamma33233 + dGammadU3113*dpk2dGamma33311 + dGammadU3123*dpk2dGamma33312 + dGammadU3133*dpk2dGamma33313 + dGammadU3213*dpk2dGamma33321 + dGammadU3223*dpk2dGamma33322 + dGammadU3233*dpk2dGamma33323 + dGammadU3313*dpk2dGamma33331 + dGammadU3323*dpk2dGamma33332 + dGammadU3333*dpk2dGamma33333 + dPhidU113*dpk2dPsi3311 + dPhidU123*dpk2dPsi3312 + dPhidU133*dpk2dPsi3313 + dPhidU213*dpk2dPsi3321 + dPhidU223*dpk2dPsi3322 + dPhidU233*dpk2dPsi3323 + dPhidU313*dpk2dPsi3331 + dPhidU323*dpk2dPsi3332 + dPhidU333*dpk2dPsi3333
    dpk2dU[I4,I3] = dCdU113*dpk2dC2311 + dCdU123*dpk2dC2312 + dCdU133*dpk2dC2313 + dCdU213*dpk2dC2321 + dCdU223*dpk2dC2322 + dCdU233*dpk2dC2323 + dCdU313*dpk2dC2331 + dCdU323*dpk2dC2332 + dCdU333*dpk2dC2333 + dGammadU1113*dpk2dGamma23111 + dGammadU1123*dpk2dGamma23112 + dGammadU1133*dpk2dGamma23113 + dGammadU1213*dpk2dGamma23121 + dGammadU1223*dpk2dGamma23122 + dGammadU1233*dpk2dGamma23123 + dGammadU1313*dpk2dGamma23131 + dGammadU1323*dpk2dGamma23132 + dGammadU1333*dpk2dGamma23133 + dGammadU2113*dpk2dGamma23211 + dGammadU2123*dpk2dGamma23212 + dGammadU2133*dpk2dGamma23213 + dGammadU2213*dpk2dGamma23221 + dGammadU2223*dpk2dGamma23222 + dGammadU2233*dpk2dGamma23223 + dGammadU2313*dpk2dGamma23231 + dGammadU2323*dpk2dGamma23232 + dGammadU2333*dpk2dGamma23233 + dGammadU3113*dpk2dGamma23311 + dGammadU3123*dpk2dGamma23312 + dGammadU3133*dpk2dGamma23313 + dGammadU3213*dpk2dGamma23321 + dGammadU3223*dpk2dGamma23322 + dGammadU3233*dpk2dGamma23323 + dGammadU3313*dpk2dGamma23331 + dGammadU3323*dpk2dGamma23332 + dGammadU3333*dpk2dGamma23333 + dPhidU113*dpk2dPsi2311 + dPhidU123*dpk2dPsi2312 + dPhidU133*dpk2dPsi2313 + dPhidU213*dpk2dPsi2321 + dPhidU223*dpk2dPsi2322 + dPhidU233*dpk2dPsi2323 + dPhidU313*dpk2dPsi2331 + dPhidU323*dpk2dPsi2332 + dPhidU333*dpk2dPsi2333
    dpk2dU[I5,I3] = dCdU113*dpk2dC1311 + dCdU123*dpk2dC1312 + dCdU133*dpk2dC1313 + dCdU213*dpk2dC1321 + dCdU223*dpk2dC1322 + dCdU233*dpk2dC1323 + dCdU313*dpk2dC1331 + dCdU323*dpk2dC1332 + dCdU333*dpk2dC1333 + dGammadU1113*dpk2dGamma13111 + dGammadU1123*dpk2dGamma13112 + dGammadU1133*dpk2dGamma13113 + dGammadU1213*dpk2dGamma13121 + dGammadU1223*dpk2dGamma13122 + dGammadU1233*dpk2dGamma13123 + dGammadU1313*dpk2dGamma13131 + dGammadU1323*dpk2dGamma13132 + dGammadU1333*dpk2dGamma13133 + dGammadU2113*dpk2dGamma13211 + dGammadU2123*dpk2dGamma13212 + dGammadU2133*dpk2dGamma13213 + dGammadU2213*dpk2dGamma13221 + dGammadU2223*dpk2dGamma13222 + dGammadU2233*dpk2dGamma13223 + dGammadU2313*dpk2dGamma13231 + dGammadU2323*dpk2dGamma13232 + dGammadU2333*dpk2dGamma13233 + dGammadU3113*dpk2dGamma13311 + dGammadU3123*dpk2dGamma13312 + dGammadU3133*dpk2dGamma13313 + dGammadU3213*dpk2dGamma13321 + dGammadU3223*dpk2dGamma13322 + dGammadU3233*dpk2dGamma13323 + dGammadU3313*dpk2dGamma13331 + dGammadU3323*dpk2dGamma13332 + dGammadU3333*dpk2dGamma13333 + dPhidU113*dpk2dPsi1311 + dPhidU123*dpk2dPsi1312 + dPhidU133*dpk2dPsi1313 + dPhidU213*dpk2dPsi1321 + dPhidU223*dpk2dPsi1322 + dPhidU233*dpk2dPsi1323 + dPhidU313*dpk2dPsi1331 + dPhidU323*dpk2dPsi1332 + dPhidU333*dpk2dPsi1333
    dpk2dU[I6,I3] = dCdU113*dpk2dC1211 + dCdU123*dpk2dC1212 + dCdU133*dpk2dC1213 + dCdU213*dpk2dC1221 + dCdU223*dpk2dC1222 + dCdU233*dpk2dC1223 + dCdU313*dpk2dC1231 + dCdU323*dpk2dC1232 + dCdU333*dpk2dC1233 + dGammadU1113*dpk2dGamma12111 + dGammadU1123*dpk2dGamma12112 + dGammadU1133*dpk2dGamma12113 + dGammadU1213*dpk2dGamma12121 + dGammadU1223*dpk2dGamma12122 + dGammadU1233*dpk2dGamma12123 + dGammadU1313*dpk2dGamma12131 + dGammadU1323*dpk2dGamma12132 + dGammadU1333*dpk2dGamma12133 + dGammadU2113*dpk2dGamma12211 + dGammadU2123*dpk2dGamma12212 + dGammadU2133*dpk2dGamma12213 + dGammadU2213*dpk2dGamma12221 + dGammadU2223*dpk2dGamma12222 + dGammadU2233*dpk2dGamma12223 + dGammadU2313*dpk2dGamma12231 + dGammadU2323*dpk2dGamma12232 + dGammadU2333*dpk2dGamma12233 + dGammadU3113*dpk2dGamma12311 + dGammadU3123*dpk2dGamma12312 + dGammadU3133*dpk2dGamma12313 + dGammadU3213*dpk2dGamma12321 + dGammadU3223*dpk2dGamma12322 + dGammadU3233*dpk2dGamma12323 + dGammadU3313*dpk2dGamma12331 + dGammadU3323*dpk2dGamma12332 + dGammadU3333*dpk2dGamma12333 + dPhidU113*dpk2dPsi1211 + dPhidU123*dpk2dPsi1212 + dPhidU133*dpk2dPsi1213 + dPhidU213*dpk2dPsi1221 + dPhidU223*dpk2dPsi1222 + dPhidU233*dpk2dPsi1223 + dPhidU313*dpk2dPsi1231 + dPhidU323*dpk2dPsi1232 + dPhidU333*dpk2dPsi1233
    dpk2dU[I7,I3] = dCdU113*dpk2dC3211 + dCdU123*dpk2dC3212 + dCdU133*dpk2dC3213 + dCdU213*dpk2dC3221 + dCdU223*dpk2dC3222 + dCdU233*dpk2dC3223 + dCdU313*dpk2dC3231 + dCdU323*dpk2dC3232 + dCdU333*dpk2dC3233 + dGammadU1113*dpk2dGamma32111 + dGammadU1123*dpk2dGamma32112 + dGammadU1133*dpk2dGamma32113 + dGammadU1213*dpk2dGamma32121 + dGammadU1223*dpk2dGamma32122 + dGammadU1233*dpk2dGamma32123 + dGammadU1313*dpk2dGamma32131 + dGammadU1323*dpk2dGamma32132 + dGammadU1333*dpk2dGamma32133 + dGammadU2113*dpk2dGamma32211 + dGammadU2123*dpk2dGamma32212 + dGammadU2133*dpk2dGamma32213 + dGammadU2213*dpk2dGamma32221 + dGammadU2223*dpk2dGamma32222 + dGammadU2233*dpk2dGamma32223 + dGammadU2313*dpk2dGamma32231 + dGammadU2323*dpk2dGamma32232 + dGammadU2333*dpk2dGamma32233 + dGammadU3113*dpk2dGamma32311 + dGammadU3123*dpk2dGamma32312 + dGammadU3133*dpk2dGamma32313 + dGammadU3213*dpk2dGamma32321 + dGammadU3223*dpk2dGamma32322 + dGammadU3233*dpk2dGamma32323 + dGammadU3313*dpk2dGamma32331 + dGammadU3323*dpk2dGamma32332 + dGammadU3333*dpk2dGamma32333 + dPhidU113*dpk2dPsi3211 + dPhidU123*dpk2dPsi3212 + dPhidU133*dpk2dPsi3213 + dPhidU213*dpk2dPsi3221 + dPhidU223*dpk2dPsi3222 + dPhidU233*dpk2dPsi3223 + dPhidU313*dpk2dPsi3231 + dPhidU323*dpk2dPsi3232 + dPhidU333*dpk2dPsi3233
    dpk2dU[I8,I3] = dCdU113*dpk2dC3111 + dCdU123*dpk2dC3112 + dCdU133*dpk2dC3113 + dCdU213*dpk2dC3121 + dCdU223*dpk2dC3122 + dCdU233*dpk2dC3123 + dCdU313*dpk2dC3131 + dCdU323*dpk2dC3132 + dCdU333*dpk2dC3133 + dGammadU1113*dpk2dGamma31111 + dGammadU1123*dpk2dGamma31112 + dGammadU1133*dpk2dGamma31113 + dGammadU1213*dpk2dGamma31121 + dGammadU1223*dpk2dGamma31122 + dGammadU1233*dpk2dGamma31123 + dGammadU1313*dpk2dGamma31131 + dGammadU1323*dpk2dGamma31132 + dGammadU1333*dpk2dGamma31133 + dGammadU2113*dpk2dGamma31211 + dGammadU2123*dpk2dGamma31212 + dGammadU2133*dpk2dGamma31213 + dGammadU2213*dpk2dGamma31221 + dGammadU2223*dpk2dGamma31222 + dGammadU2233*dpk2dGamma31223 + dGammadU2313*dpk2dGamma31231 + dGammadU2323*dpk2dGamma31232 + dGammadU2333*dpk2dGamma31233 + dGammadU3113*dpk2dGamma31311 + dGammadU3123*dpk2dGamma31312 + dGammadU3133*dpk2dGamma31313 + dGammadU3213*dpk2dGamma31321 + dGammadU3223*dpk2dGamma31322 + dGammadU3233*dpk2dGamma31323 + dGammadU3313*dpk2dGamma31331 + dGammadU3323*dpk2dGamma31332 + dGammadU3333*dpk2dGamma31333 + dPhidU113*dpk2dPsi3111 + dPhidU123*dpk2dPsi3112 + dPhidU133*dpk2dPsi3113 + dPhidU213*dpk2dPsi3121 + dPhidU223*dpk2dPsi3122 + dPhidU233*dpk2dPsi3123 + dPhidU313*dpk2dPsi3131 + dPhidU323*dpk2dPsi3132 + dPhidU333*dpk2dPsi3133
    dpk2dU[I9,I3] = dCdU113*dpk2dC2111 + dCdU123*dpk2dC2112 + dCdU133*dpk2dC2113 + dCdU213*dpk2dC2121 + dCdU223*dpk2dC2122 + dCdU233*dpk2dC2123 + dCdU313*dpk2dC2131 + dCdU323*dpk2dC2132 + dCdU333*dpk2dC2133 + dGammadU1113*dpk2dGamma21111 + dGammadU1123*dpk2dGamma21112 + dGammadU1133*dpk2dGamma21113 + dGammadU1213*dpk2dGamma21121 + dGammadU1223*dpk2dGamma21122 + dGammadU1233*dpk2dGamma21123 + dGammadU1313*dpk2dGamma21131 + dGammadU1323*dpk2dGamma21132 + dGammadU1333*dpk2dGamma21133 + dGammadU2113*dpk2dGamma21211 + dGammadU2123*dpk2dGamma21212 + dGammadU2133*dpk2dGamma21213 + dGammadU2213*dpk2dGamma21221 + dGammadU2223*dpk2dGamma21222 + dGammadU2233*dpk2dGamma21223 + dGammadU2313*dpk2dGamma21231 + dGammadU2323*dpk2dGamma21232 + dGammadU2333*dpk2dGamma21233 + dGammadU3113*dpk2dGamma21311 + dGammadU3123*dpk2dGamma21312 + dGammadU3133*dpk2dGamma21313 + dGammadU3213*dpk2dGamma21321 + dGammadU3223*dpk2dGamma21322 + dGammadU3233*dpk2dGamma21323 + dGammadU3313*dpk2dGamma21331 + dGammadU3323*dpk2dGamma21332 + dGammadU3333*dpk2dGamma21333 + dPhidU113*dpk2dPsi2111 + dPhidU123*dpk2dPsi2112 + dPhidU133*dpk2dPsi2113 + dPhidU213*dpk2dPsi2121 + dPhidU223*dpk2dPsi2122 + dPhidU233*dpk2dPsi2123 + dPhidU313*dpk2dPsi2131 + dPhidU323*dpk2dPsi2132 + dPhidU333*dpk2dPsi2133

    #Column 4
    dpk2dU[I1,I4] = dGammadU1114*dpk2dGamma11111 + dGammadU1124*dpk2dGamma11112 + dGammadU1134*dpk2dGamma11113 + dGammadU2114*dpk2dGamma11211 + dGammadU2124*dpk2dGamma11212 + dGammadU2134*dpk2dGamma11213 + dGammadU3114*dpk2dGamma11311 + dGammadU3124*dpk2dGamma11312 + dGammadU3134*dpk2dGamma11313 + dPhidU114*dpk2dPsi1111 + dPhidU214*dpk2dPsi1121 + dPhidU314*dpk2dPsi1131
    dpk2dU[I2,I4] = dGammadU1114*dpk2dGamma22111 + dGammadU1124*dpk2dGamma22112 + dGammadU1134*dpk2dGamma22113 + dGammadU2114*dpk2dGamma22211 + dGammadU2124*dpk2dGamma22212 + dGammadU2134*dpk2dGamma22213 + dGammadU3114*dpk2dGamma22311 + dGammadU3124*dpk2dGamma22312 + dGammadU3134*dpk2dGamma22313 + dPhidU114*dpk2dPsi2211 + dPhidU214*dpk2dPsi2221 + dPhidU314*dpk2dPsi2231
    dpk2dU[I3,I4] = dGammadU1114*dpk2dGamma33111 + dGammadU1124*dpk2dGamma33112 + dGammadU1134*dpk2dGamma33113 + dGammadU2114*dpk2dGamma33211 + dGammadU2124*dpk2dGamma33212 + dGammadU2134*dpk2dGamma33213 + dGammadU3114*dpk2dGamma33311 + dGammadU3124*dpk2dGamma33312 + dGammadU3134*dpk2dGamma33313 + dPhidU114*dpk2dPsi3311 + dPhidU214*dpk2dPsi3321 + dPhidU314*dpk2dPsi3331
    dpk2dU[I4,I4] = dGammadU1114*dpk2dGamma23111 + dGammadU1124*dpk2dGamma23112 + dGammadU1134*dpk2dGamma23113 + dGammadU2114*dpk2dGamma23211 + dGammadU2124*dpk2dGamma23212 + dGammadU2134*dpk2dGamma23213 + dGammadU3114*dpk2dGamma23311 + dGammadU3124*dpk2dGamma23312 + dGammadU3134*dpk2dGamma23313 + dPhidU114*dpk2dPsi2311 + dPhidU214*dpk2dPsi2321 + dPhidU314*dpk2dPsi2331
    dpk2dU[I5,I4] = dGammadU1114*dpk2dGamma13111 + dGammadU1124*dpk2dGamma13112 + dGammadU1134*dpk2dGamma13113 + dGammadU2114*dpk2dGamma13211 + dGammadU2124*dpk2dGamma13212 + dGammadU2134*dpk2dGamma13213 + dGammadU3114*dpk2dGamma13311 + dGammadU3124*dpk2dGamma13312 + dGammadU3134*dpk2dGamma13313 + dPhidU114*dpk2dPsi1311 + dPhidU214*dpk2dPsi1321 + dPhidU314*dpk2dPsi1331
    dpk2dU[I6,I4] = dGammadU1114*dpk2dGamma12111 + dGammadU1124*dpk2dGamma12112 + dGammadU1134*dpk2dGamma12113 + dGammadU2114*dpk2dGamma12211 + dGammadU2124*dpk2dGamma12212 + dGammadU2134*dpk2dGamma12213 + dGammadU3114*dpk2dGamma12311 + dGammadU3124*dpk2dGamma12312 + dGammadU3134*dpk2dGamma12313 + dPhidU114*dpk2dPsi1211 + dPhidU214*dpk2dPsi1221 + dPhidU314*dpk2dPsi1231
    dpk2dU[I7,I4] = dGammadU1114*dpk2dGamma32111 + dGammadU1124*dpk2dGamma32112 + dGammadU1134*dpk2dGamma32113 + dGammadU2114*dpk2dGamma32211 + dGammadU2124*dpk2dGamma32212 + dGammadU2134*dpk2dGamma32213 + dGammadU3114*dpk2dGamma32311 + dGammadU3124*dpk2dGamma32312 + dGammadU3134*dpk2dGamma32313 + dPhidU114*dpk2dPsi3211 + dPhidU214*dpk2dPsi3221 + dPhidU314*dpk2dPsi3231
    dpk2dU[I8,I4] = dGammadU1114*dpk2dGamma31111 + dGammadU1124*dpk2dGamma31112 + dGammadU1134*dpk2dGamma31113 + dGammadU2114*dpk2dGamma31211 + dGammadU2124*dpk2dGamma31212 + dGammadU2134*dpk2dGamma31213 + dGammadU3114*dpk2dGamma31311 + dGammadU3124*dpk2dGamma31312 + dGammadU3134*dpk2dGamma31313 + dPhidU114*dpk2dPsi3111 + dPhidU214*dpk2dPsi3121 + dPhidU314*dpk2dPsi3131
    dpk2dU[I9,I4] = dGammadU1114*dpk2dGamma21111 + dGammadU1124*dpk2dGamma21112 + dGammadU1134*dpk2dGamma21113 + dGammadU2114*dpk2dGamma21211 + dGammadU2124*dpk2dGamma21212 + dGammadU2134*dpk2dGamma21213 + dGammadU3114*dpk2dGamma21311 + dGammadU3124*dpk2dGamma21312 + dGammadU3134*dpk2dGamma21313 + dPhidU114*dpk2dPsi2111 + dPhidU214*dpk2dPsi2121 + dPhidU314*dpk2dPsi2131

    #Column 5
    dpk2dU[I1,I5] = dGammadU1215*dpk2dGamma11121 + dGammadU1225*dpk2dGamma11122 + dGammadU1235*dpk2dGamma11123 + dGammadU2215*dpk2dGamma11221 + dGammadU2225*dpk2dGamma11222 + dGammadU2235*dpk2dGamma11223 + dGammadU3215*dpk2dGamma11321 + dGammadU3225*dpk2dGamma11322 + dGammadU3235*dpk2dGamma11323 + dPhidU125*dpk2dPsi1112 + dPhidU225*dpk2dPsi1122 + dPhidU325*dpk2dPsi1132
    dpk2dU[I2,I5] = dGammadU1215*dpk2dGamma22121 + dGammadU1225*dpk2dGamma22122 + dGammadU1235*dpk2dGamma22123 + dGammadU2215*dpk2dGamma22221 + dGammadU2225*dpk2dGamma22222 + dGammadU2235*dpk2dGamma22223 + dGammadU3215*dpk2dGamma22321 + dGammadU3225*dpk2dGamma22322 + dGammadU3235*dpk2dGamma22323 + dPhidU125*dpk2dPsi2212 + dPhidU225*dpk2dPsi2222 + dPhidU325*dpk2dPsi2232
    dpk2dU[I3,I5] = dGammadU1215*dpk2dGamma33121 + dGammadU1225*dpk2dGamma33122 + dGammadU1235*dpk2dGamma33123 + dGammadU2215*dpk2dGamma33221 + dGammadU2225*dpk2dGamma33222 + dGammadU2235*dpk2dGamma33223 + dGammadU3215*dpk2dGamma33321 + dGammadU3225*dpk2dGamma33322 + dGammadU3235*dpk2dGamma33323 + dPhidU125*dpk2dPsi3312 + dPhidU225*dpk2dPsi3322 + dPhidU325*dpk2dPsi3332
    dpk2dU[I4,I5] = dGammadU1215*dpk2dGamma23121 + dGammadU1225*dpk2dGamma23122 + dGammadU1235*dpk2dGamma23123 + dGammadU2215*dpk2dGamma23221 + dGammadU2225*dpk2dGamma23222 + dGammadU2235*dpk2dGamma23223 + dGammadU3215*dpk2dGamma23321 + dGammadU3225*dpk2dGamma23322 + dGammadU3235*dpk2dGamma23323 + dPhidU125*dpk2dPsi2312 + dPhidU225*dpk2dPsi2322 + dPhidU325*dpk2dPsi2332
    dpk2dU[I5,I5] = dGammadU1215*dpk2dGamma13121 + dGammadU1225*dpk2dGamma13122 + dGammadU1235*dpk2dGamma13123 + dGammadU2215*dpk2dGamma13221 + dGammadU2225*dpk2dGamma13222 + dGammadU2235*dpk2dGamma13223 + dGammadU3215*dpk2dGamma13321 + dGammadU3225*dpk2dGamma13322 + dGammadU3235*dpk2dGamma13323 + dPhidU125*dpk2dPsi1312 + dPhidU225*dpk2dPsi1322 + dPhidU325*dpk2dPsi1332
    dpk2dU[I6,I5] = dGammadU1215*dpk2dGamma12121 + dGammadU1225*dpk2dGamma12122 + dGammadU1235*dpk2dGamma12123 + dGammadU2215*dpk2dGamma12221 + dGammadU2225*dpk2dGamma12222 + dGammadU2235*dpk2dGamma12223 + dGammadU3215*dpk2dGamma12321 + dGammadU3225*dpk2dGamma12322 + dGammadU3235*dpk2dGamma12323 + dPhidU125*dpk2dPsi1212 + dPhidU225*dpk2dPsi1222 + dPhidU325*dpk2dPsi1232
    dpk2dU[I7,I5] = dGammadU1215*dpk2dGamma32121 + dGammadU1225*dpk2dGamma32122 + dGammadU1235*dpk2dGamma32123 + dGammadU2215*dpk2dGamma32221 + dGammadU2225*dpk2dGamma32222 + dGammadU2235*dpk2dGamma32223 + dGammadU3215*dpk2dGamma32321 + dGammadU3225*dpk2dGamma32322 + dGammadU3235*dpk2dGamma32323 + dPhidU125*dpk2dPsi3212 + dPhidU225*dpk2dPsi3222 + dPhidU325*dpk2dPsi3232
    dpk2dU[I8,I5] = dGammadU1215*dpk2dGamma31121 + dGammadU1225*dpk2dGamma31122 + dGammadU1235*dpk2dGamma31123 + dGammadU2215*dpk2dGamma31221 + dGammadU2225*dpk2dGamma31222 + dGammadU2235*dpk2dGamma31223 + dGammadU3215*dpk2dGamma31321 + dGammadU3225*dpk2dGamma31322 + dGammadU3235*dpk2dGamma31323 + dPhidU125*dpk2dPsi3112 + dPhidU225*dpk2dPsi3122 + dPhidU325*dpk2dPsi3132
    dpk2dU[I9,I5] = dGammadU1215*dpk2dGamma21121 + dGammadU1225*dpk2dGamma21122 + dGammadU1235*dpk2dGamma21123 + dGammadU2215*dpk2dGamma21221 + dGammadU2225*dpk2dGamma21222 + dGammadU2235*dpk2dGamma21223 + dGammadU3215*dpk2dGamma21321 + dGammadU3225*dpk2dGamma21322 + dGammadU3235*dpk2dGamma21323 + dPhidU125*dpk2dPsi2112 + dPhidU225*dpk2dPsi2122 + dPhidU325*dpk2dPsi2132

    #Column 6
    dpk2dU[I1,I6] = dGammadU1316*dpk2dGamma11131 + dGammadU1326*dpk2dGamma11132 + dGammadU1336*dpk2dGamma11133 + dGammadU2316*dpk2dGamma11231 + dGammadU2326*dpk2dGamma11232 + dGammadU2336*dpk2dGamma11233 + dGammadU3316*dpk2dGamma11331 + dGammadU3326*dpk2dGamma11332 + dGammadU3336*dpk2dGamma11333 + dPhidU136*dpk2dPsi1113 + dPhidU236*dpk2dPsi1123 + dPhidU336*dpk2dPsi1133
    dpk2dU[I2,I6] = dGammadU1316*dpk2dGamma22131 + dGammadU1326*dpk2dGamma22132 + dGammadU1336*dpk2dGamma22133 + dGammadU2316*dpk2dGamma22231 + dGammadU2326*dpk2dGamma22232 + dGammadU2336*dpk2dGamma22233 + dGammadU3316*dpk2dGamma22331 + dGammadU3326*dpk2dGamma22332 + dGammadU3336*dpk2dGamma22333 + dPhidU136*dpk2dPsi2213 + dPhidU236*dpk2dPsi2223 + dPhidU336*dpk2dPsi2233
    dpk2dU[I3,I6] = dGammadU1316*dpk2dGamma33131 + dGammadU1326*dpk2dGamma33132 + dGammadU1336*dpk2dGamma33133 + dGammadU2316*dpk2dGamma33231 + dGammadU2326*dpk2dGamma33232 + dGammadU2336*dpk2dGamma33233 + dGammadU3316*dpk2dGamma33331 + dGammadU3326*dpk2dGamma33332 + dGammadU3336*dpk2dGamma33333 + dPhidU136*dpk2dPsi3313 + dPhidU236*dpk2dPsi3323 + dPhidU336*dpk2dPsi3333
    dpk2dU[I4,I6] = dGammadU1316*dpk2dGamma23131 + dGammadU1326*dpk2dGamma23132 + dGammadU1336*dpk2dGamma23133 + dGammadU2316*dpk2dGamma23231 + dGammadU2326*dpk2dGamma23232 + dGammadU2336*dpk2dGamma23233 + dGammadU3316*dpk2dGamma23331 + dGammadU3326*dpk2dGamma23332 + dGammadU3336*dpk2dGamma23333 + dPhidU136*dpk2dPsi2313 + dPhidU236*dpk2dPsi2323 + dPhidU336*dpk2dPsi2333
    dpk2dU[I5,I6] = dGammadU1316*dpk2dGamma13131 + dGammadU1326*dpk2dGamma13132 + dGammadU1336*dpk2dGamma13133 + dGammadU2316*dpk2dGamma13231 + dGammadU2326*dpk2dGamma13232 + dGammadU2336*dpk2dGamma13233 + dGammadU3316*dpk2dGamma13331 + dGammadU3326*dpk2dGamma13332 + dGammadU3336*dpk2dGamma13333 + dPhidU136*dpk2dPsi1313 + dPhidU236*dpk2dPsi1323 + dPhidU336*dpk2dPsi1333
    dpk2dU[I6,I6] = dGammadU1316*dpk2dGamma12131 + dGammadU1326*dpk2dGamma12132 + dGammadU1336*dpk2dGamma12133 + dGammadU2316*dpk2dGamma12231 + dGammadU2326*dpk2dGamma12232 + dGammadU2336*dpk2dGamma12233 + dGammadU3316*dpk2dGamma12331 + dGammadU3326*dpk2dGamma12332 + dGammadU3336*dpk2dGamma12333 + dPhidU136*dpk2dPsi1213 + dPhidU236*dpk2dPsi1223 + dPhidU336*dpk2dPsi1233
    dpk2dU[I7,I6] = dGammadU1316*dpk2dGamma32131 + dGammadU1326*dpk2dGamma32132 + dGammadU1336*dpk2dGamma32133 + dGammadU2316*dpk2dGamma32231 + dGammadU2326*dpk2dGamma32232 + dGammadU2336*dpk2dGamma32233 + dGammadU3316*dpk2dGamma32331 + dGammadU3326*dpk2dGamma32332 + dGammadU3336*dpk2dGamma32333 + dPhidU136*dpk2dPsi3213 + dPhidU236*dpk2dPsi3223 + dPhidU336*dpk2dPsi3233
    dpk2dU[I8,I6] = dGammadU1316*dpk2dGamma31131 + dGammadU1326*dpk2dGamma31132 + dGammadU1336*dpk2dGamma31133 + dGammadU2316*dpk2dGamma31231 + dGammadU2326*dpk2dGamma31232 + dGammadU2336*dpk2dGamma31233 + dGammadU3316*dpk2dGamma31331 + dGammadU3326*dpk2dGamma31332 + dGammadU3336*dpk2dGamma31333 + dPhidU136*dpk2dPsi3113 + dPhidU236*dpk2dPsi3123 + dPhidU336*dpk2dPsi3133
    dpk2dU[I9,I6] = dGammadU1316*dpk2dGamma21131 + dGammadU1326*dpk2dGamma21132 + dGammadU1336*dpk2dGamma21133 + dGammadU2316*dpk2dGamma21231 + dGammadU2326*dpk2dGamma21232 + dGammadU2336*dpk2dGamma21233 + dGammadU3316*dpk2dGamma21331 + dGammadU3326*dpk2dGamma21332 + dGammadU3336*dpk2dGamma21333 + dPhidU136*dpk2dPsi2113 + dPhidU236*dpk2dPsi2123 + dPhidU336*dpk2dPsi2133

    #Column 7
    dpk2dU[I1,I7] = dGammadU1317*dpk2dGamma11131 + dGammadU1327*dpk2dGamma11132 + dGammadU1337*dpk2dGamma11133 + dGammadU2317*dpk2dGamma11231 + dGammadU2327*dpk2dGamma11232 + dGammadU2337*dpk2dGamma11233 + dGammadU3317*dpk2dGamma11331 + dGammadU3327*dpk2dGamma11332 + dGammadU3337*dpk2dGamma11333 + dPhidU137*dpk2dPsi1113 + dPhidU237*dpk2dPsi1123 + dPhidU337*dpk2dPsi1133
    dpk2dU[I2,I7] = dGammadU1317*dpk2dGamma22131 + dGammadU1327*dpk2dGamma22132 + dGammadU1337*dpk2dGamma22133 + dGammadU2317*dpk2dGamma22231 + dGammadU2327*dpk2dGamma22232 + dGammadU2337*dpk2dGamma22233 + dGammadU3317*dpk2dGamma22331 + dGammadU3327*dpk2dGamma22332 + dGammadU3337*dpk2dGamma22333 + dPhidU137*dpk2dPsi2213 + dPhidU237*dpk2dPsi2223 + dPhidU337*dpk2dPsi2233
    dpk2dU[I3,I7] = dGammadU1317*dpk2dGamma33131 + dGammadU1327*dpk2dGamma33132 + dGammadU1337*dpk2dGamma33133 + dGammadU2317*dpk2dGamma33231 + dGammadU2327*dpk2dGamma33232 + dGammadU2337*dpk2dGamma33233 + dGammadU3317*dpk2dGamma33331 + dGammadU3327*dpk2dGamma33332 + dGammadU3337*dpk2dGamma33333 + dPhidU137*dpk2dPsi3313 + dPhidU237*dpk2dPsi3323 + dPhidU337*dpk2dPsi3333
    dpk2dU[I4,I7] = dGammadU1317*dpk2dGamma23131 + dGammadU1327*dpk2dGamma23132 + dGammadU1337*dpk2dGamma23133 + dGammadU2317*dpk2dGamma23231 + dGammadU2327*dpk2dGamma23232 + dGammadU2337*dpk2dGamma23233 + dGammadU3317*dpk2dGamma23331 + dGammadU3327*dpk2dGamma23332 + dGammadU3337*dpk2dGamma23333 + dPhidU137*dpk2dPsi2313 + dPhidU237*dpk2dPsi2323 + dPhidU337*dpk2dPsi2333
    dpk2dU[I5,I7] = dGammadU1317*dpk2dGamma13131 + dGammadU1327*dpk2dGamma13132 + dGammadU1337*dpk2dGamma13133 + dGammadU2317*dpk2dGamma13231 + dGammadU2327*dpk2dGamma13232 + dGammadU2337*dpk2dGamma13233 + dGammadU3317*dpk2dGamma13331 + dGammadU3327*dpk2dGamma13332 + dGammadU3337*dpk2dGamma13333 + dPhidU137*dpk2dPsi1313 + dPhidU237*dpk2dPsi1323 + dPhidU337*dpk2dPsi1333
    dpk2dU[I6,I7] = dGammadU1317*dpk2dGamma12131 + dGammadU1327*dpk2dGamma12132 + dGammadU1337*dpk2dGamma12133 + dGammadU2317*dpk2dGamma12231 + dGammadU2327*dpk2dGamma12232 + dGammadU2337*dpk2dGamma12233 + dGammadU3317*dpk2dGamma12331 + dGammadU3327*dpk2dGamma12332 + dGammadU3337*dpk2dGamma12333 + dPhidU137*dpk2dPsi1213 + dPhidU237*dpk2dPsi1223 + dPhidU337*dpk2dPsi1233
    dpk2dU[I7,I7] = dGammadU1317*dpk2dGamma32131 + dGammadU1327*dpk2dGamma32132 + dGammadU1337*dpk2dGamma32133 + dGammadU2317*dpk2dGamma32231 + dGammadU2327*dpk2dGamma32232 + dGammadU2337*dpk2dGamma32233 + dGammadU3317*dpk2dGamma32331 + dGammadU3327*dpk2dGamma32332 + dGammadU3337*dpk2dGamma32333 + dPhidU137*dpk2dPsi3213 + dPhidU237*dpk2dPsi3223 + dPhidU337*dpk2dPsi3233
    dpk2dU[I8,I7] = dGammadU1317*dpk2dGamma31131 + dGammadU1327*dpk2dGamma31132 + dGammadU1337*dpk2dGamma31133 + dGammadU2317*dpk2dGamma31231 + dGammadU2327*dpk2dGamma31232 + dGammadU2337*dpk2dGamma31233 + dGammadU3317*dpk2dGamma31331 + dGammadU3327*dpk2dGamma31332 + dGammadU3337*dpk2dGamma31333 + dPhidU137*dpk2dPsi3113 + dPhidU237*dpk2dPsi3123 + dPhidU337*dpk2dPsi3133
    dpk2dU[I9,I7] = dGammadU1317*dpk2dGamma21131 + dGammadU1327*dpk2dGamma21132 + dGammadU1337*dpk2dGamma21133 + dGammadU2317*dpk2dGamma21231 + dGammadU2327*dpk2dGamma21232 + dGammadU2337*dpk2dGamma21233 + dGammadU3317*dpk2dGamma21331 + dGammadU3327*dpk2dGamma21332 + dGammadU3337*dpk2dGamma21333 + dPhidU137*dpk2dPsi2113 + dPhidU237*dpk2dPsi2123 + dPhidU337*dpk2dPsi2133

    #Column 8
    dpk2dU[I1,I8] = dGammadU1318*dpk2dGamma11131 + dGammadU1328*dpk2dGamma11132 + dGammadU1338*dpk2dGamma11133 + dGammadU2318*dpk2dGamma11231 + dGammadU2328*dpk2dGamma11232 + dGammadU2338*dpk2dGamma11233 + dGammadU3318*dpk2dGamma11331 + dGammadU3328*dpk2dGamma11332 + dGammadU3338*dpk2dGamma11333 + dPhidU138*dpk2dPsi1113 + dPhidU238*dpk2dPsi1123 + dPhidU338*dpk2dPsi1133
    dpk2dU[I2,I8] = dGammadU1318*dpk2dGamma22131 + dGammadU1328*dpk2dGamma22132 + dGammadU1338*dpk2dGamma22133 + dGammadU2318*dpk2dGamma22231 + dGammadU2328*dpk2dGamma22232 + dGammadU2338*dpk2dGamma22233 + dGammadU3318*dpk2dGamma22331 + dGammadU3328*dpk2dGamma22332 + dGammadU3338*dpk2dGamma22333 + dPhidU138*dpk2dPsi2213 + dPhidU238*dpk2dPsi2223 + dPhidU338*dpk2dPsi2233
    dpk2dU[I3,I8] = dGammadU1318*dpk2dGamma33131 + dGammadU1328*dpk2dGamma33132 + dGammadU1338*dpk2dGamma33133 + dGammadU2318*dpk2dGamma33231 + dGammadU2328*dpk2dGamma33232 + dGammadU2338*dpk2dGamma33233 + dGammadU3318*dpk2dGamma33331 + dGammadU3328*dpk2dGamma33332 + dGammadU3338*dpk2dGamma33333 + dPhidU138*dpk2dPsi3313 + dPhidU238*dpk2dPsi3323 + dPhidU338*dpk2dPsi3333
    dpk2dU[I4,I8] = dGammadU1318*dpk2dGamma23131 + dGammadU1328*dpk2dGamma23132 + dGammadU1338*dpk2dGamma23133 + dGammadU2318*dpk2dGamma23231 + dGammadU2328*dpk2dGamma23232 + dGammadU2338*dpk2dGamma23233 + dGammadU3318*dpk2dGamma23331 + dGammadU3328*dpk2dGamma23332 + dGammadU3338*dpk2dGamma23333 + dPhidU138*dpk2dPsi2313 + dPhidU238*dpk2dPsi2323 + dPhidU338*dpk2dPsi2333
    dpk2dU[I5,I8] = dGammadU1318*dpk2dGamma13131 + dGammadU1328*dpk2dGamma13132 + dGammadU1338*dpk2dGamma13133 + dGammadU2318*dpk2dGamma13231 + dGammadU2328*dpk2dGamma13232 + dGammadU2338*dpk2dGamma13233 + dGammadU3318*dpk2dGamma13331 + dGammadU3328*dpk2dGamma13332 + dGammadU3338*dpk2dGamma13333 + dPhidU138*dpk2dPsi1313 + dPhidU238*dpk2dPsi1323 + dPhidU338*dpk2dPsi1333
    dpk2dU[I6,I8] = dGammadU1318*dpk2dGamma12131 + dGammadU1328*dpk2dGamma12132 + dGammadU1338*dpk2dGamma12133 + dGammadU2318*dpk2dGamma12231 + dGammadU2328*dpk2dGamma12232 + dGammadU2338*dpk2dGamma12233 + dGammadU3318*dpk2dGamma12331 + dGammadU3328*dpk2dGamma12332 + dGammadU3338*dpk2dGamma12333 + dPhidU138*dpk2dPsi1213 + dPhidU238*dpk2dPsi1223 + dPhidU338*dpk2dPsi1233
    dpk2dU[I7,I8] = dGammadU1318*dpk2dGamma32131 + dGammadU1328*dpk2dGamma32132 + dGammadU1338*dpk2dGamma32133 + dGammadU2318*dpk2dGamma32231 + dGammadU2328*dpk2dGamma32232 + dGammadU2338*dpk2dGamma32233 + dGammadU3318*dpk2dGamma32331 + dGammadU3328*dpk2dGamma32332 + dGammadU3338*dpk2dGamma32333 + dPhidU138*dpk2dPsi3213 + dPhidU238*dpk2dPsi3223 + dPhidU338*dpk2dPsi3233
    dpk2dU[I8,I8] = dGammadU1318*dpk2dGamma31131 + dGammadU1328*dpk2dGamma31132 + dGammadU1338*dpk2dGamma31133 + dGammadU2318*dpk2dGamma31231 + dGammadU2328*dpk2dGamma31232 + dGammadU2338*dpk2dGamma31233 + dGammadU3318*dpk2dGamma31331 + dGammadU3328*dpk2dGamma31332 + dGammadU3338*dpk2dGamma31333 + dPhidU138*dpk2dPsi3113 + dPhidU238*dpk2dPsi3123 + dPhidU338*dpk2dPsi3133
    dpk2dU[I9,I8] = dGammadU1318*dpk2dGamma21131 + dGammadU1328*dpk2dGamma21132 + dGammadU1338*dpk2dGamma21133 + dGammadU2318*dpk2dGamma21231 + dGammadU2328*dpk2dGamma21232 + dGammadU2338*dpk2dGamma21233 + dGammadU3318*dpk2dGamma21331 + dGammadU3328*dpk2dGamma21332 + dGammadU3338*dpk2dGamma21333 + dPhidU138*dpk2dPsi2113 + dPhidU238*dpk2dPsi2123 + dPhidU338*dpk2dPsi2133

    #Column 9
    dpk2dU[I1,I9] = dGammadU1219*dpk2dGamma11121 + dGammadU1229*dpk2dGamma11122 + dGammadU1239*dpk2dGamma11123 + dGammadU2219*dpk2dGamma11221 + dGammadU2229*dpk2dGamma11222 + dGammadU2239*dpk2dGamma11223 + dGammadU3219*dpk2dGamma11321 + dGammadU3229*dpk2dGamma11322 + dGammadU3239*dpk2dGamma11323 + dPhidU129*dpk2dPsi1112 + dPhidU229*dpk2dPsi1122 + dPhidU329*dpk2dPsi1132
    dpk2dU[I2,I9] = dGammadU1219*dpk2dGamma22121 + dGammadU1229*dpk2dGamma22122 + dGammadU1239*dpk2dGamma22123 + dGammadU2219*dpk2dGamma22221 + dGammadU2229*dpk2dGamma22222 + dGammadU2239*dpk2dGamma22223 + dGammadU3219*dpk2dGamma22321 + dGammadU3229*dpk2dGamma22322 + dGammadU3239*dpk2dGamma22323 + dPhidU129*dpk2dPsi2212 + dPhidU229*dpk2dPsi2222 + dPhidU329*dpk2dPsi2232
    dpk2dU[I3,I9] = dGammadU1219*dpk2dGamma33121 + dGammadU1229*dpk2dGamma33122 + dGammadU1239*dpk2dGamma33123 + dGammadU2219*dpk2dGamma33221 + dGammadU2229*dpk2dGamma33222 + dGammadU2239*dpk2dGamma33223 + dGammadU3219*dpk2dGamma33321 + dGammadU3229*dpk2dGamma33322 + dGammadU3239*dpk2dGamma33323 + dPhidU129*dpk2dPsi3312 + dPhidU229*dpk2dPsi3322 + dPhidU329*dpk2dPsi3332
    dpk2dU[I4,I9] = dGammadU1219*dpk2dGamma23121 + dGammadU1229*dpk2dGamma23122 + dGammadU1239*dpk2dGamma23123 + dGammadU2219*dpk2dGamma23221 + dGammadU2229*dpk2dGamma23222 + dGammadU2239*dpk2dGamma23223 + dGammadU3219*dpk2dGamma23321 + dGammadU3229*dpk2dGamma23322 + dGammadU3239*dpk2dGamma23323 + dPhidU129*dpk2dPsi2312 + dPhidU229*dpk2dPsi2322 + dPhidU329*dpk2dPsi2332
    dpk2dU[I5,I9] = dGammadU1219*dpk2dGamma13121 + dGammadU1229*dpk2dGamma13122 + dGammadU1239*dpk2dGamma13123 + dGammadU2219*dpk2dGamma13221 + dGammadU2229*dpk2dGamma13222 + dGammadU2239*dpk2dGamma13223 + dGammadU3219*dpk2dGamma13321 + dGammadU3229*dpk2dGamma13322 + dGammadU3239*dpk2dGamma13323 + dPhidU129*dpk2dPsi1312 + dPhidU229*dpk2dPsi1322 + dPhidU329*dpk2dPsi1332
    dpk2dU[I6,I9] = dGammadU1219*dpk2dGamma12121 + dGammadU1229*dpk2dGamma12122 + dGammadU1239*dpk2dGamma12123 + dGammadU2219*dpk2dGamma12221 + dGammadU2229*dpk2dGamma12222 + dGammadU2239*dpk2dGamma12223 + dGammadU3219*dpk2dGamma12321 + dGammadU3229*dpk2dGamma12322 + dGammadU3239*dpk2dGamma12323 + dPhidU129*dpk2dPsi1212 + dPhidU229*dpk2dPsi1222 + dPhidU329*dpk2dPsi1232
    dpk2dU[I7,I9] = dGammadU1219*dpk2dGamma32121 + dGammadU1229*dpk2dGamma32122 + dGammadU1239*dpk2dGamma32123 + dGammadU2219*dpk2dGamma32221 + dGammadU2229*dpk2dGamma32222 + dGammadU2239*dpk2dGamma32223 + dGammadU3219*dpk2dGamma32321 + dGammadU3229*dpk2dGamma32322 + dGammadU3239*dpk2dGamma32323 + dPhidU129*dpk2dPsi3212 + dPhidU229*dpk2dPsi3222 + dPhidU329*dpk2dPsi3232
    dpk2dU[I8,I9] = dGammadU1219*dpk2dGamma31121 + dGammadU1229*dpk2dGamma31122 + dGammadU1239*dpk2dGamma31123 + dGammadU2219*dpk2dGamma31221 + dGammadU2229*dpk2dGamma31222 + dGammadU2239*dpk2dGamma31223 + dGammadU3219*dpk2dGamma31321 + dGammadU3229*dpk2dGamma31322 + dGammadU3239*dpk2dGamma31323 + dPhidU129*dpk2dPsi3112 + dPhidU229*dpk2dPsi3122 + dPhidU329*dpk2dPsi3132
    dpk2dU[I9,I9] = dGammadU1219*dpk2dGamma21121 + dGammadU1229*dpk2dGamma21122 + dGammadU1239*dpk2dGamma21123 + dGammadU2219*dpk2dGamma21221 + dGammadU2229*dpk2dGamma21222 + dGammadU2239*dpk2dGamma21223 + dGammadU3219*dpk2dGamma21321 + dGammadU3229*dpk2dGamma21322 + dGammadU3239*dpk2dGamma21323 + dPhidU129*dpk2dPsi2112 + dPhidU229*dpk2dPsi2122 + dPhidU329*dpk2dPsi2132

    #Column 10
    dpk2dU[I1,I10] = dGammadU12110*dpk2dGamma11121 + dGammadU12210*dpk2dGamma11122 + dGammadU12310*dpk2dGamma11123 + dGammadU22110*dpk2dGamma11221 + dGammadU22210*dpk2dGamma11222 + dGammadU22310*dpk2dGamma11223 + dGammadU32110*dpk2dGamma11321 + dGammadU32210*dpk2dGamma11322 + dGammadU32310*dpk2dGamma11323 + dPhidU1210*dpk2dPsi1112 + dPhidU2210*dpk2dPsi1122 + dPhidU3210*dpk2dPsi1132
    dpk2dU[I2,I10] = dGammadU12110*dpk2dGamma22121 + dGammadU12210*dpk2dGamma22122 + dGammadU12310*dpk2dGamma22123 + dGammadU22110*dpk2dGamma22221 + dGammadU22210*dpk2dGamma22222 + dGammadU22310*dpk2dGamma22223 + dGammadU32110*dpk2dGamma22321 + dGammadU32210*dpk2dGamma22322 + dGammadU32310*dpk2dGamma22323 + dPhidU1210*dpk2dPsi2212 + dPhidU2210*dpk2dPsi2222 + dPhidU3210*dpk2dPsi2232
    dpk2dU[I3,I10] = dGammadU12110*dpk2dGamma33121 + dGammadU12210*dpk2dGamma33122 + dGammadU12310*dpk2dGamma33123 + dGammadU22110*dpk2dGamma33221 + dGammadU22210*dpk2dGamma33222 + dGammadU22310*dpk2dGamma33223 + dGammadU32110*dpk2dGamma33321 + dGammadU32210*dpk2dGamma33322 + dGammadU32310*dpk2dGamma33323 + dPhidU1210*dpk2dPsi3312 + dPhidU2210*dpk2dPsi3322 + dPhidU3210*dpk2dPsi3332
    dpk2dU[I4,I10] = dGammadU12110*dpk2dGamma23121 + dGammadU12210*dpk2dGamma23122 + dGammadU12310*dpk2dGamma23123 + dGammadU22110*dpk2dGamma23221 + dGammadU22210*dpk2dGamma23222 + dGammadU22310*dpk2dGamma23223 + dGammadU32110*dpk2dGamma23321 + dGammadU32210*dpk2dGamma23322 + dGammadU32310*dpk2dGamma23323 + dPhidU1210*dpk2dPsi2312 + dPhidU2210*dpk2dPsi2322 + dPhidU3210*dpk2dPsi2332
    dpk2dU[I5,I10] = dGammadU12110*dpk2dGamma13121 + dGammadU12210*dpk2dGamma13122 + dGammadU12310*dpk2dGamma13123 + dGammadU22110*dpk2dGamma13221 + dGammadU22210*dpk2dGamma13222 + dGammadU22310*dpk2dGamma13223 + dGammadU32110*dpk2dGamma13321 + dGammadU32210*dpk2dGamma13322 + dGammadU32310*dpk2dGamma13323 + dPhidU1210*dpk2dPsi1312 + dPhidU2210*dpk2dPsi1322 + dPhidU3210*dpk2dPsi1332
    dpk2dU[I6,I10] = dGammadU12110*dpk2dGamma12121 + dGammadU12210*dpk2dGamma12122 + dGammadU12310*dpk2dGamma12123 + dGammadU22110*dpk2dGamma12221 + dGammadU22210*dpk2dGamma12222 + dGammadU22310*dpk2dGamma12223 + dGammadU32110*dpk2dGamma12321 + dGammadU32210*dpk2dGamma12322 + dGammadU32310*dpk2dGamma12323 + dPhidU1210*dpk2dPsi1212 + dPhidU2210*dpk2dPsi1222 + dPhidU3210*dpk2dPsi1232
    dpk2dU[I7,I10] = dGammadU12110*dpk2dGamma32121 + dGammadU12210*dpk2dGamma32122 + dGammadU12310*dpk2dGamma32123 + dGammadU22110*dpk2dGamma32221 + dGammadU22210*dpk2dGamma32222 + dGammadU22310*dpk2dGamma32223 + dGammadU32110*dpk2dGamma32321 + dGammadU32210*dpk2dGamma32322 + dGammadU32310*dpk2dGamma32323 + dPhidU1210*dpk2dPsi3212 + dPhidU2210*dpk2dPsi3222 + dPhidU3210*dpk2dPsi3232
    dpk2dU[I8,I10] = dGammadU12110*dpk2dGamma31121 + dGammadU12210*dpk2dGamma31122 + dGammadU12310*dpk2dGamma31123 + dGammadU22110*dpk2dGamma31221 + dGammadU22210*dpk2dGamma31222 + dGammadU22310*dpk2dGamma31223 + dGammadU32110*dpk2dGamma31321 + dGammadU32210*dpk2dGamma31322 + dGammadU32310*dpk2dGamma31323 + dPhidU1210*dpk2dPsi3112 + dPhidU2210*dpk2dPsi3122 + dPhidU3210*dpk2dPsi3132
    dpk2dU[I9,I10] = dGammadU12110*dpk2dGamma21121 + dGammadU12210*dpk2dGamma21122 + dGammadU12310*dpk2dGamma21123 + dGammadU22110*dpk2dGamma21221 + dGammadU22210*dpk2dGamma21222 + dGammadU22310*dpk2dGamma21223 + dGammadU32110*dpk2dGamma21321 + dGammadU32210*dpk2dGamma21322 + dGammadU32310*dpk2dGamma21323 + dPhidU1210*dpk2dPsi2112 + dPhidU2210*dpk2dPsi2122 + dPhidU3210*dpk2dPsi2132

    #Column 11
    dpk2dU[I1,I11] = dGammadU11111*dpk2dGamma11111 + dGammadU11211*dpk2dGamma11112 + dGammadU11311*dpk2dGamma11113 + dGammadU21111*dpk2dGamma11211 + dGammadU21211*dpk2dGamma11212 + dGammadU21311*dpk2dGamma11213 + dGammadU31111*dpk2dGamma11311 + dGammadU31211*dpk2dGamma11312 + dGammadU31311*dpk2dGamma11313 + dPhidU1111*dpk2dPsi1111 + dPhidU2111*dpk2dPsi1121 + dPhidU3111*dpk2dPsi1131
    dpk2dU[I2,I11] = dGammadU11111*dpk2dGamma22111 + dGammadU11211*dpk2dGamma22112 + dGammadU11311*dpk2dGamma22113 + dGammadU21111*dpk2dGamma22211 + dGammadU21211*dpk2dGamma22212 + dGammadU21311*dpk2dGamma22213 + dGammadU31111*dpk2dGamma22311 + dGammadU31211*dpk2dGamma22312 + dGammadU31311*dpk2dGamma22313 + dPhidU1111*dpk2dPsi2211 + dPhidU2111*dpk2dPsi2221 + dPhidU3111*dpk2dPsi2231
    dpk2dU[I3,I11] = dGammadU11111*dpk2dGamma33111 + dGammadU11211*dpk2dGamma33112 + dGammadU11311*dpk2dGamma33113 + dGammadU21111*dpk2dGamma33211 + dGammadU21211*dpk2dGamma33212 + dGammadU21311*dpk2dGamma33213 + dGammadU31111*dpk2dGamma33311 + dGammadU31211*dpk2dGamma33312 + dGammadU31311*dpk2dGamma33313 + dPhidU1111*dpk2dPsi3311 + dPhidU2111*dpk2dPsi3321 + dPhidU3111*dpk2dPsi3331
    dpk2dU[I4,I11] = dGammadU11111*dpk2dGamma23111 + dGammadU11211*dpk2dGamma23112 + dGammadU11311*dpk2dGamma23113 + dGammadU21111*dpk2dGamma23211 + dGammadU21211*dpk2dGamma23212 + dGammadU21311*dpk2dGamma23213 + dGammadU31111*dpk2dGamma23311 + dGammadU31211*dpk2dGamma23312 + dGammadU31311*dpk2dGamma23313 + dPhidU1111*dpk2dPsi2311 + dPhidU2111*dpk2dPsi2321 + dPhidU3111*dpk2dPsi2331
    dpk2dU[I5,I11] = dGammadU11111*dpk2dGamma13111 + dGammadU11211*dpk2dGamma13112 + dGammadU11311*dpk2dGamma13113 + dGammadU21111*dpk2dGamma13211 + dGammadU21211*dpk2dGamma13212 + dGammadU21311*dpk2dGamma13213 + dGammadU31111*dpk2dGamma13311 + dGammadU31211*dpk2dGamma13312 + dGammadU31311*dpk2dGamma13313 + dPhidU1111*dpk2dPsi1311 + dPhidU2111*dpk2dPsi1321 + dPhidU3111*dpk2dPsi1331
    dpk2dU[I6,I11] = dGammadU11111*dpk2dGamma12111 + dGammadU11211*dpk2dGamma12112 + dGammadU11311*dpk2dGamma12113 + dGammadU21111*dpk2dGamma12211 + dGammadU21211*dpk2dGamma12212 + dGammadU21311*dpk2dGamma12213 + dGammadU31111*dpk2dGamma12311 + dGammadU31211*dpk2dGamma12312 + dGammadU31311*dpk2dGamma12313 + dPhidU1111*dpk2dPsi1211 + dPhidU2111*dpk2dPsi1221 + dPhidU3111*dpk2dPsi1231
    dpk2dU[I7,I11] = dGammadU11111*dpk2dGamma32111 + dGammadU11211*dpk2dGamma32112 + dGammadU11311*dpk2dGamma32113 + dGammadU21111*dpk2dGamma32211 + dGammadU21211*dpk2dGamma32212 + dGammadU21311*dpk2dGamma32213 + dGammadU31111*dpk2dGamma32311 + dGammadU31211*dpk2dGamma32312 + dGammadU31311*dpk2dGamma32313 + dPhidU1111*dpk2dPsi3211 + dPhidU2111*dpk2dPsi3221 + dPhidU3111*dpk2dPsi3231
    dpk2dU[I8,I11] = dGammadU11111*dpk2dGamma31111 + dGammadU11211*dpk2dGamma31112 + dGammadU11311*dpk2dGamma31113 + dGammadU21111*dpk2dGamma31211 + dGammadU21211*dpk2dGamma31212 + dGammadU21311*dpk2dGamma31213 + dGammadU31111*dpk2dGamma31311 + dGammadU31211*dpk2dGamma31312 + dGammadU31311*dpk2dGamma31313 + dPhidU1111*dpk2dPsi3111 + dPhidU2111*dpk2dPsi3121 + dPhidU3111*dpk2dPsi3131
    dpk2dU[I9,I11] = dGammadU11111*dpk2dGamma21111 + dGammadU11211*dpk2dGamma21112 + dGammadU11311*dpk2dGamma21113 + dGammadU21111*dpk2dGamma21211 + dGammadU21211*dpk2dGamma21212 + dGammadU21311*dpk2dGamma21213 + dGammadU31111*dpk2dGamma21311 + dGammadU31211*dpk2dGamma21312 + dGammadU31311*dpk2dGamma21313 + dPhidU1111*dpk2dPsi2111 + dPhidU2111*dpk2dPsi2121 + dPhidU3111*dpk2dPsi2131

    #Column 12
    dpk2dU[I1,I12] = dGammadU11112*dpk2dGamma11111 + dGammadU11212*dpk2dGamma11112 + dGammadU11312*dpk2dGamma11113 + dGammadU21112*dpk2dGamma11211 + dGammadU21212*dpk2dGamma11212 + dGammadU21312*dpk2dGamma11213 + dGammadU31112*dpk2dGamma11311 + dGammadU31212*dpk2dGamma11312 + dGammadU31312*dpk2dGamma11313 + dPhidU1112*dpk2dPsi1111 + dPhidU2112*dpk2dPsi1121 + dPhidU3112*dpk2dPsi1131
    dpk2dU[I2,I12] = dGammadU11112*dpk2dGamma22111 + dGammadU11212*dpk2dGamma22112 + dGammadU11312*dpk2dGamma22113 + dGammadU21112*dpk2dGamma22211 + dGammadU21212*dpk2dGamma22212 + dGammadU21312*dpk2dGamma22213 + dGammadU31112*dpk2dGamma22311 + dGammadU31212*dpk2dGamma22312 + dGammadU31312*dpk2dGamma22313 + dPhidU1112*dpk2dPsi2211 + dPhidU2112*dpk2dPsi2221 + dPhidU3112*dpk2dPsi2231
    dpk2dU[I3,I12] = dGammadU11112*dpk2dGamma33111 + dGammadU11212*dpk2dGamma33112 + dGammadU11312*dpk2dGamma33113 + dGammadU21112*dpk2dGamma33211 + dGammadU21212*dpk2dGamma33212 + dGammadU21312*dpk2dGamma33213 + dGammadU31112*dpk2dGamma33311 + dGammadU31212*dpk2dGamma33312 + dGammadU31312*dpk2dGamma33313 + dPhidU1112*dpk2dPsi3311 + dPhidU2112*dpk2dPsi3321 + dPhidU3112*dpk2dPsi3331
    dpk2dU[I4,I12] = dGammadU11112*dpk2dGamma23111 + dGammadU11212*dpk2dGamma23112 + dGammadU11312*dpk2dGamma23113 + dGammadU21112*dpk2dGamma23211 + dGammadU21212*dpk2dGamma23212 + dGammadU21312*dpk2dGamma23213 + dGammadU31112*dpk2dGamma23311 + dGammadU31212*dpk2dGamma23312 + dGammadU31312*dpk2dGamma23313 + dPhidU1112*dpk2dPsi2311 + dPhidU2112*dpk2dPsi2321 + dPhidU3112*dpk2dPsi2331
    dpk2dU[I5,I12] = dGammadU11112*dpk2dGamma13111 + dGammadU11212*dpk2dGamma13112 + dGammadU11312*dpk2dGamma13113 + dGammadU21112*dpk2dGamma13211 + dGammadU21212*dpk2dGamma13212 + dGammadU21312*dpk2dGamma13213 + dGammadU31112*dpk2dGamma13311 + dGammadU31212*dpk2dGamma13312 + dGammadU31312*dpk2dGamma13313 + dPhidU1112*dpk2dPsi1311 + dPhidU2112*dpk2dPsi1321 + dPhidU3112*dpk2dPsi1331
    dpk2dU[I6,I12] = dGammadU11112*dpk2dGamma12111 + dGammadU11212*dpk2dGamma12112 + dGammadU11312*dpk2dGamma12113 + dGammadU21112*dpk2dGamma12211 + dGammadU21212*dpk2dGamma12212 + dGammadU21312*dpk2dGamma12213 + dGammadU31112*dpk2dGamma12311 + dGammadU31212*dpk2dGamma12312 + dGammadU31312*dpk2dGamma12313 + dPhidU1112*dpk2dPsi1211 + dPhidU2112*dpk2dPsi1221 + dPhidU3112*dpk2dPsi1231
    dpk2dU[I7,I12] = dGammadU11112*dpk2dGamma32111 + dGammadU11212*dpk2dGamma32112 + dGammadU11312*dpk2dGamma32113 + dGammadU21112*dpk2dGamma32211 + dGammadU21212*dpk2dGamma32212 + dGammadU21312*dpk2dGamma32213 + dGammadU31112*dpk2dGamma32311 + dGammadU31212*dpk2dGamma32312 + dGammadU31312*dpk2dGamma32313 + dPhidU1112*dpk2dPsi3211 + dPhidU2112*dpk2dPsi3221 + dPhidU3112*dpk2dPsi3231
    dpk2dU[I8,I12] = dGammadU11112*dpk2dGamma31111 + dGammadU11212*dpk2dGamma31112 + dGammadU11312*dpk2dGamma31113 + dGammadU21112*dpk2dGamma31211 + dGammadU21212*dpk2dGamma31212 + dGammadU21312*dpk2dGamma31213 + dGammadU31112*dpk2dGamma31311 + dGammadU31212*dpk2dGamma31312 + dGammadU31312*dpk2dGamma31313 + dPhidU1112*dpk2dPsi3111 + dPhidU2112*dpk2dPsi3121 + dPhidU3112*dpk2dPsi3131
    dpk2dU[I9,I12] = dGammadU11112*dpk2dGamma21111 + dGammadU11212*dpk2dGamma21112 + dGammadU11312*dpk2dGamma21113 + dGammadU21112*dpk2dGamma21211 + dGammadU21212*dpk2dGamma21212 + dGammadU21312*dpk2dGamma21213 + dGammadU31112*dpk2dGamma21311 + dGammadU31212*dpk2dGamma21312 + dGammadU31312*dpk2dGamma21313 + dPhidU1112*dpk2dPsi2111 + dPhidU2112*dpk2dPsi2121 + dPhidU3112*dpk2dPsi2131

    return dpk2dU
    
def compute_dsymmetric_stressdU(dSigmadC,dSigmadPhi,dSigmadGamma,dCdU,dPhidU,dGammadU):
    """Compute the derivative of the symmetric stress with respect 
    to the degree of freedom vector"""
    
    #The calculations for the tangent for the pk2 stress and the symmetric stress are identical
    #therefore
    
    return compute_dpk2dU(dSigmadC,dSigmadPhi,dSigmadGamma,dCdU,dPhidU,dGammadU)
    
def compute_dho_stressdU(dMdC,dMdPhi,dMdGamma,dCdU,dPhidU,dGammadU):
    """Compute the derivative of the symmetric stress with respect
    to the degree of freedom vector"""
    
    #Initialize
    dMdU = np.zeros([27,12*8])
    
    for n in range(8):
        dMdUn = compute_dho_stressdUn(dMdC[:,(n*12):((n+1)*12)],dMdPhi[:,(n*12):((n+1)*12)],dMdGamma[:,(n*12):((n+1)*12)],\
                                      dCdU[:,(n*12):((n+1)*12)],dPhidU[:,(n*12):((n+1)*12)],dGammadU[:,(n*12):((n+1)*12)])
        for i in range(9):
            for j in range(12):
                dMdU[i,j+n*12] = dMdUn[i,j]
                
    return dMdU
      
def compute_dho_stressdUn(dMdC,dMdPhi,dMdGamma,dCdU,dPhidU,dGammadU):
    """Compute a submatrix of the derivative of the symmetric stress with respect
    to the degree of freedom vector"""
    
    #Initialize the derivative
    dMdU = np.zeros([27,12])
    
    #Set the indices
    I1  = 0
    I2  = 1
    I3  = 2
    I4  = 3
    I5  = 4
    I6  = 5
    I7  = 6
    I8  = 7
    I9  = 8
    I10 = 9
    I11 = 10
    I12 = 11
    
    #Extract values
    
    
    #Extract components of dCdU
    dCdU111 = dCdU[I1,I1]
    dCdU221 = dCdU[I2,I1]
    dCdU331 = dCdU[I3,I1]
    dCdU231 = dCdU[I4,I1]
    dCdU131 = dCdU[I5,I1]
    dCdU121 = dCdU[I6,I1]
    dCdU321 = dCdU[I7,I1]
    dCdU311 = dCdU[I8,I1]
    dCdU211 = dCdU[I9,I1]
    dCdU112 = dCdU[I1,I2]
    dCdU222 = dCdU[I2,I2]
    dCdU332 = dCdU[I3,I2]
    dCdU232 = dCdU[I4,I2]
    dCdU132 = dCdU[I5,I2]
    dCdU122 = dCdU[I6,I2]
    dCdU322 = dCdU[I7,I2]
    dCdU312 = dCdU[I8,I2]
    dCdU212 = dCdU[I9,I2]
    dCdU113 = dCdU[I1,I3]
    dCdU223 = dCdU[I2,I3]
    dCdU333 = dCdU[I3,I3]
    dCdU233 = dCdU[I4,I3]
    dCdU133 = dCdU[I5,I3]
    dCdU123 = dCdU[I6,I3]
    dCdU323 = dCdU[I7,I3]
    dCdU313 = dCdU[I8,I3]
    dCdU213 = dCdU[I9,I3]

    #Extract components of dPsidU
    dPhidU111 = dPsidU[I1,I1]
    dPhidU221 = dPsidU[I2,I1]
    dPhidU331 = dPsidU[I3,I1]
    dPhidU231 = dPsidU[I4,I1]
    dPhidU131 = dPsidU[I5,I1]
    dPhidU121 = dPsidU[I6,I1]
    dPhidU321 = dPsidU[I7,I1]
    dPhidU311 = dPsidU[I8,I1]
    dPhidU211 = dPsidU[I9,I1]
    dPhidU112 = dPsidU[I1,I2]
    dPhidU222 = dPsidU[I2,I2]
    dPhidU332 = dPsidU[I3,I2]
    dPhidU232 = dPsidU[I4,I2]
    dPhidU132 = dPsidU[I5,I2]
    dPhidU122 = dPsidU[I6,I2]
    dPhidU322 = dPsidU[I7,I2]
    dPhidU312 = dPsidU[I8,I2]
    dPhidU212 = dPsidU[I9,I2]
    dPhidU113 = dPsidU[I1,I3]
    dPhidU223 = dPsidU[I2,I3]
    dPhidU333 = dPsidU[I3,I3]
    dPhidU233 = dPsidU[I4,I3]
    dPhidU133 = dPsidU[I5,I3]
    dPhidU123 = dPsidU[I6,I3]
    dPhidU323 = dPsidU[I7,I3]
    dPhidU313 = dPsidU[I8,I3]
    dPhidU213 = dPsidU[I9,I3]
    dPhidU114 = dPsidU[I1,I4]
    dPhidU314 = dPsidU[I8,I4]
    dPhidU214 = dPsidU[I9,I4]
    dPhidU225 = dPsidU[I2,I5]
    dPhidU125 = dPsidU[I6,I5]
    dPhidU325 = dPsidU[I7,I5]
    dPhidU336 = dPsidU[I3,I6]
    dPhidU236 = dPsidU[I4,I6]
    dPhidU136 = dPsidU[I5,I6]
    dPhidU337 = dPsidU[I3,I7]
    dPhidU237 = dPsidU[I4,I7]
    dPhidU137 = dPsidU[I5,I7]
    dPhidU338 = dPsidU[I3,I8]
    dPhidU238 = dPsidU[I4,I8]
    dPhidU138 = dPsidU[I5,I8]
    dPhidU229 = dPsidU[I2,I9]
    dPhidU129 = dPsidU[I6,I9]
    dPhidU329 = dPsidU[I7,I9]
    dPhidU2210 = dPsidU[I2,I10]
    dPhidU1210 = dPsidU[I6,I10]
    dPhidU3210 = dPsidU[I7,I10]
    dPhidU1111 = dPsidU[I1,I11]
    dPhidU3111 = dPsidU[I8,I11]
    dPhidU2111 = dPsidU[I9,I11]
    dPhidU1112 = dPsidU[I1,I12]
    dPhidU3112 = dPsidU[I8,I12]
    dPhidU2112 = dPsidU[I9,I12]

    #Extract components of dGammadU
    dGammadU1111 = dGammadU[I1,I1]
    dGammadU2211 = dGammadU[I2,I1]
    dGammadU3311 = dGammadU[I3,I1]
    dGammadU2311 = dGammadU[I4,I1]
    dGammadU1311 = dGammadU[I5,I1]
    dGammadU1211 = dGammadU[I6,I1]
    dGammadU3211 = dGammadU[I7,I1]
    dGammadU3111 = dGammadU[I8,I1]
    dGammadU2111 = dGammadU[I9,I1]
    dGammadU1121 = dGammadU[I1+9,I1]
    dGammadU2221 = dGammadU[I2+9,I1]
    dGammadU3321 = dGammadU[I3+9,I1]
    dGammadU2321 = dGammadU[I4+9,I1]
    dGammadU1321 = dGammadU[I5+9,I1]
    dGammadU1221 = dGammadU[I6+9,I1]
    dGammadU3221 = dGammadU[I7+9,I1]
    dGammadU3121 = dGammadU[I8+9,I1]
    dGammadU2121 = dGammadU[I9+9,I1]
    dGammadU1131 = dGammadU[I1+18,I1]
    dGammadU2231 = dGammadU[I2+18,I1]
    dGammadU3331 = dGammadU[I3+18,I1]
    dGammadU2331 = dGammadU[I4+18,I1]
    dGammadU1331 = dGammadU[I5+18,I1]
    dGammadU1231 = dGammadU[I6+18,I1]
    dGammadU3231 = dGammadU[I7+18,I1]
    dGammadU3131 = dGammadU[I8+18,I1]
    dGammadU2131 = dGammadU[I9+18,I1]
    dGammadU1112 = dGammadU[I1,I2]
    dGammadU2212 = dGammadU[I2,I2]
    dGammadU3312 = dGammadU[I3,I2]
    dGammadU2312 = dGammadU[I4,I2]
    dGammadU1312 = dGammadU[I5,I2]
    dGammadU1212 = dGammadU[I6,I2]
    dGammadU3212 = dGammadU[I7,I2]
    dGammadU3112 = dGammadU[I8,I2]
    dGammadU2112 = dGammadU[I9,I2]
    dGammadU1122 = dGammadU[I1+9,I2]
    dGammadU2222 = dGammadU[I2+9,I2]
    dGammadU3322 = dGammadU[I3+9,I2]
    dGammadU2322 = dGammadU[I4+9,I2]
    dGammadU1322 = dGammadU[I5+9,I2]
    dGammadU1222 = dGammadU[I6+9,I2]
    dGammadU3222 = dGammadU[I7+9,I2]
    dGammadU3122 = dGammadU[I8+9,I2]
    dGammadU2122 = dGammadU[I9+9,I2]
    dGammadU1132 = dGammadU[I1+18,I2]
    dGammadU2232 = dGammadU[I2+18,I2]
    dGammadU3332 = dGammadU[I3+18,I2]
    dGammadU2332 = dGammadU[I4+18,I2]
    dGammadU1332 = dGammadU[I5+18,I2]
    dGammadU1232 = dGammadU[I6+18,I2]
    dGammadU3232 = dGammadU[I7+18,I2]
    dGammadU3132 = dGammadU[I8+18,I2]
    dGammadU2132 = dGammadU[I9+18,I2]
    dGammadU1113 = dGammadU[I1,I3]
    dGammadU2213 = dGammadU[I2,I3]
    dGammadU3313 = dGammadU[I3,I3]
    dGammadU2313 = dGammadU[I4,I3]
    dGammadU1313 = dGammadU[I5,I3]
    dGammadU1213 = dGammadU[I6,I3]
    dGammadU3213 = dGammadU[I7,I3]
    dGammadU3113 = dGammadU[I8,I3]
    dGammadU2113 = dGammadU[I9,I3]
    dGammadU1123 = dGammadU[I1+9,I3]
    dGammadU2223 = dGammadU[I2+9,I3]
    dGammadU3323 = dGammadU[I3+9,I3]
    dGammadU2323 = dGammadU[I4+9,I3]
    dGammadU1323 = dGammadU[I5+9,I3]
    dGammadU1223 = dGammadU[I6+9,I3]
    dGammadU3223 = dGammadU[I7+9,I3]
    dGammadU3123 = dGammadU[I8+9,I3]
    dGammadU2123 = dGammadU[I9+9,I3]
    dGammadU1133 = dGammadU[I1+18,I3]
    dGammadU2233 = dGammadU[I2+18,I3]
    dGammadU3333 = dGammadU[I3+18,I3]
    dGammadU2333 = dGammadU[I4+18,I3]
    dGammadU1333 = dGammadU[I5+18,I3]
    dGammadU1233 = dGammadU[I6+18,I3]
    dGammadU3233 = dGammadU[I7+18,I3]
    dGammadU3133 = dGammadU[I8+18,I3]
    dGammadU2133 = dGammadU[I9+18,I3]
    dGammadU1114 = dGammadU[I1,I4]
    dGammadU3114 = dGammadU[I8,I4]
    dGammadU2114 = dGammadU[I9,I4]
    dGammadU1124 = dGammadU[I1+9,I4]
    dGammadU3124 = dGammadU[I8+9,I4]
    dGammadU2124 = dGammadU[I9+9,I4]
    dGammadU1134 = dGammadU[I1+18,I4]
    dGammadU3134 = dGammadU[I8+18,I4]
    dGammadU2134 = dGammadU[I9+18,I4]
    dGammadU2215 = dGammadU[I2,I5]
    dGammadU1215 = dGammadU[I6,I5]
    dGammadU3215 = dGammadU[I7,I5]
    dGammadU2225 = dGammadU[I2+9,I5]
    dGammadU1225 = dGammadU[I6+9,I5]
    dGammadU3225 = dGammadU[I7+9,I5]
    dGammadU2235 = dGammadU[I2+18,I5]
    dGammadU1235 = dGammadU[I6+18,I5]
    dGammadU3235 = dGammadU[I7+18,I5]
    dGammadU3316 = dGammadU[I3,I6]
    dGammadU2316 = dGammadU[I4,I6]
    dGammadU1316 = dGammadU[I5,I6]
    dGammadU3326 = dGammadU[I3+9,I6]
    dGammadU2326 = dGammadU[I4+9,I6]
    dGammadU1326 = dGammadU[I5+9,I6]
    dGammadU3336 = dGammadU[I3+18,I6]
    dGammadU2336 = dGammadU[I4+18,I6]
    dGammadU1336 = dGammadU[I5+18,I6]
    dGammadU3317 = dGammadU[I3,I7]
    dGammadU2317 = dGammadU[I4,I7]
    dGammadU1317 = dGammadU[I5,I7]
    dGammadU3327 = dGammadU[I3+9,I7]
    dGammadU2327 = dGammadU[I4+9,I7]
    dGammadU1327 = dGammadU[I5+9,I7]
    dGammadU3337 = dGammadU[I3+18,I7]
    dGammadU2337 = dGammadU[I4+18,I7]
    dGammadU1337 = dGammadU[I5+18,I7]
    dGammadU3318 = dGammadU[I3,I8]
    dGammadU2318 = dGammadU[I4,I8]
    dGammadU1318 = dGammadU[I5,I8]
    dGammadU3328 = dGammadU[I3+9,I8]
    dGammadU2328 = dGammadU[I4+9,I8]
    dGammadU1328 = dGammadU[I5+9,I8]
    dGammadU3338 = dGammadU[I3+18,I8]
    dGammadU2338 = dGammadU[I4+18,I8]
    dGammadU1338 = dGammadU[I5+18,I8]
    dGammadU2219 = dGammadU[I2,I9]
    dGammadU1219 = dGammadU[I6,I9]
    dGammadU3219 = dGammadU[I7,I9]
    dGammadU2229 = dGammadU[I2+9,I9]
    dGammadU1229 = dGammadU[I6+9,I9]
    dGammadU3229 = dGammadU[I7+9,I9]
    dGammadU2239 = dGammadU[I2+18,I9]
    dGammadU1239 = dGammadU[I6+18,I9]
    dGammadU3239 = dGammadU[I7+18,I9]
    dGammadU22110 = dGammadU[I2,I10]
    dGammadU12110 = dGammadU[I6,I10]
    dGammadU32110 = dGammadU[I7,I10]
    dGammadU22210 = dGammadU[I2+9,I10]
    dGammadU12210 = dGammadU[I6+9,I10]
    dGammadU32210 = dGammadU[I7+9,I10]
    dGammadU22310 = dGammadU[I2+18,I10]
    dGammadU12310 = dGammadU[I6+18,I10]
    dGammadU32310 = dGammadU[I7+18,I10]
    dGammadU11111 = dGammadU[I1,I11]
    dGammadU31111 = dGammadU[I8,I11]
    dGammadU21111 = dGammadU[I9,I11]
    dGammadU11211 = dGammadU[I1+9,I11]
    dGammadU31211 = dGammadU[I8+9,I11]
    dGammadU21211 = dGammadU[I9+9,I11]
    dGammadU11311 = dGammadU[I1+18,I11]
    dGammadU31311 = dGammadU[I8+18,I11]
    dGammadU21311 = dGammadU[I9+18,I11]
    dGammadU11112 = dGammadU[I1,I12]
    dGammadU31112 = dGammadU[I8,I12]
    dGammadU21112 = dGammadU[I9,I12]
    dGammadU11212 = dGammadU[I1+9,I12]
    dGammadU31212 = dGammadU[I8+9,I12]
    dGammadU21212 = dGammadU[I9+9,I12]
    dGammadU11312 = dGammadU[I1+18,I12]
    dGammadU31312 = dGammadU[I8+18,I12]
    dGammadU21312 = dGammadU[I9+18,I12]

    #Extract components of dMdC
    dMdC11111 = dMdC[I1,I1]
    dMdC22111 = dMdC[I2,I1]
    dMdC33111 = dMdC[I3,I1]
    dMdC23111 = dMdC[I4,I1]
    dMdC13111 = dMdC[I5,I1]
    dMdC12111 = dMdC[I6,I1]
    dMdC32111 = dMdC[I7,I1]
    dMdC31111 = dMdC[I8,I1]
    dMdC21111 = dMdC[I9,I1]
    dMdC11121 = dMdC[I1+9,I1]
    dMdC22121 = dMdC[I2+9,I1]
    dMdC33121 = dMdC[I3+9,I1]
    dMdC23121 = dMdC[I4+9,I1]
    dMdC13121 = dMdC[I5+9,I1]
    dMdC12121 = dMdC[I6+9,I1]
    dMdC32121 = dMdC[I7+9,I1]
    dMdC31121 = dMdC[I8+9,I1]
    dMdC21121 = dMdC[I9+9,I1]
    dMdC11131 = dMdC[I1+18,I1]
    dMdC22131 = dMdC[I2+18,I1]
    dMdC33131 = dMdC[I3+18,I1]
    dMdC23131 = dMdC[I4+18,I1]
    dMdC13131 = dMdC[I5+18,I1]
    dMdC12131 = dMdC[I6+18,I1]
    dMdC32131 = dMdC[I7+18,I1]
    dMdC31131 = dMdC[I8+18,I1]
    dMdC21131 = dMdC[I9+18,I1]
    dMdC11211 = dMdC[I1+27,I1]
    dMdC22211 = dMdC[I2+27,I1]
    dMdC33211 = dMdC[I3+27,I1]
    dMdC23211 = dMdC[I4+27,I1]
    dMdC13211 = dMdC[I5+27,I1]
    dMdC12211 = dMdC[I6+27,I1]
    dMdC32211 = dMdC[I7+27,I1]
    dMdC31211 = dMdC[I8+27,I1]
    dMdC21211 = dMdC[I9+27,I1]
    dMdC11221 = dMdC[I1+36,I1]
    dMdC22221 = dMdC[I2+36,I1]
    dMdC33221 = dMdC[I3+36,I1]
    dMdC23221 = dMdC[I4+36,I1]
    dMdC13221 = dMdC[I5+36,I1]
    dMdC12221 = dMdC[I6+36,I1]
    dMdC32221 = dMdC[I7+36,I1]
    dMdC31221 = dMdC[I8+36,I1]
    dMdC21221 = dMdC[I9+36,I1]
    dMdC11231 = dMdC[I1+45,I1]
    dMdC22231 = dMdC[I2+45,I1]
    dMdC33231 = dMdC[I3+45,I1]
    dMdC23231 = dMdC[I4+45,I1]
    dMdC13231 = dMdC[I5+45,I1]
    dMdC12231 = dMdC[I6+45,I1]
    dMdC32231 = dMdC[I7+45,I1]
    dMdC31231 = dMdC[I8+45,I1]
    dMdC21231 = dMdC[I9+45,I1]
    dMdC11311 = dMdC[I1+54,I1]
    dMdC22311 = dMdC[I2+54,I1]
    dMdC33311 = dMdC[I3+54,I1]
    dMdC23311 = dMdC[I4+54,I1]
    dMdC13311 = dMdC[I5+54,I1]
    dMdC12311 = dMdC[I6+54,I1]
    dMdC32311 = dMdC[I7+54,I1]
    dMdC31311 = dMdC[I8+54,I1]
    dMdC21311 = dMdC[I9+54,I1]
    dMdC11321 = dMdC[I1+63,I1]
    dMdC22321 = dMdC[I2+63,I1]
    dMdC33321 = dMdC[I3+63,I1]
    dMdC23321 = dMdC[I4+63,I1]
    dMdC13321 = dMdC[I5+63,I1]
    dMdC12321 = dMdC[I6+63,I1]
    dMdC32321 = dMdC[I7+63,I1]
    dMdC31321 = dMdC[I8+63,I1]
    dMdC21321 = dMdC[I9+63,I1]
    dMdC11331 = dMdC[I1+72,I1]
    dMdC22331 = dMdC[I2+72,I1]
    dMdC33331 = dMdC[I3+72,I1]
    dMdC23331 = dMdC[I4+72,I1]
    dMdC13331 = dMdC[I5+72,I1]
    dMdC12331 = dMdC[I6+72,I1]
    dMdC32331 = dMdC[I7+72,I1]
    dMdC31331 = dMdC[I8+72,I1]
    dMdC21331 = dMdC[I9+72,I1]
    dMdC11112 = dMdC[I1,I2]
    dMdC22112 = dMdC[I2,I2]
    dMdC33112 = dMdC[I3,I2]
    dMdC23112 = dMdC[I4,I2]
    dMdC13112 = dMdC[I5,I2]
    dMdC12112 = dMdC[I6,I2]
    dMdC32112 = dMdC[I7,I2]
    dMdC31112 = dMdC[I8,I2]
    dMdC21112 = dMdC[I9,I2]
    dMdC11122 = dMdC[I1+9,I2]
    dMdC22122 = dMdC[I2+9,I2]
    dMdC33122 = dMdC[I3+9,I2]
    dMdC23122 = dMdC[I4+9,I2]
    dMdC13122 = dMdC[I5+9,I2]
    dMdC12122 = dMdC[I6+9,I2]
    dMdC32122 = dMdC[I7+9,I2]
    dMdC31122 = dMdC[I8+9,I2]
    dMdC21122 = dMdC[I9+9,I2]
    dMdC11132 = dMdC[I1+18,I2]
    dMdC22132 = dMdC[I2+18,I2]
    dMdC33132 = dMdC[I3+18,I2]
    dMdC23132 = dMdC[I4+18,I2]
    dMdC13132 = dMdC[I5+18,I2]
    dMdC12132 = dMdC[I6+18,I2]
    dMdC32132 = dMdC[I7+18,I2]
    dMdC31132 = dMdC[I8+18,I2]
    dMdC21132 = dMdC[I9+18,I2]
    dMdC11212 = dMdC[I1+27,I2]
    dMdC22212 = dMdC[I2+27,I2]
    dMdC33212 = dMdC[I3+27,I2]
    dMdC23212 = dMdC[I4+27,I2]
    dMdC13212 = dMdC[I5+27,I2]
    dMdC12212 = dMdC[I6+27,I2]
    dMdC32212 = dMdC[I7+27,I2]
    dMdC31212 = dMdC[I8+27,I2]
    dMdC21212 = dMdC[I9+27,I2]
    dMdC11222 = dMdC[I1+36,I2]
    dMdC22222 = dMdC[I2+36,I2]
    dMdC33222 = dMdC[I3+36,I2]
    dMdC23222 = dMdC[I4+36,I2]
    dMdC13222 = dMdC[I5+36,I2]
    dMdC12222 = dMdC[I6+36,I2]
    dMdC32222 = dMdC[I7+36,I2]
    dMdC31222 = dMdC[I8+36,I2]
    dMdC21222 = dMdC[I9+36,I2]
    dMdC11232 = dMdC[I1+45,I2]
    dMdC22232 = dMdC[I2+45,I2]
    dMdC33232 = dMdC[I3+45,I2]
    dMdC23232 = dMdC[I4+45,I2]
    dMdC13232 = dMdC[I5+45,I2]
    dMdC12232 = dMdC[I6+45,I2]
    dMdC32232 = dMdC[I7+45,I2]
    dMdC31232 = dMdC[I8+45,I2]
    dMdC21232 = dMdC[I9+45,I2]
    dMdC11312 = dMdC[I1+54,I2]
    dMdC22312 = dMdC[I2+54,I2]
    dMdC33312 = dMdC[I3+54,I2]
    dMdC23312 = dMdC[I4+54,I2]
    dMdC13312 = dMdC[I5+54,I2]
    dMdC12312 = dMdC[I6+54,I2]
    dMdC32312 = dMdC[I7+54,I2]
    dMdC31312 = dMdC[I8+54,I2]
    dMdC21312 = dMdC[I9+54,I2]
    dMdC11322 = dMdC[I1+63,I2]
    dMdC22322 = dMdC[I2+63,I2]
    dMdC33322 = dMdC[I3+63,I2]
    dMdC23322 = dMdC[I4+63,I2]
    dMdC13322 = dMdC[I5+63,I2]
    dMdC12322 = dMdC[I6+63,I2]
    dMdC32322 = dMdC[I7+63,I2]
    dMdC31322 = dMdC[I8+63,I2]
    dMdC21322 = dMdC[I9+63,I2]
    dMdC11332 = dMdC[I1+72,I2]
    dMdC22332 = dMdC[I2+72,I2]
    dMdC33332 = dMdC[I3+72,I2]
    dMdC23332 = dMdC[I4+72,I2]
    dMdC13332 = dMdC[I5+72,I2]
    dMdC12332 = dMdC[I6+72,I2]
    dMdC32332 = dMdC[I7+72,I2]
    dMdC31332 = dMdC[I8+72,I2]
    dMdC21332 = dMdC[I9+72,I2]
    dMdC11113 = dMdC[I1,I3]
    dMdC22113 = dMdC[I2,I3]
    dMdC33113 = dMdC[I3,I3]
    dMdC23113 = dMdC[I4,I3]
    dMdC13113 = dMdC[I5,I3]
    dMdC12113 = dMdC[I6,I3]
    dMdC32113 = dMdC[I7,I3]
    dMdC31113 = dMdC[I8,I3]
    dMdC21113 = dMdC[I9,I3]
    dMdC11123 = dMdC[I1+9,I3]
    dMdC22123 = dMdC[I2+9,I3]
    dMdC33123 = dMdC[I3+9,I3]
    dMdC23123 = dMdC[I4+9,I3]
    dMdC13123 = dMdC[I5+9,I3]
    dMdC12123 = dMdC[I6+9,I3]
    dMdC32123 = dMdC[I7+9,I3]
    dMdC31123 = dMdC[I8+9,I3]
    dMdC21123 = dMdC[I9+9,I3]
    dMdC11133 = dMdC[I1+18,I3]
    dMdC22133 = dMdC[I2+18,I3]
    dMdC33133 = dMdC[I3+18,I3]
    dMdC23133 = dMdC[I4+18,I3]
    dMdC13133 = dMdC[I5+18,I3]
    dMdC12133 = dMdC[I6+18,I3]
    dMdC32133 = dMdC[I7+18,I3]
    dMdC31133 = dMdC[I8+18,I3]
    dMdC21133 = dMdC[I9+18,I3]
    dMdC11213 = dMdC[I1+27,I3]
    dMdC22213 = dMdC[I2+27,I3]
    dMdC33213 = dMdC[I3+27,I3]
    dMdC23213 = dMdC[I4+27,I3]
    dMdC13213 = dMdC[I5+27,I3]
    dMdC12213 = dMdC[I6+27,I3]
    dMdC32213 = dMdC[I7+27,I3]
    dMdC31213 = dMdC[I8+27,I3]
    dMdC21213 = dMdC[I9+27,I3]
    dMdC11223 = dMdC[I1+36,I3]
    dMdC22223 = dMdC[I2+36,I3]
    dMdC33223 = dMdC[I3+36,I3]
    dMdC23223 = dMdC[I4+36,I3]
    dMdC13223 = dMdC[I5+36,I3]
    dMdC12223 = dMdC[I6+36,I3]
    dMdC32223 = dMdC[I7+36,I3]
    dMdC31223 = dMdC[I8+36,I3]
    dMdC21223 = dMdC[I9+36,I3]
    dMdC11233 = dMdC[I1+45,I3]
    dMdC22233 = dMdC[I2+45,I3]
    dMdC33233 = dMdC[I3+45,I3]
    dMdC23233 = dMdC[I4+45,I3]
    dMdC13233 = dMdC[I5+45,I3]
    dMdC12233 = dMdC[I6+45,I3]
    dMdC32233 = dMdC[I7+45,I3]
    dMdC31233 = dMdC[I8+45,I3]
    dMdC21233 = dMdC[I9+45,I3]
    dMdC11313 = dMdC[I1+54,I3]
    dMdC22313 = dMdC[I2+54,I3]
    dMdC33313 = dMdC[I3+54,I3]
    dMdC23313 = dMdC[I4+54,I3]
    dMdC13313 = dMdC[I5+54,I3]
    dMdC12313 = dMdC[I6+54,I3]
    dMdC32313 = dMdC[I7+54,I3]
    dMdC31313 = dMdC[I8+54,I3]
    dMdC21313 = dMdC[I9+54,I3]
    dMdC11323 = dMdC[I1+63,I3]
    dMdC22323 = dMdC[I2+63,I3]
    dMdC33323 = dMdC[I3+63,I3]
    dMdC23323 = dMdC[I4+63,I3]
    dMdC13323 = dMdC[I5+63,I3]
    dMdC12323 = dMdC[I6+63,I3]
    dMdC32323 = dMdC[I7+63,I3]
    dMdC31323 = dMdC[I8+63,I3]
    dMdC21323 = dMdC[I9+63,I3]
    dMdC11333 = dMdC[I1+72,I3]
    dMdC22333 = dMdC[I2+72,I3]
    dMdC33333 = dMdC[I3+72,I3]
    dMdC23333 = dMdC[I4+72,I3]
    dMdC13333 = dMdC[I5+72,I3]
    dMdC12333 = dMdC[I6+72,I3]
    dMdC32333 = dMdC[I7+72,I3]
    dMdC31333 = dMdC[I8+72,I3]
    dMdC21333 = dMdC[I9+72,I3]

    #Extract components of dMdPsi
    dMdPsi11111 = dMdPsi[I1,I1]
    dMdPsi22111 = dMdPsi[I2,I1]
    dMdPsi33111 = dMdPsi[I3,I1]
    dMdPsi23111 = dMdPsi[I4,I1]
    dMdPsi13111 = dMdPsi[I5,I1]
    dMdPsi12111 = dMdPsi[I6,I1]
    dMdPsi32111 = dMdPsi[I7,I1]
    dMdPsi31111 = dMdPsi[I8,I1]
    dMdPsi21111 = dMdPsi[I9,I1]
    dMdPsi11121 = dMdPsi[I1+9,I1]
    dMdPsi22121 = dMdPsi[I2+9,I1]
    dMdPsi33121 = dMdPsi[I3+9,I1]
    dMdPsi23121 = dMdPsi[I4+9,I1]
    dMdPsi13121 = dMdPsi[I5+9,I1]
    dMdPsi12121 = dMdPsi[I6+9,I1]
    dMdPsi32121 = dMdPsi[I7+9,I1]
    dMdPsi31121 = dMdPsi[I8+9,I1]
    dMdPsi21121 = dMdPsi[I9+9,I1]
    dMdPsi11131 = dMdPsi[I1+18,I1]
    dMdPsi22131 = dMdPsi[I2+18,I1]
    dMdPsi33131 = dMdPsi[I3+18,I1]
    dMdPsi23131 = dMdPsi[I4+18,I1]
    dMdPsi13131 = dMdPsi[I5+18,I1]
    dMdPsi12131 = dMdPsi[I6+18,I1]
    dMdPsi32131 = dMdPsi[I7+18,I1]
    dMdPsi31131 = dMdPsi[I8+18,I1]
    dMdPsi21131 = dMdPsi[I9+18,I1]
    dMdPsi11211 = dMdPsi[I1+27,I1]
    dMdPsi22211 = dMdPsi[I2+27,I1]
    dMdPsi33211 = dMdPsi[I3+27,I1]
    dMdPsi23211 = dMdPsi[I4+27,I1]
    dMdPsi13211 = dMdPsi[I5+27,I1]
    dMdPsi12211 = dMdPsi[I6+27,I1]
    dMdPsi32211 = dMdPsi[I7+27,I1]
    dMdPsi31211 = dMdPsi[I8+27,I1]
    dMdPsi21211 = dMdPsi[I9+27,I1]
    dMdPsi11221 = dMdPsi[I1+36,I1]
    dMdPsi22221 = dMdPsi[I2+36,I1]
    dMdPsi33221 = dMdPsi[I3+36,I1]
    dMdPsi23221 = dMdPsi[I4+36,I1]
    dMdPsi13221 = dMdPsi[I5+36,I1]
    dMdPsi12221 = dMdPsi[I6+36,I1]
    dMdPsi32221 = dMdPsi[I7+36,I1]
    dMdPsi31221 = dMdPsi[I8+36,I1]
    dMdPsi21221 = dMdPsi[I9+36,I1]
    dMdPsi11231 = dMdPsi[I1+45,I1]
    dMdPsi22231 = dMdPsi[I2+45,I1]
    dMdPsi33231 = dMdPsi[I3+45,I1]
    dMdPsi23231 = dMdPsi[I4+45,I1]
    dMdPsi13231 = dMdPsi[I5+45,I1]
    dMdPsi12231 = dMdPsi[I6+45,I1]
    dMdPsi32231 = dMdPsi[I7+45,I1]
    dMdPsi31231 = dMdPsi[I8+45,I1]
    dMdPsi21231 = dMdPsi[I9+45,I1]
    dMdPsi11311 = dMdPsi[I1+54,I1]
    dMdPsi22311 = dMdPsi[I2+54,I1]
    dMdPsi33311 = dMdPsi[I3+54,I1]
    dMdPsi23311 = dMdPsi[I4+54,I1]
    dMdPsi13311 = dMdPsi[I5+54,I1]
    dMdPsi12311 = dMdPsi[I6+54,I1]
    dMdPsi32311 = dMdPsi[I7+54,I1]
    dMdPsi31311 = dMdPsi[I8+54,I1]
    dMdPsi21311 = dMdPsi[I9+54,I1]
    dMdPsi11321 = dMdPsi[I1+63,I1]
    dMdPsi22321 = dMdPsi[I2+63,I1]
    dMdPsi33321 = dMdPsi[I3+63,I1]
    dMdPsi23321 = dMdPsi[I4+63,I1]
    dMdPsi13321 = dMdPsi[I5+63,I1]
    dMdPsi12321 = dMdPsi[I6+63,I1]
    dMdPsi32321 = dMdPsi[I7+63,I1]
    dMdPsi31321 = dMdPsi[I8+63,I1]
    dMdPsi21321 = dMdPsi[I9+63,I1]
    dMdPsi11331 = dMdPsi[I1+72,I1]
    dMdPsi22331 = dMdPsi[I2+72,I1]
    dMdPsi33331 = dMdPsi[I3+72,I1]
    dMdPsi23331 = dMdPsi[I4+72,I1]
    dMdPsi13331 = dMdPsi[I5+72,I1]
    dMdPsi12331 = dMdPsi[I6+72,I1]
    dMdPsi32331 = dMdPsi[I7+72,I1]
    dMdPsi31331 = dMdPsi[I8+72,I1]
    dMdPsi21331 = dMdPsi[I9+72,I1]
    dMdPsi11112 = dMdPsi[I1,I2]
    dMdPsi22112 = dMdPsi[I2,I2]
    dMdPsi33112 = dMdPsi[I3,I2]
    dMdPsi23112 = dMdPsi[I4,I2]
    dMdPsi13112 = dMdPsi[I5,I2]
    dMdPsi12112 = dMdPsi[I6,I2]
    dMdPsi32112 = dMdPsi[I7,I2]
    dMdPsi31112 = dMdPsi[I8,I2]
    dMdPsi21112 = dMdPsi[I9,I2]
    dMdPsi11122 = dMdPsi[I1+9,I2]
    dMdPsi22122 = dMdPsi[I2+9,I2]
    dMdPsi33122 = dMdPsi[I3+9,I2]
    dMdPsi23122 = dMdPsi[I4+9,I2]
    dMdPsi13122 = dMdPsi[I5+9,I2]
    dMdPsi12122 = dMdPsi[I6+9,I2]
    dMdPsi32122 = dMdPsi[I7+9,I2]
    dMdPsi31122 = dMdPsi[I8+9,I2]
    dMdPsi21122 = dMdPsi[I9+9,I2]
    dMdPsi11132 = dMdPsi[I1+18,I2]
    dMdPsi22132 = dMdPsi[I2+18,I2]
    dMdPsi33132 = dMdPsi[I3+18,I2]
    dMdPsi23132 = dMdPsi[I4+18,I2]
    dMdPsi13132 = dMdPsi[I5+18,I2]
    dMdPsi12132 = dMdPsi[I6+18,I2]
    dMdPsi32132 = dMdPsi[I7+18,I2]
    dMdPsi31132 = dMdPsi[I8+18,I2]
    dMdPsi21132 = dMdPsi[I9+18,I2]
    dMdPsi11212 = dMdPsi[I1+27,I2]
    dMdPsi22212 = dMdPsi[I2+27,I2]
    dMdPsi33212 = dMdPsi[I3+27,I2]
    dMdPsi23212 = dMdPsi[I4+27,I2]
    dMdPsi13212 = dMdPsi[I5+27,I2]
    dMdPsi12212 = dMdPsi[I6+27,I2]
    dMdPsi32212 = dMdPsi[I7+27,I2]
    dMdPsi31212 = dMdPsi[I8+27,I2]
    dMdPsi21212 = dMdPsi[I9+27,I2]
    dMdPsi11222 = dMdPsi[I1+36,I2]
    dMdPsi22222 = dMdPsi[I2+36,I2]
    dMdPsi33222 = dMdPsi[I3+36,I2]
    dMdPsi23222 = dMdPsi[I4+36,I2]
    dMdPsi13222 = dMdPsi[I5+36,I2]
    dMdPsi12222 = dMdPsi[I6+36,I2]
    dMdPsi32222 = dMdPsi[I7+36,I2]
    dMdPsi31222 = dMdPsi[I8+36,I2]
    dMdPsi21222 = dMdPsi[I9+36,I2]
    dMdPsi11232 = dMdPsi[I1+45,I2]
    dMdPsi22232 = dMdPsi[I2+45,I2]
    dMdPsi33232 = dMdPsi[I3+45,I2]
    dMdPsi23232 = dMdPsi[I4+45,I2]
    dMdPsi13232 = dMdPsi[I5+45,I2]
    dMdPsi12232 = dMdPsi[I6+45,I2]
    dMdPsi32232 = dMdPsi[I7+45,I2]
    dMdPsi31232 = dMdPsi[I8+45,I2]
    dMdPsi21232 = dMdPsi[I9+45,I2]
    dMdPsi11312 = dMdPsi[I1+54,I2]
    dMdPsi22312 = dMdPsi[I2+54,I2]
    dMdPsi33312 = dMdPsi[I3+54,I2]
    dMdPsi23312 = dMdPsi[I4+54,I2]
    dMdPsi13312 = dMdPsi[I5+54,I2]
    dMdPsi12312 = dMdPsi[I6+54,I2]
    dMdPsi32312 = dMdPsi[I7+54,I2]
    dMdPsi31312 = dMdPsi[I8+54,I2]
    dMdPsi21312 = dMdPsi[I9+54,I2]
    dMdPsi11322 = dMdPsi[I1+63,I2]
    dMdPsi22322 = dMdPsi[I2+63,I2]
    dMdPsi33322 = dMdPsi[I3+63,I2]
    dMdPsi23322 = dMdPsi[I4+63,I2]
    dMdPsi13322 = dMdPsi[I5+63,I2]
    dMdPsi12322 = dMdPsi[I6+63,I2]
    dMdPsi32322 = dMdPsi[I7+63,I2]
    dMdPsi31322 = dMdPsi[I8+63,I2]
    dMdPsi21322 = dMdPsi[I9+63,I2]
    dMdPsi11332 = dMdPsi[I1+72,I2]
    dMdPsi22332 = dMdPsi[I2+72,I2]
    dMdPsi33332 = dMdPsi[I3+72,I2]
    dMdPsi23332 = dMdPsi[I4+72,I2]
    dMdPsi13332 = dMdPsi[I5+72,I2]
    dMdPsi12332 = dMdPsi[I6+72,I2]
    dMdPsi32332 = dMdPsi[I7+72,I2]
    dMdPsi31332 = dMdPsi[I8+72,I2]
    dMdPsi21332 = dMdPsi[I9+72,I2]
    dMdPsi11113 = dMdPsi[I1,I3]
    dMdPsi22113 = dMdPsi[I2,I3]
    dMdPsi33113 = dMdPsi[I3,I3]
    dMdPsi23113 = dMdPsi[I4,I3]
    dMdPsi13113 = dMdPsi[I5,I3]
    dMdPsi12113 = dMdPsi[I6,I3]
    dMdPsi32113 = dMdPsi[I7,I3]
    dMdPsi31113 = dMdPsi[I8,I3]
    dMdPsi21113 = dMdPsi[I9,I3]
    dMdPsi11123 = dMdPsi[I1+9,I3]
    dMdPsi22123 = dMdPsi[I2+9,I3]
    dMdPsi33123 = dMdPsi[I3+9,I3]
    dMdPsi23123 = dMdPsi[I4+9,I3]
    dMdPsi13123 = dMdPsi[I5+9,I3]
    dMdPsi12123 = dMdPsi[I6+9,I3]
    dMdPsi32123 = dMdPsi[I7+9,I3]
    dMdPsi31123 = dMdPsi[I8+9,I3]
    dMdPsi21123 = dMdPsi[I9+9,I3]
    dMdPsi11133 = dMdPsi[I1+18,I3]
    dMdPsi22133 = dMdPsi[I2+18,I3]
    dMdPsi33133 = dMdPsi[I3+18,I3]
    dMdPsi23133 = dMdPsi[I4+18,I3]
    dMdPsi13133 = dMdPsi[I5+18,I3]
    dMdPsi12133 = dMdPsi[I6+18,I3]
    dMdPsi32133 = dMdPsi[I7+18,I3]
    dMdPsi31133 = dMdPsi[I8+18,I3]
    dMdPsi21133 = dMdPsi[I9+18,I3]
    dMdPsi11213 = dMdPsi[I1+27,I3]
    dMdPsi22213 = dMdPsi[I2+27,I3]
    dMdPsi33213 = dMdPsi[I3+27,I3]
    dMdPsi23213 = dMdPsi[I4+27,I3]
    dMdPsi13213 = dMdPsi[I5+27,I3]
    dMdPsi12213 = dMdPsi[I6+27,I3]
    dMdPsi32213 = dMdPsi[I7+27,I3]
    dMdPsi31213 = dMdPsi[I8+27,I3]
    dMdPsi21213 = dMdPsi[I9+27,I3]
    dMdPsi11223 = dMdPsi[I1+36,I3]
    dMdPsi22223 = dMdPsi[I2+36,I3]
    dMdPsi33223 = dMdPsi[I3+36,I3]
    dMdPsi23223 = dMdPsi[I4+36,I3]
    dMdPsi13223 = dMdPsi[I5+36,I3]
    dMdPsi12223 = dMdPsi[I6+36,I3]
    dMdPsi32223 = dMdPsi[I7+36,I3]
    dMdPsi31223 = dMdPsi[I8+36,I3]
    dMdPsi21223 = dMdPsi[I9+36,I3]
    dMdPsi11233 = dMdPsi[I1+45,I3]
    dMdPsi22233 = dMdPsi[I2+45,I3]
    dMdPsi33233 = dMdPsi[I3+45,I3]
    dMdPsi23233 = dMdPsi[I4+45,I3]
    dMdPsi13233 = dMdPsi[I5+45,I3]
    dMdPsi12233 = dMdPsi[I6+45,I3]
    dMdPsi32233 = dMdPsi[I7+45,I3]
    dMdPsi31233 = dMdPsi[I8+45,I3]
    dMdPsi21233 = dMdPsi[I9+45,I3]
    dMdPsi11313 = dMdPsi[I1+54,I3]
    dMdPsi22313 = dMdPsi[I2+54,I3]
    dMdPsi33313 = dMdPsi[I3+54,I3]
    dMdPsi23313 = dMdPsi[I4+54,I3]
    dMdPsi13313 = dMdPsi[I5+54,I3]
    dMdPsi12313 = dMdPsi[I6+54,I3]
    dMdPsi32313 = dMdPsi[I7+54,I3]
    dMdPsi31313 = dMdPsi[I8+54,I3]
    dMdPsi21313 = dMdPsi[I9+54,I3]
    dMdPsi11323 = dMdPsi[I1+63,I3]
    dMdPsi22323 = dMdPsi[I2+63,I3]
    dMdPsi33323 = dMdPsi[I3+63,I3]
    dMdPsi23323 = dMdPsi[I4+63,I3]
    dMdPsi13323 = dMdPsi[I5+63,I3]
    dMdPsi12323 = dMdPsi[I6+63,I3]
    dMdPsi32323 = dMdPsi[I7+63,I3]
    dMdPsi31323 = dMdPsi[I8+63,I3]
    dMdPsi21323 = dMdPsi[I9+63,I3]
    dMdPsi11333 = dMdPsi[I1+72,I3]
    dMdPsi22333 = dMdPsi[I2+72,I3]
    dMdPsi33333 = dMdPsi[I3+72,I3]
    dMdPsi23333 = dMdPsi[I4+72,I3]
    dMdPsi13333 = dMdPsi[I5+72,I3]
    dMdPsi12333 = dMdPsi[I6+72,I3]
    dMdPsi32333 = dMdPsi[I7+72,I3]
    dMdPsi31333 = dMdPsi[I8+72,I3]
    dMdPsi21333 = dMdPsi[I9+72,I3]

    #Extract components of dMdGamma
    dMdGamma111111 = dMdGamma[I1,I1]
    dMdGamma221111 = dMdGamma[I2,I1]
    dMdGamma331111 = dMdGamma[I3,I1]
    dMdGamma231111 = dMdGamma[I4,I1]
    dMdGamma131111 = dMdGamma[I5,I1]
    dMdGamma121111 = dMdGamma[I6,I1]
    dMdGamma321111 = dMdGamma[I7,I1]
    dMdGamma311111 = dMdGamma[I8,I1]
    dMdGamma211111 = dMdGamma[I9,I1]
    dMdGamma111121 = dMdGamma[I1+9,I1]
    dMdGamma221121 = dMdGamma[I2+9,I1]
    dMdGamma331121 = dMdGamma[I3+9,I1]
    dMdGamma231121 = dMdGamma[I4+9,I1]
    dMdGamma131121 = dMdGamma[I5+9,I1]
    dMdGamma121121 = dMdGamma[I6+9,I1]
    dMdGamma321121 = dMdGamma[I7+9,I1]
    dMdGamma311121 = dMdGamma[I8+9,I1]
    dMdGamma211121 = dMdGamma[I9+9,I1]
    dMdGamma111131 = dMdGamma[I1+18,I1]
    dMdGamma221131 = dMdGamma[I2+18,I1]
    dMdGamma331131 = dMdGamma[I3+18,I1]
    dMdGamma231131 = dMdGamma[I4+18,I1]
    dMdGamma131131 = dMdGamma[I5+18,I1]
    dMdGamma121131 = dMdGamma[I6+18,I1]
    dMdGamma321131 = dMdGamma[I7+18,I1]
    dMdGamma311131 = dMdGamma[I8+18,I1]
    dMdGamma211131 = dMdGamma[I9+18,I1]
    dMdGamma111211 = dMdGamma[I1+27,I1]
    dMdGamma221211 = dMdGamma[I2+27,I1]
    dMdGamma331211 = dMdGamma[I3+27,I1]
    dMdGamma231211 = dMdGamma[I4+27,I1]
    dMdGamma131211 = dMdGamma[I5+27,I1]
    dMdGamma121211 = dMdGamma[I6+27,I1]
    dMdGamma321211 = dMdGamma[I7+27,I1]
    dMdGamma311211 = dMdGamma[I8+27,I1]
    dMdGamma211211 = dMdGamma[I9+27,I1]
    dMdGamma111221 = dMdGamma[I1+36,I1]
    dMdGamma221221 = dMdGamma[I2+36,I1]
    dMdGamma331221 = dMdGamma[I3+36,I1]
    dMdGamma231221 = dMdGamma[I4+36,I1]
    dMdGamma131221 = dMdGamma[I5+36,I1]
    dMdGamma121221 = dMdGamma[I6+36,I1]
    dMdGamma321221 = dMdGamma[I7+36,I1]
    dMdGamma311221 = dMdGamma[I8+36,I1]
    dMdGamma211221 = dMdGamma[I9+36,I1]
    dMdGamma111231 = dMdGamma[I1+45,I1]
    dMdGamma221231 = dMdGamma[I2+45,I1]
    dMdGamma331231 = dMdGamma[I3+45,I1]
    dMdGamma231231 = dMdGamma[I4+45,I1]
    dMdGamma131231 = dMdGamma[I5+45,I1]
    dMdGamma121231 = dMdGamma[I6+45,I1]
    dMdGamma321231 = dMdGamma[I7+45,I1]
    dMdGamma311231 = dMdGamma[I8+45,I1]
    dMdGamma211231 = dMdGamma[I9+45,I1]
    dMdGamma111311 = dMdGamma[I1+54,I1]
    dMdGamma221311 = dMdGamma[I2+54,I1]
    dMdGamma331311 = dMdGamma[I3+54,I1]
    dMdGamma231311 = dMdGamma[I4+54,I1]
    dMdGamma131311 = dMdGamma[I5+54,I1]
    dMdGamma121311 = dMdGamma[I6+54,I1]
    dMdGamma321311 = dMdGamma[I7+54,I1]
    dMdGamma311311 = dMdGamma[I8+54,I1]
    dMdGamma211311 = dMdGamma[I9+54,I1]
    dMdGamma111321 = dMdGamma[I1+63,I1]
    dMdGamma221321 = dMdGamma[I2+63,I1]
    dMdGamma331321 = dMdGamma[I3+63,I1]
    dMdGamma231321 = dMdGamma[I4+63,I1]
    dMdGamma131321 = dMdGamma[I5+63,I1]
    dMdGamma121321 = dMdGamma[I6+63,I1]
    dMdGamma321321 = dMdGamma[I7+63,I1]
    dMdGamma311321 = dMdGamma[I8+63,I1]
    dMdGamma211321 = dMdGamma[I9+63,I1]
    dMdGamma111331 = dMdGamma[I1+72,I1]
    dMdGamma221331 = dMdGamma[I2+72,I1]
    dMdGamma331331 = dMdGamma[I3+72,I1]
    dMdGamma231331 = dMdGamma[I4+72,I1]
    dMdGamma131331 = dMdGamma[I5+72,I1]
    dMdGamma121331 = dMdGamma[I6+72,I1]
    dMdGamma321331 = dMdGamma[I7+72,I1]
    dMdGamma311331 = dMdGamma[I8+72,I1]
    dMdGamma211331 = dMdGamma[I9+72,I1]
    dMdGamma112111 = dMdGamma[I1+81,I1]
    dMdGamma222111 = dMdGamma[I2+81,I1]
    dMdGamma332111 = dMdGamma[I3+81,I1]
    dMdGamma232111 = dMdGamma[I4+81,I1]
    dMdGamma132111 = dMdGamma[I5+81,I1]
    dMdGamma122111 = dMdGamma[I6+81,I1]
    dMdGamma322111 = dMdGamma[I7+81,I1]
    dMdGamma312111 = dMdGamma[I8+81,I1]
    dMdGamma212111 = dMdGamma[I9+81,I1]
    dMdGamma112121 = dMdGamma[I1+90,I1]
    dMdGamma222121 = dMdGamma[I2+90,I1]
    dMdGamma332121 = dMdGamma[I3+90,I1]
    dMdGamma232121 = dMdGamma[I4+90,I1]
    dMdGamma132121 = dMdGamma[I5+90,I1]
    dMdGamma122121 = dMdGamma[I6+90,I1]
    dMdGamma322121 = dMdGamma[I7+90,I1]
    dMdGamma312121 = dMdGamma[I8+90,I1]
    dMdGamma212121 = dMdGamma[I9+90,I1]
    dMdGamma112131 = dMdGamma[I1+99,I1]
    dMdGamma222131 = dMdGamma[I2+99,I1]
    dMdGamma332131 = dMdGamma[I3+99,I1]
    dMdGamma232131 = dMdGamma[I4+99,I1]
    dMdGamma132131 = dMdGamma[I5+99,I1]
    dMdGamma122131 = dMdGamma[I6+99,I1]
    dMdGamma322131 = dMdGamma[I7+99,I1]
    dMdGamma312131 = dMdGamma[I8+99,I1]
    dMdGamma212131 = dMdGamma[I9+99,I1]
    dMdGamma112211 = dMdGamma[I1+108,I1]
    dMdGamma222211 = dMdGamma[I2+108,I1]
    dMdGamma332211 = dMdGamma[I3+108,I1]
    dMdGamma232211 = dMdGamma[I4+108,I1]
    dMdGamma132211 = dMdGamma[I5+108,I1]
    dMdGamma122211 = dMdGamma[I6+108,I1]
    dMdGamma322211 = dMdGamma[I7+108,I1]
    dMdGamma312211 = dMdGamma[I8+108,I1]
    dMdGamma212211 = dMdGamma[I9+108,I1]
    dMdGamma112221 = dMdGamma[I1+117,I1]
    dMdGamma222221 = dMdGamma[I2+117,I1]
    dMdGamma332221 = dMdGamma[I3+117,I1]
    dMdGamma232221 = dMdGamma[I4+117,I1]
    dMdGamma132221 = dMdGamma[I5+117,I1]
    dMdGamma122221 = dMdGamma[I6+117,I1]
    dMdGamma322221 = dMdGamma[I7+117,I1]
    dMdGamma312221 = dMdGamma[I8+117,I1]
    dMdGamma212221 = dMdGamma[I9+117,I1]
    dMdGamma112231 = dMdGamma[I1+126,I1]
    dMdGamma222231 = dMdGamma[I2+126,I1]
    dMdGamma332231 = dMdGamma[I3+126,I1]
    dMdGamma232231 = dMdGamma[I4+126,I1]
    dMdGamma132231 = dMdGamma[I5+126,I1]
    dMdGamma122231 = dMdGamma[I6+126,I1]
    dMdGamma322231 = dMdGamma[I7+126,I1]
    dMdGamma312231 = dMdGamma[I8+126,I1]
    dMdGamma212231 = dMdGamma[I9+126,I1]
    dMdGamma112311 = dMdGamma[I1+135,I1]
    dMdGamma222311 = dMdGamma[I2+135,I1]
    dMdGamma332311 = dMdGamma[I3+135,I1]
    dMdGamma232311 = dMdGamma[I4+135,I1]
    dMdGamma132311 = dMdGamma[I5+135,I1]
    dMdGamma122311 = dMdGamma[I6+135,I1]
    dMdGamma322311 = dMdGamma[I7+135,I1]
    dMdGamma312311 = dMdGamma[I8+135,I1]
    dMdGamma212311 = dMdGamma[I9+135,I1]
    dMdGamma112321 = dMdGamma[I1+144,I1]
    dMdGamma222321 = dMdGamma[I2+144,I1]
    dMdGamma332321 = dMdGamma[I3+144,I1]
    dMdGamma232321 = dMdGamma[I4+144,I1]
    dMdGamma132321 = dMdGamma[I5+144,I1]
    dMdGamma122321 = dMdGamma[I6+144,I1]
    dMdGamma322321 = dMdGamma[I7+144,I1]
    dMdGamma312321 = dMdGamma[I8+144,I1]
    dMdGamma212321 = dMdGamma[I9+144,I1]
    dMdGamma112331 = dMdGamma[I1+153,I1]
    dMdGamma222331 = dMdGamma[I2+153,I1]
    dMdGamma332331 = dMdGamma[I3+153,I1]
    dMdGamma232331 = dMdGamma[I4+153,I1]
    dMdGamma132331 = dMdGamma[I5+153,I1]
    dMdGamma122331 = dMdGamma[I6+153,I1]
    dMdGamma322331 = dMdGamma[I7+153,I1]
    dMdGamma312331 = dMdGamma[I8+153,I1]
    dMdGamma212331 = dMdGamma[I9+153,I1]
    dMdGamma113111 = dMdGamma[I1+162,I1]
    dMdGamma223111 = dMdGamma[I2+162,I1]
    dMdGamma333111 = dMdGamma[I3+162,I1]
    dMdGamma233111 = dMdGamma[I4+162,I1]
    dMdGamma133111 = dMdGamma[I5+162,I1]
    dMdGamma123111 = dMdGamma[I6+162,I1]
    dMdGamma323111 = dMdGamma[I7+162,I1]
    dMdGamma313111 = dMdGamma[I8+162,I1]
    dMdGamma213111 = dMdGamma[I9+162,I1]
    dMdGamma113121 = dMdGamma[I1+171,I1]
    dMdGamma223121 = dMdGamma[I2+171,I1]
    dMdGamma333121 = dMdGamma[I3+171,I1]
    dMdGamma233121 = dMdGamma[I4+171,I1]
    dMdGamma133121 = dMdGamma[I5+171,I1]
    dMdGamma123121 = dMdGamma[I6+171,I1]
    dMdGamma323121 = dMdGamma[I7+171,I1]
    dMdGamma313121 = dMdGamma[I8+171,I1]
    dMdGamma213121 = dMdGamma[I9+171,I1]
    dMdGamma113131 = dMdGamma[I1+180,I1]
    dMdGamma223131 = dMdGamma[I2+180,I1]
    dMdGamma333131 = dMdGamma[I3+180,I1]
    dMdGamma233131 = dMdGamma[I4+180,I1]
    dMdGamma133131 = dMdGamma[I5+180,I1]
    dMdGamma123131 = dMdGamma[I6+180,I1]
    dMdGamma323131 = dMdGamma[I7+180,I1]
    dMdGamma313131 = dMdGamma[I8+180,I1]
    dMdGamma213131 = dMdGamma[I9+180,I1]
    dMdGamma113211 = dMdGamma[I1+189,I1]
    dMdGamma223211 = dMdGamma[I2+189,I1]
    dMdGamma333211 = dMdGamma[I3+189,I1]
    dMdGamma233211 = dMdGamma[I4+189,I1]
    dMdGamma133211 = dMdGamma[I5+189,I1]
    dMdGamma123211 = dMdGamma[I6+189,I1]
    dMdGamma323211 = dMdGamma[I7+189,I1]
    dMdGamma313211 = dMdGamma[I8+189,I1]
    dMdGamma213211 = dMdGamma[I9+189,I1]
    dMdGamma113221 = dMdGamma[I1+198,I1]
    dMdGamma223221 = dMdGamma[I2+198,I1]
    dMdGamma333221 = dMdGamma[I3+198,I1]
    dMdGamma233221 = dMdGamma[I4+198,I1]
    dMdGamma133221 = dMdGamma[I5+198,I1]
    dMdGamma123221 = dMdGamma[I6+198,I1]
    dMdGamma323221 = dMdGamma[I7+198,I1]
    dMdGamma313221 = dMdGamma[I8+198,I1]
    dMdGamma213221 = dMdGamma[I9+198,I1]
    dMdGamma113231 = dMdGamma[I1+207,I1]
    dMdGamma223231 = dMdGamma[I2+207,I1]
    dMdGamma333231 = dMdGamma[I3+207,I1]
    dMdGamma233231 = dMdGamma[I4+207,I1]
    dMdGamma133231 = dMdGamma[I5+207,I1]
    dMdGamma123231 = dMdGamma[I6+207,I1]
    dMdGamma323231 = dMdGamma[I7+207,I1]
    dMdGamma313231 = dMdGamma[I8+207,I1]
    dMdGamma213231 = dMdGamma[I9+207,I1]
    dMdGamma113311 = dMdGamma[I1+216,I1]
    dMdGamma223311 = dMdGamma[I2+216,I1]
    dMdGamma333311 = dMdGamma[I3+216,I1]
    dMdGamma233311 = dMdGamma[I4+216,I1]
    dMdGamma133311 = dMdGamma[I5+216,I1]
    dMdGamma123311 = dMdGamma[I6+216,I1]
    dMdGamma323311 = dMdGamma[I7+216,I1]
    dMdGamma313311 = dMdGamma[I8+216,I1]
    dMdGamma213311 = dMdGamma[I9+216,I1]
    dMdGamma113321 = dMdGamma[I1+225,I1]
    dMdGamma223321 = dMdGamma[I2+225,I1]
    dMdGamma333321 = dMdGamma[I3+225,I1]
    dMdGamma233321 = dMdGamma[I4+225,I1]
    dMdGamma133321 = dMdGamma[I5+225,I1]
    dMdGamma123321 = dMdGamma[I6+225,I1]
    dMdGamma323321 = dMdGamma[I7+225,I1]
    dMdGamma313321 = dMdGamma[I8+225,I1]
    dMdGamma213321 = dMdGamma[I9+225,I1]
    dMdGamma113331 = dMdGamma[I1+234,I1]
    dMdGamma223331 = dMdGamma[I2+234,I1]
    dMdGamma333331 = dMdGamma[I3+234,I1]
    dMdGamma233331 = dMdGamma[I4+234,I1]
    dMdGamma133331 = dMdGamma[I5+234,I1]
    dMdGamma123331 = dMdGamma[I6+234,I1]
    dMdGamma323331 = dMdGamma[I7+234,I1]
    dMdGamma313331 = dMdGamma[I8+234,I1]
    dMdGamma213331 = dMdGamma[I9+234,I1]
    dMdGamma111112 = dMdGamma[I1,I2]
    dMdGamma221112 = dMdGamma[I2,I2]
    dMdGamma331112 = dMdGamma[I3,I2]
    dMdGamma231112 = dMdGamma[I4,I2]
    dMdGamma131112 = dMdGamma[I5,I2]
    dMdGamma121112 = dMdGamma[I6,I2]
    dMdGamma321112 = dMdGamma[I7,I2]
    dMdGamma311112 = dMdGamma[I8,I2]
    dMdGamma211112 = dMdGamma[I9,I2]
    dMdGamma111122 = dMdGamma[I1+9,I2]
    dMdGamma221122 = dMdGamma[I2+9,I2]
    dMdGamma331122 = dMdGamma[I3+9,I2]
    dMdGamma231122 = dMdGamma[I4+9,I2]
    dMdGamma131122 = dMdGamma[I5+9,I2]
    dMdGamma121122 = dMdGamma[I6+9,I2]
    dMdGamma321122 = dMdGamma[I7+9,I2]
    dMdGamma311122 = dMdGamma[I8+9,I2]
    dMdGamma211122 = dMdGamma[I9+9,I2]
    dMdGamma111132 = dMdGamma[I1+18,I2]
    dMdGamma221132 = dMdGamma[I2+18,I2]
    dMdGamma331132 = dMdGamma[I3+18,I2]
    dMdGamma231132 = dMdGamma[I4+18,I2]
    dMdGamma131132 = dMdGamma[I5+18,I2]
    dMdGamma121132 = dMdGamma[I6+18,I2]
    dMdGamma321132 = dMdGamma[I7+18,I2]
    dMdGamma311132 = dMdGamma[I8+18,I2]
    dMdGamma211132 = dMdGamma[I9+18,I2]
    dMdGamma111212 = dMdGamma[I1+27,I2]
    dMdGamma221212 = dMdGamma[I2+27,I2]
    dMdGamma331212 = dMdGamma[I3+27,I2]
    dMdGamma231212 = dMdGamma[I4+27,I2]
    dMdGamma131212 = dMdGamma[I5+27,I2]
    dMdGamma121212 = dMdGamma[I6+27,I2]
    dMdGamma321212 = dMdGamma[I7+27,I2]
    dMdGamma311212 = dMdGamma[I8+27,I2]
    dMdGamma211212 = dMdGamma[I9+27,I2]
    dMdGamma111222 = dMdGamma[I1+36,I2]
    dMdGamma221222 = dMdGamma[I2+36,I2]
    dMdGamma331222 = dMdGamma[I3+36,I2]
    dMdGamma231222 = dMdGamma[I4+36,I2]
    dMdGamma131222 = dMdGamma[I5+36,I2]
    dMdGamma121222 = dMdGamma[I6+36,I2]
    dMdGamma321222 = dMdGamma[I7+36,I2]
    dMdGamma311222 = dMdGamma[I8+36,I2]
    dMdGamma211222 = dMdGamma[I9+36,I2]
    dMdGamma111232 = dMdGamma[I1+45,I2]
    dMdGamma221232 = dMdGamma[I2+45,I2]
    dMdGamma331232 = dMdGamma[I3+45,I2]
    dMdGamma231232 = dMdGamma[I4+45,I2]
    dMdGamma131232 = dMdGamma[I5+45,I2]
    dMdGamma121232 = dMdGamma[I6+45,I2]
    dMdGamma321232 = dMdGamma[I7+45,I2]
    dMdGamma311232 = dMdGamma[I8+45,I2]
    dMdGamma211232 = dMdGamma[I9+45,I2]
    dMdGamma111312 = dMdGamma[I1+54,I2]
    dMdGamma221312 = dMdGamma[I2+54,I2]
    dMdGamma331312 = dMdGamma[I3+54,I2]
    dMdGamma231312 = dMdGamma[I4+54,I2]
    dMdGamma131312 = dMdGamma[I5+54,I2]
    dMdGamma121312 = dMdGamma[I6+54,I2]
    dMdGamma321312 = dMdGamma[I7+54,I2]
    dMdGamma311312 = dMdGamma[I8+54,I2]
    dMdGamma211312 = dMdGamma[I9+54,I2]
    dMdGamma111322 = dMdGamma[I1+63,I2]
    dMdGamma221322 = dMdGamma[I2+63,I2]
    dMdGamma331322 = dMdGamma[I3+63,I2]
    dMdGamma231322 = dMdGamma[I4+63,I2]
    dMdGamma131322 = dMdGamma[I5+63,I2]
    dMdGamma121322 = dMdGamma[I6+63,I2]
    dMdGamma321322 = dMdGamma[I7+63,I2]
    dMdGamma311322 = dMdGamma[I8+63,I2]
    dMdGamma211322 = dMdGamma[I9+63,I2]
    dMdGamma111332 = dMdGamma[I1+72,I2]
    dMdGamma221332 = dMdGamma[I2+72,I2]
    dMdGamma331332 = dMdGamma[I3+72,I2]
    dMdGamma231332 = dMdGamma[I4+72,I2]
    dMdGamma131332 = dMdGamma[I5+72,I2]
    dMdGamma121332 = dMdGamma[I6+72,I2]
    dMdGamma321332 = dMdGamma[I7+72,I2]
    dMdGamma311332 = dMdGamma[I8+72,I2]
    dMdGamma211332 = dMdGamma[I9+72,I2]
    dMdGamma112112 = dMdGamma[I1+81,I2]
    dMdGamma222112 = dMdGamma[I2+81,I2]
    dMdGamma332112 = dMdGamma[I3+81,I2]
    dMdGamma232112 = dMdGamma[I4+81,I2]
    dMdGamma132112 = dMdGamma[I5+81,I2]
    dMdGamma122112 = dMdGamma[I6+81,I2]
    dMdGamma322112 = dMdGamma[I7+81,I2]
    dMdGamma312112 = dMdGamma[I8+81,I2]
    dMdGamma212112 = dMdGamma[I9+81,I2]
    dMdGamma112122 = dMdGamma[I1+90,I2]
    dMdGamma222122 = dMdGamma[I2+90,I2]
    dMdGamma332122 = dMdGamma[I3+90,I2]
    dMdGamma232122 = dMdGamma[I4+90,I2]
    dMdGamma132122 = dMdGamma[I5+90,I2]
    dMdGamma122122 = dMdGamma[I6+90,I2]
    dMdGamma322122 = dMdGamma[I7+90,I2]
    dMdGamma312122 = dMdGamma[I8+90,I2]
    dMdGamma212122 = dMdGamma[I9+90,I2]
    dMdGamma112132 = dMdGamma[I1+99,I2]
    dMdGamma222132 = dMdGamma[I2+99,I2]
    dMdGamma332132 = dMdGamma[I3+99,I2]
    dMdGamma232132 = dMdGamma[I4+99,I2]
    dMdGamma132132 = dMdGamma[I5+99,I2]
    dMdGamma122132 = dMdGamma[I6+99,I2]
    dMdGamma322132 = dMdGamma[I7+99,I2]
    dMdGamma312132 = dMdGamma[I8+99,I2]
    dMdGamma212132 = dMdGamma[I9+99,I2]
    dMdGamma112212 = dMdGamma[I1+108,I2]
    dMdGamma222212 = dMdGamma[I2+108,I2]
    dMdGamma332212 = dMdGamma[I3+108,I2]
    dMdGamma232212 = dMdGamma[I4+108,I2]
    dMdGamma132212 = dMdGamma[I5+108,I2]
    dMdGamma122212 = dMdGamma[I6+108,I2]
    dMdGamma322212 = dMdGamma[I7+108,I2]
    dMdGamma312212 = dMdGamma[I8+108,I2]
    dMdGamma212212 = dMdGamma[I9+108,I2]
    dMdGamma112222 = dMdGamma[I1+117,I2]
    dMdGamma222222 = dMdGamma[I2+117,I2]
    dMdGamma332222 = dMdGamma[I3+117,I2]
    dMdGamma232222 = dMdGamma[I4+117,I2]
    dMdGamma132222 = dMdGamma[I5+117,I2]
    dMdGamma122222 = dMdGamma[I6+117,I2]
    dMdGamma322222 = dMdGamma[I7+117,I2]
    dMdGamma312222 = dMdGamma[I8+117,I2]
    dMdGamma212222 = dMdGamma[I9+117,I2]
    dMdGamma112232 = dMdGamma[I1+126,I2]
    dMdGamma222232 = dMdGamma[I2+126,I2]
    dMdGamma332232 = dMdGamma[I3+126,I2]
    dMdGamma232232 = dMdGamma[I4+126,I2]
    dMdGamma132232 = dMdGamma[I5+126,I2]
    dMdGamma122232 = dMdGamma[I6+126,I2]
    dMdGamma322232 = dMdGamma[I7+126,I2]
    dMdGamma312232 = dMdGamma[I8+126,I2]
    dMdGamma212232 = dMdGamma[I9+126,I2]
    dMdGamma112312 = dMdGamma[I1+135,I2]
    dMdGamma222312 = dMdGamma[I2+135,I2]
    dMdGamma332312 = dMdGamma[I3+135,I2]
    dMdGamma232312 = dMdGamma[I4+135,I2]
    dMdGamma132312 = dMdGamma[I5+135,I2]
    dMdGamma122312 = dMdGamma[I6+135,I2]
    dMdGamma322312 = dMdGamma[I7+135,I2]
    dMdGamma312312 = dMdGamma[I8+135,I2]
    dMdGamma212312 = dMdGamma[I9+135,I2]
    dMdGamma112322 = dMdGamma[I1+144,I2]
    dMdGamma222322 = dMdGamma[I2+144,I2]
    dMdGamma332322 = dMdGamma[I3+144,I2]
    dMdGamma232322 = dMdGamma[I4+144,I2]
    dMdGamma132322 = dMdGamma[I5+144,I2]
    dMdGamma122322 = dMdGamma[I6+144,I2]
    dMdGamma322322 = dMdGamma[I7+144,I2]
    dMdGamma312322 = dMdGamma[I8+144,I2]
    dMdGamma212322 = dMdGamma[I9+144,I2]
    dMdGamma112332 = dMdGamma[I1+153,I2]
    dMdGamma222332 = dMdGamma[I2+153,I2]
    dMdGamma332332 = dMdGamma[I3+153,I2]
    dMdGamma232332 = dMdGamma[I4+153,I2]
    dMdGamma132332 = dMdGamma[I5+153,I2]
    dMdGamma122332 = dMdGamma[I6+153,I2]
    dMdGamma322332 = dMdGamma[I7+153,I2]
    dMdGamma312332 = dMdGamma[I8+153,I2]
    dMdGamma212332 = dMdGamma[I9+153,I2]
    dMdGamma113112 = dMdGamma[I1+162,I2]
    dMdGamma223112 = dMdGamma[I2+162,I2]
    dMdGamma333112 = dMdGamma[I3+162,I2]
    dMdGamma233112 = dMdGamma[I4+162,I2]
    dMdGamma133112 = dMdGamma[I5+162,I2]
    dMdGamma123112 = dMdGamma[I6+162,I2]
    dMdGamma323112 = dMdGamma[I7+162,I2]
    dMdGamma313112 = dMdGamma[I8+162,I2]
    dMdGamma213112 = dMdGamma[I9+162,I2]
    dMdGamma113122 = dMdGamma[I1+171,I2]
    dMdGamma223122 = dMdGamma[I2+171,I2]
    dMdGamma333122 = dMdGamma[I3+171,I2]
    dMdGamma233122 = dMdGamma[I4+171,I2]
    dMdGamma133122 = dMdGamma[I5+171,I2]
    dMdGamma123122 = dMdGamma[I6+171,I2]
    dMdGamma323122 = dMdGamma[I7+171,I2]
    dMdGamma313122 = dMdGamma[I8+171,I2]
    dMdGamma213122 = dMdGamma[I9+171,I2]
    dMdGamma113132 = dMdGamma[I1+180,I2]
    dMdGamma223132 = dMdGamma[I2+180,I2]
    dMdGamma333132 = dMdGamma[I3+180,I2]
    dMdGamma233132 = dMdGamma[I4+180,I2]
    dMdGamma133132 = dMdGamma[I5+180,I2]
    dMdGamma123132 = dMdGamma[I6+180,I2]
    dMdGamma323132 = dMdGamma[I7+180,I2]
    dMdGamma313132 = dMdGamma[I8+180,I2]
    dMdGamma213132 = dMdGamma[I9+180,I2]
    dMdGamma113212 = dMdGamma[I1+189,I2]
    dMdGamma223212 = dMdGamma[I2+189,I2]
    dMdGamma333212 = dMdGamma[I3+189,I2]
    dMdGamma233212 = dMdGamma[I4+189,I2]
    dMdGamma133212 = dMdGamma[I5+189,I2]
    dMdGamma123212 = dMdGamma[I6+189,I2]
    dMdGamma323212 = dMdGamma[I7+189,I2]
    dMdGamma313212 = dMdGamma[I8+189,I2]
    dMdGamma213212 = dMdGamma[I9+189,I2]
    dMdGamma113222 = dMdGamma[I1+198,I2]
    dMdGamma223222 = dMdGamma[I2+198,I2]
    dMdGamma333222 = dMdGamma[I3+198,I2]
    dMdGamma233222 = dMdGamma[I4+198,I2]
    dMdGamma133222 = dMdGamma[I5+198,I2]
    dMdGamma123222 = dMdGamma[I6+198,I2]
    dMdGamma323222 = dMdGamma[I7+198,I2]
    dMdGamma313222 = dMdGamma[I8+198,I2]
    dMdGamma213222 = dMdGamma[I9+198,I2]
    dMdGamma113232 = dMdGamma[I1+207,I2]
    dMdGamma223232 = dMdGamma[I2+207,I2]
    dMdGamma333232 = dMdGamma[I3+207,I2]
    dMdGamma233232 = dMdGamma[I4+207,I2]
    dMdGamma133232 = dMdGamma[I5+207,I2]
    dMdGamma123232 = dMdGamma[I6+207,I2]
    dMdGamma323232 = dMdGamma[I7+207,I2]
    dMdGamma313232 = dMdGamma[I8+207,I2]
    dMdGamma213232 = dMdGamma[I9+207,I2]
    dMdGamma113312 = dMdGamma[I1+216,I2]
    dMdGamma223312 = dMdGamma[I2+216,I2]
    dMdGamma333312 = dMdGamma[I3+216,I2]
    dMdGamma233312 = dMdGamma[I4+216,I2]
    dMdGamma133312 = dMdGamma[I5+216,I2]
    dMdGamma123312 = dMdGamma[I6+216,I2]
    dMdGamma323312 = dMdGamma[I7+216,I2]
    dMdGamma313312 = dMdGamma[I8+216,I2]
    dMdGamma213312 = dMdGamma[I9+216,I2]
    dMdGamma113322 = dMdGamma[I1+225,I2]
    dMdGamma223322 = dMdGamma[I2+225,I2]
    dMdGamma333322 = dMdGamma[I3+225,I2]
    dMdGamma233322 = dMdGamma[I4+225,I2]
    dMdGamma133322 = dMdGamma[I5+225,I2]
    dMdGamma123322 = dMdGamma[I6+225,I2]
    dMdGamma323322 = dMdGamma[I7+225,I2]
    dMdGamma313322 = dMdGamma[I8+225,I2]
    dMdGamma213322 = dMdGamma[I9+225,I2]
    dMdGamma113332 = dMdGamma[I1+234,I2]
    dMdGamma223332 = dMdGamma[I2+234,I2]
    dMdGamma333332 = dMdGamma[I3+234,I2]
    dMdGamma233332 = dMdGamma[I4+234,I2]
    dMdGamma133332 = dMdGamma[I5+234,I2]
    dMdGamma123332 = dMdGamma[I6+234,I2]
    dMdGamma323332 = dMdGamma[I7+234,I2]
    dMdGamma313332 = dMdGamma[I8+234,I2]
    dMdGamma213332 = dMdGamma[I9+234,I2]
    dMdGamma111113 = dMdGamma[I1,I3]
    dMdGamma221113 = dMdGamma[I2,I3]
    dMdGamma331113 = dMdGamma[I3,I3]
    dMdGamma231113 = dMdGamma[I4,I3]
    dMdGamma131113 = dMdGamma[I5,I3]
    dMdGamma121113 = dMdGamma[I6,I3]
    dMdGamma321113 = dMdGamma[I7,I3]
    dMdGamma311113 = dMdGamma[I8,I3]
    dMdGamma211113 = dMdGamma[I9,I3]
    dMdGamma111123 = dMdGamma[I1+9,I3]
    dMdGamma221123 = dMdGamma[I2+9,I3]
    dMdGamma331123 = dMdGamma[I3+9,I3]
    dMdGamma231123 = dMdGamma[I4+9,I3]
    dMdGamma131123 = dMdGamma[I5+9,I3]
    dMdGamma121123 = dMdGamma[I6+9,I3]
    dMdGamma321123 = dMdGamma[I7+9,I3]
    dMdGamma311123 = dMdGamma[I8+9,I3]
    dMdGamma211123 = dMdGamma[I9+9,I3]
    dMdGamma111133 = dMdGamma[I1+18,I3]
    dMdGamma221133 = dMdGamma[I2+18,I3]
    dMdGamma331133 = dMdGamma[I3+18,I3]
    dMdGamma231133 = dMdGamma[I4+18,I3]
    dMdGamma131133 = dMdGamma[I5+18,I3]
    dMdGamma121133 = dMdGamma[I6+18,I3]
    dMdGamma321133 = dMdGamma[I7+18,I3]
    dMdGamma311133 = dMdGamma[I8+18,I3]
    dMdGamma211133 = dMdGamma[I9+18,I3]
    dMdGamma111213 = dMdGamma[I1+27,I3]
    dMdGamma221213 = dMdGamma[I2+27,I3]
    dMdGamma331213 = dMdGamma[I3+27,I3]
    dMdGamma231213 = dMdGamma[I4+27,I3]
    dMdGamma131213 = dMdGamma[I5+27,I3]
    dMdGamma121213 = dMdGamma[I6+27,I3]
    dMdGamma321213 = dMdGamma[I7+27,I3]
    dMdGamma311213 = dMdGamma[I8+27,I3]
    dMdGamma211213 = dMdGamma[I9+27,I3]
    dMdGamma111223 = dMdGamma[I1+36,I3]
    dMdGamma221223 = dMdGamma[I2+36,I3]
    dMdGamma331223 = dMdGamma[I3+36,I3]
    dMdGamma231223 = dMdGamma[I4+36,I3]
    dMdGamma131223 = dMdGamma[I5+36,I3]
    dMdGamma121223 = dMdGamma[I6+36,I3]
    dMdGamma321223 = dMdGamma[I7+36,I3]
    dMdGamma311223 = dMdGamma[I8+36,I3]
    dMdGamma211223 = dMdGamma[I9+36,I3]
    dMdGamma111233 = dMdGamma[I1+45,I3]
    dMdGamma221233 = dMdGamma[I2+45,I3]
    dMdGamma331233 = dMdGamma[I3+45,I3]
    dMdGamma231233 = dMdGamma[I4+45,I3]
    dMdGamma131233 = dMdGamma[I5+45,I3]
    dMdGamma121233 = dMdGamma[I6+45,I3]
    dMdGamma321233 = dMdGamma[I7+45,I3]
    dMdGamma311233 = dMdGamma[I8+45,I3]
    dMdGamma211233 = dMdGamma[I9+45,I3]
    dMdGamma111313 = dMdGamma[I1+54,I3]
    dMdGamma221313 = dMdGamma[I2+54,I3]
    dMdGamma331313 = dMdGamma[I3+54,I3]
    dMdGamma231313 = dMdGamma[I4+54,I3]
    dMdGamma131313 = dMdGamma[I5+54,I3]
    dMdGamma121313 = dMdGamma[I6+54,I3]
    dMdGamma321313 = dMdGamma[I7+54,I3]
    dMdGamma311313 = dMdGamma[I8+54,I3]
    dMdGamma211313 = dMdGamma[I9+54,I3]
    dMdGamma111323 = dMdGamma[I1+63,I3]
    dMdGamma221323 = dMdGamma[I2+63,I3]
    dMdGamma331323 = dMdGamma[I3+63,I3]
    dMdGamma231323 = dMdGamma[I4+63,I3]
    dMdGamma131323 = dMdGamma[I5+63,I3]
    dMdGamma121323 = dMdGamma[I6+63,I3]
    dMdGamma321323 = dMdGamma[I7+63,I3]
    dMdGamma311323 = dMdGamma[I8+63,I3]
    dMdGamma211323 = dMdGamma[I9+63,I3]
    dMdGamma111333 = dMdGamma[I1+72,I3]
    dMdGamma221333 = dMdGamma[I2+72,I3]
    dMdGamma331333 = dMdGamma[I3+72,I3]
    dMdGamma231333 = dMdGamma[I4+72,I3]
    dMdGamma131333 = dMdGamma[I5+72,I3]
    dMdGamma121333 = dMdGamma[I6+72,I3]
    dMdGamma321333 = dMdGamma[I7+72,I3]
    dMdGamma311333 = dMdGamma[I8+72,I3]
    dMdGamma211333 = dMdGamma[I9+72,I3]
    dMdGamma112113 = dMdGamma[I1+81,I3]
    dMdGamma222113 = dMdGamma[I2+81,I3]
    dMdGamma332113 = dMdGamma[I3+81,I3]
    dMdGamma232113 = dMdGamma[I4+81,I3]
    dMdGamma132113 = dMdGamma[I5+81,I3]
    dMdGamma122113 = dMdGamma[I6+81,I3]
    dMdGamma322113 = dMdGamma[I7+81,I3]
    dMdGamma312113 = dMdGamma[I8+81,I3]
    dMdGamma212113 = dMdGamma[I9+81,I3]
    dMdGamma112123 = dMdGamma[I1+90,I3]
    dMdGamma222123 = dMdGamma[I2+90,I3]
    dMdGamma332123 = dMdGamma[I3+90,I3]
    dMdGamma232123 = dMdGamma[I4+90,I3]
    dMdGamma132123 = dMdGamma[I5+90,I3]
    dMdGamma122123 = dMdGamma[I6+90,I3]
    dMdGamma322123 = dMdGamma[I7+90,I3]
    dMdGamma312123 = dMdGamma[I8+90,I3]
    dMdGamma212123 = dMdGamma[I9+90,I3]
    dMdGamma112133 = dMdGamma[I1+99,I3]
    dMdGamma222133 = dMdGamma[I2+99,I3]
    dMdGamma332133 = dMdGamma[I3+99,I3]
    dMdGamma232133 = dMdGamma[I4+99,I3]
    dMdGamma132133 = dMdGamma[I5+99,I3]
    dMdGamma122133 = dMdGamma[I6+99,I3]
    dMdGamma322133 = dMdGamma[I7+99,I3]
    dMdGamma312133 = dMdGamma[I8+99,I3]
    dMdGamma212133 = dMdGamma[I9+99,I3]
    dMdGamma112213 = dMdGamma[I1+108,I3]
    dMdGamma222213 = dMdGamma[I2+108,I3]
    dMdGamma332213 = dMdGamma[I3+108,I3]
    dMdGamma232213 = dMdGamma[I4+108,I3]
    dMdGamma132213 = dMdGamma[I5+108,I3]
    dMdGamma122213 = dMdGamma[I6+108,I3]
    dMdGamma322213 = dMdGamma[I7+108,I3]
    dMdGamma312213 = dMdGamma[I8+108,I3]
    dMdGamma212213 = dMdGamma[I9+108,I3]
    dMdGamma112223 = dMdGamma[I1+117,I3]
    dMdGamma222223 = dMdGamma[I2+117,I3]
    dMdGamma332223 = dMdGamma[I3+117,I3]
    dMdGamma232223 = dMdGamma[I4+117,I3]
    dMdGamma132223 = dMdGamma[I5+117,I3]
    dMdGamma122223 = dMdGamma[I6+117,I3]
    dMdGamma322223 = dMdGamma[I7+117,I3]
    dMdGamma312223 = dMdGamma[I8+117,I3]
    dMdGamma212223 = dMdGamma[I9+117,I3]
    dMdGamma112233 = dMdGamma[I1+126,I3]
    dMdGamma222233 = dMdGamma[I2+126,I3]
    dMdGamma332233 = dMdGamma[I3+126,I3]
    dMdGamma232233 = dMdGamma[I4+126,I3]
    dMdGamma132233 = dMdGamma[I5+126,I3]
    dMdGamma122233 = dMdGamma[I6+126,I3]
    dMdGamma322233 = dMdGamma[I7+126,I3]
    dMdGamma312233 = dMdGamma[I8+126,I3]
    dMdGamma212233 = dMdGamma[I9+126,I3]
    dMdGamma112313 = dMdGamma[I1+135,I3]
    dMdGamma222313 = dMdGamma[I2+135,I3]
    dMdGamma332313 = dMdGamma[I3+135,I3]
    dMdGamma232313 = dMdGamma[I4+135,I3]
    dMdGamma132313 = dMdGamma[I5+135,I3]
    dMdGamma122313 = dMdGamma[I6+135,I3]
    dMdGamma322313 = dMdGamma[I7+135,I3]
    dMdGamma312313 = dMdGamma[I8+135,I3]
    dMdGamma212313 = dMdGamma[I9+135,I3]
    dMdGamma112323 = dMdGamma[I1+144,I3]
    dMdGamma222323 = dMdGamma[I2+144,I3]
    dMdGamma332323 = dMdGamma[I3+144,I3]
    dMdGamma232323 = dMdGamma[I4+144,I3]
    dMdGamma132323 = dMdGamma[I5+144,I3]
    dMdGamma122323 = dMdGamma[I6+144,I3]
    dMdGamma322323 = dMdGamma[I7+144,I3]
    dMdGamma312323 = dMdGamma[I8+144,I3]
    dMdGamma212323 = dMdGamma[I9+144,I3]
    dMdGamma112333 = dMdGamma[I1+153,I3]
    dMdGamma222333 = dMdGamma[I2+153,I3]
    dMdGamma332333 = dMdGamma[I3+153,I3]
    dMdGamma232333 = dMdGamma[I4+153,I3]
    dMdGamma132333 = dMdGamma[I5+153,I3]
    dMdGamma122333 = dMdGamma[I6+153,I3]
    dMdGamma322333 = dMdGamma[I7+153,I3]
    dMdGamma312333 = dMdGamma[I8+153,I3]
    dMdGamma212333 = dMdGamma[I9+153,I3]
    dMdGamma113113 = dMdGamma[I1+162,I3]
    dMdGamma223113 = dMdGamma[I2+162,I3]
    dMdGamma333113 = dMdGamma[I3+162,I3]
    dMdGamma233113 = dMdGamma[I4+162,I3]
    dMdGamma133113 = dMdGamma[I5+162,I3]
    dMdGamma123113 = dMdGamma[I6+162,I3]
    dMdGamma323113 = dMdGamma[I7+162,I3]
    dMdGamma313113 = dMdGamma[I8+162,I3]
    dMdGamma213113 = dMdGamma[I9+162,I3]
    dMdGamma113123 = dMdGamma[I1+171,I3]
    dMdGamma223123 = dMdGamma[I2+171,I3]
    dMdGamma333123 = dMdGamma[I3+171,I3]
    dMdGamma233123 = dMdGamma[I4+171,I3]
    dMdGamma133123 = dMdGamma[I5+171,I3]
    dMdGamma123123 = dMdGamma[I6+171,I3]
    dMdGamma323123 = dMdGamma[I7+171,I3]
    dMdGamma313123 = dMdGamma[I8+171,I3]
    dMdGamma213123 = dMdGamma[I9+171,I3]
    dMdGamma113133 = dMdGamma[I1+180,I3]
    dMdGamma223133 = dMdGamma[I2+180,I3]
    dMdGamma333133 = dMdGamma[I3+180,I3]
    dMdGamma233133 = dMdGamma[I4+180,I3]
    dMdGamma133133 = dMdGamma[I5+180,I3]
    dMdGamma123133 = dMdGamma[I6+180,I3]
    dMdGamma323133 = dMdGamma[I7+180,I3]
    dMdGamma313133 = dMdGamma[I8+180,I3]
    dMdGamma213133 = dMdGamma[I9+180,I3]
    dMdGamma113213 = dMdGamma[I1+189,I3]
    dMdGamma223213 = dMdGamma[I2+189,I3]
    dMdGamma333213 = dMdGamma[I3+189,I3]
    dMdGamma233213 = dMdGamma[I4+189,I3]
    dMdGamma133213 = dMdGamma[I5+189,I3]
    dMdGamma123213 = dMdGamma[I6+189,I3]
    dMdGamma323213 = dMdGamma[I7+189,I3]
    dMdGamma313213 = dMdGamma[I8+189,I3]
    dMdGamma213213 = dMdGamma[I9+189,I3]
    dMdGamma113223 = dMdGamma[I1+198,I3]
    dMdGamma223223 = dMdGamma[I2+198,I3]
    dMdGamma333223 = dMdGamma[I3+198,I3]
    dMdGamma233223 = dMdGamma[I4+198,I3]
    dMdGamma133223 = dMdGamma[I5+198,I3]
    dMdGamma123223 = dMdGamma[I6+198,I3]
    dMdGamma323223 = dMdGamma[I7+198,I3]
    dMdGamma313223 = dMdGamma[I8+198,I3]
    dMdGamma213223 = dMdGamma[I9+198,I3]
    dMdGamma113233 = dMdGamma[I1+207,I3]
    dMdGamma223233 = dMdGamma[I2+207,I3]
    dMdGamma333233 = dMdGamma[I3+207,I3]
    dMdGamma233233 = dMdGamma[I4+207,I3]
    dMdGamma133233 = dMdGamma[I5+207,I3]
    dMdGamma123233 = dMdGamma[I6+207,I3]
    dMdGamma323233 = dMdGamma[I7+207,I3]
    dMdGamma313233 = dMdGamma[I8+207,I3]
    dMdGamma213233 = dMdGamma[I9+207,I3]
    dMdGamma113313 = dMdGamma[I1+216,I3]
    dMdGamma223313 = dMdGamma[I2+216,I3]
    dMdGamma333313 = dMdGamma[I3+216,I3]
    dMdGamma233313 = dMdGamma[I4+216,I3]
    dMdGamma133313 = dMdGamma[I5+216,I3]
    dMdGamma123313 = dMdGamma[I6+216,I3]
    dMdGamma323313 = dMdGamma[I7+216,I3]
    dMdGamma313313 = dMdGamma[I8+216,I3]
    dMdGamma213313 = dMdGamma[I9+216,I3]
    dMdGamma113323 = dMdGamma[I1+225,I3]
    dMdGamma223323 = dMdGamma[I2+225,I3]
    dMdGamma333323 = dMdGamma[I3+225,I3]
    dMdGamma233323 = dMdGamma[I4+225,I3]
    dMdGamma133323 = dMdGamma[I5+225,I3]
    dMdGamma123323 = dMdGamma[I6+225,I3]
    dMdGamma323323 = dMdGamma[I7+225,I3]
    dMdGamma313323 = dMdGamma[I8+225,I3]
    dMdGamma213323 = dMdGamma[I9+225,I3]
    dMdGamma113333 = dMdGamma[I1+234,I3]
    dMdGamma223333 = dMdGamma[I2+234,I3]
    dMdGamma333333 = dMdGamma[I3+234,I3]
    dMdGamma233333 = dMdGamma[I4+234,I3]
    dMdGamma133333 = dMdGamma[I5+234,I3]
    dMdGamma123333 = dMdGamma[I6+234,I3]
    dMdGamma323333 = dMdGamma[I7+234,I3]
    dMdGamma313333 = dMdGamma[I8+234,I3]
    dMdGamma213333 = dMdGamma[I9+234,I3]

    #Compute Tangent
    
    #Column 1
    dMdU[I1,I1] = dCdU111*dMdC11111 + dCdU121*dMdC11112 + dCdU131*dMdC11113 + dCdU211*dMdC11121 + dCdU221*dMdC11122 + dCdU231*dMdC11123 + dCdU311*dMdC11131 + dCdU321*dMdC11132 + dCdU331*dMdC11133 + dGammadU1111*dMdGamma111111 + dGammadU1121*dMdGamma111112 + dGammadU1131*dMdGamma111113 + dGammadU1211*dMdGamma111121 + dGammadU1221*dMdGamma111122 + dGammadU1231*dMdGamma111123 + dGammadU1311*dMdGamma111131 + dGammadU1321*dMdGamma111132 + dGammadU1331*dMdGamma111133 + dGammadU2111*dMdGamma111211 + dGammadU2121*dMdGamma111212 + dGammadU2131*dMdGamma111213 + dGammadU2211*dMdGamma111221 + dGammadU2221*dMdGamma111222 + dGammadU2231*dMdGamma111223 + dGammadU2311*dMdGamma111231 + dGammadU2321*dMdGamma111232 + dGammadU2331*dMdGamma111233 + dGammadU3111*dMdGamma111311 + dGammadU3121*dMdGamma111312 + dGammadU3131*dMdGamma111313 + dGammadU3211*dMdGamma111321 + dGammadU3221*dMdGamma111322 + dGammadU3231*dMdGamma111323 + dGammadU3311*dMdGamma111331 + dGammadU3321*dMdGamma111332 + dGammadU3331*dMdGamma111333 + dMdPsi11111*dPhidU111 + dMdPsi11112*dPhidU121 + dMdPsi11113*dPhidU131 + dMdPsi11121*dPhidU211 + dMdPsi11122*dPhidU221 + dMdPsi11123*dPhidU231 + dMdPsi11131*dPhidU311 + dMdPsi11132*dPhidU321 + dMdPsi11133*dPhidU331
    dMdU[I2,I1] = dCdU111*dMdC22111 + dCdU121*dMdC22112 + dCdU131*dMdC22113 + dCdU211*dMdC22121 + dCdU221*dMdC22122 + dCdU231*dMdC22123 + dCdU311*dMdC22131 + dCdU321*dMdC22132 + dCdU331*dMdC22133 + dGammadU1111*dMdGamma221111 + dGammadU1121*dMdGamma221112 + dGammadU1131*dMdGamma221113 + dGammadU1211*dMdGamma221121 + dGammadU1221*dMdGamma221122 + dGammadU1231*dMdGamma221123 + dGammadU1311*dMdGamma221131 + dGammadU1321*dMdGamma221132 + dGammadU1331*dMdGamma221133 + dGammadU2111*dMdGamma221211 + dGammadU2121*dMdGamma221212 + dGammadU2131*dMdGamma221213 + dGammadU2211*dMdGamma221221 + dGammadU2221*dMdGamma221222 + dGammadU2231*dMdGamma221223 + dGammadU2311*dMdGamma221231 + dGammadU2321*dMdGamma221232 + dGammadU2331*dMdGamma221233 + dGammadU3111*dMdGamma221311 + dGammadU3121*dMdGamma221312 + dGammadU3131*dMdGamma221313 + dGammadU3211*dMdGamma221321 + dGammadU3221*dMdGamma221322 + dGammadU3231*dMdGamma221323 + dGammadU3311*dMdGamma221331 + dGammadU3321*dMdGamma221332 + dGammadU3331*dMdGamma221333 + dMdPsi22111*dPhidU111 + dMdPsi22112*dPhidU121 + dMdPsi22113*dPhidU131 + dMdPsi22121*dPhidU211 + dMdPsi22122*dPhidU221 + dMdPsi22123*dPhidU231 + dMdPsi22131*dPhidU311 + dMdPsi22132*dPhidU321 + dMdPsi22133*dPhidU331
    dMdU[I3,I1] = dCdU111*dMdC33111 + dCdU121*dMdC33112 + dCdU131*dMdC33113 + dCdU211*dMdC33121 + dCdU221*dMdC33122 + dCdU231*dMdC33123 + dCdU311*dMdC33131 + dCdU321*dMdC33132 + dCdU331*dMdC33133 + dGammadU1111*dMdGamma331111 + dGammadU1121*dMdGamma331112 + dGammadU1131*dMdGamma331113 + dGammadU1211*dMdGamma331121 + dGammadU1221*dMdGamma331122 + dGammadU1231*dMdGamma331123 + dGammadU1311*dMdGamma331131 + dGammadU1321*dMdGamma331132 + dGammadU1331*dMdGamma331133 + dGammadU2111*dMdGamma331211 + dGammadU2121*dMdGamma331212 + dGammadU2131*dMdGamma331213 + dGammadU2211*dMdGamma331221 + dGammadU2221*dMdGamma331222 + dGammadU2231*dMdGamma331223 + dGammadU2311*dMdGamma331231 + dGammadU2321*dMdGamma331232 + dGammadU2331*dMdGamma331233 + dGammadU3111*dMdGamma331311 + dGammadU3121*dMdGamma331312 + dGammadU3131*dMdGamma331313 + dGammadU3211*dMdGamma331321 + dGammadU3221*dMdGamma331322 + dGammadU3231*dMdGamma331323 + dGammadU3311*dMdGamma331331 + dGammadU3321*dMdGamma331332 + dGammadU3331*dMdGamma331333 + dMdPsi33111*dPhidU111 + dMdPsi33112*dPhidU121 + dMdPsi33113*dPhidU131 + dMdPsi33121*dPhidU211 + dMdPsi33122*dPhidU221 + dMdPsi33123*dPhidU231 + dMdPsi33131*dPhidU311 + dMdPsi33132*dPhidU321 + dMdPsi33133*dPhidU331
    dMdU[I4,I1] = dCdU111*dMdC23111 + dCdU121*dMdC23112 + dCdU131*dMdC23113 + dCdU211*dMdC23121 + dCdU221*dMdC23122 + dCdU231*dMdC23123 + dCdU311*dMdC23131 + dCdU321*dMdC23132 + dCdU331*dMdC23133 + dGammadU1111*dMdGamma231111 + dGammadU1121*dMdGamma231112 + dGammadU1131*dMdGamma231113 + dGammadU1211*dMdGamma231121 + dGammadU1221*dMdGamma231122 + dGammadU1231*dMdGamma231123 + dGammadU1311*dMdGamma231131 + dGammadU1321*dMdGamma231132 + dGammadU1331*dMdGamma231133 + dGammadU2111*dMdGamma231211 + dGammadU2121*dMdGamma231212 + dGammadU2131*dMdGamma231213 + dGammadU2211*dMdGamma231221 + dGammadU2221*dMdGamma231222 + dGammadU2231*dMdGamma231223 + dGammadU2311*dMdGamma231231 + dGammadU2321*dMdGamma231232 + dGammadU2331*dMdGamma231233 + dGammadU3111*dMdGamma231311 + dGammadU3121*dMdGamma231312 + dGammadU3131*dMdGamma231313 + dGammadU3211*dMdGamma231321 + dGammadU3221*dMdGamma231322 + dGammadU3231*dMdGamma231323 + dGammadU3311*dMdGamma231331 + dGammadU3321*dMdGamma231332 + dGammadU3331*dMdGamma231333 + dMdPsi23111*dPhidU111 + dMdPsi23112*dPhidU121 + dMdPsi23113*dPhidU131 + dMdPsi23121*dPhidU211 + dMdPsi23122*dPhidU221 + dMdPsi23123*dPhidU231 + dMdPsi23131*dPhidU311 + dMdPsi23132*dPhidU321 + dMdPsi23133*dPhidU331
    dMdU[I5,I1] = dCdU111*dMdC13111 + dCdU121*dMdC13112 + dCdU131*dMdC13113 + dCdU211*dMdC13121 + dCdU221*dMdC13122 + dCdU231*dMdC13123 + dCdU311*dMdC13131 + dCdU321*dMdC13132 + dCdU331*dMdC13133 + dGammadU1111*dMdGamma131111 + dGammadU1121*dMdGamma131112 + dGammadU1131*dMdGamma131113 + dGammadU1211*dMdGamma131121 + dGammadU1221*dMdGamma131122 + dGammadU1231*dMdGamma131123 + dGammadU1311*dMdGamma131131 + dGammadU1321*dMdGamma131132 + dGammadU1331*dMdGamma131133 + dGammadU2111*dMdGamma131211 + dGammadU2121*dMdGamma131212 + dGammadU2131*dMdGamma131213 + dGammadU2211*dMdGamma131221 + dGammadU2221*dMdGamma131222 + dGammadU2231*dMdGamma131223 + dGammadU2311*dMdGamma131231 + dGammadU2321*dMdGamma131232 + dGammadU2331*dMdGamma131233 + dGammadU3111*dMdGamma131311 + dGammadU3121*dMdGamma131312 + dGammadU3131*dMdGamma131313 + dGammadU3211*dMdGamma131321 + dGammadU3221*dMdGamma131322 + dGammadU3231*dMdGamma131323 + dGammadU3311*dMdGamma131331 + dGammadU3321*dMdGamma131332 + dGammadU3331*dMdGamma131333 + dMdPsi13111*dPhidU111 + dMdPsi13112*dPhidU121 + dMdPsi13113*dPhidU131 + dMdPsi13121*dPhidU211 + dMdPsi13122*dPhidU221 + dMdPsi13123*dPhidU231 + dMdPsi13131*dPhidU311 + dMdPsi13132*dPhidU321 + dMdPsi13133*dPhidU331
    dMdU[I6,I1] = dCdU111*dMdC12111 + dCdU121*dMdC12112 + dCdU131*dMdC12113 + dCdU211*dMdC12121 + dCdU221*dMdC12122 + dCdU231*dMdC12123 + dCdU311*dMdC12131 + dCdU321*dMdC12132 + dCdU331*dMdC12133 + dGammadU1111*dMdGamma121111 + dGammadU1121*dMdGamma121112 + dGammadU1131*dMdGamma121113 + dGammadU1211*dMdGamma121121 + dGammadU1221*dMdGamma121122 + dGammadU1231*dMdGamma121123 + dGammadU1311*dMdGamma121131 + dGammadU1321*dMdGamma121132 + dGammadU1331*dMdGamma121133 + dGammadU2111*dMdGamma121211 + dGammadU2121*dMdGamma121212 + dGammadU2131*dMdGamma121213 + dGammadU2211*dMdGamma121221 + dGammadU2221*dMdGamma121222 + dGammadU2231*dMdGamma121223 + dGammadU2311*dMdGamma121231 + dGammadU2321*dMdGamma121232 + dGammadU2331*dMdGamma121233 + dGammadU3111*dMdGamma121311 + dGammadU3121*dMdGamma121312 + dGammadU3131*dMdGamma121313 + dGammadU3211*dMdGamma121321 + dGammadU3221*dMdGamma121322 + dGammadU3231*dMdGamma121323 + dGammadU3311*dMdGamma121331 + dGammadU3321*dMdGamma121332 + dGammadU3331*dMdGamma121333 + dMdPsi12111*dPhidU111 + dMdPsi12112*dPhidU121 + dMdPsi12113*dPhidU131 + dMdPsi12121*dPhidU211 + dMdPsi12122*dPhidU221 + dMdPsi12123*dPhidU231 + dMdPsi12131*dPhidU311 + dMdPsi12132*dPhidU321 + dMdPsi12133*dPhidU331
    dMdU[I7,I1] = dCdU111*dMdC32111 + dCdU121*dMdC32112 + dCdU131*dMdC32113 + dCdU211*dMdC32121 + dCdU221*dMdC32122 + dCdU231*dMdC32123 + dCdU311*dMdC32131 + dCdU321*dMdC32132 + dCdU331*dMdC32133 + dGammadU1111*dMdGamma321111 + dGammadU1121*dMdGamma321112 + dGammadU1131*dMdGamma321113 + dGammadU1211*dMdGamma321121 + dGammadU1221*dMdGamma321122 + dGammadU1231*dMdGamma321123 + dGammadU1311*dMdGamma321131 + dGammadU1321*dMdGamma321132 + dGammadU1331*dMdGamma321133 + dGammadU2111*dMdGamma321211 + dGammadU2121*dMdGamma321212 + dGammadU2131*dMdGamma321213 + dGammadU2211*dMdGamma321221 + dGammadU2221*dMdGamma321222 + dGammadU2231*dMdGamma321223 + dGammadU2311*dMdGamma321231 + dGammadU2321*dMdGamma321232 + dGammadU2331*dMdGamma321233 + dGammadU3111*dMdGamma321311 + dGammadU3121*dMdGamma321312 + dGammadU3131*dMdGamma321313 + dGammadU3211*dMdGamma321321 + dGammadU3221*dMdGamma321322 + dGammadU3231*dMdGamma321323 + dGammadU3311*dMdGamma321331 + dGammadU3321*dMdGamma321332 + dGammadU3331*dMdGamma321333 + dMdPsi32111*dPhidU111 + dMdPsi32112*dPhidU121 + dMdPsi32113*dPhidU131 + dMdPsi32121*dPhidU211 + dMdPsi32122*dPhidU221 + dMdPsi32123*dPhidU231 + dMdPsi32131*dPhidU311 + dMdPsi32132*dPhidU321 + dMdPsi32133*dPhidU331
    dMdU[I8,I1] = dCdU111*dMdC31111 + dCdU121*dMdC31112 + dCdU131*dMdC31113 + dCdU211*dMdC31121 + dCdU221*dMdC31122 + dCdU231*dMdC31123 + dCdU311*dMdC31131 + dCdU321*dMdC31132 + dCdU331*dMdC31133 + dGammadU1111*dMdGamma311111 + dGammadU1121*dMdGamma311112 + dGammadU1131*dMdGamma311113 + dGammadU1211*dMdGamma311121 + dGammadU1221*dMdGamma311122 + dGammadU1231*dMdGamma311123 + dGammadU1311*dMdGamma311131 + dGammadU1321*dMdGamma311132 + dGammadU1331*dMdGamma311133 + dGammadU2111*dMdGamma311211 + dGammadU2121*dMdGamma311212 + dGammadU2131*dMdGamma311213 + dGammadU2211*dMdGamma311221 + dGammadU2221*dMdGamma311222 + dGammadU2231*dMdGamma311223 + dGammadU2311*dMdGamma311231 + dGammadU2321*dMdGamma311232 + dGammadU2331*dMdGamma311233 + dGammadU3111*dMdGamma311311 + dGammadU3121*dMdGamma311312 + dGammadU3131*dMdGamma311313 + dGammadU3211*dMdGamma311321 + dGammadU3221*dMdGamma311322 + dGammadU3231*dMdGamma311323 + dGammadU3311*dMdGamma311331 + dGammadU3321*dMdGamma311332 + dGammadU3331*dMdGamma311333 + dMdPsi31111*dPhidU111 + dMdPsi31112*dPhidU121 + dMdPsi31113*dPhidU131 + dMdPsi31121*dPhidU211 + dMdPsi31122*dPhidU221 + dMdPsi31123*dPhidU231 + dMdPsi31131*dPhidU311 + dMdPsi31132*dPhidU321 + dMdPsi31133*dPhidU331
    dMdU[I9,I1] = dCdU111*dMdC21111 + dCdU121*dMdC21112 + dCdU131*dMdC21113 + dCdU211*dMdC21121 + dCdU221*dMdC21122 + dCdU231*dMdC21123 + dCdU311*dMdC21131 + dCdU321*dMdC21132 + dCdU331*dMdC21133 + dGammadU1111*dMdGamma211111 + dGammadU1121*dMdGamma211112 + dGammadU1131*dMdGamma211113 + dGammadU1211*dMdGamma211121 + dGammadU1221*dMdGamma211122 + dGammadU1231*dMdGamma211123 + dGammadU1311*dMdGamma211131 + dGammadU1321*dMdGamma211132 + dGammadU1331*dMdGamma211133 + dGammadU2111*dMdGamma211211 + dGammadU2121*dMdGamma211212 + dGammadU2131*dMdGamma211213 + dGammadU2211*dMdGamma211221 + dGammadU2221*dMdGamma211222 + dGammadU2231*dMdGamma211223 + dGammadU2311*dMdGamma211231 + dGammadU2321*dMdGamma211232 + dGammadU2331*dMdGamma211233 + dGammadU3111*dMdGamma211311 + dGammadU3121*dMdGamma211312 + dGammadU3131*dMdGamma211313 + dGammadU3211*dMdGamma211321 + dGammadU3221*dMdGamma211322 + dGammadU3231*dMdGamma211323 + dGammadU3311*dMdGamma211331 + dGammadU3321*dMdGamma211332 + dGammadU3331*dMdGamma211333 + dMdPsi21111*dPhidU111 + dMdPsi21112*dPhidU121 + dMdPsi21113*dPhidU131 + dMdPsi21121*dPhidU211 + dMdPsi21122*dPhidU221 + dMdPsi21123*dPhidU231 + dMdPsi21131*dPhidU311 + dMdPsi21132*dPhidU321 + dMdPsi21133*dPhidU331
    dMdU[I1+9,I1] = dCdU111*dMdC11211 + dCdU121*dMdC11212 + dCdU131*dMdC11213 + dCdU211*dMdC11221 + dCdU221*dMdC11222 + dCdU231*dMdC11223 + dCdU311*dMdC11231 + dCdU321*dMdC11232 + dCdU331*dMdC11233 + dGammadU1111*dMdGamma112111 + dGammadU1121*dMdGamma112112 + dGammadU1131*dMdGamma112113 + dGammadU1211*dMdGamma112121 + dGammadU1221*dMdGamma112122 + dGammadU1231*dMdGamma112123 + dGammadU1311*dMdGamma112131 + dGammadU1321*dMdGamma112132 + dGammadU1331*dMdGamma112133 + dGammadU2111*dMdGamma112211 + dGammadU2121*dMdGamma112212 + dGammadU2131*dMdGamma112213 + dGammadU2211*dMdGamma112221 + dGammadU2221*dMdGamma112222 + dGammadU2231*dMdGamma112223 + dGammadU2311*dMdGamma112231 + dGammadU2321*dMdGamma112232 + dGammadU2331*dMdGamma112233 + dGammadU3111*dMdGamma112311 + dGammadU3121*dMdGamma112312 + dGammadU3131*dMdGamma112313 + dGammadU3211*dMdGamma112321 + dGammadU3221*dMdGamma112322 + dGammadU3231*dMdGamma112323 + dGammadU3311*dMdGamma112331 + dGammadU3321*dMdGamma112332 + dGammadU3331*dMdGamma112333 + dMdPsi11211*dPhidU111 + dMdPsi11212*dPhidU121 + dMdPsi11213*dPhidU131 + dMdPsi11221*dPhidU211 + dMdPsi11222*dPhidU221 + dMdPsi11223*dPhidU231 + dMdPsi11231*dPhidU311 + dMdPsi11232*dPhidU321 + dMdPsi11233*dPhidU331
    dMdU[I2+9,I1] = dCdU111*dMdC22211 + dCdU121*dMdC22212 + dCdU131*dMdC22213 + dCdU211*dMdC22221 + dCdU221*dMdC22222 + dCdU231*dMdC22223 + dCdU311*dMdC22231 + dCdU321*dMdC22232 + dCdU331*dMdC22233 + dGammadU1111*dMdGamma222111 + dGammadU1121*dMdGamma222112 + dGammadU1131*dMdGamma222113 + dGammadU1211*dMdGamma222121 + dGammadU1221*dMdGamma222122 + dGammadU1231*dMdGamma222123 + dGammadU1311*dMdGamma222131 + dGammadU1321*dMdGamma222132 + dGammadU1331*dMdGamma222133 + dGammadU2111*dMdGamma222211 + dGammadU2121*dMdGamma222212 + dGammadU2131*dMdGamma222213 + dGammadU2211*dMdGamma222221 + dGammadU2221*dMdGamma222222 + dGammadU2231*dMdGamma222223 + dGammadU2311*dMdGamma222231 + dGammadU2321*dMdGamma222232 + dGammadU2331*dMdGamma222233 + dGammadU3111*dMdGamma222311 + dGammadU3121*dMdGamma222312 + dGammadU3131*dMdGamma222313 + dGammadU3211*dMdGamma222321 + dGammadU3221*dMdGamma222322 + dGammadU3231*dMdGamma222323 + dGammadU3311*dMdGamma222331 + dGammadU3321*dMdGamma222332 + dGammadU3331*dMdGamma222333 + dMdPsi22211*dPhidU111 + dMdPsi22212*dPhidU121 + dMdPsi22213*dPhidU131 + dMdPsi22221*dPhidU211 + dMdPsi22222*dPhidU221 + dMdPsi22223*dPhidU231 + dMdPsi22231*dPhidU311 + dMdPsi22232*dPhidU321 + dMdPsi22233*dPhidU331
    dMdU[I3+9,I1] = dCdU111*dMdC33211 + dCdU121*dMdC33212 + dCdU131*dMdC33213 + dCdU211*dMdC33221 + dCdU221*dMdC33222 + dCdU231*dMdC33223 + dCdU311*dMdC33231 + dCdU321*dMdC33232 + dCdU331*dMdC33233 + dGammadU1111*dMdGamma332111 + dGammadU1121*dMdGamma332112 + dGammadU1131*dMdGamma332113 + dGammadU1211*dMdGamma332121 + dGammadU1221*dMdGamma332122 + dGammadU1231*dMdGamma332123 + dGammadU1311*dMdGamma332131 + dGammadU1321*dMdGamma332132 + dGammadU1331*dMdGamma332133 + dGammadU2111*dMdGamma332211 + dGammadU2121*dMdGamma332212 + dGammadU2131*dMdGamma332213 + dGammadU2211*dMdGamma332221 + dGammadU2221*dMdGamma332222 + dGammadU2231*dMdGamma332223 + dGammadU2311*dMdGamma332231 + dGammadU2321*dMdGamma332232 + dGammadU2331*dMdGamma332233 + dGammadU3111*dMdGamma332311 + dGammadU3121*dMdGamma332312 + dGammadU3131*dMdGamma332313 + dGammadU3211*dMdGamma332321 + dGammadU3221*dMdGamma332322 + dGammadU3231*dMdGamma332323 + dGammadU3311*dMdGamma332331 + dGammadU3321*dMdGamma332332 + dGammadU3331*dMdGamma332333 + dMdPsi33211*dPhidU111 + dMdPsi33212*dPhidU121 + dMdPsi33213*dPhidU131 + dMdPsi33221*dPhidU211 + dMdPsi33222*dPhidU221 + dMdPsi33223*dPhidU231 + dMdPsi33231*dPhidU311 + dMdPsi33232*dPhidU321 + dMdPsi33233*dPhidU331
    dMdU[I4+9,I1] = dCdU111*dMdC23211 + dCdU121*dMdC23212 + dCdU131*dMdC23213 + dCdU211*dMdC23221 + dCdU221*dMdC23222 + dCdU231*dMdC23223 + dCdU311*dMdC23231 + dCdU321*dMdC23232 + dCdU331*dMdC23233 + dGammadU1111*dMdGamma232111 + dGammadU1121*dMdGamma232112 + dGammadU1131*dMdGamma232113 + dGammadU1211*dMdGamma232121 + dGammadU1221*dMdGamma232122 + dGammadU1231*dMdGamma232123 + dGammadU1311*dMdGamma232131 + dGammadU1321*dMdGamma232132 + dGammadU1331*dMdGamma232133 + dGammadU2111*dMdGamma232211 + dGammadU2121*dMdGamma232212 + dGammadU2131*dMdGamma232213 + dGammadU2211*dMdGamma232221 + dGammadU2221*dMdGamma232222 + dGammadU2231*dMdGamma232223 + dGammadU2311*dMdGamma232231 + dGammadU2321*dMdGamma232232 + dGammadU2331*dMdGamma232233 + dGammadU3111*dMdGamma232311 + dGammadU3121*dMdGamma232312 + dGammadU3131*dMdGamma232313 + dGammadU3211*dMdGamma232321 + dGammadU3221*dMdGamma232322 + dGammadU3231*dMdGamma232323 + dGammadU3311*dMdGamma232331 + dGammadU3321*dMdGamma232332 + dGammadU3331*dMdGamma232333 + dMdPsi23211*dPhidU111 + dMdPsi23212*dPhidU121 + dMdPsi23213*dPhidU131 + dMdPsi23221*dPhidU211 + dMdPsi23222*dPhidU221 + dMdPsi23223*dPhidU231 + dMdPsi23231*dPhidU311 + dMdPsi23232*dPhidU321 + dMdPsi23233*dPhidU331
    dMdU[I5+9,I1] = dCdU111*dMdC13211 + dCdU121*dMdC13212 + dCdU131*dMdC13213 + dCdU211*dMdC13221 + dCdU221*dMdC13222 + dCdU231*dMdC13223 + dCdU311*dMdC13231 + dCdU321*dMdC13232 + dCdU331*dMdC13233 + dGammadU1111*dMdGamma132111 + dGammadU1121*dMdGamma132112 + dGammadU1131*dMdGamma132113 + dGammadU1211*dMdGamma132121 + dGammadU1221*dMdGamma132122 + dGammadU1231*dMdGamma132123 + dGammadU1311*dMdGamma132131 + dGammadU1321*dMdGamma132132 + dGammadU1331*dMdGamma132133 + dGammadU2111*dMdGamma132211 + dGammadU2121*dMdGamma132212 + dGammadU2131*dMdGamma132213 + dGammadU2211*dMdGamma132221 + dGammadU2221*dMdGamma132222 + dGammadU2231*dMdGamma132223 + dGammadU2311*dMdGamma132231 + dGammadU2321*dMdGamma132232 + dGammadU2331*dMdGamma132233 + dGammadU3111*dMdGamma132311 + dGammadU3121*dMdGamma132312 + dGammadU3131*dMdGamma132313 + dGammadU3211*dMdGamma132321 + dGammadU3221*dMdGamma132322 + dGammadU3231*dMdGamma132323 + dGammadU3311*dMdGamma132331 + dGammadU3321*dMdGamma132332 + dGammadU3331*dMdGamma132333 + dMdPsi13211*dPhidU111 + dMdPsi13212*dPhidU121 + dMdPsi13213*dPhidU131 + dMdPsi13221*dPhidU211 + dMdPsi13222*dPhidU221 + dMdPsi13223*dPhidU231 + dMdPsi13231*dPhidU311 + dMdPsi13232*dPhidU321 + dMdPsi13233*dPhidU331
    dMdU[I6+9,I1] = dCdU111*dMdC12211 + dCdU121*dMdC12212 + dCdU131*dMdC12213 + dCdU211*dMdC12221 + dCdU221*dMdC12222 + dCdU231*dMdC12223 + dCdU311*dMdC12231 + dCdU321*dMdC12232 + dCdU331*dMdC12233 + dGammadU1111*dMdGamma122111 + dGammadU1121*dMdGamma122112 + dGammadU1131*dMdGamma122113 + dGammadU1211*dMdGamma122121 + dGammadU1221*dMdGamma122122 + dGammadU1231*dMdGamma122123 + dGammadU1311*dMdGamma122131 + dGammadU1321*dMdGamma122132 + dGammadU1331*dMdGamma122133 + dGammadU2111*dMdGamma122211 + dGammadU2121*dMdGamma122212 + dGammadU2131*dMdGamma122213 + dGammadU2211*dMdGamma122221 + dGammadU2221*dMdGamma122222 + dGammadU2231*dMdGamma122223 + dGammadU2311*dMdGamma122231 + dGammadU2321*dMdGamma122232 + dGammadU2331*dMdGamma122233 + dGammadU3111*dMdGamma122311 + dGammadU3121*dMdGamma122312 + dGammadU3131*dMdGamma122313 + dGammadU3211*dMdGamma122321 + dGammadU3221*dMdGamma122322 + dGammadU3231*dMdGamma122323 + dGammadU3311*dMdGamma122331 + dGammadU3321*dMdGamma122332 + dGammadU3331*dMdGamma122333 + dMdPsi12211*dPhidU111 + dMdPsi12212*dPhidU121 + dMdPsi12213*dPhidU131 + dMdPsi12221*dPhidU211 + dMdPsi12222*dPhidU221 + dMdPsi12223*dPhidU231 + dMdPsi12231*dPhidU311 + dMdPsi12232*dPhidU321 + dMdPsi12233*dPhidU331
    dMdU[I7+9,I1] = dCdU111*dMdC32211 + dCdU121*dMdC32212 + dCdU131*dMdC32213 + dCdU211*dMdC32221 + dCdU221*dMdC32222 + dCdU231*dMdC32223 + dCdU311*dMdC32231 + dCdU321*dMdC32232 + dCdU331*dMdC32233 + dGammadU1111*dMdGamma322111 + dGammadU1121*dMdGamma322112 + dGammadU1131*dMdGamma322113 + dGammadU1211*dMdGamma322121 + dGammadU1221*dMdGamma322122 + dGammadU1231*dMdGamma322123 + dGammadU1311*dMdGamma322131 + dGammadU1321*dMdGamma322132 + dGammadU1331*dMdGamma322133 + dGammadU2111*dMdGamma322211 + dGammadU2121*dMdGamma322212 + dGammadU2131*dMdGamma322213 + dGammadU2211*dMdGamma322221 + dGammadU2221*dMdGamma322222 + dGammadU2231*dMdGamma322223 + dGammadU2311*dMdGamma322231 + dGammadU2321*dMdGamma322232 + dGammadU2331*dMdGamma322233 + dGammadU3111*dMdGamma322311 + dGammadU3121*dMdGamma322312 + dGammadU3131*dMdGamma322313 + dGammadU3211*dMdGamma322321 + dGammadU3221*dMdGamma322322 + dGammadU3231*dMdGamma322323 + dGammadU3311*dMdGamma322331 + dGammadU3321*dMdGamma322332 + dGammadU3331*dMdGamma322333 + dMdPsi32211*dPhidU111 + dMdPsi32212*dPhidU121 + dMdPsi32213*dPhidU131 + dMdPsi32221*dPhidU211 + dMdPsi32222*dPhidU221 + dMdPsi32223*dPhidU231 + dMdPsi32231*dPhidU311 + dMdPsi32232*dPhidU321 + dMdPsi32233*dPhidU331
    dMdU[I8+9,I1] = dCdU111*dMdC31211 + dCdU121*dMdC31212 + dCdU131*dMdC31213 + dCdU211*dMdC31221 + dCdU221*dMdC31222 + dCdU231*dMdC31223 + dCdU311*dMdC31231 + dCdU321*dMdC31232 + dCdU331*dMdC31233 + dGammadU1111*dMdGamma312111 + dGammadU1121*dMdGamma312112 + dGammadU1131*dMdGamma312113 + dGammadU1211*dMdGamma312121 + dGammadU1221*dMdGamma312122 + dGammadU1231*dMdGamma312123 + dGammadU1311*dMdGamma312131 + dGammadU1321*dMdGamma312132 + dGammadU1331*dMdGamma312133 + dGammadU2111*dMdGamma312211 + dGammadU2121*dMdGamma312212 + dGammadU2131*dMdGamma312213 + dGammadU2211*dMdGamma312221 + dGammadU2221*dMdGamma312222 + dGammadU2231*dMdGamma312223 + dGammadU2311*dMdGamma312231 + dGammadU2321*dMdGamma312232 + dGammadU2331*dMdGamma312233 + dGammadU3111*dMdGamma312311 + dGammadU3121*dMdGamma312312 + dGammadU3131*dMdGamma312313 + dGammadU3211*dMdGamma312321 + dGammadU3221*dMdGamma312322 + dGammadU3231*dMdGamma312323 + dGammadU3311*dMdGamma312331 + dGammadU3321*dMdGamma312332 + dGammadU3331*dMdGamma312333 + dMdPsi31211*dPhidU111 + dMdPsi31212*dPhidU121 + dMdPsi31213*dPhidU131 + dMdPsi31221*dPhidU211 + dMdPsi31222*dPhidU221 + dMdPsi31223*dPhidU231 + dMdPsi31231*dPhidU311 + dMdPsi31232*dPhidU321 + dMdPsi31233*dPhidU331
    dMdU[I9+9,I1] = dCdU111*dMdC21211 + dCdU121*dMdC21212 + dCdU131*dMdC21213 + dCdU211*dMdC21221 + dCdU221*dMdC21222 + dCdU231*dMdC21223 + dCdU311*dMdC21231 + dCdU321*dMdC21232 + dCdU331*dMdC21233 + dGammadU1111*dMdGamma212111 + dGammadU1121*dMdGamma212112 + dGammadU1131*dMdGamma212113 + dGammadU1211*dMdGamma212121 + dGammadU1221*dMdGamma212122 + dGammadU1231*dMdGamma212123 + dGammadU1311*dMdGamma212131 + dGammadU1321*dMdGamma212132 + dGammadU1331*dMdGamma212133 + dGammadU2111*dMdGamma212211 + dGammadU2121*dMdGamma212212 + dGammadU2131*dMdGamma212213 + dGammadU2211*dMdGamma212221 + dGammadU2221*dMdGamma212222 + dGammadU2231*dMdGamma212223 + dGammadU2311*dMdGamma212231 + dGammadU2321*dMdGamma212232 + dGammadU2331*dMdGamma212233 + dGammadU3111*dMdGamma212311 + dGammadU3121*dMdGamma212312 + dGammadU3131*dMdGamma212313 + dGammadU3211*dMdGamma212321 + dGammadU3221*dMdGamma212322 + dGammadU3231*dMdGamma212323 + dGammadU3311*dMdGamma212331 + dGammadU3321*dMdGamma212332 + dGammadU3331*dMdGamma212333 + dMdPsi21211*dPhidU111 + dMdPsi21212*dPhidU121 + dMdPsi21213*dPhidU131 + dMdPsi21221*dPhidU211 + dMdPsi21222*dPhidU221 + dMdPsi21223*dPhidU231 + dMdPsi21231*dPhidU311 + dMdPsi21232*dPhidU321 + dMdPsi21233*dPhidU331
    dMdU[I1+18,I1] = dCdU111*dMdC11311 + dCdU121*dMdC11312 + dCdU131*dMdC11313 + dCdU211*dMdC11321 + dCdU221*dMdC11322 + dCdU231*dMdC11323 + dCdU311*dMdC11331 + dCdU321*dMdC11332 + dCdU331*dMdC11333 + dGammadU1111*dMdGamma113111 + dGammadU1121*dMdGamma113112 + dGammadU1131*dMdGamma113113 + dGammadU1211*dMdGamma113121 + dGammadU1221*dMdGamma113122 + dGammadU1231*dMdGamma113123 + dGammadU1311*dMdGamma113131 + dGammadU1321*dMdGamma113132 + dGammadU1331*dMdGamma113133 + dGammadU2111*dMdGamma113211 + dGammadU2121*dMdGamma113212 + dGammadU2131*dMdGamma113213 + dGammadU2211*dMdGamma113221 + dGammadU2221*dMdGamma113222 + dGammadU2231*dMdGamma113223 + dGammadU2311*dMdGamma113231 + dGammadU2321*dMdGamma113232 + dGammadU2331*dMdGamma113233 + dGammadU3111*dMdGamma113311 + dGammadU3121*dMdGamma113312 + dGammadU3131*dMdGamma113313 + dGammadU3211*dMdGamma113321 + dGammadU3221*dMdGamma113322 + dGammadU3231*dMdGamma113323 + dGammadU3311*dMdGamma113331 + dGammadU3321*dMdGamma113332 + dGammadU3331*dMdGamma113333 + dMdPsi11311*dPhidU111 + dMdPsi11312*dPhidU121 + dMdPsi11313*dPhidU131 + dMdPsi11321*dPhidU211 + dMdPsi11322*dPhidU221 + dMdPsi11323*dPhidU231 + dMdPsi11331*dPhidU311 + dMdPsi11332*dPhidU321 + dMdPsi11333*dPhidU331
    dMdU[I2+18,I1] = dCdU111*dMdC22311 + dCdU121*dMdC22312 + dCdU131*dMdC22313 + dCdU211*dMdC22321 + dCdU221*dMdC22322 + dCdU231*dMdC22323 + dCdU311*dMdC22331 + dCdU321*dMdC22332 + dCdU331*dMdC22333 + dGammadU1111*dMdGamma223111 + dGammadU1121*dMdGamma223112 + dGammadU1131*dMdGamma223113 + dGammadU1211*dMdGamma223121 + dGammadU1221*dMdGamma223122 + dGammadU1231*dMdGamma223123 + dGammadU1311*dMdGamma223131 + dGammadU1321*dMdGamma223132 + dGammadU1331*dMdGamma223133 + dGammadU2111*dMdGamma223211 + dGammadU2121*dMdGamma223212 + dGammadU2131*dMdGamma223213 + dGammadU2211*dMdGamma223221 + dGammadU2221*dMdGamma223222 + dGammadU2231*dMdGamma223223 + dGammadU2311*dMdGamma223231 + dGammadU2321*dMdGamma223232 + dGammadU2331*dMdGamma223233 + dGammadU3111*dMdGamma223311 + dGammadU3121*dMdGamma223312 + dGammadU3131*dMdGamma223313 + dGammadU3211*dMdGamma223321 + dGammadU3221*dMdGamma223322 + dGammadU3231*dMdGamma223323 + dGammadU3311*dMdGamma223331 + dGammadU3321*dMdGamma223332 + dGammadU3331*dMdGamma223333 + dMdPsi22311*dPhidU111 + dMdPsi22312*dPhidU121 + dMdPsi22313*dPhidU131 + dMdPsi22321*dPhidU211 + dMdPsi22322*dPhidU221 + dMdPsi22323*dPhidU231 + dMdPsi22331*dPhidU311 + dMdPsi22332*dPhidU321 + dMdPsi22333*dPhidU331
    dMdU[I3+18,I1] = dCdU111*dMdC33311 + dCdU121*dMdC33312 + dCdU131*dMdC33313 + dCdU211*dMdC33321 + dCdU221*dMdC33322 + dCdU231*dMdC33323 + dCdU311*dMdC33331 + dCdU321*dMdC33332 + dCdU331*dMdC33333 + dGammadU1111*dMdGamma333111 + dGammadU1121*dMdGamma333112 + dGammadU1131*dMdGamma333113 + dGammadU1211*dMdGamma333121 + dGammadU1221*dMdGamma333122 + dGammadU1231*dMdGamma333123 + dGammadU1311*dMdGamma333131 + dGammadU1321*dMdGamma333132 + dGammadU1331*dMdGamma333133 + dGammadU2111*dMdGamma333211 + dGammadU2121*dMdGamma333212 + dGammadU2131*dMdGamma333213 + dGammadU2211*dMdGamma333221 + dGammadU2221*dMdGamma333222 + dGammadU2231*dMdGamma333223 + dGammadU2311*dMdGamma333231 + dGammadU2321*dMdGamma333232 + dGammadU2331*dMdGamma333233 + dGammadU3111*dMdGamma333311 + dGammadU3121*dMdGamma333312 + dGammadU3131*dMdGamma333313 + dGammadU3211*dMdGamma333321 + dGammadU3221*dMdGamma333322 + dGammadU3231*dMdGamma333323 + dGammadU3311*dMdGamma333331 + dGammadU3321*dMdGamma333332 + dGammadU3331*dMdGamma333333 + dMdPsi33311*dPhidU111 + dMdPsi33312*dPhidU121 + dMdPsi33313*dPhidU131 + dMdPsi33321*dPhidU211 + dMdPsi33322*dPhidU221 + dMdPsi33323*dPhidU231 + dMdPsi33331*dPhidU311 + dMdPsi33332*dPhidU321 + dMdPsi33333*dPhidU331
    dMdU[I4+18,I1] = dCdU111*dMdC23311 + dCdU121*dMdC23312 + dCdU131*dMdC23313 + dCdU211*dMdC23321 + dCdU221*dMdC23322 + dCdU231*dMdC23323 + dCdU311*dMdC23331 + dCdU321*dMdC23332 + dCdU331*dMdC23333 + dGammadU1111*dMdGamma233111 + dGammadU1121*dMdGamma233112 + dGammadU1131*dMdGamma233113 + dGammadU1211*dMdGamma233121 + dGammadU1221*dMdGamma233122 + dGammadU1231*dMdGamma233123 + dGammadU1311*dMdGamma233131 + dGammadU1321*dMdGamma233132 + dGammadU1331*dMdGamma233133 + dGammadU2111*dMdGamma233211 + dGammadU2121*dMdGamma233212 + dGammadU2131*dMdGamma233213 + dGammadU2211*dMdGamma233221 + dGammadU2221*dMdGamma233222 + dGammadU2231*dMdGamma233223 + dGammadU2311*dMdGamma233231 + dGammadU2321*dMdGamma233232 + dGammadU2331*dMdGamma233233 + dGammadU3111*dMdGamma233311 + dGammadU3121*dMdGamma233312 + dGammadU3131*dMdGamma233313 + dGammadU3211*dMdGamma233321 + dGammadU3221*dMdGamma233322 + dGammadU3231*dMdGamma233323 + dGammadU3311*dMdGamma233331 + dGammadU3321*dMdGamma233332 + dGammadU3331*dMdGamma233333 + dMdPsi23311*dPhidU111 + dMdPsi23312*dPhidU121 + dMdPsi23313*dPhidU131 + dMdPsi23321*dPhidU211 + dMdPsi23322*dPhidU221 + dMdPsi23323*dPhidU231 + dMdPsi23331*dPhidU311 + dMdPsi23332*dPhidU321 + dMdPsi23333*dPhidU331
    dMdU[I5+18,I1] = dCdU111*dMdC13311 + dCdU121*dMdC13312 + dCdU131*dMdC13313 + dCdU211*dMdC13321 + dCdU221*dMdC13322 + dCdU231*dMdC13323 + dCdU311*dMdC13331 + dCdU321*dMdC13332 + dCdU331*dMdC13333 + dGammadU1111*dMdGamma133111 + dGammadU1121*dMdGamma133112 + dGammadU1131*dMdGamma133113 + dGammadU1211*dMdGamma133121 + dGammadU1221*dMdGamma133122 + dGammadU1231*dMdGamma133123 + dGammadU1311*dMdGamma133131 + dGammadU1321*dMdGamma133132 + dGammadU1331*dMdGamma133133 + dGammadU2111*dMdGamma133211 + dGammadU2121*dMdGamma133212 + dGammadU2131*dMdGamma133213 + dGammadU2211*dMdGamma133221 + dGammadU2221*dMdGamma133222 + dGammadU2231*dMdGamma133223 + dGammadU2311*dMdGamma133231 + dGammadU2321*dMdGamma133232 + dGammadU2331*dMdGamma133233 + dGammadU3111*dMdGamma133311 + dGammadU3121*dMdGamma133312 + dGammadU3131*dMdGamma133313 + dGammadU3211*dMdGamma133321 + dGammadU3221*dMdGamma133322 + dGammadU3231*dMdGamma133323 + dGammadU3311*dMdGamma133331 + dGammadU3321*dMdGamma133332 + dGammadU3331*dMdGamma133333 + dMdPsi13311*dPhidU111 + dMdPsi13312*dPhidU121 + dMdPsi13313*dPhidU131 + dMdPsi13321*dPhidU211 + dMdPsi13322*dPhidU221 + dMdPsi13323*dPhidU231 + dMdPsi13331*dPhidU311 + dMdPsi13332*dPhidU321 + dMdPsi13333*dPhidU331
    dMdU[I6+18,I1] = dCdU111*dMdC12311 + dCdU121*dMdC12312 + dCdU131*dMdC12313 + dCdU211*dMdC12321 + dCdU221*dMdC12322 + dCdU231*dMdC12323 + dCdU311*dMdC12331 + dCdU321*dMdC12332 + dCdU331*dMdC12333 + dGammadU1111*dMdGamma123111 + dGammadU1121*dMdGamma123112 + dGammadU1131*dMdGamma123113 + dGammadU1211*dMdGamma123121 + dGammadU1221*dMdGamma123122 + dGammadU1231*dMdGamma123123 + dGammadU1311*dMdGamma123131 + dGammadU1321*dMdGamma123132 + dGammadU1331*dMdGamma123133 + dGammadU2111*dMdGamma123211 + dGammadU2121*dMdGamma123212 + dGammadU2131*dMdGamma123213 + dGammadU2211*dMdGamma123221 + dGammadU2221*dMdGamma123222 + dGammadU2231*dMdGamma123223 + dGammadU2311*dMdGamma123231 + dGammadU2321*dMdGamma123232 + dGammadU2331*dMdGamma123233 + dGammadU3111*dMdGamma123311 + dGammadU3121*dMdGamma123312 + dGammadU3131*dMdGamma123313 + dGammadU3211*dMdGamma123321 + dGammadU3221*dMdGamma123322 + dGammadU3231*dMdGamma123323 + dGammadU3311*dMdGamma123331 + dGammadU3321*dMdGamma123332 + dGammadU3331*dMdGamma123333 + dMdPsi12311*dPhidU111 + dMdPsi12312*dPhidU121 + dMdPsi12313*dPhidU131 + dMdPsi12321*dPhidU211 + dMdPsi12322*dPhidU221 + dMdPsi12323*dPhidU231 + dMdPsi12331*dPhidU311 + dMdPsi12332*dPhidU321 + dMdPsi12333*dPhidU331
    dMdU[I7+18,I1] = dCdU111*dMdC32311 + dCdU121*dMdC32312 + dCdU131*dMdC32313 + dCdU211*dMdC32321 + dCdU221*dMdC32322 + dCdU231*dMdC32323 + dCdU311*dMdC32331 + dCdU321*dMdC32332 + dCdU331*dMdC32333 + dGammadU1111*dMdGamma323111 + dGammadU1121*dMdGamma323112 + dGammadU1131*dMdGamma323113 + dGammadU1211*dMdGamma323121 + dGammadU1221*dMdGamma323122 + dGammadU1231*dMdGamma323123 + dGammadU1311*dMdGamma323131 + dGammadU1321*dMdGamma323132 + dGammadU1331*dMdGamma323133 + dGammadU2111*dMdGamma323211 + dGammadU2121*dMdGamma323212 + dGammadU2131*dMdGamma323213 + dGammadU2211*dMdGamma323221 + dGammadU2221*dMdGamma323222 + dGammadU2231*dMdGamma323223 + dGammadU2311*dMdGamma323231 + dGammadU2321*dMdGamma323232 + dGammadU2331*dMdGamma323233 + dGammadU3111*dMdGamma323311 + dGammadU3121*dMdGamma323312 + dGammadU3131*dMdGamma323313 + dGammadU3211*dMdGamma323321 + dGammadU3221*dMdGamma323322 + dGammadU3231*dMdGamma323323 + dGammadU3311*dMdGamma323331 + dGammadU3321*dMdGamma323332 + dGammadU3331*dMdGamma323333 + dMdPsi32311*dPhidU111 + dMdPsi32312*dPhidU121 + dMdPsi32313*dPhidU131 + dMdPsi32321*dPhidU211 + dMdPsi32322*dPhidU221 + dMdPsi32323*dPhidU231 + dMdPsi32331*dPhidU311 + dMdPsi32332*dPhidU321 + dMdPsi32333*dPhidU331
    dMdU[I8+18,I1] = dCdU111*dMdC31311 + dCdU121*dMdC31312 + dCdU131*dMdC31313 + dCdU211*dMdC31321 + dCdU221*dMdC31322 + dCdU231*dMdC31323 + dCdU311*dMdC31331 + dCdU321*dMdC31332 + dCdU331*dMdC31333 + dGammadU1111*dMdGamma313111 + dGammadU1121*dMdGamma313112 + dGammadU1131*dMdGamma313113 + dGammadU1211*dMdGamma313121 + dGammadU1221*dMdGamma313122 + dGammadU1231*dMdGamma313123 + dGammadU1311*dMdGamma313131 + dGammadU1321*dMdGamma313132 + dGammadU1331*dMdGamma313133 + dGammadU2111*dMdGamma313211 + dGammadU2121*dMdGamma313212 + dGammadU2131*dMdGamma313213 + dGammadU2211*dMdGamma313221 + dGammadU2221*dMdGamma313222 + dGammadU2231*dMdGamma313223 + dGammadU2311*dMdGamma313231 + dGammadU2321*dMdGamma313232 + dGammadU2331*dMdGamma313233 + dGammadU3111*dMdGamma313311 + dGammadU3121*dMdGamma313312 + dGammadU3131*dMdGamma313313 + dGammadU3211*dMdGamma313321 + dGammadU3221*dMdGamma313322 + dGammadU3231*dMdGamma313323 + dGammadU3311*dMdGamma313331 + dGammadU3321*dMdGamma313332 + dGammadU3331*dMdGamma313333 + dMdPsi31311*dPhidU111 + dMdPsi31312*dPhidU121 + dMdPsi31313*dPhidU131 + dMdPsi31321*dPhidU211 + dMdPsi31322*dPhidU221 + dMdPsi31323*dPhidU231 + dMdPsi31331*dPhidU311 + dMdPsi31332*dPhidU321 + dMdPsi31333*dPhidU331
    dMdU[I9+18,I1] = dCdU111*dMdC21311 + dCdU121*dMdC21312 + dCdU131*dMdC21313 + dCdU211*dMdC21321 + dCdU221*dMdC21322 + dCdU231*dMdC21323 + dCdU311*dMdC21331 + dCdU321*dMdC21332 + dCdU331*dMdC21333 + dGammadU1111*dMdGamma213111 + dGammadU1121*dMdGamma213112 + dGammadU1131*dMdGamma213113 + dGammadU1211*dMdGamma213121 + dGammadU1221*dMdGamma213122 + dGammadU1231*dMdGamma213123 + dGammadU1311*dMdGamma213131 + dGammadU1321*dMdGamma213132 + dGammadU1331*dMdGamma213133 + dGammadU2111*dMdGamma213211 + dGammadU2121*dMdGamma213212 + dGammadU2131*dMdGamma213213 + dGammadU2211*dMdGamma213221 + dGammadU2221*dMdGamma213222 + dGammadU2231*dMdGamma213223 + dGammadU2311*dMdGamma213231 + dGammadU2321*dMdGamma213232 + dGammadU2331*dMdGamma213233 + dGammadU3111*dMdGamma213311 + dGammadU3121*dMdGamma213312 + dGammadU3131*dMdGamma213313 + dGammadU3211*dMdGamma213321 + dGammadU3221*dMdGamma213322 + dGammadU3231*dMdGamma213323 + dGammadU3311*dMdGamma213331 + dGammadU3321*dMdGamma213332 + dGammadU3331*dMdGamma213333 + dMdPsi21311*dPhidU111 + dMdPsi21312*dPhidU121 + dMdPsi21313*dPhidU131 + dMdPsi21321*dPhidU211 + dMdPsi21322*dPhidU221 + dMdPsi21323*dPhidU231 + dMdPsi21331*dPhidU311 + dMdPsi21332*dPhidU321 + dMdPsi21333*dPhidU331

    #Column 2
    dMdU[I1,I2] = dCdU112*dMdC11111 + dCdU122*dMdC11112 + dCdU132*dMdC11113 + dCdU212*dMdC11121 + dCdU222*dMdC11122 + dCdU232*dMdC11123 + dCdU312*dMdC11131 + dCdU322*dMdC11132 + dCdU332*dMdC11133 + dGammadU1112*dMdGamma111111 + dGammadU1122*dMdGamma111112 + dGammadU1132*dMdGamma111113 + dGammadU1212*dMdGamma111121 + dGammadU1222*dMdGamma111122 + dGammadU1232*dMdGamma111123 + dGammadU1312*dMdGamma111131 + dGammadU1322*dMdGamma111132 + dGammadU1332*dMdGamma111133 + dGammadU2112*dMdGamma111211 + dGammadU2122*dMdGamma111212 + dGammadU2132*dMdGamma111213 + dGammadU2212*dMdGamma111221 + dGammadU2222*dMdGamma111222 + dGammadU2232*dMdGamma111223 + dGammadU2312*dMdGamma111231 + dGammadU2322*dMdGamma111232 + dGammadU2332*dMdGamma111233 + dGammadU3112*dMdGamma111311 + dGammadU3122*dMdGamma111312 + dGammadU3132*dMdGamma111313 + dGammadU3212*dMdGamma111321 + dGammadU3222*dMdGamma111322 + dGammadU3232*dMdGamma111323 + dGammadU3312*dMdGamma111331 + dGammadU3322*dMdGamma111332 + dGammadU3332*dMdGamma111333 + dMdPsi11111*dPhidU112 + dMdPsi11112*dPhidU122 + dMdPsi11113*dPhidU132 + dMdPsi11121*dPhidU212 + dMdPsi11122*dPhidU222 + dMdPsi11123*dPhidU232 + dMdPsi11131*dPhidU312 + dMdPsi11132*dPhidU322 + dMdPsi11133*dPhidU332
    dMdU[I2,I2] = dCdU112*dMdC22111 + dCdU122*dMdC22112 + dCdU132*dMdC22113 + dCdU212*dMdC22121 + dCdU222*dMdC22122 + dCdU232*dMdC22123 + dCdU312*dMdC22131 + dCdU322*dMdC22132 + dCdU332*dMdC22133 + dGammadU1112*dMdGamma221111 + dGammadU1122*dMdGamma221112 + dGammadU1132*dMdGamma221113 + dGammadU1212*dMdGamma221121 + dGammadU1222*dMdGamma221122 + dGammadU1232*dMdGamma221123 + dGammadU1312*dMdGamma221131 + dGammadU1322*dMdGamma221132 + dGammadU1332*dMdGamma221133 + dGammadU2112*dMdGamma221211 + dGammadU2122*dMdGamma221212 + dGammadU2132*dMdGamma221213 + dGammadU2212*dMdGamma221221 + dGammadU2222*dMdGamma221222 + dGammadU2232*dMdGamma221223 + dGammadU2312*dMdGamma221231 + dGammadU2322*dMdGamma221232 + dGammadU2332*dMdGamma221233 + dGammadU3112*dMdGamma221311 + dGammadU3122*dMdGamma221312 + dGammadU3132*dMdGamma221313 + dGammadU3212*dMdGamma221321 + dGammadU3222*dMdGamma221322 + dGammadU3232*dMdGamma221323 + dGammadU3312*dMdGamma221331 + dGammadU3322*dMdGamma221332 + dGammadU3332*dMdGamma221333 + dMdPsi22111*dPhidU112 + dMdPsi22112*dPhidU122 + dMdPsi22113*dPhidU132 + dMdPsi22121*dPhidU212 + dMdPsi22122*dPhidU222 + dMdPsi22123*dPhidU232 + dMdPsi22131*dPhidU312 + dMdPsi22132*dPhidU322 + dMdPsi22133*dPhidU332
    dMdU[I3,I2] = dCdU112*dMdC33111 + dCdU122*dMdC33112 + dCdU132*dMdC33113 + dCdU212*dMdC33121 + dCdU222*dMdC33122 + dCdU232*dMdC33123 + dCdU312*dMdC33131 + dCdU322*dMdC33132 + dCdU332*dMdC33133 + dGammadU1112*dMdGamma331111 + dGammadU1122*dMdGamma331112 + dGammadU1132*dMdGamma331113 + dGammadU1212*dMdGamma331121 + dGammadU1222*dMdGamma331122 + dGammadU1232*dMdGamma331123 + dGammadU1312*dMdGamma331131 + dGammadU1322*dMdGamma331132 + dGammadU1332*dMdGamma331133 + dGammadU2112*dMdGamma331211 + dGammadU2122*dMdGamma331212 + dGammadU2132*dMdGamma331213 + dGammadU2212*dMdGamma331221 + dGammadU2222*dMdGamma331222 + dGammadU2232*dMdGamma331223 + dGammadU2312*dMdGamma331231 + dGammadU2322*dMdGamma331232 + dGammadU2332*dMdGamma331233 + dGammadU3112*dMdGamma331311 + dGammadU3122*dMdGamma331312 + dGammadU3132*dMdGamma331313 + dGammadU3212*dMdGamma331321 + dGammadU3222*dMdGamma331322 + dGammadU3232*dMdGamma331323 + dGammadU3312*dMdGamma331331 + dGammadU3322*dMdGamma331332 + dGammadU3332*dMdGamma331333 + dMdPsi33111*dPhidU112 + dMdPsi33112*dPhidU122 + dMdPsi33113*dPhidU132 + dMdPsi33121*dPhidU212 + dMdPsi33122*dPhidU222 + dMdPsi33123*dPhidU232 + dMdPsi33131*dPhidU312 + dMdPsi33132*dPhidU322 + dMdPsi33133*dPhidU332
    dMdU[I4,I2] = dCdU112*dMdC23111 + dCdU122*dMdC23112 + dCdU132*dMdC23113 + dCdU212*dMdC23121 + dCdU222*dMdC23122 + dCdU232*dMdC23123 + dCdU312*dMdC23131 + dCdU322*dMdC23132 + dCdU332*dMdC23133 + dGammadU1112*dMdGamma231111 + dGammadU1122*dMdGamma231112 + dGammadU1132*dMdGamma231113 + dGammadU1212*dMdGamma231121 + dGammadU1222*dMdGamma231122 + dGammadU1232*dMdGamma231123 + dGammadU1312*dMdGamma231131 + dGammadU1322*dMdGamma231132 + dGammadU1332*dMdGamma231133 + dGammadU2112*dMdGamma231211 + dGammadU2122*dMdGamma231212 + dGammadU2132*dMdGamma231213 + dGammadU2212*dMdGamma231221 + dGammadU2222*dMdGamma231222 + dGammadU2232*dMdGamma231223 + dGammadU2312*dMdGamma231231 + dGammadU2322*dMdGamma231232 + dGammadU2332*dMdGamma231233 + dGammadU3112*dMdGamma231311 + dGammadU3122*dMdGamma231312 + dGammadU3132*dMdGamma231313 + dGammadU3212*dMdGamma231321 + dGammadU3222*dMdGamma231322 + dGammadU3232*dMdGamma231323 + dGammadU3312*dMdGamma231331 + dGammadU3322*dMdGamma231332 + dGammadU3332*dMdGamma231333 + dMdPsi23111*dPhidU112 + dMdPsi23112*dPhidU122 + dMdPsi23113*dPhidU132 + dMdPsi23121*dPhidU212 + dMdPsi23122*dPhidU222 + dMdPsi23123*dPhidU232 + dMdPsi23131*dPhidU312 + dMdPsi23132*dPhidU322 + dMdPsi23133*dPhidU332
    dMdU[I5,I2] = dCdU112*dMdC13111 + dCdU122*dMdC13112 + dCdU132*dMdC13113 + dCdU212*dMdC13121 + dCdU222*dMdC13122 + dCdU232*dMdC13123 + dCdU312*dMdC13131 + dCdU322*dMdC13132 + dCdU332*dMdC13133 + dGammadU1112*dMdGamma131111 + dGammadU1122*dMdGamma131112 + dGammadU1132*dMdGamma131113 + dGammadU1212*dMdGamma131121 + dGammadU1222*dMdGamma131122 + dGammadU1232*dMdGamma131123 + dGammadU1312*dMdGamma131131 + dGammadU1322*dMdGamma131132 + dGammadU1332*dMdGamma131133 + dGammadU2112*dMdGamma131211 + dGammadU2122*dMdGamma131212 + dGammadU2132*dMdGamma131213 + dGammadU2212*dMdGamma131221 + dGammadU2222*dMdGamma131222 + dGammadU2232*dMdGamma131223 + dGammadU2312*dMdGamma131231 + dGammadU2322*dMdGamma131232 + dGammadU2332*dMdGamma131233 + dGammadU3112*dMdGamma131311 + dGammadU3122*dMdGamma131312 + dGammadU3132*dMdGamma131313 + dGammadU3212*dMdGamma131321 + dGammadU3222*dMdGamma131322 + dGammadU3232*dMdGamma131323 + dGammadU3312*dMdGamma131331 + dGammadU3322*dMdGamma131332 + dGammadU3332*dMdGamma131333 + dMdPsi13111*dPhidU112 + dMdPsi13112*dPhidU122 + dMdPsi13113*dPhidU132 + dMdPsi13121*dPhidU212 + dMdPsi13122*dPhidU222 + dMdPsi13123*dPhidU232 + dMdPsi13131*dPhidU312 + dMdPsi13132*dPhidU322 + dMdPsi13133*dPhidU332
    dMdU[I6,I2] = dCdU112*dMdC12111 + dCdU122*dMdC12112 + dCdU132*dMdC12113 + dCdU212*dMdC12121 + dCdU222*dMdC12122 + dCdU232*dMdC12123 + dCdU312*dMdC12131 + dCdU322*dMdC12132 + dCdU332*dMdC12133 + dGammadU1112*dMdGamma121111 + dGammadU1122*dMdGamma121112 + dGammadU1132*dMdGamma121113 + dGammadU1212*dMdGamma121121 + dGammadU1222*dMdGamma121122 + dGammadU1232*dMdGamma121123 + dGammadU1312*dMdGamma121131 + dGammadU1322*dMdGamma121132 + dGammadU1332*dMdGamma121133 + dGammadU2112*dMdGamma121211 + dGammadU2122*dMdGamma121212 + dGammadU2132*dMdGamma121213 + dGammadU2212*dMdGamma121221 + dGammadU2222*dMdGamma121222 + dGammadU2232*dMdGamma121223 + dGammadU2312*dMdGamma121231 + dGammadU2322*dMdGamma121232 + dGammadU2332*dMdGamma121233 + dGammadU3112*dMdGamma121311 + dGammadU3122*dMdGamma121312 + dGammadU3132*dMdGamma121313 + dGammadU3212*dMdGamma121321 + dGammadU3222*dMdGamma121322 + dGammadU3232*dMdGamma121323 + dGammadU3312*dMdGamma121331 + dGammadU3322*dMdGamma121332 + dGammadU3332*dMdGamma121333 + dMdPsi12111*dPhidU112 + dMdPsi12112*dPhidU122 + dMdPsi12113*dPhidU132 + dMdPsi12121*dPhidU212 + dMdPsi12122*dPhidU222 + dMdPsi12123*dPhidU232 + dMdPsi12131*dPhidU312 + dMdPsi12132*dPhidU322 + dMdPsi12133*dPhidU332
    dMdU[I7,I2] = dCdU112*dMdC32111 + dCdU122*dMdC32112 + dCdU132*dMdC32113 + dCdU212*dMdC32121 + dCdU222*dMdC32122 + dCdU232*dMdC32123 + dCdU312*dMdC32131 + dCdU322*dMdC32132 + dCdU332*dMdC32133 + dGammadU1112*dMdGamma321111 + dGammadU1122*dMdGamma321112 + dGammadU1132*dMdGamma321113 + dGammadU1212*dMdGamma321121 + dGammadU1222*dMdGamma321122 + dGammadU1232*dMdGamma321123 + dGammadU1312*dMdGamma321131 + dGammadU1322*dMdGamma321132 + dGammadU1332*dMdGamma321133 + dGammadU2112*dMdGamma321211 + dGammadU2122*dMdGamma321212 + dGammadU2132*dMdGamma321213 + dGammadU2212*dMdGamma321221 + dGammadU2222*dMdGamma321222 + dGammadU2232*dMdGamma321223 + dGammadU2312*dMdGamma321231 + dGammadU2322*dMdGamma321232 + dGammadU2332*dMdGamma321233 + dGammadU3112*dMdGamma321311 + dGammadU3122*dMdGamma321312 + dGammadU3132*dMdGamma321313 + dGammadU3212*dMdGamma321321 + dGammadU3222*dMdGamma321322 + dGammadU3232*dMdGamma321323 + dGammadU3312*dMdGamma321331 + dGammadU3322*dMdGamma321332 + dGammadU3332*dMdGamma321333 + dMdPsi32111*dPhidU112 + dMdPsi32112*dPhidU122 + dMdPsi32113*dPhidU132 + dMdPsi32121*dPhidU212 + dMdPsi32122*dPhidU222 + dMdPsi32123*dPhidU232 + dMdPsi32131*dPhidU312 + dMdPsi32132*dPhidU322 + dMdPsi32133*dPhidU332
    dMdU[I8,I2] = dCdU112*dMdC31111 + dCdU122*dMdC31112 + dCdU132*dMdC31113 + dCdU212*dMdC31121 + dCdU222*dMdC31122 + dCdU232*dMdC31123 + dCdU312*dMdC31131 + dCdU322*dMdC31132 + dCdU332*dMdC31133 + dGammadU1112*dMdGamma311111 + dGammadU1122*dMdGamma311112 + dGammadU1132*dMdGamma311113 + dGammadU1212*dMdGamma311121 + dGammadU1222*dMdGamma311122 + dGammadU1232*dMdGamma311123 + dGammadU1312*dMdGamma311131 + dGammadU1322*dMdGamma311132 + dGammadU1332*dMdGamma311133 + dGammadU2112*dMdGamma311211 + dGammadU2122*dMdGamma311212 + dGammadU2132*dMdGamma311213 + dGammadU2212*dMdGamma311221 + dGammadU2222*dMdGamma311222 + dGammadU2232*dMdGamma311223 + dGammadU2312*dMdGamma311231 + dGammadU2322*dMdGamma311232 + dGammadU2332*dMdGamma311233 + dGammadU3112*dMdGamma311311 + dGammadU3122*dMdGamma311312 + dGammadU3132*dMdGamma311313 + dGammadU3212*dMdGamma311321 + dGammadU3222*dMdGamma311322 + dGammadU3232*dMdGamma311323 + dGammadU3312*dMdGamma311331 + dGammadU3322*dMdGamma311332 + dGammadU3332*dMdGamma311333 + dMdPsi31111*dPhidU112 + dMdPsi31112*dPhidU122 + dMdPsi31113*dPhidU132 + dMdPsi31121*dPhidU212 + dMdPsi31122*dPhidU222 + dMdPsi31123*dPhidU232 + dMdPsi31131*dPhidU312 + dMdPsi31132*dPhidU322 + dMdPsi31133*dPhidU332
    dMdU[I9,I2] = dCdU112*dMdC21111 + dCdU122*dMdC21112 + dCdU132*dMdC21113 + dCdU212*dMdC21121 + dCdU222*dMdC21122 + dCdU232*dMdC21123 + dCdU312*dMdC21131 + dCdU322*dMdC21132 + dCdU332*dMdC21133 + dGammadU1112*dMdGamma211111 + dGammadU1122*dMdGamma211112 + dGammadU1132*dMdGamma211113 + dGammadU1212*dMdGamma211121 + dGammadU1222*dMdGamma211122 + dGammadU1232*dMdGamma211123 + dGammadU1312*dMdGamma211131 + dGammadU1322*dMdGamma211132 + dGammadU1332*dMdGamma211133 + dGammadU2112*dMdGamma211211 + dGammadU2122*dMdGamma211212 + dGammadU2132*dMdGamma211213 + dGammadU2212*dMdGamma211221 + dGammadU2222*dMdGamma211222 + dGammadU2232*dMdGamma211223 + dGammadU2312*dMdGamma211231 + dGammadU2322*dMdGamma211232 + dGammadU2332*dMdGamma211233 + dGammadU3112*dMdGamma211311 + dGammadU3122*dMdGamma211312 + dGammadU3132*dMdGamma211313 + dGammadU3212*dMdGamma211321 + dGammadU3222*dMdGamma211322 + dGammadU3232*dMdGamma211323 + dGammadU3312*dMdGamma211331 + dGammadU3322*dMdGamma211332 + dGammadU3332*dMdGamma211333 + dMdPsi21111*dPhidU112 + dMdPsi21112*dPhidU122 + dMdPsi21113*dPhidU132 + dMdPsi21121*dPhidU212 + dMdPsi21122*dPhidU222 + dMdPsi21123*dPhidU232 + dMdPsi21131*dPhidU312 + dMdPsi21132*dPhidU322 + dMdPsi21133*dPhidU332
    dMdU[I1+9,I2] = dCdU112*dMdC11211 + dCdU122*dMdC11212 + dCdU132*dMdC11213 + dCdU212*dMdC11221 + dCdU222*dMdC11222 + dCdU232*dMdC11223 + dCdU312*dMdC11231 + dCdU322*dMdC11232 + dCdU332*dMdC11233 + dGammadU1112*dMdGamma112111 + dGammadU1122*dMdGamma112112 + dGammadU1132*dMdGamma112113 + dGammadU1212*dMdGamma112121 + dGammadU1222*dMdGamma112122 + dGammadU1232*dMdGamma112123 + dGammadU1312*dMdGamma112131 + dGammadU1322*dMdGamma112132 + dGammadU1332*dMdGamma112133 + dGammadU2112*dMdGamma112211 + dGammadU2122*dMdGamma112212 + dGammadU2132*dMdGamma112213 + dGammadU2212*dMdGamma112221 + dGammadU2222*dMdGamma112222 + dGammadU2232*dMdGamma112223 + dGammadU2312*dMdGamma112231 + dGammadU2322*dMdGamma112232 + dGammadU2332*dMdGamma112233 + dGammadU3112*dMdGamma112311 + dGammadU3122*dMdGamma112312 + dGammadU3132*dMdGamma112313 + dGammadU3212*dMdGamma112321 + dGammadU3222*dMdGamma112322 + dGammadU3232*dMdGamma112323 + dGammadU3312*dMdGamma112331 + dGammadU3322*dMdGamma112332 + dGammadU3332*dMdGamma112333 + dMdPsi11211*dPhidU112 + dMdPsi11212*dPhidU122 + dMdPsi11213*dPhidU132 + dMdPsi11221*dPhidU212 + dMdPsi11222*dPhidU222 + dMdPsi11223*dPhidU232 + dMdPsi11231*dPhidU312 + dMdPsi11232*dPhidU322 + dMdPsi11233*dPhidU332
    dMdU[I2+9,I2] = dCdU112*dMdC22211 + dCdU122*dMdC22212 + dCdU132*dMdC22213 + dCdU212*dMdC22221 + dCdU222*dMdC22222 + dCdU232*dMdC22223 + dCdU312*dMdC22231 + dCdU322*dMdC22232 + dCdU332*dMdC22233 + dGammadU1112*dMdGamma222111 + dGammadU1122*dMdGamma222112 + dGammadU1132*dMdGamma222113 + dGammadU1212*dMdGamma222121 + dGammadU1222*dMdGamma222122 + dGammadU1232*dMdGamma222123 + dGammadU1312*dMdGamma222131 + dGammadU1322*dMdGamma222132 + dGammadU1332*dMdGamma222133 + dGammadU2112*dMdGamma222211 + dGammadU2122*dMdGamma222212 + dGammadU2132*dMdGamma222213 + dGammadU2212*dMdGamma222221 + dGammadU2222*dMdGamma222222 + dGammadU2232*dMdGamma222223 + dGammadU2312*dMdGamma222231 + dGammadU2322*dMdGamma222232 + dGammadU2332*dMdGamma222233 + dGammadU3112*dMdGamma222311 + dGammadU3122*dMdGamma222312 + dGammadU3132*dMdGamma222313 + dGammadU3212*dMdGamma222321 + dGammadU3222*dMdGamma222322 + dGammadU3232*dMdGamma222323 + dGammadU3312*dMdGamma222331 + dGammadU3322*dMdGamma222332 + dGammadU3332*dMdGamma222333 + dMdPsi22211*dPhidU112 + dMdPsi22212*dPhidU122 + dMdPsi22213*dPhidU132 + dMdPsi22221*dPhidU212 + dMdPsi22222*dPhidU222 + dMdPsi22223*dPhidU232 + dMdPsi22231*dPhidU312 + dMdPsi22232*dPhidU322 + dMdPsi22233*dPhidU332
    dMdU[I3+9,I2] = dCdU112*dMdC33211 + dCdU122*dMdC33212 + dCdU132*dMdC33213 + dCdU212*dMdC33221 + dCdU222*dMdC33222 + dCdU232*dMdC33223 + dCdU312*dMdC33231 + dCdU322*dMdC33232 + dCdU332*dMdC33233 + dGammadU1112*dMdGamma332111 + dGammadU1122*dMdGamma332112 + dGammadU1132*dMdGamma332113 + dGammadU1212*dMdGamma332121 + dGammadU1222*dMdGamma332122 + dGammadU1232*dMdGamma332123 + dGammadU1312*dMdGamma332131 + dGammadU1322*dMdGamma332132 + dGammadU1332*dMdGamma332133 + dGammadU2112*dMdGamma332211 + dGammadU2122*dMdGamma332212 + dGammadU2132*dMdGamma332213 + dGammadU2212*dMdGamma332221 + dGammadU2222*dMdGamma332222 + dGammadU2232*dMdGamma332223 + dGammadU2312*dMdGamma332231 + dGammadU2322*dMdGamma332232 + dGammadU2332*dMdGamma332233 + dGammadU3112*dMdGamma332311 + dGammadU3122*dMdGamma332312 + dGammadU3132*dMdGamma332313 + dGammadU3212*dMdGamma332321 + dGammadU3222*dMdGamma332322 + dGammadU3232*dMdGamma332323 + dGammadU3312*dMdGamma332331 + dGammadU3322*dMdGamma332332 + dGammadU3332*dMdGamma332333 + dMdPsi33211*dPhidU112 + dMdPsi33212*dPhidU122 + dMdPsi33213*dPhidU132 + dMdPsi33221*dPhidU212 + dMdPsi33222*dPhidU222 + dMdPsi33223*dPhidU232 + dMdPsi33231*dPhidU312 + dMdPsi33232*dPhidU322 + dMdPsi33233*dPhidU332
    dMdU[I4+9,I2] = dCdU112*dMdC23211 + dCdU122*dMdC23212 + dCdU132*dMdC23213 + dCdU212*dMdC23221 + dCdU222*dMdC23222 + dCdU232*dMdC23223 + dCdU312*dMdC23231 + dCdU322*dMdC23232 + dCdU332*dMdC23233 + dGammadU1112*dMdGamma232111 + dGammadU1122*dMdGamma232112 + dGammadU1132*dMdGamma232113 + dGammadU1212*dMdGamma232121 + dGammadU1222*dMdGamma232122 + dGammadU1232*dMdGamma232123 + dGammadU1312*dMdGamma232131 + dGammadU1322*dMdGamma232132 + dGammadU1332*dMdGamma232133 + dGammadU2112*dMdGamma232211 + dGammadU2122*dMdGamma232212 + dGammadU2132*dMdGamma232213 + dGammadU2212*dMdGamma232221 + dGammadU2222*dMdGamma232222 + dGammadU2232*dMdGamma232223 + dGammadU2312*dMdGamma232231 + dGammadU2322*dMdGamma232232 + dGammadU2332*dMdGamma232233 + dGammadU3112*dMdGamma232311 + dGammadU3122*dMdGamma232312 + dGammadU3132*dMdGamma232313 + dGammadU3212*dMdGamma232321 + dGammadU3222*dMdGamma232322 + dGammadU3232*dMdGamma232323 + dGammadU3312*dMdGamma232331 + dGammadU3322*dMdGamma232332 + dGammadU3332*dMdGamma232333 + dMdPsi23211*dPhidU112 + dMdPsi23212*dPhidU122 + dMdPsi23213*dPhidU132 + dMdPsi23221*dPhidU212 + dMdPsi23222*dPhidU222 + dMdPsi23223*dPhidU232 + dMdPsi23231*dPhidU312 + dMdPsi23232*dPhidU322 + dMdPsi23233*dPhidU332
    dMdU[I5+9,I2] = dCdU112*dMdC13211 + dCdU122*dMdC13212 + dCdU132*dMdC13213 + dCdU212*dMdC13221 + dCdU222*dMdC13222 + dCdU232*dMdC13223 + dCdU312*dMdC13231 + dCdU322*dMdC13232 + dCdU332*dMdC13233 + dGammadU1112*dMdGamma132111 + dGammadU1122*dMdGamma132112 + dGammadU1132*dMdGamma132113 + dGammadU1212*dMdGamma132121 + dGammadU1222*dMdGamma132122 + dGammadU1232*dMdGamma132123 + dGammadU1312*dMdGamma132131 + dGammadU1322*dMdGamma132132 + dGammadU1332*dMdGamma132133 + dGammadU2112*dMdGamma132211 + dGammadU2122*dMdGamma132212 + dGammadU2132*dMdGamma132213 + dGammadU2212*dMdGamma132221 + dGammadU2222*dMdGamma132222 + dGammadU2232*dMdGamma132223 + dGammadU2312*dMdGamma132231 + dGammadU2322*dMdGamma132232 + dGammadU2332*dMdGamma132233 + dGammadU3112*dMdGamma132311 + dGammadU3122*dMdGamma132312 + dGammadU3132*dMdGamma132313 + dGammadU3212*dMdGamma132321 + dGammadU3222*dMdGamma132322 + dGammadU3232*dMdGamma132323 + dGammadU3312*dMdGamma132331 + dGammadU3322*dMdGamma132332 + dGammadU3332*dMdGamma132333 + dMdPsi13211*dPhidU112 + dMdPsi13212*dPhidU122 + dMdPsi13213*dPhidU132 + dMdPsi13221*dPhidU212 + dMdPsi13222*dPhidU222 + dMdPsi13223*dPhidU232 + dMdPsi13231*dPhidU312 + dMdPsi13232*dPhidU322 + dMdPsi13233*dPhidU332
    dMdU[I6+9,I2] = dCdU112*dMdC12211 + dCdU122*dMdC12212 + dCdU132*dMdC12213 + dCdU212*dMdC12221 + dCdU222*dMdC12222 + dCdU232*dMdC12223 + dCdU312*dMdC12231 + dCdU322*dMdC12232 + dCdU332*dMdC12233 + dGammadU1112*dMdGamma122111 + dGammadU1122*dMdGamma122112 + dGammadU1132*dMdGamma122113 + dGammadU1212*dMdGamma122121 + dGammadU1222*dMdGamma122122 + dGammadU1232*dMdGamma122123 + dGammadU1312*dMdGamma122131 + dGammadU1322*dMdGamma122132 + dGammadU1332*dMdGamma122133 + dGammadU2112*dMdGamma122211 + dGammadU2122*dMdGamma122212 + dGammadU2132*dMdGamma122213 + dGammadU2212*dMdGamma122221 + dGammadU2222*dMdGamma122222 + dGammadU2232*dMdGamma122223 + dGammadU2312*dMdGamma122231 + dGammadU2322*dMdGamma122232 + dGammadU2332*dMdGamma122233 + dGammadU3112*dMdGamma122311 + dGammadU3122*dMdGamma122312 + dGammadU3132*dMdGamma122313 + dGammadU3212*dMdGamma122321 + dGammadU3222*dMdGamma122322 + dGammadU3232*dMdGamma122323 + dGammadU3312*dMdGamma122331 + dGammadU3322*dMdGamma122332 + dGammadU3332*dMdGamma122333 + dMdPsi12211*dPhidU112 + dMdPsi12212*dPhidU122 + dMdPsi12213*dPhidU132 + dMdPsi12221*dPhidU212 + dMdPsi12222*dPhidU222 + dMdPsi12223*dPhidU232 + dMdPsi12231*dPhidU312 + dMdPsi12232*dPhidU322 + dMdPsi12233*dPhidU332
    dMdU[I7+9,I2] = dCdU112*dMdC32211 + dCdU122*dMdC32212 + dCdU132*dMdC32213 + dCdU212*dMdC32221 + dCdU222*dMdC32222 + dCdU232*dMdC32223 + dCdU312*dMdC32231 + dCdU322*dMdC32232 + dCdU332*dMdC32233 + dGammadU1112*dMdGamma322111 + dGammadU1122*dMdGamma322112 + dGammadU1132*dMdGamma322113 + dGammadU1212*dMdGamma322121 + dGammadU1222*dMdGamma322122 + dGammadU1232*dMdGamma322123 + dGammadU1312*dMdGamma322131 + dGammadU1322*dMdGamma322132 + dGammadU1332*dMdGamma322133 + dGammadU2112*dMdGamma322211 + dGammadU2122*dMdGamma322212 + dGammadU2132*dMdGamma322213 + dGammadU2212*dMdGamma322221 + dGammadU2222*dMdGamma322222 + dGammadU2232*dMdGamma322223 + dGammadU2312*dMdGamma322231 + dGammadU2322*dMdGamma322232 + dGammadU2332*dMdGamma322233 + dGammadU3112*dMdGamma322311 + dGammadU3122*dMdGamma322312 + dGammadU3132*dMdGamma322313 + dGammadU3212*dMdGamma322321 + dGammadU3222*dMdGamma322322 + dGammadU3232*dMdGamma322323 + dGammadU3312*dMdGamma322331 + dGammadU3322*dMdGamma322332 + dGammadU3332*dMdGamma322333 + dMdPsi32211*dPhidU112 + dMdPsi32212*dPhidU122 + dMdPsi32213*dPhidU132 + dMdPsi32221*dPhidU212 + dMdPsi32222*dPhidU222 + dMdPsi32223*dPhidU232 + dMdPsi32231*dPhidU312 + dMdPsi32232*dPhidU322 + dMdPsi32233*dPhidU332
    dMdU[I8+9,I2] = dCdU112*dMdC31211 + dCdU122*dMdC31212 + dCdU132*dMdC31213 + dCdU212*dMdC31221 + dCdU222*dMdC31222 + dCdU232*dMdC31223 + dCdU312*dMdC31231 + dCdU322*dMdC31232 + dCdU332*dMdC31233 + dGammadU1112*dMdGamma312111 + dGammadU1122*dMdGamma312112 + dGammadU1132*dMdGamma312113 + dGammadU1212*dMdGamma312121 + dGammadU1222*dMdGamma312122 + dGammadU1232*dMdGamma312123 + dGammadU1312*dMdGamma312131 + dGammadU1322*dMdGamma312132 + dGammadU1332*dMdGamma312133 + dGammadU2112*dMdGamma312211 + dGammadU2122*dMdGamma312212 + dGammadU2132*dMdGamma312213 + dGammadU2212*dMdGamma312221 + dGammadU2222*dMdGamma312222 + dGammadU2232*dMdGamma312223 + dGammadU2312*dMdGamma312231 + dGammadU2322*dMdGamma312232 + dGammadU2332*dMdGamma312233 + dGammadU3112*dMdGamma312311 + dGammadU3122*dMdGamma312312 + dGammadU3132*dMdGamma312313 + dGammadU3212*dMdGamma312321 + dGammadU3222*dMdGamma312322 + dGammadU3232*dMdGamma312323 + dGammadU3312*dMdGamma312331 + dGammadU3322*dMdGamma312332 + dGammadU3332*dMdGamma312333 + dMdPsi31211*dPhidU112 + dMdPsi31212*dPhidU122 + dMdPsi31213*dPhidU132 + dMdPsi31221*dPhidU212 + dMdPsi31222*dPhidU222 + dMdPsi31223*dPhidU232 + dMdPsi31231*dPhidU312 + dMdPsi31232*dPhidU322 + dMdPsi31233*dPhidU332
    dMdU[I9+9,I2] = dCdU112*dMdC21211 + dCdU122*dMdC21212 + dCdU132*dMdC21213 + dCdU212*dMdC21221 + dCdU222*dMdC21222 + dCdU232*dMdC21223 + dCdU312*dMdC21231 + dCdU322*dMdC21232 + dCdU332*dMdC21233 + dGammadU1112*dMdGamma212111 + dGammadU1122*dMdGamma212112 + dGammadU1132*dMdGamma212113 + dGammadU1212*dMdGamma212121 + dGammadU1222*dMdGamma212122 + dGammadU1232*dMdGamma212123 + dGammadU1312*dMdGamma212131 + dGammadU1322*dMdGamma212132 + dGammadU1332*dMdGamma212133 + dGammadU2112*dMdGamma212211 + dGammadU2122*dMdGamma212212 + dGammadU2132*dMdGamma212213 + dGammadU2212*dMdGamma212221 + dGammadU2222*dMdGamma212222 + dGammadU2232*dMdGamma212223 + dGammadU2312*dMdGamma212231 + dGammadU2322*dMdGamma212232 + dGammadU2332*dMdGamma212233 + dGammadU3112*dMdGamma212311 + dGammadU3122*dMdGamma212312 + dGammadU3132*dMdGamma212313 + dGammadU3212*dMdGamma212321 + dGammadU3222*dMdGamma212322 + dGammadU3232*dMdGamma212323 + dGammadU3312*dMdGamma212331 + dGammadU3322*dMdGamma212332 + dGammadU3332*dMdGamma212333 + dMdPsi21211*dPhidU112 + dMdPsi21212*dPhidU122 + dMdPsi21213*dPhidU132 + dMdPsi21221*dPhidU212 + dMdPsi21222*dPhidU222 + dMdPsi21223*dPhidU232 + dMdPsi21231*dPhidU312 + dMdPsi21232*dPhidU322 + dMdPsi21233*dPhidU332
    dMdU[I1+18,I2] = dCdU112*dMdC11311 + dCdU122*dMdC11312 + dCdU132*dMdC11313 + dCdU212*dMdC11321 + dCdU222*dMdC11322 + dCdU232*dMdC11323 + dCdU312*dMdC11331 + dCdU322*dMdC11332 + dCdU332*dMdC11333 + dGammadU1112*dMdGamma113111 + dGammadU1122*dMdGamma113112 + dGammadU1132*dMdGamma113113 + dGammadU1212*dMdGamma113121 + dGammadU1222*dMdGamma113122 + dGammadU1232*dMdGamma113123 + dGammadU1312*dMdGamma113131 + dGammadU1322*dMdGamma113132 + dGammadU1332*dMdGamma113133 + dGammadU2112*dMdGamma113211 + dGammadU2122*dMdGamma113212 + dGammadU2132*dMdGamma113213 + dGammadU2212*dMdGamma113221 + dGammadU2222*dMdGamma113222 + dGammadU2232*dMdGamma113223 + dGammadU2312*dMdGamma113231 + dGammadU2322*dMdGamma113232 + dGammadU2332*dMdGamma113233 + dGammadU3112*dMdGamma113311 + dGammadU3122*dMdGamma113312 + dGammadU3132*dMdGamma113313 + dGammadU3212*dMdGamma113321 + dGammadU3222*dMdGamma113322 + dGammadU3232*dMdGamma113323 + dGammadU3312*dMdGamma113331 + dGammadU3322*dMdGamma113332 + dGammadU3332*dMdGamma113333 + dMdPsi11311*dPhidU112 + dMdPsi11312*dPhidU122 + dMdPsi11313*dPhidU132 + dMdPsi11321*dPhidU212 + dMdPsi11322*dPhidU222 + dMdPsi11323*dPhidU232 + dMdPsi11331*dPhidU312 + dMdPsi11332*dPhidU322 + dMdPsi11333*dPhidU332
    dMdU[I2+18,I2] = dCdU112*dMdC22311 + dCdU122*dMdC22312 + dCdU132*dMdC22313 + dCdU212*dMdC22321 + dCdU222*dMdC22322 + dCdU232*dMdC22323 + dCdU312*dMdC22331 + dCdU322*dMdC22332 + dCdU332*dMdC22333 + dGammadU1112*dMdGamma223111 + dGammadU1122*dMdGamma223112 + dGammadU1132*dMdGamma223113 + dGammadU1212*dMdGamma223121 + dGammadU1222*dMdGamma223122 + dGammadU1232*dMdGamma223123 + dGammadU1312*dMdGamma223131 + dGammadU1322*dMdGamma223132 + dGammadU1332*dMdGamma223133 + dGammadU2112*dMdGamma223211 + dGammadU2122*dMdGamma223212 + dGammadU2132*dMdGamma223213 + dGammadU2212*dMdGamma223221 + dGammadU2222*dMdGamma223222 + dGammadU2232*dMdGamma223223 + dGammadU2312*dMdGamma223231 + dGammadU2322*dMdGamma223232 + dGammadU2332*dMdGamma223233 + dGammadU3112*dMdGamma223311 + dGammadU3122*dMdGamma223312 + dGammadU3132*dMdGamma223313 + dGammadU3212*dMdGamma223321 + dGammadU3222*dMdGamma223322 + dGammadU3232*dMdGamma223323 + dGammadU3312*dMdGamma223331 + dGammadU3322*dMdGamma223332 + dGammadU3332*dMdGamma223333 + dMdPsi22311*dPhidU112 + dMdPsi22312*dPhidU122 + dMdPsi22313*dPhidU132 + dMdPsi22321*dPhidU212 + dMdPsi22322*dPhidU222 + dMdPsi22323*dPhidU232 + dMdPsi22331*dPhidU312 + dMdPsi22332*dPhidU322 + dMdPsi22333*dPhidU332
    dMdU[I3+18,I2] = dCdU112*dMdC33311 + dCdU122*dMdC33312 + dCdU132*dMdC33313 + dCdU212*dMdC33321 + dCdU222*dMdC33322 + dCdU232*dMdC33323 + dCdU312*dMdC33331 + dCdU322*dMdC33332 + dCdU332*dMdC33333 + dGammadU1112*dMdGamma333111 + dGammadU1122*dMdGamma333112 + dGammadU1132*dMdGamma333113 + dGammadU1212*dMdGamma333121 + dGammadU1222*dMdGamma333122 + dGammadU1232*dMdGamma333123 + dGammadU1312*dMdGamma333131 + dGammadU1322*dMdGamma333132 + dGammadU1332*dMdGamma333133 + dGammadU2112*dMdGamma333211 + dGammadU2122*dMdGamma333212 + dGammadU2132*dMdGamma333213 + dGammadU2212*dMdGamma333221 + dGammadU2222*dMdGamma333222 + dGammadU2232*dMdGamma333223 + dGammadU2312*dMdGamma333231 + dGammadU2322*dMdGamma333232 + dGammadU2332*dMdGamma333233 + dGammadU3112*dMdGamma333311 + dGammadU3122*dMdGamma333312 + dGammadU3132*dMdGamma333313 + dGammadU3212*dMdGamma333321 + dGammadU3222*dMdGamma333322 + dGammadU3232*dMdGamma333323 + dGammadU3312*dMdGamma333331 + dGammadU3322*dMdGamma333332 + dGammadU3332*dMdGamma333333 + dMdPsi33311*dPhidU112 + dMdPsi33312*dPhidU122 + dMdPsi33313*dPhidU132 + dMdPsi33321*dPhidU212 + dMdPsi33322*dPhidU222 + dMdPsi33323*dPhidU232 + dMdPsi33331*dPhidU312 + dMdPsi33332*dPhidU322 + dMdPsi33333*dPhidU332
    dMdU[I4+18,I2] = dCdU112*dMdC23311 + dCdU122*dMdC23312 + dCdU132*dMdC23313 + dCdU212*dMdC23321 + dCdU222*dMdC23322 + dCdU232*dMdC23323 + dCdU312*dMdC23331 + dCdU322*dMdC23332 + dCdU332*dMdC23333 + dGammadU1112*dMdGamma233111 + dGammadU1122*dMdGamma233112 + dGammadU1132*dMdGamma233113 + dGammadU1212*dMdGamma233121 + dGammadU1222*dMdGamma233122 + dGammadU1232*dMdGamma233123 + dGammadU1312*dMdGamma233131 + dGammadU1322*dMdGamma233132 + dGammadU1332*dMdGamma233133 + dGammadU2112*dMdGamma233211 + dGammadU2122*dMdGamma233212 + dGammadU2132*dMdGamma233213 + dGammadU2212*dMdGamma233221 + dGammadU2222*dMdGamma233222 + dGammadU2232*dMdGamma233223 + dGammadU2312*dMdGamma233231 + dGammadU2322*dMdGamma233232 + dGammadU2332*dMdGamma233233 + dGammadU3112*dMdGamma233311 + dGammadU3122*dMdGamma233312 + dGammadU3132*dMdGamma233313 + dGammadU3212*dMdGamma233321 + dGammadU3222*dMdGamma233322 + dGammadU3232*dMdGamma233323 + dGammadU3312*dMdGamma233331 + dGammadU3322*dMdGamma233332 + dGammadU3332*dMdGamma233333 + dMdPsi23311*dPhidU112 + dMdPsi23312*dPhidU122 + dMdPsi23313*dPhidU132 + dMdPsi23321*dPhidU212 + dMdPsi23322*dPhidU222 + dMdPsi23323*dPhidU232 + dMdPsi23331*dPhidU312 + dMdPsi23332*dPhidU322 + dMdPsi23333*dPhidU332
    dMdU[I5+18,I2] = dCdU112*dMdC13311 + dCdU122*dMdC13312 + dCdU132*dMdC13313 + dCdU212*dMdC13321 + dCdU222*dMdC13322 + dCdU232*dMdC13323 + dCdU312*dMdC13331 + dCdU322*dMdC13332 + dCdU332*dMdC13333 + dGammadU1112*dMdGamma133111 + dGammadU1122*dMdGamma133112 + dGammadU1132*dMdGamma133113 + dGammadU1212*dMdGamma133121 + dGammadU1222*dMdGamma133122 + dGammadU1232*dMdGamma133123 + dGammadU1312*dMdGamma133131 + dGammadU1322*dMdGamma133132 + dGammadU1332*dMdGamma133133 + dGammadU2112*dMdGamma133211 + dGammadU2122*dMdGamma133212 + dGammadU2132*dMdGamma133213 + dGammadU2212*dMdGamma133221 + dGammadU2222*dMdGamma133222 + dGammadU2232*dMdGamma133223 + dGammadU2312*dMdGamma133231 + dGammadU2322*dMdGamma133232 + dGammadU2332*dMdGamma133233 + dGammadU3112*dMdGamma133311 + dGammadU3122*dMdGamma133312 + dGammadU3132*dMdGamma133313 + dGammadU3212*dMdGamma133321 + dGammadU3222*dMdGamma133322 + dGammadU3232*dMdGamma133323 + dGammadU3312*dMdGamma133331 + dGammadU3322*dMdGamma133332 + dGammadU3332*dMdGamma133333 + dMdPsi13311*dPhidU112 + dMdPsi13312*dPhidU122 + dMdPsi13313*dPhidU132 + dMdPsi13321*dPhidU212 + dMdPsi13322*dPhidU222 + dMdPsi13323*dPhidU232 + dMdPsi13331*dPhidU312 + dMdPsi13332*dPhidU322 + dMdPsi13333*dPhidU332
    dMdU[I6+18,I2] = dCdU112*dMdC12311 + dCdU122*dMdC12312 + dCdU132*dMdC12313 + dCdU212*dMdC12321 + dCdU222*dMdC12322 + dCdU232*dMdC12323 + dCdU312*dMdC12331 + dCdU322*dMdC12332 + dCdU332*dMdC12333 + dGammadU1112*dMdGamma123111 + dGammadU1122*dMdGamma123112 + dGammadU1132*dMdGamma123113 + dGammadU1212*dMdGamma123121 + dGammadU1222*dMdGamma123122 + dGammadU1232*dMdGamma123123 + dGammadU1312*dMdGamma123131 + dGammadU1322*dMdGamma123132 + dGammadU1332*dMdGamma123133 + dGammadU2112*dMdGamma123211 + dGammadU2122*dMdGamma123212 + dGammadU2132*dMdGamma123213 + dGammadU2212*dMdGamma123221 + dGammadU2222*dMdGamma123222 + dGammadU2232*dMdGamma123223 + dGammadU2312*dMdGamma123231 + dGammadU2322*dMdGamma123232 + dGammadU2332*dMdGamma123233 + dGammadU3112*dMdGamma123311 + dGammadU3122*dMdGamma123312 + dGammadU3132*dMdGamma123313 + dGammadU3212*dMdGamma123321 + dGammadU3222*dMdGamma123322 + dGammadU3232*dMdGamma123323 + dGammadU3312*dMdGamma123331 + dGammadU3322*dMdGamma123332 + dGammadU3332*dMdGamma123333 + dMdPsi12311*dPhidU112 + dMdPsi12312*dPhidU122 + dMdPsi12313*dPhidU132 + dMdPsi12321*dPhidU212 + dMdPsi12322*dPhidU222 + dMdPsi12323*dPhidU232 + dMdPsi12331*dPhidU312 + dMdPsi12332*dPhidU322 + dMdPsi12333*dPhidU332
    dMdU[I7+18,I2] = dCdU112*dMdC32311 + dCdU122*dMdC32312 + dCdU132*dMdC32313 + dCdU212*dMdC32321 + dCdU222*dMdC32322 + dCdU232*dMdC32323 + dCdU312*dMdC32331 + dCdU322*dMdC32332 + dCdU332*dMdC32333 + dGammadU1112*dMdGamma323111 + dGammadU1122*dMdGamma323112 + dGammadU1132*dMdGamma323113 + dGammadU1212*dMdGamma323121 + dGammadU1222*dMdGamma323122 + dGammadU1232*dMdGamma323123 + dGammadU1312*dMdGamma323131 + dGammadU1322*dMdGamma323132 + dGammadU1332*dMdGamma323133 + dGammadU2112*dMdGamma323211 + dGammadU2122*dMdGamma323212 + dGammadU2132*dMdGamma323213 + dGammadU2212*dMdGamma323221 + dGammadU2222*dMdGamma323222 + dGammadU2232*dMdGamma323223 + dGammadU2312*dMdGamma323231 + dGammadU2322*dMdGamma323232 + dGammadU2332*dMdGamma323233 + dGammadU3112*dMdGamma323311 + dGammadU3122*dMdGamma323312 + dGammadU3132*dMdGamma323313 + dGammadU3212*dMdGamma323321 + dGammadU3222*dMdGamma323322 + dGammadU3232*dMdGamma323323 + dGammadU3312*dMdGamma323331 + dGammadU3322*dMdGamma323332 + dGammadU3332*dMdGamma323333 + dMdPsi32311*dPhidU112 + dMdPsi32312*dPhidU122 + dMdPsi32313*dPhidU132 + dMdPsi32321*dPhidU212 + dMdPsi32322*dPhidU222 + dMdPsi32323*dPhidU232 + dMdPsi32331*dPhidU312 + dMdPsi32332*dPhidU322 + dMdPsi32333*dPhidU332
    dMdU[I8+18,I2] = dCdU112*dMdC31311 + dCdU122*dMdC31312 + dCdU132*dMdC31313 + dCdU212*dMdC31321 + dCdU222*dMdC31322 + dCdU232*dMdC31323 + dCdU312*dMdC31331 + dCdU322*dMdC31332 + dCdU332*dMdC31333 + dGammadU1112*dMdGamma313111 + dGammadU1122*dMdGamma313112 + dGammadU1132*dMdGamma313113 + dGammadU1212*dMdGamma313121 + dGammadU1222*dMdGamma313122 + dGammadU1232*dMdGamma313123 + dGammadU1312*dMdGamma313131 + dGammadU1322*dMdGamma313132 + dGammadU1332*dMdGamma313133 + dGammadU2112*dMdGamma313211 + dGammadU2122*dMdGamma313212 + dGammadU2132*dMdGamma313213 + dGammadU2212*dMdGamma313221 + dGammadU2222*dMdGamma313222 + dGammadU2232*dMdGamma313223 + dGammadU2312*dMdGamma313231 + dGammadU2322*dMdGamma313232 + dGammadU2332*dMdGamma313233 + dGammadU3112*dMdGamma313311 + dGammadU3122*dMdGamma313312 + dGammadU3132*dMdGamma313313 + dGammadU3212*dMdGamma313321 + dGammadU3222*dMdGamma313322 + dGammadU3232*dMdGamma313323 + dGammadU3312*dMdGamma313331 + dGammadU3322*dMdGamma313332 + dGammadU3332*dMdGamma313333 + dMdPsi31311*dPhidU112 + dMdPsi31312*dPhidU122 + dMdPsi31313*dPhidU132 + dMdPsi31321*dPhidU212 + dMdPsi31322*dPhidU222 + dMdPsi31323*dPhidU232 + dMdPsi31331*dPhidU312 + dMdPsi31332*dPhidU322 + dMdPsi31333*dPhidU332
    dMdU[I9+18,I2] = dCdU112*dMdC21311 + dCdU122*dMdC21312 + dCdU132*dMdC21313 + dCdU212*dMdC21321 + dCdU222*dMdC21322 + dCdU232*dMdC21323 + dCdU312*dMdC21331 + dCdU322*dMdC21332 + dCdU332*dMdC21333 + dGammadU1112*dMdGamma213111 + dGammadU1122*dMdGamma213112 + dGammadU1132*dMdGamma213113 + dGammadU1212*dMdGamma213121 + dGammadU1222*dMdGamma213122 + dGammadU1232*dMdGamma213123 + dGammadU1312*dMdGamma213131 + dGammadU1322*dMdGamma213132 + dGammadU1332*dMdGamma213133 + dGammadU2112*dMdGamma213211 + dGammadU2122*dMdGamma213212 + dGammadU2132*dMdGamma213213 + dGammadU2212*dMdGamma213221 + dGammadU2222*dMdGamma213222 + dGammadU2232*dMdGamma213223 + dGammadU2312*dMdGamma213231 + dGammadU2322*dMdGamma213232 + dGammadU2332*dMdGamma213233 + dGammadU3112*dMdGamma213311 + dGammadU3122*dMdGamma213312 + dGammadU3132*dMdGamma213313 + dGammadU3212*dMdGamma213321 + dGammadU3222*dMdGamma213322 + dGammadU3232*dMdGamma213323 + dGammadU3312*dMdGamma213331 + dGammadU3322*dMdGamma213332 + dGammadU3332*dMdGamma213333 + dMdPsi21311*dPhidU112 + dMdPsi21312*dPhidU122 + dMdPsi21313*dPhidU132 + dMdPsi21321*dPhidU212 + dMdPsi21322*dPhidU222 + dMdPsi21323*dPhidU232 + dMdPsi21331*dPhidU312 + dMdPsi21332*dPhidU322 + dMdPsi21333*dPhidU332

    #Column 3
    dMdU[I1,I3] = dCdU113*dMdC11111 + dCdU123*dMdC11112 + dCdU133*dMdC11113 + dCdU213*dMdC11121 + dCdU223*dMdC11122 + dCdU233*dMdC11123 + dCdU313*dMdC11131 + dCdU323*dMdC11132 + dCdU333*dMdC11133 + dGammadU1113*dMdGamma111111 + dGammadU1123*dMdGamma111112 + dGammadU1133*dMdGamma111113 + dGammadU1213*dMdGamma111121 + dGammadU1223*dMdGamma111122 + dGammadU1233*dMdGamma111123 + dGammadU1313*dMdGamma111131 + dGammadU1323*dMdGamma111132 + dGammadU1333*dMdGamma111133 + dGammadU2113*dMdGamma111211 + dGammadU2123*dMdGamma111212 + dGammadU2133*dMdGamma111213 + dGammadU2213*dMdGamma111221 + dGammadU2223*dMdGamma111222 + dGammadU2233*dMdGamma111223 + dGammadU2313*dMdGamma111231 + dGammadU2323*dMdGamma111232 + dGammadU2333*dMdGamma111233 + dGammadU3113*dMdGamma111311 + dGammadU3123*dMdGamma111312 + dGammadU3133*dMdGamma111313 + dGammadU3213*dMdGamma111321 + dGammadU3223*dMdGamma111322 + dGammadU3233*dMdGamma111323 + dGammadU3313*dMdGamma111331 + dGammadU3323*dMdGamma111332 + dGammadU3333*dMdGamma111333 + dMdPsi11111*dPhidU113 + dMdPsi11112*dPhidU123 + dMdPsi11113*dPhidU133 + dMdPsi11121*dPhidU213 + dMdPsi11122*dPhidU223 + dMdPsi11123*dPhidU233 + dMdPsi11131*dPhidU313 + dMdPsi11132*dPhidU323 + dMdPsi11133*dPhidU333
    dMdU[I2,I3] = dCdU113*dMdC22111 + dCdU123*dMdC22112 + dCdU133*dMdC22113 + dCdU213*dMdC22121 + dCdU223*dMdC22122 + dCdU233*dMdC22123 + dCdU313*dMdC22131 + dCdU323*dMdC22132 + dCdU333*dMdC22133 + dGammadU1113*dMdGamma221111 + dGammadU1123*dMdGamma221112 + dGammadU1133*dMdGamma221113 + dGammadU1213*dMdGamma221121 + dGammadU1223*dMdGamma221122 + dGammadU1233*dMdGamma221123 + dGammadU1313*dMdGamma221131 + dGammadU1323*dMdGamma221132 + dGammadU1333*dMdGamma221133 + dGammadU2113*dMdGamma221211 + dGammadU2123*dMdGamma221212 + dGammadU2133*dMdGamma221213 + dGammadU2213*dMdGamma221221 + dGammadU2223*dMdGamma221222 + dGammadU2233*dMdGamma221223 + dGammadU2313*dMdGamma221231 + dGammadU2323*dMdGamma221232 + dGammadU2333*dMdGamma221233 + dGammadU3113*dMdGamma221311 + dGammadU3123*dMdGamma221312 + dGammadU3133*dMdGamma221313 + dGammadU3213*dMdGamma221321 + dGammadU3223*dMdGamma221322 + dGammadU3233*dMdGamma221323 + dGammadU3313*dMdGamma221331 + dGammadU3323*dMdGamma221332 + dGammadU3333*dMdGamma221333 + dMdPsi22111*dPhidU113 + dMdPsi22112*dPhidU123 + dMdPsi22113*dPhidU133 + dMdPsi22121*dPhidU213 + dMdPsi22122*dPhidU223 + dMdPsi22123*dPhidU233 + dMdPsi22131*dPhidU313 + dMdPsi22132*dPhidU323 + dMdPsi22133*dPhidU333
    dMdU[I3,I3] = dCdU113*dMdC33111 + dCdU123*dMdC33112 + dCdU133*dMdC33113 + dCdU213*dMdC33121 + dCdU223*dMdC33122 + dCdU233*dMdC33123 + dCdU313*dMdC33131 + dCdU323*dMdC33132 + dCdU333*dMdC33133 + dGammadU1113*dMdGamma331111 + dGammadU1123*dMdGamma331112 + dGammadU1133*dMdGamma331113 + dGammadU1213*dMdGamma331121 + dGammadU1223*dMdGamma331122 + dGammadU1233*dMdGamma331123 + dGammadU1313*dMdGamma331131 + dGammadU1323*dMdGamma331132 + dGammadU1333*dMdGamma331133 + dGammadU2113*dMdGamma331211 + dGammadU2123*dMdGamma331212 + dGammadU2133*dMdGamma331213 + dGammadU2213*dMdGamma331221 + dGammadU2223*dMdGamma331222 + dGammadU2233*dMdGamma331223 + dGammadU2313*dMdGamma331231 + dGammadU2323*dMdGamma331232 + dGammadU2333*dMdGamma331233 + dGammadU3113*dMdGamma331311 + dGammadU3123*dMdGamma331312 + dGammadU3133*dMdGamma331313 + dGammadU3213*dMdGamma331321 + dGammadU3223*dMdGamma331322 + dGammadU3233*dMdGamma331323 + dGammadU3313*dMdGamma331331 + dGammadU3323*dMdGamma331332 + dGammadU3333*dMdGamma331333 + dMdPsi33111*dPhidU113 + dMdPsi33112*dPhidU123 + dMdPsi33113*dPhidU133 + dMdPsi33121*dPhidU213 + dMdPsi33122*dPhidU223 + dMdPsi33123*dPhidU233 + dMdPsi33131*dPhidU313 + dMdPsi33132*dPhidU323 + dMdPsi33133*dPhidU333
    dMdU[I4,I3] = dCdU113*dMdC23111 + dCdU123*dMdC23112 + dCdU133*dMdC23113 + dCdU213*dMdC23121 + dCdU223*dMdC23122 + dCdU233*dMdC23123 + dCdU313*dMdC23131 + dCdU323*dMdC23132 + dCdU333*dMdC23133 + dGammadU1113*dMdGamma231111 + dGammadU1123*dMdGamma231112 + dGammadU1133*dMdGamma231113 + dGammadU1213*dMdGamma231121 + dGammadU1223*dMdGamma231122 + dGammadU1233*dMdGamma231123 + dGammadU1313*dMdGamma231131 + dGammadU1323*dMdGamma231132 + dGammadU1333*dMdGamma231133 + dGammadU2113*dMdGamma231211 + dGammadU2123*dMdGamma231212 + dGammadU2133*dMdGamma231213 + dGammadU2213*dMdGamma231221 + dGammadU2223*dMdGamma231222 + dGammadU2233*dMdGamma231223 + dGammadU2313*dMdGamma231231 + dGammadU2323*dMdGamma231232 + dGammadU2333*dMdGamma231233 + dGammadU3113*dMdGamma231311 + dGammadU3123*dMdGamma231312 + dGammadU3133*dMdGamma231313 + dGammadU3213*dMdGamma231321 + dGammadU3223*dMdGamma231322 + dGammadU3233*dMdGamma231323 + dGammadU3313*dMdGamma231331 + dGammadU3323*dMdGamma231332 + dGammadU3333*dMdGamma231333 + dMdPsi23111*dPhidU113 + dMdPsi23112*dPhidU123 + dMdPsi23113*dPhidU133 + dMdPsi23121*dPhidU213 + dMdPsi23122*dPhidU223 + dMdPsi23123*dPhidU233 + dMdPsi23131*dPhidU313 + dMdPsi23132*dPhidU323 + dMdPsi23133*dPhidU333
    dMdU[I5,I3] = dCdU113*dMdC13111 + dCdU123*dMdC13112 + dCdU133*dMdC13113 + dCdU213*dMdC13121 + dCdU223*dMdC13122 + dCdU233*dMdC13123 + dCdU313*dMdC13131 + dCdU323*dMdC13132 + dCdU333*dMdC13133 + dGammadU1113*dMdGamma131111 + dGammadU1123*dMdGamma131112 + dGammadU1133*dMdGamma131113 + dGammadU1213*dMdGamma131121 + dGammadU1223*dMdGamma131122 + dGammadU1233*dMdGamma131123 + dGammadU1313*dMdGamma131131 + dGammadU1323*dMdGamma131132 + dGammadU1333*dMdGamma131133 + dGammadU2113*dMdGamma131211 + dGammadU2123*dMdGamma131212 + dGammadU2133*dMdGamma131213 + dGammadU2213*dMdGamma131221 + dGammadU2223*dMdGamma131222 + dGammadU2233*dMdGamma131223 + dGammadU2313*dMdGamma131231 + dGammadU2323*dMdGamma131232 + dGammadU2333*dMdGamma131233 + dGammadU3113*dMdGamma131311 + dGammadU3123*dMdGamma131312 + dGammadU3133*dMdGamma131313 + dGammadU3213*dMdGamma131321 + dGammadU3223*dMdGamma131322 + dGammadU3233*dMdGamma131323 + dGammadU3313*dMdGamma131331 + dGammadU3323*dMdGamma131332 + dGammadU3333*dMdGamma131333 + dMdPsi13111*dPhidU113 + dMdPsi13112*dPhidU123 + dMdPsi13113*dPhidU133 + dMdPsi13121*dPhidU213 + dMdPsi13122*dPhidU223 + dMdPsi13123*dPhidU233 + dMdPsi13131*dPhidU313 + dMdPsi13132*dPhidU323 + dMdPsi13133*dPhidU333
    dMdU[I6,I3] = dCdU113*dMdC12111 + dCdU123*dMdC12112 + dCdU133*dMdC12113 + dCdU213*dMdC12121 + dCdU223*dMdC12122 + dCdU233*dMdC12123 + dCdU313*dMdC12131 + dCdU323*dMdC12132 + dCdU333*dMdC12133 + dGammadU1113*dMdGamma121111 + dGammadU1123*dMdGamma121112 + dGammadU1133*dMdGamma121113 + dGammadU1213*dMdGamma121121 + dGammadU1223*dMdGamma121122 + dGammadU1233*dMdGamma121123 + dGammadU1313*dMdGamma121131 + dGammadU1323*dMdGamma121132 + dGammadU1333*dMdGamma121133 + dGammadU2113*dMdGamma121211 + dGammadU2123*dMdGamma121212 + dGammadU2133*dMdGamma121213 + dGammadU2213*dMdGamma121221 + dGammadU2223*dMdGamma121222 + dGammadU2233*dMdGamma121223 + dGammadU2313*dMdGamma121231 + dGammadU2323*dMdGamma121232 + dGammadU2333*dMdGamma121233 + dGammadU3113*dMdGamma121311 + dGammadU3123*dMdGamma121312 + dGammadU3133*dMdGamma121313 + dGammadU3213*dMdGamma121321 + dGammadU3223*dMdGamma121322 + dGammadU3233*dMdGamma121323 + dGammadU3313*dMdGamma121331 + dGammadU3323*dMdGamma121332 + dGammadU3333*dMdGamma121333 + dMdPsi12111*dPhidU113 + dMdPsi12112*dPhidU123 + dMdPsi12113*dPhidU133 + dMdPsi12121*dPhidU213 + dMdPsi12122*dPhidU223 + dMdPsi12123*dPhidU233 + dMdPsi12131*dPhidU313 + dMdPsi12132*dPhidU323 + dMdPsi12133*dPhidU333
    dMdU[I7,I3] = dCdU113*dMdC32111 + dCdU123*dMdC32112 + dCdU133*dMdC32113 + dCdU213*dMdC32121 + dCdU223*dMdC32122 + dCdU233*dMdC32123 + dCdU313*dMdC32131 + dCdU323*dMdC32132 + dCdU333*dMdC32133 + dGammadU1113*dMdGamma321111 + dGammadU1123*dMdGamma321112 + dGammadU1133*dMdGamma321113 + dGammadU1213*dMdGamma321121 + dGammadU1223*dMdGamma321122 + dGammadU1233*dMdGamma321123 + dGammadU1313*dMdGamma321131 + dGammadU1323*dMdGamma321132 + dGammadU1333*dMdGamma321133 + dGammadU2113*dMdGamma321211 + dGammadU2123*dMdGamma321212 + dGammadU2133*dMdGamma321213 + dGammadU2213*dMdGamma321221 + dGammadU2223*dMdGamma321222 + dGammadU2233*dMdGamma321223 + dGammadU2313*dMdGamma321231 + dGammadU2323*dMdGamma321232 + dGammadU2333*dMdGamma321233 + dGammadU3113*dMdGamma321311 + dGammadU3123*dMdGamma321312 + dGammadU3133*dMdGamma321313 + dGammadU3213*dMdGamma321321 + dGammadU3223*dMdGamma321322 + dGammadU3233*dMdGamma321323 + dGammadU3313*dMdGamma321331 + dGammadU3323*dMdGamma321332 + dGammadU3333*dMdGamma321333 + dMdPsi32111*dPhidU113 + dMdPsi32112*dPhidU123 + dMdPsi32113*dPhidU133 + dMdPsi32121*dPhidU213 + dMdPsi32122*dPhidU223 + dMdPsi32123*dPhidU233 + dMdPsi32131*dPhidU313 + dMdPsi32132*dPhidU323 + dMdPsi32133*dPhidU333
    dMdU[I8,I3] = dCdU113*dMdC31111 + dCdU123*dMdC31112 + dCdU133*dMdC31113 + dCdU213*dMdC31121 + dCdU223*dMdC31122 + dCdU233*dMdC31123 + dCdU313*dMdC31131 + dCdU323*dMdC31132 + dCdU333*dMdC31133 + dGammadU1113*dMdGamma311111 + dGammadU1123*dMdGamma311112 + dGammadU1133*dMdGamma311113 + dGammadU1213*dMdGamma311121 + dGammadU1223*dMdGamma311122 + dGammadU1233*dMdGamma311123 + dGammadU1313*dMdGamma311131 + dGammadU1323*dMdGamma311132 + dGammadU1333*dMdGamma311133 + dGammadU2113*dMdGamma311211 + dGammadU2123*dMdGamma311212 + dGammadU2133*dMdGamma311213 + dGammadU2213*dMdGamma311221 + dGammadU2223*dMdGamma311222 + dGammadU2233*dMdGamma311223 + dGammadU2313*dMdGamma311231 + dGammadU2323*dMdGamma311232 + dGammadU2333*dMdGamma311233 + dGammadU3113*dMdGamma311311 + dGammadU3123*dMdGamma311312 + dGammadU3133*dMdGamma311313 + dGammadU3213*dMdGamma311321 + dGammadU3223*dMdGamma311322 + dGammadU3233*dMdGamma311323 + dGammadU3313*dMdGamma311331 + dGammadU3323*dMdGamma311332 + dGammadU3333*dMdGamma311333 + dMdPsi31111*dPhidU113 + dMdPsi31112*dPhidU123 + dMdPsi31113*dPhidU133 + dMdPsi31121*dPhidU213 + dMdPsi31122*dPhidU223 + dMdPsi31123*dPhidU233 + dMdPsi31131*dPhidU313 + dMdPsi31132*dPhidU323 + dMdPsi31133*dPhidU333
    dMdU[I9,I3] = dCdU113*dMdC21111 + dCdU123*dMdC21112 + dCdU133*dMdC21113 + dCdU213*dMdC21121 + dCdU223*dMdC21122 + dCdU233*dMdC21123 + dCdU313*dMdC21131 + dCdU323*dMdC21132 + dCdU333*dMdC21133 + dGammadU1113*dMdGamma211111 + dGammadU1123*dMdGamma211112 + dGammadU1133*dMdGamma211113 + dGammadU1213*dMdGamma211121 + dGammadU1223*dMdGamma211122 + dGammadU1233*dMdGamma211123 + dGammadU1313*dMdGamma211131 + dGammadU1323*dMdGamma211132 + dGammadU1333*dMdGamma211133 + dGammadU2113*dMdGamma211211 + dGammadU2123*dMdGamma211212 + dGammadU2133*dMdGamma211213 + dGammadU2213*dMdGamma211221 + dGammadU2223*dMdGamma211222 + dGammadU2233*dMdGamma211223 + dGammadU2313*dMdGamma211231 + dGammadU2323*dMdGamma211232 + dGammadU2333*dMdGamma211233 + dGammadU3113*dMdGamma211311 + dGammadU3123*dMdGamma211312 + dGammadU3133*dMdGamma211313 + dGammadU3213*dMdGamma211321 + dGammadU3223*dMdGamma211322 + dGammadU3233*dMdGamma211323 + dGammadU3313*dMdGamma211331 + dGammadU3323*dMdGamma211332 + dGammadU3333*dMdGamma211333 + dMdPsi21111*dPhidU113 + dMdPsi21112*dPhidU123 + dMdPsi21113*dPhidU133 + dMdPsi21121*dPhidU213 + dMdPsi21122*dPhidU223 + dMdPsi21123*dPhidU233 + dMdPsi21131*dPhidU313 + dMdPsi21132*dPhidU323 + dMdPsi21133*dPhidU333
    dMdU[I1+9,I3] = dCdU113*dMdC11211 + dCdU123*dMdC11212 + dCdU133*dMdC11213 + dCdU213*dMdC11221 + dCdU223*dMdC11222 + dCdU233*dMdC11223 + dCdU313*dMdC11231 + dCdU323*dMdC11232 + dCdU333*dMdC11233 + dGammadU1113*dMdGamma112111 + dGammadU1123*dMdGamma112112 + dGammadU1133*dMdGamma112113 + dGammadU1213*dMdGamma112121 + dGammadU1223*dMdGamma112122 + dGammadU1233*dMdGamma112123 + dGammadU1313*dMdGamma112131 + dGammadU1323*dMdGamma112132 + dGammadU1333*dMdGamma112133 + dGammadU2113*dMdGamma112211 + dGammadU2123*dMdGamma112212 + dGammadU2133*dMdGamma112213 + dGammadU2213*dMdGamma112221 + dGammadU2223*dMdGamma112222 + dGammadU2233*dMdGamma112223 + dGammadU2313*dMdGamma112231 + dGammadU2323*dMdGamma112232 + dGammadU2333*dMdGamma112233 + dGammadU3113*dMdGamma112311 + dGammadU3123*dMdGamma112312 + dGammadU3133*dMdGamma112313 + dGammadU3213*dMdGamma112321 + dGammadU3223*dMdGamma112322 + dGammadU3233*dMdGamma112323 + dGammadU3313*dMdGamma112331 + dGammadU3323*dMdGamma112332 + dGammadU3333*dMdGamma112333 + dMdPsi11211*dPhidU113 + dMdPsi11212*dPhidU123 + dMdPsi11213*dPhidU133 + dMdPsi11221*dPhidU213 + dMdPsi11222*dPhidU223 + dMdPsi11223*dPhidU233 + dMdPsi11231*dPhidU313 + dMdPsi11232*dPhidU323 + dMdPsi11233*dPhidU333
    dMdU[I2+9,I3] = dCdU113*dMdC22211 + dCdU123*dMdC22212 + dCdU133*dMdC22213 + dCdU213*dMdC22221 + dCdU223*dMdC22222 + dCdU233*dMdC22223 + dCdU313*dMdC22231 + dCdU323*dMdC22232 + dCdU333*dMdC22233 + dGammadU1113*dMdGamma222111 + dGammadU1123*dMdGamma222112 + dGammadU1133*dMdGamma222113 + dGammadU1213*dMdGamma222121 + dGammadU1223*dMdGamma222122 + dGammadU1233*dMdGamma222123 + dGammadU1313*dMdGamma222131 + dGammadU1323*dMdGamma222132 + dGammadU1333*dMdGamma222133 + dGammadU2113*dMdGamma222211 + dGammadU2123*dMdGamma222212 + dGammadU2133*dMdGamma222213 + dGammadU2213*dMdGamma222221 + dGammadU2223*dMdGamma222222 + dGammadU2233*dMdGamma222223 + dGammadU2313*dMdGamma222231 + dGammadU2323*dMdGamma222232 + dGammadU2333*dMdGamma222233 + dGammadU3113*dMdGamma222311 + dGammadU3123*dMdGamma222312 + dGammadU3133*dMdGamma222313 + dGammadU3213*dMdGamma222321 + dGammadU3223*dMdGamma222322 + dGammadU3233*dMdGamma222323 + dGammadU3313*dMdGamma222331 + dGammadU3323*dMdGamma222332 + dGammadU3333*dMdGamma222333 + dMdPsi22211*dPhidU113 + dMdPsi22212*dPhidU123 + dMdPsi22213*dPhidU133 + dMdPsi22221*dPhidU213 + dMdPsi22222*dPhidU223 + dMdPsi22223*dPhidU233 + dMdPsi22231*dPhidU313 + dMdPsi22232*dPhidU323 + dMdPsi22233*dPhidU333
    dMdU[I3+9,I3] = dCdU113*dMdC33211 + dCdU123*dMdC33212 + dCdU133*dMdC33213 + dCdU213*dMdC33221 + dCdU223*dMdC33222 + dCdU233*dMdC33223 + dCdU313*dMdC33231 + dCdU323*dMdC33232 + dCdU333*dMdC33233 + dGammadU1113*dMdGamma332111 + dGammadU1123*dMdGamma332112 + dGammadU1133*dMdGamma332113 + dGammadU1213*dMdGamma332121 + dGammadU1223*dMdGamma332122 + dGammadU1233*dMdGamma332123 + dGammadU1313*dMdGamma332131 + dGammadU1323*dMdGamma332132 + dGammadU1333*dMdGamma332133 + dGammadU2113*dMdGamma332211 + dGammadU2123*dMdGamma332212 + dGammadU2133*dMdGamma332213 + dGammadU2213*dMdGamma332221 + dGammadU2223*dMdGamma332222 + dGammadU2233*dMdGamma332223 + dGammadU2313*dMdGamma332231 + dGammadU2323*dMdGamma332232 + dGammadU2333*dMdGamma332233 + dGammadU3113*dMdGamma332311 + dGammadU3123*dMdGamma332312 + dGammadU3133*dMdGamma332313 + dGammadU3213*dMdGamma332321 + dGammadU3223*dMdGamma332322 + dGammadU3233*dMdGamma332323 + dGammadU3313*dMdGamma332331 + dGammadU3323*dMdGamma332332 + dGammadU3333*dMdGamma332333 + dMdPsi33211*dPhidU113 + dMdPsi33212*dPhidU123 + dMdPsi33213*dPhidU133 + dMdPsi33221*dPhidU213 + dMdPsi33222*dPhidU223 + dMdPsi33223*dPhidU233 + dMdPsi33231*dPhidU313 + dMdPsi33232*dPhidU323 + dMdPsi33233*dPhidU333
    dMdU[I4+9,I3] = dCdU113*dMdC23211 + dCdU123*dMdC23212 + dCdU133*dMdC23213 + dCdU213*dMdC23221 + dCdU223*dMdC23222 + dCdU233*dMdC23223 + dCdU313*dMdC23231 + dCdU323*dMdC23232 + dCdU333*dMdC23233 + dGammadU1113*dMdGamma232111 + dGammadU1123*dMdGamma232112 + dGammadU1133*dMdGamma232113 + dGammadU1213*dMdGamma232121 + dGammadU1223*dMdGamma232122 + dGammadU1233*dMdGamma232123 + dGammadU1313*dMdGamma232131 + dGammadU1323*dMdGamma232132 + dGammadU1333*dMdGamma232133 + dGammadU2113*dMdGamma232211 + dGammadU2123*dMdGamma232212 + dGammadU2133*dMdGamma232213 + dGammadU2213*dMdGamma232221 + dGammadU2223*dMdGamma232222 + dGammadU2233*dMdGamma232223 + dGammadU2313*dMdGamma232231 + dGammadU2323*dMdGamma232232 + dGammadU2333*dMdGamma232233 + dGammadU3113*dMdGamma232311 + dGammadU3123*dMdGamma232312 + dGammadU3133*dMdGamma232313 + dGammadU3213*dMdGamma232321 + dGammadU3223*dMdGamma232322 + dGammadU3233*dMdGamma232323 + dGammadU3313*dMdGamma232331 + dGammadU3323*dMdGamma232332 + dGammadU3333*dMdGamma232333 + dMdPsi23211*dPhidU113 + dMdPsi23212*dPhidU123 + dMdPsi23213*dPhidU133 + dMdPsi23221*dPhidU213 + dMdPsi23222*dPhidU223 + dMdPsi23223*dPhidU233 + dMdPsi23231*dPhidU313 + dMdPsi23232*dPhidU323 + dMdPsi23233*dPhidU333
    dMdU[I5+9,I3] = dCdU113*dMdC13211 + dCdU123*dMdC13212 + dCdU133*dMdC13213 + dCdU213*dMdC13221 + dCdU223*dMdC13222 + dCdU233*dMdC13223 + dCdU313*dMdC13231 + dCdU323*dMdC13232 + dCdU333*dMdC13233 + dGammadU1113*dMdGamma132111 + dGammadU1123*dMdGamma132112 + dGammadU1133*dMdGamma132113 + dGammadU1213*dMdGamma132121 + dGammadU1223*dMdGamma132122 + dGammadU1233*dMdGamma132123 + dGammadU1313*dMdGamma132131 + dGammadU1323*dMdGamma132132 + dGammadU1333*dMdGamma132133 + dGammadU2113*dMdGamma132211 + dGammadU2123*dMdGamma132212 + dGammadU2133*dMdGamma132213 + dGammadU2213*dMdGamma132221 + dGammadU2223*dMdGamma132222 + dGammadU2233*dMdGamma132223 + dGammadU2313*dMdGamma132231 + dGammadU2323*dMdGamma132232 + dGammadU2333*dMdGamma132233 + dGammadU3113*dMdGamma132311 + dGammadU3123*dMdGamma132312 + dGammadU3133*dMdGamma132313 + dGammadU3213*dMdGamma132321 + dGammadU3223*dMdGamma132322 + dGammadU3233*dMdGamma132323 + dGammadU3313*dMdGamma132331 + dGammadU3323*dMdGamma132332 + dGammadU3333*dMdGamma132333 + dMdPsi13211*dPhidU113 + dMdPsi13212*dPhidU123 + dMdPsi13213*dPhidU133 + dMdPsi13221*dPhidU213 + dMdPsi13222*dPhidU223 + dMdPsi13223*dPhidU233 + dMdPsi13231*dPhidU313 + dMdPsi13232*dPhidU323 + dMdPsi13233*dPhidU333
    dMdU[I6+9,I3] = dCdU113*dMdC12211 + dCdU123*dMdC12212 + dCdU133*dMdC12213 + dCdU213*dMdC12221 + dCdU223*dMdC12222 + dCdU233*dMdC12223 + dCdU313*dMdC12231 + dCdU323*dMdC12232 + dCdU333*dMdC12233 + dGammadU1113*dMdGamma122111 + dGammadU1123*dMdGamma122112 + dGammadU1133*dMdGamma122113 + dGammadU1213*dMdGamma122121 + dGammadU1223*dMdGamma122122 + dGammadU1233*dMdGamma122123 + dGammadU1313*dMdGamma122131 + dGammadU1323*dMdGamma122132 + dGammadU1333*dMdGamma122133 + dGammadU2113*dMdGamma122211 + dGammadU2123*dMdGamma122212 + dGammadU2133*dMdGamma122213 + dGammadU2213*dMdGamma122221 + dGammadU2223*dMdGamma122222 + dGammadU2233*dMdGamma122223 + dGammadU2313*dMdGamma122231 + dGammadU2323*dMdGamma122232 + dGammadU2333*dMdGamma122233 + dGammadU3113*dMdGamma122311 + dGammadU3123*dMdGamma122312 + dGammadU3133*dMdGamma122313 + dGammadU3213*dMdGamma122321 + dGammadU3223*dMdGamma122322 + dGammadU3233*dMdGamma122323 + dGammadU3313*dMdGamma122331 + dGammadU3323*dMdGamma122332 + dGammadU3333*dMdGamma122333 + dMdPsi12211*dPhidU113 + dMdPsi12212*dPhidU123 + dMdPsi12213*dPhidU133 + dMdPsi12221*dPhidU213 + dMdPsi12222*dPhidU223 + dMdPsi12223*dPhidU233 + dMdPsi12231*dPhidU313 + dMdPsi12232*dPhidU323 + dMdPsi12233*dPhidU333
    dMdU[I7+9,I3] = dCdU113*dMdC32211 + dCdU123*dMdC32212 + dCdU133*dMdC32213 + dCdU213*dMdC32221 + dCdU223*dMdC32222 + dCdU233*dMdC32223 + dCdU313*dMdC32231 + dCdU323*dMdC32232 + dCdU333*dMdC32233 + dGammadU1113*dMdGamma322111 + dGammadU1123*dMdGamma322112 + dGammadU1133*dMdGamma322113 + dGammadU1213*dMdGamma322121 + dGammadU1223*dMdGamma322122 + dGammadU1233*dMdGamma322123 + dGammadU1313*dMdGamma322131 + dGammadU1323*dMdGamma322132 + dGammadU1333*dMdGamma322133 + dGammadU2113*dMdGamma322211 + dGammadU2123*dMdGamma322212 + dGammadU2133*dMdGamma322213 + dGammadU2213*dMdGamma322221 + dGammadU2223*dMdGamma322222 + dGammadU2233*dMdGamma322223 + dGammadU2313*dMdGamma322231 + dGammadU2323*dMdGamma322232 + dGammadU2333*dMdGamma322233 + dGammadU3113*dMdGamma322311 + dGammadU3123*dMdGamma322312 + dGammadU3133*dMdGamma322313 + dGammadU3213*dMdGamma322321 + dGammadU3223*dMdGamma322322 + dGammadU3233*dMdGamma322323 + dGammadU3313*dMdGamma322331 + dGammadU3323*dMdGamma322332 + dGammadU3333*dMdGamma322333 + dMdPsi32211*dPhidU113 + dMdPsi32212*dPhidU123 + dMdPsi32213*dPhidU133 + dMdPsi32221*dPhidU213 + dMdPsi32222*dPhidU223 + dMdPsi32223*dPhidU233 + dMdPsi32231*dPhidU313 + dMdPsi32232*dPhidU323 + dMdPsi32233*dPhidU333
    dMdU[I8+9,I3] = dCdU113*dMdC31211 + dCdU123*dMdC31212 + dCdU133*dMdC31213 + dCdU213*dMdC31221 + dCdU223*dMdC31222 + dCdU233*dMdC31223 + dCdU313*dMdC31231 + dCdU323*dMdC31232 + dCdU333*dMdC31233 + dGammadU1113*dMdGamma312111 + dGammadU1123*dMdGamma312112 + dGammadU1133*dMdGamma312113 + dGammadU1213*dMdGamma312121 + dGammadU1223*dMdGamma312122 + dGammadU1233*dMdGamma312123 + dGammadU1313*dMdGamma312131 + dGammadU1323*dMdGamma312132 + dGammadU1333*dMdGamma312133 + dGammadU2113*dMdGamma312211 + dGammadU2123*dMdGamma312212 + dGammadU2133*dMdGamma312213 + dGammadU2213*dMdGamma312221 + dGammadU2223*dMdGamma312222 + dGammadU2233*dMdGamma312223 + dGammadU2313*dMdGamma312231 + dGammadU2323*dMdGamma312232 + dGammadU2333*dMdGamma312233 + dGammadU3113*dMdGamma312311 + dGammadU3123*dMdGamma312312 + dGammadU3133*dMdGamma312313 + dGammadU3213*dMdGamma312321 + dGammadU3223*dMdGamma312322 + dGammadU3233*dMdGamma312323 + dGammadU3313*dMdGamma312331 + dGammadU3323*dMdGamma312332 + dGammadU3333*dMdGamma312333 + dMdPsi31211*dPhidU113 + dMdPsi31212*dPhidU123 + dMdPsi31213*dPhidU133 + dMdPsi31221*dPhidU213 + dMdPsi31222*dPhidU223 + dMdPsi31223*dPhidU233 + dMdPsi31231*dPhidU313 + dMdPsi31232*dPhidU323 + dMdPsi31233*dPhidU333
    dMdU[I9+9,I3] = dCdU113*dMdC21211 + dCdU123*dMdC21212 + dCdU133*dMdC21213 + dCdU213*dMdC21221 + dCdU223*dMdC21222 + dCdU233*dMdC21223 + dCdU313*dMdC21231 + dCdU323*dMdC21232 + dCdU333*dMdC21233 + dGammadU1113*dMdGamma212111 + dGammadU1123*dMdGamma212112 + dGammadU1133*dMdGamma212113 + dGammadU1213*dMdGamma212121 + dGammadU1223*dMdGamma212122 + dGammadU1233*dMdGamma212123 + dGammadU1313*dMdGamma212131 + dGammadU1323*dMdGamma212132 + dGammadU1333*dMdGamma212133 + dGammadU2113*dMdGamma212211 + dGammadU2123*dMdGamma212212 + dGammadU2133*dMdGamma212213 + dGammadU2213*dMdGamma212221 + dGammadU2223*dMdGamma212222 + dGammadU2233*dMdGamma212223 + dGammadU2313*dMdGamma212231 + dGammadU2323*dMdGamma212232 + dGammadU2333*dMdGamma212233 + dGammadU3113*dMdGamma212311 + dGammadU3123*dMdGamma212312 + dGammadU3133*dMdGamma212313 + dGammadU3213*dMdGamma212321 + dGammadU3223*dMdGamma212322 + dGammadU3233*dMdGamma212323 + dGammadU3313*dMdGamma212331 + dGammadU3323*dMdGamma212332 + dGammadU3333*dMdGamma212333 + dMdPsi21211*dPhidU113 + dMdPsi21212*dPhidU123 + dMdPsi21213*dPhidU133 + dMdPsi21221*dPhidU213 + dMdPsi21222*dPhidU223 + dMdPsi21223*dPhidU233 + dMdPsi21231*dPhidU313 + dMdPsi21232*dPhidU323 + dMdPsi21233*dPhidU333
    dMdU[I1+18,I3] = dCdU113*dMdC11311 + dCdU123*dMdC11312 + dCdU133*dMdC11313 + dCdU213*dMdC11321 + dCdU223*dMdC11322 + dCdU233*dMdC11323 + dCdU313*dMdC11331 + dCdU323*dMdC11332 + dCdU333*dMdC11333 + dGammadU1113*dMdGamma113111 + dGammadU1123*dMdGamma113112 + dGammadU1133*dMdGamma113113 + dGammadU1213*dMdGamma113121 + dGammadU1223*dMdGamma113122 + dGammadU1233*dMdGamma113123 + dGammadU1313*dMdGamma113131 + dGammadU1323*dMdGamma113132 + dGammadU1333*dMdGamma113133 + dGammadU2113*dMdGamma113211 + dGammadU2123*dMdGamma113212 + dGammadU2133*dMdGamma113213 + dGammadU2213*dMdGamma113221 + dGammadU2223*dMdGamma113222 + dGammadU2233*dMdGamma113223 + dGammadU2313*dMdGamma113231 + dGammadU2323*dMdGamma113232 + dGammadU2333*dMdGamma113233 + dGammadU3113*dMdGamma113311 + dGammadU3123*dMdGamma113312 + dGammadU3133*dMdGamma113313 + dGammadU3213*dMdGamma113321 + dGammadU3223*dMdGamma113322 + dGammadU3233*dMdGamma113323 + dGammadU3313*dMdGamma113331 + dGammadU3323*dMdGamma113332 + dGammadU3333*dMdGamma113333 + dMdPsi11311*dPhidU113 + dMdPsi11312*dPhidU123 + dMdPsi11313*dPhidU133 + dMdPsi11321*dPhidU213 + dMdPsi11322*dPhidU223 + dMdPsi11323*dPhidU233 + dMdPsi11331*dPhidU313 + dMdPsi11332*dPhidU323 + dMdPsi11333*dPhidU333
    dMdU[I2+18,I3] = dCdU113*dMdC22311 + dCdU123*dMdC22312 + dCdU133*dMdC22313 + dCdU213*dMdC22321 + dCdU223*dMdC22322 + dCdU233*dMdC22323 + dCdU313*dMdC22331 + dCdU323*dMdC22332 + dCdU333*dMdC22333 + dGammadU1113*dMdGamma223111 + dGammadU1123*dMdGamma223112 + dGammadU1133*dMdGamma223113 + dGammadU1213*dMdGamma223121 + dGammadU1223*dMdGamma223122 + dGammadU1233*dMdGamma223123 + dGammadU1313*dMdGamma223131 + dGammadU1323*dMdGamma223132 + dGammadU1333*dMdGamma223133 + dGammadU2113*dMdGamma223211 + dGammadU2123*dMdGamma223212 + dGammadU2133*dMdGamma223213 + dGammadU2213*dMdGamma223221 + dGammadU2223*dMdGamma223222 + dGammadU2233*dMdGamma223223 + dGammadU2313*dMdGamma223231 + dGammadU2323*dMdGamma223232 + dGammadU2333*dMdGamma223233 + dGammadU3113*dMdGamma223311 + dGammadU3123*dMdGamma223312 + dGammadU3133*dMdGamma223313 + dGammadU3213*dMdGamma223321 + dGammadU3223*dMdGamma223322 + dGammadU3233*dMdGamma223323 + dGammadU3313*dMdGamma223331 + dGammadU3323*dMdGamma223332 + dGammadU3333*dMdGamma223333 + dMdPsi22311*dPhidU113 + dMdPsi22312*dPhidU123 + dMdPsi22313*dPhidU133 + dMdPsi22321*dPhidU213 + dMdPsi22322*dPhidU223 + dMdPsi22323*dPhidU233 + dMdPsi22331*dPhidU313 + dMdPsi22332*dPhidU323 + dMdPsi22333*dPhidU333
    dMdU[I3+18,I3] = dCdU113*dMdC33311 + dCdU123*dMdC33312 + dCdU133*dMdC33313 + dCdU213*dMdC33321 + dCdU223*dMdC33322 + dCdU233*dMdC33323 + dCdU313*dMdC33331 + dCdU323*dMdC33332 + dCdU333*dMdC33333 + dGammadU1113*dMdGamma333111 + dGammadU1123*dMdGamma333112 + dGammadU1133*dMdGamma333113 + dGammadU1213*dMdGamma333121 + dGammadU1223*dMdGamma333122 + dGammadU1233*dMdGamma333123 + dGammadU1313*dMdGamma333131 + dGammadU1323*dMdGamma333132 + dGammadU1333*dMdGamma333133 + dGammadU2113*dMdGamma333211 + dGammadU2123*dMdGamma333212 + dGammadU2133*dMdGamma333213 + dGammadU2213*dMdGamma333221 + dGammadU2223*dMdGamma333222 + dGammadU2233*dMdGamma333223 + dGammadU2313*dMdGamma333231 + dGammadU2323*dMdGamma333232 + dGammadU2333*dMdGamma333233 + dGammadU3113*dMdGamma333311 + dGammadU3123*dMdGamma333312 + dGammadU3133*dMdGamma333313 + dGammadU3213*dMdGamma333321 + dGammadU3223*dMdGamma333322 + dGammadU3233*dMdGamma333323 + dGammadU3313*dMdGamma333331 + dGammadU3323*dMdGamma333332 + dGammadU3333*dMdGamma333333 + dMdPsi33311*dPhidU113 + dMdPsi33312*dPhidU123 + dMdPsi33313*dPhidU133 + dMdPsi33321*dPhidU213 + dMdPsi33322*dPhidU223 + dMdPsi33323*dPhidU233 + dMdPsi33331*dPhidU313 + dMdPsi33332*dPhidU323 + dMdPsi33333*dPhidU333
    dMdU[I4+18,I3] = dCdU113*dMdC23311 + dCdU123*dMdC23312 + dCdU133*dMdC23313 + dCdU213*dMdC23321 + dCdU223*dMdC23322 + dCdU233*dMdC23323 + dCdU313*dMdC23331 + dCdU323*dMdC23332 + dCdU333*dMdC23333 + dGammadU1113*dMdGamma233111 + dGammadU1123*dMdGamma233112 + dGammadU1133*dMdGamma233113 + dGammadU1213*dMdGamma233121 + dGammadU1223*dMdGamma233122 + dGammadU1233*dMdGamma233123 + dGammadU1313*dMdGamma233131 + dGammadU1323*dMdGamma233132 + dGammadU1333*dMdGamma233133 + dGammadU2113*dMdGamma233211 + dGammadU2123*dMdGamma233212 + dGammadU2133*dMdGamma233213 + dGammadU2213*dMdGamma233221 + dGammadU2223*dMdGamma233222 + dGammadU2233*dMdGamma233223 + dGammadU2313*dMdGamma233231 + dGammadU2323*dMdGamma233232 + dGammadU2333*dMdGamma233233 + dGammadU3113*dMdGamma233311 + dGammadU3123*dMdGamma233312 + dGammadU3133*dMdGamma233313 + dGammadU3213*dMdGamma233321 + dGammadU3223*dMdGamma233322 + dGammadU3233*dMdGamma233323 + dGammadU3313*dMdGamma233331 + dGammadU3323*dMdGamma233332 + dGammadU3333*dMdGamma233333 + dMdPsi23311*dPhidU113 + dMdPsi23312*dPhidU123 + dMdPsi23313*dPhidU133 + dMdPsi23321*dPhidU213 + dMdPsi23322*dPhidU223 + dMdPsi23323*dPhidU233 + dMdPsi23331*dPhidU313 + dMdPsi23332*dPhidU323 + dMdPsi23333*dPhidU333
    dMdU[I5+18,I3] = dCdU113*dMdC13311 + dCdU123*dMdC13312 + dCdU133*dMdC13313 + dCdU213*dMdC13321 + dCdU223*dMdC13322 + dCdU233*dMdC13323 + dCdU313*dMdC13331 + dCdU323*dMdC13332 + dCdU333*dMdC13333 + dGammadU1113*dMdGamma133111 + dGammadU1123*dMdGamma133112 + dGammadU1133*dMdGamma133113 + dGammadU1213*dMdGamma133121 + dGammadU1223*dMdGamma133122 + dGammadU1233*dMdGamma133123 + dGammadU1313*dMdGamma133131 + dGammadU1323*dMdGamma133132 + dGammadU1333*dMdGamma133133 + dGammadU2113*dMdGamma133211 + dGammadU2123*dMdGamma133212 + dGammadU2133*dMdGamma133213 + dGammadU2213*dMdGamma133221 + dGammadU2223*dMdGamma133222 + dGammadU2233*dMdGamma133223 + dGammadU2313*dMdGamma133231 + dGammadU2323*dMdGamma133232 + dGammadU2333*dMdGamma133233 + dGammadU3113*dMdGamma133311 + dGammadU3123*dMdGamma133312 + dGammadU3133*dMdGamma133313 + dGammadU3213*dMdGamma133321 + dGammadU3223*dMdGamma133322 + dGammadU3233*dMdGamma133323 + dGammadU3313*dMdGamma133331 + dGammadU3323*dMdGamma133332 + dGammadU3333*dMdGamma133333 + dMdPsi13311*dPhidU113 + dMdPsi13312*dPhidU123 + dMdPsi13313*dPhidU133 + dMdPsi13321*dPhidU213 + dMdPsi13322*dPhidU223 + dMdPsi13323*dPhidU233 + dMdPsi13331*dPhidU313 + dMdPsi13332*dPhidU323 + dMdPsi13333*dPhidU333
    dMdU[I6+18,I3] = dCdU113*dMdC12311 + dCdU123*dMdC12312 + dCdU133*dMdC12313 + dCdU213*dMdC12321 + dCdU223*dMdC12322 + dCdU233*dMdC12323 + dCdU313*dMdC12331 + dCdU323*dMdC12332 + dCdU333*dMdC12333 + dGammadU1113*dMdGamma123111 + dGammadU1123*dMdGamma123112 + dGammadU1133*dMdGamma123113 + dGammadU1213*dMdGamma123121 + dGammadU1223*dMdGamma123122 + dGammadU1233*dMdGamma123123 + dGammadU1313*dMdGamma123131 + dGammadU1323*dMdGamma123132 + dGammadU1333*dMdGamma123133 + dGammadU2113*dMdGamma123211 + dGammadU2123*dMdGamma123212 + dGammadU2133*dMdGamma123213 + dGammadU2213*dMdGamma123221 + dGammadU2223*dMdGamma123222 + dGammadU2233*dMdGamma123223 + dGammadU2313*dMdGamma123231 + dGammadU2323*dMdGamma123232 + dGammadU2333*dMdGamma123233 + dGammadU3113*dMdGamma123311 + dGammadU3123*dMdGamma123312 + dGammadU3133*dMdGamma123313 + dGammadU3213*dMdGamma123321 + dGammadU3223*dMdGamma123322 + dGammadU3233*dMdGamma123323 + dGammadU3313*dMdGamma123331 + dGammadU3323*dMdGamma123332 + dGammadU3333*dMdGamma123333 + dMdPsi12311*dPhidU113 + dMdPsi12312*dPhidU123 + dMdPsi12313*dPhidU133 + dMdPsi12321*dPhidU213 + dMdPsi12322*dPhidU223 + dMdPsi12323*dPhidU233 + dMdPsi12331*dPhidU313 + dMdPsi12332*dPhidU323 + dMdPsi12333*dPhidU333
    dMdU[I7+18,I3] = dCdU113*dMdC32311 + dCdU123*dMdC32312 + dCdU133*dMdC32313 + dCdU213*dMdC32321 + dCdU223*dMdC32322 + dCdU233*dMdC32323 + dCdU313*dMdC32331 + dCdU323*dMdC32332 + dCdU333*dMdC32333 + dGammadU1113*dMdGamma323111 + dGammadU1123*dMdGamma323112 + dGammadU1133*dMdGamma323113 + dGammadU1213*dMdGamma323121 + dGammadU1223*dMdGamma323122 + dGammadU1233*dMdGamma323123 + dGammadU1313*dMdGamma323131 + dGammadU1323*dMdGamma323132 + dGammadU1333*dMdGamma323133 + dGammadU2113*dMdGamma323211 + dGammadU2123*dMdGamma323212 + dGammadU2133*dMdGamma323213 + dGammadU2213*dMdGamma323221 + dGammadU2223*dMdGamma323222 + dGammadU2233*dMdGamma323223 + dGammadU2313*dMdGamma323231 + dGammadU2323*dMdGamma323232 + dGammadU2333*dMdGamma323233 + dGammadU3113*dMdGamma323311 + dGammadU3123*dMdGamma323312 + dGammadU3133*dMdGamma323313 + dGammadU3213*dMdGamma323321 + dGammadU3223*dMdGamma323322 + dGammadU3233*dMdGamma323323 + dGammadU3313*dMdGamma323331 + dGammadU3323*dMdGamma323332 + dGammadU3333*dMdGamma323333 + dMdPsi32311*dPhidU113 + dMdPsi32312*dPhidU123 + dMdPsi32313*dPhidU133 + dMdPsi32321*dPhidU213 + dMdPsi32322*dPhidU223 + dMdPsi32323*dPhidU233 + dMdPsi32331*dPhidU313 + dMdPsi32332*dPhidU323 + dMdPsi32333*dPhidU333
    dMdU[I8+18,I3] = dCdU113*dMdC31311 + dCdU123*dMdC31312 + dCdU133*dMdC31313 + dCdU213*dMdC31321 + dCdU223*dMdC31322 + dCdU233*dMdC31323 + dCdU313*dMdC31331 + dCdU323*dMdC31332 + dCdU333*dMdC31333 + dGammadU1113*dMdGamma313111 + dGammadU1123*dMdGamma313112 + dGammadU1133*dMdGamma313113 + dGammadU1213*dMdGamma313121 + dGammadU1223*dMdGamma313122 + dGammadU1233*dMdGamma313123 + dGammadU1313*dMdGamma313131 + dGammadU1323*dMdGamma313132 + dGammadU1333*dMdGamma313133 + dGammadU2113*dMdGamma313211 + dGammadU2123*dMdGamma313212 + dGammadU2133*dMdGamma313213 + dGammadU2213*dMdGamma313221 + dGammadU2223*dMdGamma313222 + dGammadU2233*dMdGamma313223 + dGammadU2313*dMdGamma313231 + dGammadU2323*dMdGamma313232 + dGammadU2333*dMdGamma313233 + dGammadU3113*dMdGamma313311 + dGammadU3123*dMdGamma313312 + dGammadU3133*dMdGamma313313 + dGammadU3213*dMdGamma313321 + dGammadU3223*dMdGamma313322 + dGammadU3233*dMdGamma313323 + dGammadU3313*dMdGamma313331 + dGammadU3323*dMdGamma313332 + dGammadU3333*dMdGamma313333 + dMdPsi31311*dPhidU113 + dMdPsi31312*dPhidU123 + dMdPsi31313*dPhidU133 + dMdPsi31321*dPhidU213 + dMdPsi31322*dPhidU223 + dMdPsi31323*dPhidU233 + dMdPsi31331*dPhidU313 + dMdPsi31332*dPhidU323 + dMdPsi31333*dPhidU333
    dMdU[I9+18,I3] = dCdU113*dMdC21311 + dCdU123*dMdC21312 + dCdU133*dMdC21313 + dCdU213*dMdC21321 + dCdU223*dMdC21322 + dCdU233*dMdC21323 + dCdU313*dMdC21331 + dCdU323*dMdC21332 + dCdU333*dMdC21333 + dGammadU1113*dMdGamma213111 + dGammadU1123*dMdGamma213112 + dGammadU1133*dMdGamma213113 + dGammadU1213*dMdGamma213121 + dGammadU1223*dMdGamma213122 + dGammadU1233*dMdGamma213123 + dGammadU1313*dMdGamma213131 + dGammadU1323*dMdGamma213132 + dGammadU1333*dMdGamma213133 + dGammadU2113*dMdGamma213211 + dGammadU2123*dMdGamma213212 + dGammadU2133*dMdGamma213213 + dGammadU2213*dMdGamma213221 + dGammadU2223*dMdGamma213222 + dGammadU2233*dMdGamma213223 + dGammadU2313*dMdGamma213231 + dGammadU2323*dMdGamma213232 + dGammadU2333*dMdGamma213233 + dGammadU3113*dMdGamma213311 + dGammadU3123*dMdGamma213312 + dGammadU3133*dMdGamma213313 + dGammadU3213*dMdGamma213321 + dGammadU3223*dMdGamma213322 + dGammadU3233*dMdGamma213323 + dGammadU3313*dMdGamma213331 + dGammadU3323*dMdGamma213332 + dGammadU3333*dMdGamma213333 + dMdPsi21311*dPhidU113 + dMdPsi21312*dPhidU123 + dMdPsi21313*dPhidU133 + dMdPsi21321*dPhidU213 + dMdPsi21322*dPhidU223 + dMdPsi21323*dPhidU233 + dMdPsi21331*dPhidU313 + dMdPsi21332*dPhidU323 + dMdPsi21333*dPhidU333

    #Column 4
    dMdU[I1,I4] = dGammadU1114*dMdGamma111111 + dGammadU1124*dMdGamma111112 + dGammadU1134*dMdGamma111113 + dGammadU2114*dMdGamma111211 + dGammadU2124*dMdGamma111212 + dGammadU2134*dMdGamma111213 + dGammadU3114*dMdGamma111311 + dGammadU3124*dMdGamma111312 + dGammadU3134*dMdGamma111313 + dMdPsi11111*dPhidU114 + dMdPsi11121*dPhidU214 + dMdPsi11131*dPhidU314
    dMdU[I2,I4] = dGammadU1114*dMdGamma221111 + dGammadU1124*dMdGamma221112 + dGammadU1134*dMdGamma221113 + dGammadU2114*dMdGamma221211 + dGammadU2124*dMdGamma221212 + dGammadU2134*dMdGamma221213 + dGammadU3114*dMdGamma221311 + dGammadU3124*dMdGamma221312 + dGammadU3134*dMdGamma221313 + dMdPsi22111*dPhidU114 + dMdPsi22121*dPhidU214 + dMdPsi22131*dPhidU314
    dMdU[I3,I4] = dGammadU1114*dMdGamma331111 + dGammadU1124*dMdGamma331112 + dGammadU1134*dMdGamma331113 + dGammadU2114*dMdGamma331211 + dGammadU2124*dMdGamma331212 + dGammadU2134*dMdGamma331213 + dGammadU3114*dMdGamma331311 + dGammadU3124*dMdGamma331312 + dGammadU3134*dMdGamma331313 + dMdPsi33111*dPhidU114 + dMdPsi33121*dPhidU214 + dMdPsi33131*dPhidU314
    dMdU[I4,I4] = dGammadU1114*dMdGamma231111 + dGammadU1124*dMdGamma231112 + dGammadU1134*dMdGamma231113 + dGammadU2114*dMdGamma231211 + dGammadU2124*dMdGamma231212 + dGammadU2134*dMdGamma231213 + dGammadU3114*dMdGamma231311 + dGammadU3124*dMdGamma231312 + dGammadU3134*dMdGamma231313 + dMdPsi23111*dPhidU114 + dMdPsi23121*dPhidU214 + dMdPsi23131*dPhidU314
    dMdU[I5,I4] = dGammadU1114*dMdGamma131111 + dGammadU1124*dMdGamma131112 + dGammadU1134*dMdGamma131113 + dGammadU2114*dMdGamma131211 + dGammadU2124*dMdGamma131212 + dGammadU2134*dMdGamma131213 + dGammadU3114*dMdGamma131311 + dGammadU3124*dMdGamma131312 + dGammadU3134*dMdGamma131313 + dMdPsi13111*dPhidU114 + dMdPsi13121*dPhidU214 + dMdPsi13131*dPhidU314
    dMdU[I6,I4] = dGammadU1114*dMdGamma121111 + dGammadU1124*dMdGamma121112 + dGammadU1134*dMdGamma121113 + dGammadU2114*dMdGamma121211 + dGammadU2124*dMdGamma121212 + dGammadU2134*dMdGamma121213 + dGammadU3114*dMdGamma121311 + dGammadU3124*dMdGamma121312 + dGammadU3134*dMdGamma121313 + dMdPsi12111*dPhidU114 + dMdPsi12121*dPhidU214 + dMdPsi12131*dPhidU314
    dMdU[I7,I4] = dGammadU1114*dMdGamma321111 + dGammadU1124*dMdGamma321112 + dGammadU1134*dMdGamma321113 + dGammadU2114*dMdGamma321211 + dGammadU2124*dMdGamma321212 + dGammadU2134*dMdGamma321213 + dGammadU3114*dMdGamma321311 + dGammadU3124*dMdGamma321312 + dGammadU3134*dMdGamma321313 + dMdPsi32111*dPhidU114 + dMdPsi32121*dPhidU214 + dMdPsi32131*dPhidU314
    dMdU[I8,I4] = dGammadU1114*dMdGamma311111 + dGammadU1124*dMdGamma311112 + dGammadU1134*dMdGamma311113 + dGammadU2114*dMdGamma311211 + dGammadU2124*dMdGamma311212 + dGammadU2134*dMdGamma311213 + dGammadU3114*dMdGamma311311 + dGammadU3124*dMdGamma311312 + dGammadU3134*dMdGamma311313 + dMdPsi31111*dPhidU114 + dMdPsi31121*dPhidU214 + dMdPsi31131*dPhidU314
    dMdU[I9,I4] = dGammadU1114*dMdGamma211111 + dGammadU1124*dMdGamma211112 + dGammadU1134*dMdGamma211113 + dGammadU2114*dMdGamma211211 + dGammadU2124*dMdGamma211212 + dGammadU2134*dMdGamma211213 + dGammadU3114*dMdGamma211311 + dGammadU3124*dMdGamma211312 + dGammadU3134*dMdGamma211313 + dMdPsi21111*dPhidU114 + dMdPsi21121*dPhidU214 + dMdPsi21131*dPhidU314
    dMdU[I1+9,I4] = dGammadU1114*dMdGamma112111 + dGammadU1124*dMdGamma112112 + dGammadU1134*dMdGamma112113 + dGammadU2114*dMdGamma112211 + dGammadU2124*dMdGamma112212 + dGammadU2134*dMdGamma112213 + dGammadU3114*dMdGamma112311 + dGammadU3124*dMdGamma112312 + dGammadU3134*dMdGamma112313 + dMdPsi11211*dPhidU114 + dMdPsi11221*dPhidU214 + dMdPsi11231*dPhidU314
    dMdU[I2+9,I4] = dGammadU1114*dMdGamma222111 + dGammadU1124*dMdGamma222112 + dGammadU1134*dMdGamma222113 + dGammadU2114*dMdGamma222211 + dGammadU2124*dMdGamma222212 + dGammadU2134*dMdGamma222213 + dGammadU3114*dMdGamma222311 + dGammadU3124*dMdGamma222312 + dGammadU3134*dMdGamma222313 + dMdPsi22211*dPhidU114 + dMdPsi22221*dPhidU214 + dMdPsi22231*dPhidU314
    dMdU[I3+9,I4] = dGammadU1114*dMdGamma332111 + dGammadU1124*dMdGamma332112 + dGammadU1134*dMdGamma332113 + dGammadU2114*dMdGamma332211 + dGammadU2124*dMdGamma332212 + dGammadU2134*dMdGamma332213 + dGammadU3114*dMdGamma332311 + dGammadU3124*dMdGamma332312 + dGammadU3134*dMdGamma332313 + dMdPsi33211*dPhidU114 + dMdPsi33221*dPhidU214 + dMdPsi33231*dPhidU314
    dMdU[I4+9,I4] = dGammadU1114*dMdGamma232111 + dGammadU1124*dMdGamma232112 + dGammadU1134*dMdGamma232113 + dGammadU2114*dMdGamma232211 + dGammadU2124*dMdGamma232212 + dGammadU2134*dMdGamma232213 + dGammadU3114*dMdGamma232311 + dGammadU3124*dMdGamma232312 + dGammadU3134*dMdGamma232313 + dMdPsi23211*dPhidU114 + dMdPsi23221*dPhidU214 + dMdPsi23231*dPhidU314
    dMdU[I5+9,I4] = dGammadU1114*dMdGamma132111 + dGammadU1124*dMdGamma132112 + dGammadU1134*dMdGamma132113 + dGammadU2114*dMdGamma132211 + dGammadU2124*dMdGamma132212 + dGammadU2134*dMdGamma132213 + dGammadU3114*dMdGamma132311 + dGammadU3124*dMdGamma132312 + dGammadU3134*dMdGamma132313 + dMdPsi13211*dPhidU114 + dMdPsi13221*dPhidU214 + dMdPsi13231*dPhidU314
    dMdU[I6+9,I4] = dGammadU1114*dMdGamma122111 + dGammadU1124*dMdGamma122112 + dGammadU1134*dMdGamma122113 + dGammadU2114*dMdGamma122211 + dGammadU2124*dMdGamma122212 + dGammadU2134*dMdGamma122213 + dGammadU3114*dMdGamma122311 + dGammadU3124*dMdGamma122312 + dGammadU3134*dMdGamma122313 + dMdPsi12211*dPhidU114 + dMdPsi12221*dPhidU214 + dMdPsi12231*dPhidU314
    dMdU[I7+9,I4] = dGammadU1114*dMdGamma322111 + dGammadU1124*dMdGamma322112 + dGammadU1134*dMdGamma322113 + dGammadU2114*dMdGamma322211 + dGammadU2124*dMdGamma322212 + dGammadU2134*dMdGamma322213 + dGammadU3114*dMdGamma322311 + dGammadU3124*dMdGamma322312 + dGammadU3134*dMdGamma322313 + dMdPsi32211*dPhidU114 + dMdPsi32221*dPhidU214 + dMdPsi32231*dPhidU314
    dMdU[I8+9,I4] = dGammadU1114*dMdGamma312111 + dGammadU1124*dMdGamma312112 + dGammadU1134*dMdGamma312113 + dGammadU2114*dMdGamma312211 + dGammadU2124*dMdGamma312212 + dGammadU2134*dMdGamma312213 + dGammadU3114*dMdGamma312311 + dGammadU3124*dMdGamma312312 + dGammadU3134*dMdGamma312313 + dMdPsi31211*dPhidU114 + dMdPsi31221*dPhidU214 + dMdPsi31231*dPhidU314
    dMdU[I9+9,I4] = dGammadU1114*dMdGamma212111 + dGammadU1124*dMdGamma212112 + dGammadU1134*dMdGamma212113 + dGammadU2114*dMdGamma212211 + dGammadU2124*dMdGamma212212 + dGammadU2134*dMdGamma212213 + dGammadU3114*dMdGamma212311 + dGammadU3124*dMdGamma212312 + dGammadU3134*dMdGamma212313 + dMdPsi21211*dPhidU114 + dMdPsi21221*dPhidU214 + dMdPsi21231*dPhidU314
    dMdU[I1+18,I4] = dGammadU1114*dMdGamma113111 + dGammadU1124*dMdGamma113112 + dGammadU1134*dMdGamma113113 + dGammadU2114*dMdGamma113211 + dGammadU2124*dMdGamma113212 + dGammadU2134*dMdGamma113213 + dGammadU3114*dMdGamma113311 + dGammadU3124*dMdGamma113312 + dGammadU3134*dMdGamma113313 + dMdPsi11311*dPhidU114 + dMdPsi11321*dPhidU214 + dMdPsi11331*dPhidU314
    dMdU[I2+18,I4] = dGammadU1114*dMdGamma223111 + dGammadU1124*dMdGamma223112 + dGammadU1134*dMdGamma223113 + dGammadU2114*dMdGamma223211 + dGammadU2124*dMdGamma223212 + dGammadU2134*dMdGamma223213 + dGammadU3114*dMdGamma223311 + dGammadU3124*dMdGamma223312 + dGammadU3134*dMdGamma223313 + dMdPsi22311*dPhidU114 + dMdPsi22321*dPhidU214 + dMdPsi22331*dPhidU314
    dMdU[I3+18,I4] = dGammadU1114*dMdGamma333111 + dGammadU1124*dMdGamma333112 + dGammadU1134*dMdGamma333113 + dGammadU2114*dMdGamma333211 + dGammadU2124*dMdGamma333212 + dGammadU2134*dMdGamma333213 + dGammadU3114*dMdGamma333311 + dGammadU3124*dMdGamma333312 + dGammadU3134*dMdGamma333313 + dMdPsi33311*dPhidU114 + dMdPsi33321*dPhidU214 + dMdPsi33331*dPhidU314
    dMdU[I4+18,I4] = dGammadU1114*dMdGamma233111 + dGammadU1124*dMdGamma233112 + dGammadU1134*dMdGamma233113 + dGammadU2114*dMdGamma233211 + dGammadU2124*dMdGamma233212 + dGammadU2134*dMdGamma233213 + dGammadU3114*dMdGamma233311 + dGammadU3124*dMdGamma233312 + dGammadU3134*dMdGamma233313 + dMdPsi23311*dPhidU114 + dMdPsi23321*dPhidU214 + dMdPsi23331*dPhidU314
    dMdU[I5+18,I4] = dGammadU1114*dMdGamma133111 + dGammadU1124*dMdGamma133112 + dGammadU1134*dMdGamma133113 + dGammadU2114*dMdGamma133211 + dGammadU2124*dMdGamma133212 + dGammadU2134*dMdGamma133213 + dGammadU3114*dMdGamma133311 + dGammadU3124*dMdGamma133312 + dGammadU3134*dMdGamma133313 + dMdPsi13311*dPhidU114 + dMdPsi13321*dPhidU214 + dMdPsi13331*dPhidU314
    dMdU[I6+18,I4] = dGammadU1114*dMdGamma123111 + dGammadU1124*dMdGamma123112 + dGammadU1134*dMdGamma123113 + dGammadU2114*dMdGamma123211 + dGammadU2124*dMdGamma123212 + dGammadU2134*dMdGamma123213 + dGammadU3114*dMdGamma123311 + dGammadU3124*dMdGamma123312 + dGammadU3134*dMdGamma123313 + dMdPsi12311*dPhidU114 + dMdPsi12321*dPhidU214 + dMdPsi12331*dPhidU314
    dMdU[I7+18,I4] = dGammadU1114*dMdGamma323111 + dGammadU1124*dMdGamma323112 + dGammadU1134*dMdGamma323113 + dGammadU2114*dMdGamma323211 + dGammadU2124*dMdGamma323212 + dGammadU2134*dMdGamma323213 + dGammadU3114*dMdGamma323311 + dGammadU3124*dMdGamma323312 + dGammadU3134*dMdGamma323313 + dMdPsi32311*dPhidU114 + dMdPsi32321*dPhidU214 + dMdPsi32331*dPhidU314
    dMdU[I8+18,I4] = dGammadU1114*dMdGamma313111 + dGammadU1124*dMdGamma313112 + dGammadU1134*dMdGamma313113 + dGammadU2114*dMdGamma313211 + dGammadU2124*dMdGamma313212 + dGammadU2134*dMdGamma313213 + dGammadU3114*dMdGamma313311 + dGammadU3124*dMdGamma313312 + dGammadU3134*dMdGamma313313 + dMdPsi31311*dPhidU114 + dMdPsi31321*dPhidU214 + dMdPsi31331*dPhidU314
    dMdU[I9+18,I4] = dGammadU1114*dMdGamma213111 + dGammadU1124*dMdGamma213112 + dGammadU1134*dMdGamma213113 + dGammadU2114*dMdGamma213211 + dGammadU2124*dMdGamma213212 + dGammadU2134*dMdGamma213213 + dGammadU3114*dMdGamma213311 + dGammadU3124*dMdGamma213312 + dGammadU3134*dMdGamma213313 + dMdPsi21311*dPhidU114 + dMdPsi21321*dPhidU214 + dMdPsi21331*dPhidU314

    #Column 5
    dMdU[I1,I5] = dGammadU1215*dMdGamma111121 + dGammadU1225*dMdGamma111122 + dGammadU1235*dMdGamma111123 + dGammadU2215*dMdGamma111221 + dGammadU2225*dMdGamma111222 + dGammadU2235*dMdGamma111223 + dGammadU3215*dMdGamma111321 + dGammadU3225*dMdGamma111322 + dGammadU3235*dMdGamma111323 + dMdPsi11112*dPhidU125 + dMdPsi11122*dPhidU225 + dMdPsi11132*dPhidU325
    dMdU[I2,I5] = dGammadU1215*dMdGamma221121 + dGammadU1225*dMdGamma221122 + dGammadU1235*dMdGamma221123 + dGammadU2215*dMdGamma221221 + dGammadU2225*dMdGamma221222 + dGammadU2235*dMdGamma221223 + dGammadU3215*dMdGamma221321 + dGammadU3225*dMdGamma221322 + dGammadU3235*dMdGamma221323 + dMdPsi22112*dPhidU125 + dMdPsi22122*dPhidU225 + dMdPsi22132*dPhidU325
    dMdU[I3,I5] = dGammadU1215*dMdGamma331121 + dGammadU1225*dMdGamma331122 + dGammadU1235*dMdGamma331123 + dGammadU2215*dMdGamma331221 + dGammadU2225*dMdGamma331222 + dGammadU2235*dMdGamma331223 + dGammadU3215*dMdGamma331321 + dGammadU3225*dMdGamma331322 + dGammadU3235*dMdGamma331323 + dMdPsi33112*dPhidU125 + dMdPsi33122*dPhidU225 + dMdPsi33132*dPhidU325
    dMdU[I4,I5] = dGammadU1215*dMdGamma231121 + dGammadU1225*dMdGamma231122 + dGammadU1235*dMdGamma231123 + dGammadU2215*dMdGamma231221 + dGammadU2225*dMdGamma231222 + dGammadU2235*dMdGamma231223 + dGammadU3215*dMdGamma231321 + dGammadU3225*dMdGamma231322 + dGammadU3235*dMdGamma231323 + dMdPsi23112*dPhidU125 + dMdPsi23122*dPhidU225 + dMdPsi23132*dPhidU325
    dMdU[I5,I5] = dGammadU1215*dMdGamma131121 + dGammadU1225*dMdGamma131122 + dGammadU1235*dMdGamma131123 + dGammadU2215*dMdGamma131221 + dGammadU2225*dMdGamma131222 + dGammadU2235*dMdGamma131223 + dGammadU3215*dMdGamma131321 + dGammadU3225*dMdGamma131322 + dGammadU3235*dMdGamma131323 + dMdPsi13112*dPhidU125 + dMdPsi13122*dPhidU225 + dMdPsi13132*dPhidU325
    dMdU[I6,I5] = dGammadU1215*dMdGamma121121 + dGammadU1225*dMdGamma121122 + dGammadU1235*dMdGamma121123 + dGammadU2215*dMdGamma121221 + dGammadU2225*dMdGamma121222 + dGammadU2235*dMdGamma121223 + dGammadU3215*dMdGamma121321 + dGammadU3225*dMdGamma121322 + dGammadU3235*dMdGamma121323 + dMdPsi12112*dPhidU125 + dMdPsi12122*dPhidU225 + dMdPsi12132*dPhidU325
    dMdU[I7,I5] = dGammadU1215*dMdGamma321121 + dGammadU1225*dMdGamma321122 + dGammadU1235*dMdGamma321123 + dGammadU2215*dMdGamma321221 + dGammadU2225*dMdGamma321222 + dGammadU2235*dMdGamma321223 + dGammadU3215*dMdGamma321321 + dGammadU3225*dMdGamma321322 + dGammadU3235*dMdGamma321323 + dMdPsi32112*dPhidU125 + dMdPsi32122*dPhidU225 + dMdPsi32132*dPhidU325
    dMdU[I8,I5] = dGammadU1215*dMdGamma311121 + dGammadU1225*dMdGamma311122 + dGammadU1235*dMdGamma311123 + dGammadU2215*dMdGamma311221 + dGammadU2225*dMdGamma311222 + dGammadU2235*dMdGamma311223 + dGammadU3215*dMdGamma311321 + dGammadU3225*dMdGamma311322 + dGammadU3235*dMdGamma311323 + dMdPsi31112*dPhidU125 + dMdPsi31122*dPhidU225 + dMdPsi31132*dPhidU325
    dMdU[I9,I5] = dGammadU1215*dMdGamma211121 + dGammadU1225*dMdGamma211122 + dGammadU1235*dMdGamma211123 + dGammadU2215*dMdGamma211221 + dGammadU2225*dMdGamma211222 + dGammadU2235*dMdGamma211223 + dGammadU3215*dMdGamma211321 + dGammadU3225*dMdGamma211322 + dGammadU3235*dMdGamma211323 + dMdPsi21112*dPhidU125 + dMdPsi21122*dPhidU225 + dMdPsi21132*dPhidU325
    dMdU[I1+9,I5] = dGammadU1215*dMdGamma112121 + dGammadU1225*dMdGamma112122 + dGammadU1235*dMdGamma112123 + dGammadU2215*dMdGamma112221 + dGammadU2225*dMdGamma112222 + dGammadU2235*dMdGamma112223 + dGammadU3215*dMdGamma112321 + dGammadU3225*dMdGamma112322 + dGammadU3235*dMdGamma112323 + dMdPsi11212*dPhidU125 + dMdPsi11222*dPhidU225 + dMdPsi11232*dPhidU325
    dMdU[I2+9,I5] = dGammadU1215*dMdGamma222121 + dGammadU1225*dMdGamma222122 + dGammadU1235*dMdGamma222123 + dGammadU2215*dMdGamma222221 + dGammadU2225*dMdGamma222222 + dGammadU2235*dMdGamma222223 + dGammadU3215*dMdGamma222321 + dGammadU3225*dMdGamma222322 + dGammadU3235*dMdGamma222323 + dMdPsi22212*dPhidU125 + dMdPsi22222*dPhidU225 + dMdPsi22232*dPhidU325
    dMdU[I3+9,I5] = dGammadU1215*dMdGamma332121 + dGammadU1225*dMdGamma332122 + dGammadU1235*dMdGamma332123 + dGammadU2215*dMdGamma332221 + dGammadU2225*dMdGamma332222 + dGammadU2235*dMdGamma332223 + dGammadU3215*dMdGamma332321 + dGammadU3225*dMdGamma332322 + dGammadU3235*dMdGamma332323 + dMdPsi33212*dPhidU125 + dMdPsi33222*dPhidU225 + dMdPsi33232*dPhidU325
    dMdU[I4+9,I5] = dGammadU1215*dMdGamma232121 + dGammadU1225*dMdGamma232122 + dGammadU1235*dMdGamma232123 + dGammadU2215*dMdGamma232221 + dGammadU2225*dMdGamma232222 + dGammadU2235*dMdGamma232223 + dGammadU3215*dMdGamma232321 + dGammadU3225*dMdGamma232322 + dGammadU3235*dMdGamma232323 + dMdPsi23212*dPhidU125 + dMdPsi23222*dPhidU225 + dMdPsi23232*dPhidU325
    dMdU[I5+9,I5] = dGammadU1215*dMdGamma132121 + dGammadU1225*dMdGamma132122 + dGammadU1235*dMdGamma132123 + dGammadU2215*dMdGamma132221 + dGammadU2225*dMdGamma132222 + dGammadU2235*dMdGamma132223 + dGammadU3215*dMdGamma132321 + dGammadU3225*dMdGamma132322 + dGammadU3235*dMdGamma132323 + dMdPsi13212*dPhidU125 + dMdPsi13222*dPhidU225 + dMdPsi13232*dPhidU325
    dMdU[I6+9,I5] = dGammadU1215*dMdGamma122121 + dGammadU1225*dMdGamma122122 + dGammadU1235*dMdGamma122123 + dGammadU2215*dMdGamma122221 + dGammadU2225*dMdGamma122222 + dGammadU2235*dMdGamma122223 + dGammadU3215*dMdGamma122321 + dGammadU3225*dMdGamma122322 + dGammadU3235*dMdGamma122323 + dMdPsi12212*dPhidU125 + dMdPsi12222*dPhidU225 + dMdPsi12232*dPhidU325
    dMdU[I7+9,I5] = dGammadU1215*dMdGamma322121 + dGammadU1225*dMdGamma322122 + dGammadU1235*dMdGamma322123 + dGammadU2215*dMdGamma322221 + dGammadU2225*dMdGamma322222 + dGammadU2235*dMdGamma322223 + dGammadU3215*dMdGamma322321 + dGammadU3225*dMdGamma322322 + dGammadU3235*dMdGamma322323 + dMdPsi32212*dPhidU125 + dMdPsi32222*dPhidU225 + dMdPsi32232*dPhidU325
    dMdU[I8+9,I5] = dGammadU1215*dMdGamma312121 + dGammadU1225*dMdGamma312122 + dGammadU1235*dMdGamma312123 + dGammadU2215*dMdGamma312221 + dGammadU2225*dMdGamma312222 + dGammadU2235*dMdGamma312223 + dGammadU3215*dMdGamma312321 + dGammadU3225*dMdGamma312322 + dGammadU3235*dMdGamma312323 + dMdPsi31212*dPhidU125 + dMdPsi31222*dPhidU225 + dMdPsi31232*dPhidU325
    dMdU[I9+9,I5] = dGammadU1215*dMdGamma212121 + dGammadU1225*dMdGamma212122 + dGammadU1235*dMdGamma212123 + dGammadU2215*dMdGamma212221 + dGammadU2225*dMdGamma212222 + dGammadU2235*dMdGamma212223 + dGammadU3215*dMdGamma212321 + dGammadU3225*dMdGamma212322 + dGammadU3235*dMdGamma212323 + dMdPsi21212*dPhidU125 + dMdPsi21222*dPhidU225 + dMdPsi21232*dPhidU325
    dMdU[I1+18,I5] = dGammadU1215*dMdGamma113121 + dGammadU1225*dMdGamma113122 + dGammadU1235*dMdGamma113123 + dGammadU2215*dMdGamma113221 + dGammadU2225*dMdGamma113222 + dGammadU2235*dMdGamma113223 + dGammadU3215*dMdGamma113321 + dGammadU3225*dMdGamma113322 + dGammadU3235*dMdGamma113323 + dMdPsi11312*dPhidU125 + dMdPsi11322*dPhidU225 + dMdPsi11332*dPhidU325
    dMdU[I2+18,I5] = dGammadU1215*dMdGamma223121 + dGammadU1225*dMdGamma223122 + dGammadU1235*dMdGamma223123 + dGammadU2215*dMdGamma223221 + dGammadU2225*dMdGamma223222 + dGammadU2235*dMdGamma223223 + dGammadU3215*dMdGamma223321 + dGammadU3225*dMdGamma223322 + dGammadU3235*dMdGamma223323 + dMdPsi22312*dPhidU125 + dMdPsi22322*dPhidU225 + dMdPsi22332*dPhidU325
    dMdU[I3+18,I5] = dGammadU1215*dMdGamma333121 + dGammadU1225*dMdGamma333122 + dGammadU1235*dMdGamma333123 + dGammadU2215*dMdGamma333221 + dGammadU2225*dMdGamma333222 + dGammadU2235*dMdGamma333223 + dGammadU3215*dMdGamma333321 + dGammadU3225*dMdGamma333322 + dGammadU3235*dMdGamma333323 + dMdPsi33312*dPhidU125 + dMdPsi33322*dPhidU225 + dMdPsi33332*dPhidU325
    dMdU[I4+18,I5] = dGammadU1215*dMdGamma233121 + dGammadU1225*dMdGamma233122 + dGammadU1235*dMdGamma233123 + dGammadU2215*dMdGamma233221 + dGammadU2225*dMdGamma233222 + dGammadU2235*dMdGamma233223 + dGammadU3215*dMdGamma233321 + dGammadU3225*dMdGamma233322 + dGammadU3235*dMdGamma233323 + dMdPsi23312*dPhidU125 + dMdPsi23322*dPhidU225 + dMdPsi23332*dPhidU325
    dMdU[I5+18,I5] = dGammadU1215*dMdGamma133121 + dGammadU1225*dMdGamma133122 + dGammadU1235*dMdGamma133123 + dGammadU2215*dMdGamma133221 + dGammadU2225*dMdGamma133222 + dGammadU2235*dMdGamma133223 + dGammadU3215*dMdGamma133321 + dGammadU3225*dMdGamma133322 + dGammadU3235*dMdGamma133323 + dMdPsi13312*dPhidU125 + dMdPsi13322*dPhidU225 + dMdPsi13332*dPhidU325
    dMdU[I6+18,I5] = dGammadU1215*dMdGamma123121 + dGammadU1225*dMdGamma123122 + dGammadU1235*dMdGamma123123 + dGammadU2215*dMdGamma123221 + dGammadU2225*dMdGamma123222 + dGammadU2235*dMdGamma123223 + dGammadU3215*dMdGamma123321 + dGammadU3225*dMdGamma123322 + dGammadU3235*dMdGamma123323 + dMdPsi12312*dPhidU125 + dMdPsi12322*dPhidU225 + dMdPsi12332*dPhidU325
    dMdU[I7+18,I5] = dGammadU1215*dMdGamma323121 + dGammadU1225*dMdGamma323122 + dGammadU1235*dMdGamma323123 + dGammadU2215*dMdGamma323221 + dGammadU2225*dMdGamma323222 + dGammadU2235*dMdGamma323223 + dGammadU3215*dMdGamma323321 + dGammadU3225*dMdGamma323322 + dGammadU3235*dMdGamma323323 + dMdPsi32312*dPhidU125 + dMdPsi32322*dPhidU225 + dMdPsi32332*dPhidU325
    dMdU[I8+18,I5] = dGammadU1215*dMdGamma313121 + dGammadU1225*dMdGamma313122 + dGammadU1235*dMdGamma313123 + dGammadU2215*dMdGamma313221 + dGammadU2225*dMdGamma313222 + dGammadU2235*dMdGamma313223 + dGammadU3215*dMdGamma313321 + dGammadU3225*dMdGamma313322 + dGammadU3235*dMdGamma313323 + dMdPsi31312*dPhidU125 + dMdPsi31322*dPhidU225 + dMdPsi31332*dPhidU325
    dMdU[I9+18,I5] = dGammadU1215*dMdGamma213121 + dGammadU1225*dMdGamma213122 + dGammadU1235*dMdGamma213123 + dGammadU2215*dMdGamma213221 + dGammadU2225*dMdGamma213222 + dGammadU2235*dMdGamma213223 + dGammadU3215*dMdGamma213321 + dGammadU3225*dMdGamma213322 + dGammadU3235*dMdGamma213323 + dMdPsi21312*dPhidU125 + dMdPsi21322*dPhidU225 + dMdPsi21332*dPhidU325

    #Column 6
    dMdU[I1,I6] = dGammadU1316*dMdGamma111131 + dGammadU1326*dMdGamma111132 + dGammadU1336*dMdGamma111133 + dGammadU2316*dMdGamma111231 + dGammadU2326*dMdGamma111232 + dGammadU2336*dMdGamma111233 + dGammadU3316*dMdGamma111331 + dGammadU3326*dMdGamma111332 + dGammadU3336*dMdGamma111333 + dMdPsi11113*dPhidU136 + dMdPsi11123*dPhidU236 + dMdPsi11133*dPhidU336
    dMdU[I2,I6] = dGammadU1316*dMdGamma221131 + dGammadU1326*dMdGamma221132 + dGammadU1336*dMdGamma221133 + dGammadU2316*dMdGamma221231 + dGammadU2326*dMdGamma221232 + dGammadU2336*dMdGamma221233 + dGammadU3316*dMdGamma221331 + dGammadU3326*dMdGamma221332 + dGammadU3336*dMdGamma221333 + dMdPsi22113*dPhidU136 + dMdPsi22123*dPhidU236 + dMdPsi22133*dPhidU336
    dMdU[I3,I6] = dGammadU1316*dMdGamma331131 + dGammadU1326*dMdGamma331132 + dGammadU1336*dMdGamma331133 + dGammadU2316*dMdGamma331231 + dGammadU2326*dMdGamma331232 + dGammadU2336*dMdGamma331233 + dGammadU3316*dMdGamma331331 + dGammadU3326*dMdGamma331332 + dGammadU3336*dMdGamma331333 + dMdPsi33113*dPhidU136 + dMdPsi33123*dPhidU236 + dMdPsi33133*dPhidU336
    dMdU[I4,I6] = dGammadU1316*dMdGamma231131 + dGammadU1326*dMdGamma231132 + dGammadU1336*dMdGamma231133 + dGammadU2316*dMdGamma231231 + dGammadU2326*dMdGamma231232 + dGammadU2336*dMdGamma231233 + dGammadU3316*dMdGamma231331 + dGammadU3326*dMdGamma231332 + dGammadU3336*dMdGamma231333 + dMdPsi23113*dPhidU136 + dMdPsi23123*dPhidU236 + dMdPsi23133*dPhidU336
    dMdU[I5,I6] = dGammadU1316*dMdGamma131131 + dGammadU1326*dMdGamma131132 + dGammadU1336*dMdGamma131133 + dGammadU2316*dMdGamma131231 + dGammadU2326*dMdGamma131232 + dGammadU2336*dMdGamma131233 + dGammadU3316*dMdGamma131331 + dGammadU3326*dMdGamma131332 + dGammadU3336*dMdGamma131333 + dMdPsi13113*dPhidU136 + dMdPsi13123*dPhidU236 + dMdPsi13133*dPhidU336
    dMdU[I6,I6] = dGammadU1316*dMdGamma121131 + dGammadU1326*dMdGamma121132 + dGammadU1336*dMdGamma121133 + dGammadU2316*dMdGamma121231 + dGammadU2326*dMdGamma121232 + dGammadU2336*dMdGamma121233 + dGammadU3316*dMdGamma121331 + dGammadU3326*dMdGamma121332 + dGammadU3336*dMdGamma121333 + dMdPsi12113*dPhidU136 + dMdPsi12123*dPhidU236 + dMdPsi12133*dPhidU336
    dMdU[I7,I6] = dGammadU1316*dMdGamma321131 + dGammadU1326*dMdGamma321132 + dGammadU1336*dMdGamma321133 + dGammadU2316*dMdGamma321231 + dGammadU2326*dMdGamma321232 + dGammadU2336*dMdGamma321233 + dGammadU3316*dMdGamma321331 + dGammadU3326*dMdGamma321332 + dGammadU3336*dMdGamma321333 + dMdPsi32113*dPhidU136 + dMdPsi32123*dPhidU236 + dMdPsi32133*dPhidU336
    dMdU[I8,I6] = dGammadU1316*dMdGamma311131 + dGammadU1326*dMdGamma311132 + dGammadU1336*dMdGamma311133 + dGammadU2316*dMdGamma311231 + dGammadU2326*dMdGamma311232 + dGammadU2336*dMdGamma311233 + dGammadU3316*dMdGamma311331 + dGammadU3326*dMdGamma311332 + dGammadU3336*dMdGamma311333 + dMdPsi31113*dPhidU136 + dMdPsi31123*dPhidU236 + dMdPsi31133*dPhidU336
    dMdU[I9,I6] = dGammadU1316*dMdGamma211131 + dGammadU1326*dMdGamma211132 + dGammadU1336*dMdGamma211133 + dGammadU2316*dMdGamma211231 + dGammadU2326*dMdGamma211232 + dGammadU2336*dMdGamma211233 + dGammadU3316*dMdGamma211331 + dGammadU3326*dMdGamma211332 + dGammadU3336*dMdGamma211333 + dMdPsi21113*dPhidU136 + dMdPsi21123*dPhidU236 + dMdPsi21133*dPhidU336
    dMdU[I1+9,I6] = dGammadU1316*dMdGamma112131 + dGammadU1326*dMdGamma112132 + dGammadU1336*dMdGamma112133 + dGammadU2316*dMdGamma112231 + dGammadU2326*dMdGamma112232 + dGammadU2336*dMdGamma112233 + dGammadU3316*dMdGamma112331 + dGammadU3326*dMdGamma112332 + dGammadU3336*dMdGamma112333 + dMdPsi11213*dPhidU136 + dMdPsi11223*dPhidU236 + dMdPsi11233*dPhidU336
    dMdU[I2+9,I6] = dGammadU1316*dMdGamma222131 + dGammadU1326*dMdGamma222132 + dGammadU1336*dMdGamma222133 + dGammadU2316*dMdGamma222231 + dGammadU2326*dMdGamma222232 + dGammadU2336*dMdGamma222233 + dGammadU3316*dMdGamma222331 + dGammadU3326*dMdGamma222332 + dGammadU3336*dMdGamma222333 + dMdPsi22213*dPhidU136 + dMdPsi22223*dPhidU236 + dMdPsi22233*dPhidU336
    dMdU[I3+9,I6] = dGammadU1316*dMdGamma332131 + dGammadU1326*dMdGamma332132 + dGammadU1336*dMdGamma332133 + dGammadU2316*dMdGamma332231 + dGammadU2326*dMdGamma332232 + dGammadU2336*dMdGamma332233 + dGammadU3316*dMdGamma332331 + dGammadU3326*dMdGamma332332 + dGammadU3336*dMdGamma332333 + dMdPsi33213*dPhidU136 + dMdPsi33223*dPhidU236 + dMdPsi33233*dPhidU336
    dMdU[I4+9,I6] = dGammadU1316*dMdGamma232131 + dGammadU1326*dMdGamma232132 + dGammadU1336*dMdGamma232133 + dGammadU2316*dMdGamma232231 + dGammadU2326*dMdGamma232232 + dGammadU2336*dMdGamma232233 + dGammadU3316*dMdGamma232331 + dGammadU3326*dMdGamma232332 + dGammadU3336*dMdGamma232333 + dMdPsi23213*dPhidU136 + dMdPsi23223*dPhidU236 + dMdPsi23233*dPhidU336
    dMdU[I5+9,I6] = dGammadU1316*dMdGamma132131 + dGammadU1326*dMdGamma132132 + dGammadU1336*dMdGamma132133 + dGammadU2316*dMdGamma132231 + dGammadU2326*dMdGamma132232 + dGammadU2336*dMdGamma132233 + dGammadU3316*dMdGamma132331 + dGammadU3326*dMdGamma132332 + dGammadU3336*dMdGamma132333 + dMdPsi13213*dPhidU136 + dMdPsi13223*dPhidU236 + dMdPsi13233*dPhidU336
    dMdU[I6+9,I6] = dGammadU1316*dMdGamma122131 + dGammadU1326*dMdGamma122132 + dGammadU1336*dMdGamma122133 + dGammadU2316*dMdGamma122231 + dGammadU2326*dMdGamma122232 + dGammadU2336*dMdGamma122233 + dGammadU3316*dMdGamma122331 + dGammadU3326*dMdGamma122332 + dGammadU3336*dMdGamma122333 + dMdPsi12213*dPhidU136 + dMdPsi12223*dPhidU236 + dMdPsi12233*dPhidU336
    dMdU[I7+9,I6] = dGammadU1316*dMdGamma322131 + dGammadU1326*dMdGamma322132 + dGammadU1336*dMdGamma322133 + dGammadU2316*dMdGamma322231 + dGammadU2326*dMdGamma322232 + dGammadU2336*dMdGamma322233 + dGammadU3316*dMdGamma322331 + dGammadU3326*dMdGamma322332 + dGammadU3336*dMdGamma322333 + dMdPsi32213*dPhidU136 + dMdPsi32223*dPhidU236 + dMdPsi32233*dPhidU336
    dMdU[I8+9,I6] = dGammadU1316*dMdGamma312131 + dGammadU1326*dMdGamma312132 + dGammadU1336*dMdGamma312133 + dGammadU2316*dMdGamma312231 + dGammadU2326*dMdGamma312232 + dGammadU2336*dMdGamma312233 + dGammadU3316*dMdGamma312331 + dGammadU3326*dMdGamma312332 + dGammadU3336*dMdGamma312333 + dMdPsi31213*dPhidU136 + dMdPsi31223*dPhidU236 + dMdPsi31233*dPhidU336
    dMdU[I9+9,I6] = dGammadU1316*dMdGamma212131 + dGammadU1326*dMdGamma212132 + dGammadU1336*dMdGamma212133 + dGammadU2316*dMdGamma212231 + dGammadU2326*dMdGamma212232 + dGammadU2336*dMdGamma212233 + dGammadU3316*dMdGamma212331 + dGammadU3326*dMdGamma212332 + dGammadU3336*dMdGamma212333 + dMdPsi21213*dPhidU136 + dMdPsi21223*dPhidU236 + dMdPsi21233*dPhidU336
    dMdU[I1+18,I6] = dGammadU1316*dMdGamma113131 + dGammadU1326*dMdGamma113132 + dGammadU1336*dMdGamma113133 + dGammadU2316*dMdGamma113231 + dGammadU2326*dMdGamma113232 + dGammadU2336*dMdGamma113233 + dGammadU3316*dMdGamma113331 + dGammadU3326*dMdGamma113332 + dGammadU3336*dMdGamma113333 + dMdPsi11313*dPhidU136 + dMdPsi11323*dPhidU236 + dMdPsi11333*dPhidU336
    dMdU[I2+18,I6] = dGammadU1316*dMdGamma223131 + dGammadU1326*dMdGamma223132 + dGammadU1336*dMdGamma223133 + dGammadU2316*dMdGamma223231 + dGammadU2326*dMdGamma223232 + dGammadU2336*dMdGamma223233 + dGammadU3316*dMdGamma223331 + dGammadU3326*dMdGamma223332 + dGammadU3336*dMdGamma223333 + dMdPsi22313*dPhidU136 + dMdPsi22323*dPhidU236 + dMdPsi22333*dPhidU336
    dMdU[I3+18,I6] = dGammadU1316*dMdGamma333131 + dGammadU1326*dMdGamma333132 + dGammadU1336*dMdGamma333133 + dGammadU2316*dMdGamma333231 + dGammadU2326*dMdGamma333232 + dGammadU2336*dMdGamma333233 + dGammadU3316*dMdGamma333331 + dGammadU3326*dMdGamma333332 + dGammadU3336*dMdGamma333333 + dMdPsi33313*dPhidU136 + dMdPsi33323*dPhidU236 + dMdPsi33333*dPhidU336
    dMdU[I4+18,I6] = dGammadU1316*dMdGamma233131 + dGammadU1326*dMdGamma233132 + dGammadU1336*dMdGamma233133 + dGammadU2316*dMdGamma233231 + dGammadU2326*dMdGamma233232 + dGammadU2336*dMdGamma233233 + dGammadU3316*dMdGamma233331 + dGammadU3326*dMdGamma233332 + dGammadU3336*dMdGamma233333 + dMdPsi23313*dPhidU136 + dMdPsi23323*dPhidU236 + dMdPsi23333*dPhidU336
    dMdU[I5+18,I6] = dGammadU1316*dMdGamma133131 + dGammadU1326*dMdGamma133132 + dGammadU1336*dMdGamma133133 + dGammadU2316*dMdGamma133231 + dGammadU2326*dMdGamma133232 + dGammadU2336*dMdGamma133233 + dGammadU3316*dMdGamma133331 + dGammadU3326*dMdGamma133332 + dGammadU3336*dMdGamma133333 + dMdPsi13313*dPhidU136 + dMdPsi13323*dPhidU236 + dMdPsi13333*dPhidU336
    dMdU[I6+18,I6] = dGammadU1316*dMdGamma123131 + dGammadU1326*dMdGamma123132 + dGammadU1336*dMdGamma123133 + dGammadU2316*dMdGamma123231 + dGammadU2326*dMdGamma123232 + dGammadU2336*dMdGamma123233 + dGammadU3316*dMdGamma123331 + dGammadU3326*dMdGamma123332 + dGammadU3336*dMdGamma123333 + dMdPsi12313*dPhidU136 + dMdPsi12323*dPhidU236 + dMdPsi12333*dPhidU336
    dMdU[I7+18,I6] = dGammadU1316*dMdGamma323131 + dGammadU1326*dMdGamma323132 + dGammadU1336*dMdGamma323133 + dGammadU2316*dMdGamma323231 + dGammadU2326*dMdGamma323232 + dGammadU2336*dMdGamma323233 + dGammadU3316*dMdGamma323331 + dGammadU3326*dMdGamma323332 + dGammadU3336*dMdGamma323333 + dMdPsi32313*dPhidU136 + dMdPsi32323*dPhidU236 + dMdPsi32333*dPhidU336
    dMdU[I8+18,I6] = dGammadU1316*dMdGamma313131 + dGammadU1326*dMdGamma313132 + dGammadU1336*dMdGamma313133 + dGammadU2316*dMdGamma313231 + dGammadU2326*dMdGamma313232 + dGammadU2336*dMdGamma313233 + dGammadU3316*dMdGamma313331 + dGammadU3326*dMdGamma313332 + dGammadU3336*dMdGamma313333 + dMdPsi31313*dPhidU136 + dMdPsi31323*dPhidU236 + dMdPsi31333*dPhidU336
    dMdU[I9+18,I6] = dGammadU1316*dMdGamma213131 + dGammadU1326*dMdGamma213132 + dGammadU1336*dMdGamma213133 + dGammadU2316*dMdGamma213231 + dGammadU2326*dMdGamma213232 + dGammadU2336*dMdGamma213233 + dGammadU3316*dMdGamma213331 + dGammadU3326*dMdGamma213332 + dGammadU3336*dMdGamma213333 + dMdPsi21313*dPhidU136 + dMdPsi21323*dPhidU236 + dMdPsi21333*dPhidU336

    #Column 7
    dMdU[I1,I7] = dGammadU1317*dMdGamma111131 + dGammadU1327*dMdGamma111132 + dGammadU1337*dMdGamma111133 + dGammadU2317*dMdGamma111231 + dGammadU2327*dMdGamma111232 + dGammadU2337*dMdGamma111233 + dGammadU3317*dMdGamma111331 + dGammadU3327*dMdGamma111332 + dGammadU3337*dMdGamma111333 + dMdPsi11113*dPhidU137 + dMdPsi11123*dPhidU237 + dMdPsi11133*dPhidU337
    dMdU[I2,I7] = dGammadU1317*dMdGamma221131 + dGammadU1327*dMdGamma221132 + dGammadU1337*dMdGamma221133 + dGammadU2317*dMdGamma221231 + dGammadU2327*dMdGamma221232 + dGammadU2337*dMdGamma221233 + dGammadU3317*dMdGamma221331 + dGammadU3327*dMdGamma221332 + dGammadU3337*dMdGamma221333 + dMdPsi22113*dPhidU137 + dMdPsi22123*dPhidU237 + dMdPsi22133*dPhidU337
    dMdU[I3,I7] = dGammadU1317*dMdGamma331131 + dGammadU1327*dMdGamma331132 + dGammadU1337*dMdGamma331133 + dGammadU2317*dMdGamma331231 + dGammadU2327*dMdGamma331232 + dGammadU2337*dMdGamma331233 + dGammadU3317*dMdGamma331331 + dGammadU3327*dMdGamma331332 + dGammadU3337*dMdGamma331333 + dMdPsi33113*dPhidU137 + dMdPsi33123*dPhidU237 + dMdPsi33133*dPhidU337
    dMdU[I4,I7] = dGammadU1317*dMdGamma231131 + dGammadU1327*dMdGamma231132 + dGammadU1337*dMdGamma231133 + dGammadU2317*dMdGamma231231 + dGammadU2327*dMdGamma231232 + dGammadU2337*dMdGamma231233 + dGammadU3317*dMdGamma231331 + dGammadU3327*dMdGamma231332 + dGammadU3337*dMdGamma231333 + dMdPsi23113*dPhidU137 + dMdPsi23123*dPhidU237 + dMdPsi23133*dPhidU337
    dMdU[I5,I7] = dGammadU1317*dMdGamma131131 + dGammadU1327*dMdGamma131132 + dGammadU1337*dMdGamma131133 + dGammadU2317*dMdGamma131231 + dGammadU2327*dMdGamma131232 + dGammadU2337*dMdGamma131233 + dGammadU3317*dMdGamma131331 + dGammadU3327*dMdGamma131332 + dGammadU3337*dMdGamma131333 + dMdPsi13113*dPhidU137 + dMdPsi13123*dPhidU237 + dMdPsi13133*dPhidU337
    dMdU[I6,I7] = dGammadU1317*dMdGamma121131 + dGammadU1327*dMdGamma121132 + dGammadU1337*dMdGamma121133 + dGammadU2317*dMdGamma121231 + dGammadU2327*dMdGamma121232 + dGammadU2337*dMdGamma121233 + dGammadU3317*dMdGamma121331 + dGammadU3327*dMdGamma121332 + dGammadU3337*dMdGamma121333 + dMdPsi12113*dPhidU137 + dMdPsi12123*dPhidU237 + dMdPsi12133*dPhidU337
    dMdU[I7,I7] = dGammadU1317*dMdGamma321131 + dGammadU1327*dMdGamma321132 + dGammadU1337*dMdGamma321133 + dGammadU2317*dMdGamma321231 + dGammadU2327*dMdGamma321232 + dGammadU2337*dMdGamma321233 + dGammadU3317*dMdGamma321331 + dGammadU3327*dMdGamma321332 + dGammadU3337*dMdGamma321333 + dMdPsi32113*dPhidU137 + dMdPsi32123*dPhidU237 + dMdPsi32133*dPhidU337
    dMdU[I8,I7] = dGammadU1317*dMdGamma311131 + dGammadU1327*dMdGamma311132 + dGammadU1337*dMdGamma311133 + dGammadU2317*dMdGamma311231 + dGammadU2327*dMdGamma311232 + dGammadU2337*dMdGamma311233 + dGammadU3317*dMdGamma311331 + dGammadU3327*dMdGamma311332 + dGammadU3337*dMdGamma311333 + dMdPsi31113*dPhidU137 + dMdPsi31123*dPhidU237 + dMdPsi31133*dPhidU337
    dMdU[I9,I7] = dGammadU1317*dMdGamma211131 + dGammadU1327*dMdGamma211132 + dGammadU1337*dMdGamma211133 + dGammadU2317*dMdGamma211231 + dGammadU2327*dMdGamma211232 + dGammadU2337*dMdGamma211233 + dGammadU3317*dMdGamma211331 + dGammadU3327*dMdGamma211332 + dGammadU3337*dMdGamma211333 + dMdPsi21113*dPhidU137 + dMdPsi21123*dPhidU237 + dMdPsi21133*dPhidU337
    dMdU[I1+9,I7] = dGammadU1317*dMdGamma112131 + dGammadU1327*dMdGamma112132 + dGammadU1337*dMdGamma112133 + dGammadU2317*dMdGamma112231 + dGammadU2327*dMdGamma112232 + dGammadU2337*dMdGamma112233 + dGammadU3317*dMdGamma112331 + dGammadU3327*dMdGamma112332 + dGammadU3337*dMdGamma112333 + dMdPsi11213*dPhidU137 + dMdPsi11223*dPhidU237 + dMdPsi11233*dPhidU337
    dMdU[I2+9,I7] = dGammadU1317*dMdGamma222131 + dGammadU1327*dMdGamma222132 + dGammadU1337*dMdGamma222133 + dGammadU2317*dMdGamma222231 + dGammadU2327*dMdGamma222232 + dGammadU2337*dMdGamma222233 + dGammadU3317*dMdGamma222331 + dGammadU3327*dMdGamma222332 + dGammadU3337*dMdGamma222333 + dMdPsi22213*dPhidU137 + dMdPsi22223*dPhidU237 + dMdPsi22233*dPhidU337
    dMdU[I3+9,I7] = dGammadU1317*dMdGamma332131 + dGammadU1327*dMdGamma332132 + dGammadU1337*dMdGamma332133 + dGammadU2317*dMdGamma332231 + dGammadU2327*dMdGamma332232 + dGammadU2337*dMdGamma332233 + dGammadU3317*dMdGamma332331 + dGammadU3327*dMdGamma332332 + dGammadU3337*dMdGamma332333 + dMdPsi33213*dPhidU137 + dMdPsi33223*dPhidU237 + dMdPsi33233*dPhidU337
    dMdU[I4+9,I7] = dGammadU1317*dMdGamma232131 + dGammadU1327*dMdGamma232132 + dGammadU1337*dMdGamma232133 + dGammadU2317*dMdGamma232231 + dGammadU2327*dMdGamma232232 + dGammadU2337*dMdGamma232233 + dGammadU3317*dMdGamma232331 + dGammadU3327*dMdGamma232332 + dGammadU3337*dMdGamma232333 + dMdPsi23213*dPhidU137 + dMdPsi23223*dPhidU237 + dMdPsi23233*dPhidU337
    dMdU[I5+9,I7] = dGammadU1317*dMdGamma132131 + dGammadU1327*dMdGamma132132 + dGammadU1337*dMdGamma132133 + dGammadU2317*dMdGamma132231 + dGammadU2327*dMdGamma132232 + dGammadU2337*dMdGamma132233 + dGammadU3317*dMdGamma132331 + dGammadU3327*dMdGamma132332 + dGammadU3337*dMdGamma132333 + dMdPsi13213*dPhidU137 + dMdPsi13223*dPhidU237 + dMdPsi13233*dPhidU337
    dMdU[I6+9,I7] = dGammadU1317*dMdGamma122131 + dGammadU1327*dMdGamma122132 + dGammadU1337*dMdGamma122133 + dGammadU2317*dMdGamma122231 + dGammadU2327*dMdGamma122232 + dGammadU2337*dMdGamma122233 + dGammadU3317*dMdGamma122331 + dGammadU3327*dMdGamma122332 + dGammadU3337*dMdGamma122333 + dMdPsi12213*dPhidU137 + dMdPsi12223*dPhidU237 + dMdPsi12233*dPhidU337
    dMdU[I7+9,I7] = dGammadU1317*dMdGamma322131 + dGammadU1327*dMdGamma322132 + dGammadU1337*dMdGamma322133 + dGammadU2317*dMdGamma322231 + dGammadU2327*dMdGamma322232 + dGammadU2337*dMdGamma322233 + dGammadU3317*dMdGamma322331 + dGammadU3327*dMdGamma322332 + dGammadU3337*dMdGamma322333 + dMdPsi32213*dPhidU137 + dMdPsi32223*dPhidU237 + dMdPsi32233*dPhidU337
    dMdU[I8+9,I7] = dGammadU1317*dMdGamma312131 + dGammadU1327*dMdGamma312132 + dGammadU1337*dMdGamma312133 + dGammadU2317*dMdGamma312231 + dGammadU2327*dMdGamma312232 + dGammadU2337*dMdGamma312233 + dGammadU3317*dMdGamma312331 + dGammadU3327*dMdGamma312332 + dGammadU3337*dMdGamma312333 + dMdPsi31213*dPhidU137 + dMdPsi31223*dPhidU237 + dMdPsi31233*dPhidU337
    dMdU[I9+9,I7] = dGammadU1317*dMdGamma212131 + dGammadU1327*dMdGamma212132 + dGammadU1337*dMdGamma212133 + dGammadU2317*dMdGamma212231 + dGammadU2327*dMdGamma212232 + dGammadU2337*dMdGamma212233 + dGammadU3317*dMdGamma212331 + dGammadU3327*dMdGamma212332 + dGammadU3337*dMdGamma212333 + dMdPsi21213*dPhidU137 + dMdPsi21223*dPhidU237 + dMdPsi21233*dPhidU337
    dMdU[I1+18,I7] = dGammadU1317*dMdGamma113131 + dGammadU1327*dMdGamma113132 + dGammadU1337*dMdGamma113133 + dGammadU2317*dMdGamma113231 + dGammadU2327*dMdGamma113232 + dGammadU2337*dMdGamma113233 + dGammadU3317*dMdGamma113331 + dGammadU3327*dMdGamma113332 + dGammadU3337*dMdGamma113333 + dMdPsi11313*dPhidU137 + dMdPsi11323*dPhidU237 + dMdPsi11333*dPhidU337
    dMdU[I2+18,I7] = dGammadU1317*dMdGamma223131 + dGammadU1327*dMdGamma223132 + dGammadU1337*dMdGamma223133 + dGammadU2317*dMdGamma223231 + dGammadU2327*dMdGamma223232 + dGammadU2337*dMdGamma223233 + dGammadU3317*dMdGamma223331 + dGammadU3327*dMdGamma223332 + dGammadU3337*dMdGamma223333 + dMdPsi22313*dPhidU137 + dMdPsi22323*dPhidU237 + dMdPsi22333*dPhidU337
    dMdU[I3+18,I7] = dGammadU1317*dMdGamma333131 + dGammadU1327*dMdGamma333132 + dGammadU1337*dMdGamma333133 + dGammadU2317*dMdGamma333231 + dGammadU2327*dMdGamma333232 + dGammadU2337*dMdGamma333233 + dGammadU3317*dMdGamma333331 + dGammadU3327*dMdGamma333332 + dGammadU3337*dMdGamma333333 + dMdPsi33313*dPhidU137 + dMdPsi33323*dPhidU237 + dMdPsi33333*dPhidU337
    dMdU[I4+18,I7] = dGammadU1317*dMdGamma233131 + dGammadU1327*dMdGamma233132 + dGammadU1337*dMdGamma233133 + dGammadU2317*dMdGamma233231 + dGammadU2327*dMdGamma233232 + dGammadU2337*dMdGamma233233 + dGammadU3317*dMdGamma233331 + dGammadU3327*dMdGamma233332 + dGammadU3337*dMdGamma233333 + dMdPsi23313*dPhidU137 + dMdPsi23323*dPhidU237 + dMdPsi23333*dPhidU337
    dMdU[I5+18,I7] = dGammadU1317*dMdGamma133131 + dGammadU1327*dMdGamma133132 + dGammadU1337*dMdGamma133133 + dGammadU2317*dMdGamma133231 + dGammadU2327*dMdGamma133232 + dGammadU2337*dMdGamma133233 + dGammadU3317*dMdGamma133331 + dGammadU3327*dMdGamma133332 + dGammadU3337*dMdGamma133333 + dMdPsi13313*dPhidU137 + dMdPsi13323*dPhidU237 + dMdPsi13333*dPhidU337
    dMdU[I6+18,I7] = dGammadU1317*dMdGamma123131 + dGammadU1327*dMdGamma123132 + dGammadU1337*dMdGamma123133 + dGammadU2317*dMdGamma123231 + dGammadU2327*dMdGamma123232 + dGammadU2337*dMdGamma123233 + dGammadU3317*dMdGamma123331 + dGammadU3327*dMdGamma123332 + dGammadU3337*dMdGamma123333 + dMdPsi12313*dPhidU137 + dMdPsi12323*dPhidU237 + dMdPsi12333*dPhidU337
    dMdU[I7+18,I7] = dGammadU1317*dMdGamma323131 + dGammadU1327*dMdGamma323132 + dGammadU1337*dMdGamma323133 + dGammadU2317*dMdGamma323231 + dGammadU2327*dMdGamma323232 + dGammadU2337*dMdGamma323233 + dGammadU3317*dMdGamma323331 + dGammadU3327*dMdGamma323332 + dGammadU3337*dMdGamma323333 + dMdPsi32313*dPhidU137 + dMdPsi32323*dPhidU237 + dMdPsi32333*dPhidU337
    dMdU[I8+18,I7] = dGammadU1317*dMdGamma313131 + dGammadU1327*dMdGamma313132 + dGammadU1337*dMdGamma313133 + dGammadU2317*dMdGamma313231 + dGammadU2327*dMdGamma313232 + dGammadU2337*dMdGamma313233 + dGammadU3317*dMdGamma313331 + dGammadU3327*dMdGamma313332 + dGammadU3337*dMdGamma313333 + dMdPsi31313*dPhidU137 + dMdPsi31323*dPhidU237 + dMdPsi31333*dPhidU337
    dMdU[I9+18,I7] = dGammadU1317*dMdGamma213131 + dGammadU1327*dMdGamma213132 + dGammadU1337*dMdGamma213133 + dGammadU2317*dMdGamma213231 + dGammadU2327*dMdGamma213232 + dGammadU2337*dMdGamma213233 + dGammadU3317*dMdGamma213331 + dGammadU3327*dMdGamma213332 + dGammadU3337*dMdGamma213333 + dMdPsi21313*dPhidU137 + dMdPsi21323*dPhidU237 + dMdPsi21333*dPhidU337

    #Column 8
    dMdU[I1,I8] = dGammadU1318*dMdGamma111131 + dGammadU1328*dMdGamma111132 + dGammadU1338*dMdGamma111133 + dGammadU2318*dMdGamma111231 + dGammadU2328*dMdGamma111232 + dGammadU2338*dMdGamma111233 + dGammadU3318*dMdGamma111331 + dGammadU3328*dMdGamma111332 + dGammadU3338*dMdGamma111333 + dMdPsi11113*dPhidU138 + dMdPsi11123*dPhidU238 + dMdPsi11133*dPhidU338
    dMdU[I2,I8] = dGammadU1318*dMdGamma221131 + dGammadU1328*dMdGamma221132 + dGammadU1338*dMdGamma221133 + dGammadU2318*dMdGamma221231 + dGammadU2328*dMdGamma221232 + dGammadU2338*dMdGamma221233 + dGammadU3318*dMdGamma221331 + dGammadU3328*dMdGamma221332 + dGammadU3338*dMdGamma221333 + dMdPsi22113*dPhidU138 + dMdPsi22123*dPhidU238 + dMdPsi22133*dPhidU338
    dMdU[I3,I8] = dGammadU1318*dMdGamma331131 + dGammadU1328*dMdGamma331132 + dGammadU1338*dMdGamma331133 + dGammadU2318*dMdGamma331231 + dGammadU2328*dMdGamma331232 + dGammadU2338*dMdGamma331233 + dGammadU3318*dMdGamma331331 + dGammadU3328*dMdGamma331332 + dGammadU3338*dMdGamma331333 + dMdPsi33113*dPhidU138 + dMdPsi33123*dPhidU238 + dMdPsi33133*dPhidU338
    dMdU[I4,I8] = dGammadU1318*dMdGamma231131 + dGammadU1328*dMdGamma231132 + dGammadU1338*dMdGamma231133 + dGammadU2318*dMdGamma231231 + dGammadU2328*dMdGamma231232 + dGammadU2338*dMdGamma231233 + dGammadU3318*dMdGamma231331 + dGammadU3328*dMdGamma231332 + dGammadU3338*dMdGamma231333 + dMdPsi23113*dPhidU138 + dMdPsi23123*dPhidU238 + dMdPsi23133*dPhidU338
    dMdU[I5,I8] = dGammadU1318*dMdGamma131131 + dGammadU1328*dMdGamma131132 + dGammadU1338*dMdGamma131133 + dGammadU2318*dMdGamma131231 + dGammadU2328*dMdGamma131232 + dGammadU2338*dMdGamma131233 + dGammadU3318*dMdGamma131331 + dGammadU3328*dMdGamma131332 + dGammadU3338*dMdGamma131333 + dMdPsi13113*dPhidU138 + dMdPsi13123*dPhidU238 + dMdPsi13133*dPhidU338
    dMdU[I6,I8] = dGammadU1318*dMdGamma121131 + dGammadU1328*dMdGamma121132 + dGammadU1338*dMdGamma121133 + dGammadU2318*dMdGamma121231 + dGammadU2328*dMdGamma121232 + dGammadU2338*dMdGamma121233 + dGammadU3318*dMdGamma121331 + dGammadU3328*dMdGamma121332 + dGammadU3338*dMdGamma121333 + dMdPsi12113*dPhidU138 + dMdPsi12123*dPhidU238 + dMdPsi12133*dPhidU338
    dMdU[I7,I8] = dGammadU1318*dMdGamma321131 + dGammadU1328*dMdGamma321132 + dGammadU1338*dMdGamma321133 + dGammadU2318*dMdGamma321231 + dGammadU2328*dMdGamma321232 + dGammadU2338*dMdGamma321233 + dGammadU3318*dMdGamma321331 + dGammadU3328*dMdGamma321332 + dGammadU3338*dMdGamma321333 + dMdPsi32113*dPhidU138 + dMdPsi32123*dPhidU238 + dMdPsi32133*dPhidU338
    dMdU[I8,I8] = dGammadU1318*dMdGamma311131 + dGammadU1328*dMdGamma311132 + dGammadU1338*dMdGamma311133 + dGammadU2318*dMdGamma311231 + dGammadU2328*dMdGamma311232 + dGammadU2338*dMdGamma311233 + dGammadU3318*dMdGamma311331 + dGammadU3328*dMdGamma311332 + dGammadU3338*dMdGamma311333 + dMdPsi31113*dPhidU138 + dMdPsi31123*dPhidU238 + dMdPsi31133*dPhidU338
    dMdU[I9,I8] = dGammadU1318*dMdGamma211131 + dGammadU1328*dMdGamma211132 + dGammadU1338*dMdGamma211133 + dGammadU2318*dMdGamma211231 + dGammadU2328*dMdGamma211232 + dGammadU2338*dMdGamma211233 + dGammadU3318*dMdGamma211331 + dGammadU3328*dMdGamma211332 + dGammadU3338*dMdGamma211333 + dMdPsi21113*dPhidU138 + dMdPsi21123*dPhidU238 + dMdPsi21133*dPhidU338
    dMdU[I1+9,I8] = dGammadU1318*dMdGamma112131 + dGammadU1328*dMdGamma112132 + dGammadU1338*dMdGamma112133 + dGammadU2318*dMdGamma112231 + dGammadU2328*dMdGamma112232 + dGammadU2338*dMdGamma112233 + dGammadU3318*dMdGamma112331 + dGammadU3328*dMdGamma112332 + dGammadU3338*dMdGamma112333 + dMdPsi11213*dPhidU138 + dMdPsi11223*dPhidU238 + dMdPsi11233*dPhidU338
    dMdU[I2+9,I8] = dGammadU1318*dMdGamma222131 + dGammadU1328*dMdGamma222132 + dGammadU1338*dMdGamma222133 + dGammadU2318*dMdGamma222231 + dGammadU2328*dMdGamma222232 + dGammadU2338*dMdGamma222233 + dGammadU3318*dMdGamma222331 + dGammadU3328*dMdGamma222332 + dGammadU3338*dMdGamma222333 + dMdPsi22213*dPhidU138 + dMdPsi22223*dPhidU238 + dMdPsi22233*dPhidU338
    dMdU[I3+9,I8] = dGammadU1318*dMdGamma332131 + dGammadU1328*dMdGamma332132 + dGammadU1338*dMdGamma332133 + dGammadU2318*dMdGamma332231 + dGammadU2328*dMdGamma332232 + dGammadU2338*dMdGamma332233 + dGammadU3318*dMdGamma332331 + dGammadU3328*dMdGamma332332 + dGammadU3338*dMdGamma332333 + dMdPsi33213*dPhidU138 + dMdPsi33223*dPhidU238 + dMdPsi33233*dPhidU338
    dMdU[I4+9,I8] = dGammadU1318*dMdGamma232131 + dGammadU1328*dMdGamma232132 + dGammadU1338*dMdGamma232133 + dGammadU2318*dMdGamma232231 + dGammadU2328*dMdGamma232232 + dGammadU2338*dMdGamma232233 + dGammadU3318*dMdGamma232331 + dGammadU3328*dMdGamma232332 + dGammadU3338*dMdGamma232333 + dMdPsi23213*dPhidU138 + dMdPsi23223*dPhidU238 + dMdPsi23233*dPhidU338
    dMdU[I5+9,I8] = dGammadU1318*dMdGamma132131 + dGammadU1328*dMdGamma132132 + dGammadU1338*dMdGamma132133 + dGammadU2318*dMdGamma132231 + dGammadU2328*dMdGamma132232 + dGammadU2338*dMdGamma132233 + dGammadU3318*dMdGamma132331 + dGammadU3328*dMdGamma132332 + dGammadU3338*dMdGamma132333 + dMdPsi13213*dPhidU138 + dMdPsi13223*dPhidU238 + dMdPsi13233*dPhidU338
    dMdU[I6+9,I8] = dGammadU1318*dMdGamma122131 + dGammadU1328*dMdGamma122132 + dGammadU1338*dMdGamma122133 + dGammadU2318*dMdGamma122231 + dGammadU2328*dMdGamma122232 + dGammadU2338*dMdGamma122233 + dGammadU3318*dMdGamma122331 + dGammadU3328*dMdGamma122332 + dGammadU3338*dMdGamma122333 + dMdPsi12213*dPhidU138 + dMdPsi12223*dPhidU238 + dMdPsi12233*dPhidU338
    dMdU[I7+9,I8] = dGammadU1318*dMdGamma322131 + dGammadU1328*dMdGamma322132 + dGammadU1338*dMdGamma322133 + dGammadU2318*dMdGamma322231 + dGammadU2328*dMdGamma322232 + dGammadU2338*dMdGamma322233 + dGammadU3318*dMdGamma322331 + dGammadU3328*dMdGamma322332 + dGammadU3338*dMdGamma322333 + dMdPsi32213*dPhidU138 + dMdPsi32223*dPhidU238 + dMdPsi32233*dPhidU338
    dMdU[I8+9,I8] = dGammadU1318*dMdGamma312131 + dGammadU1328*dMdGamma312132 + dGammadU1338*dMdGamma312133 + dGammadU2318*dMdGamma312231 + dGammadU2328*dMdGamma312232 + dGammadU2338*dMdGamma312233 + dGammadU3318*dMdGamma312331 + dGammadU3328*dMdGamma312332 + dGammadU3338*dMdGamma312333 + dMdPsi31213*dPhidU138 + dMdPsi31223*dPhidU238 + dMdPsi31233*dPhidU338
    dMdU[I9+9,I8] = dGammadU1318*dMdGamma212131 + dGammadU1328*dMdGamma212132 + dGammadU1338*dMdGamma212133 + dGammadU2318*dMdGamma212231 + dGammadU2328*dMdGamma212232 + dGammadU2338*dMdGamma212233 + dGammadU3318*dMdGamma212331 + dGammadU3328*dMdGamma212332 + dGammadU3338*dMdGamma212333 + dMdPsi21213*dPhidU138 + dMdPsi21223*dPhidU238 + dMdPsi21233*dPhidU338
    dMdU[I1+18,I8] = dGammadU1318*dMdGamma113131 + dGammadU1328*dMdGamma113132 + dGammadU1338*dMdGamma113133 + dGammadU2318*dMdGamma113231 + dGammadU2328*dMdGamma113232 + dGammadU2338*dMdGamma113233 + dGammadU3318*dMdGamma113331 + dGammadU3328*dMdGamma113332 + dGammadU3338*dMdGamma113333 + dMdPsi11313*dPhidU138 + dMdPsi11323*dPhidU238 + dMdPsi11333*dPhidU338
    dMdU[I2+18,I8] = dGammadU1318*dMdGamma223131 + dGammadU1328*dMdGamma223132 + dGammadU1338*dMdGamma223133 + dGammadU2318*dMdGamma223231 + dGammadU2328*dMdGamma223232 + dGammadU2338*dMdGamma223233 + dGammadU3318*dMdGamma223331 + dGammadU3328*dMdGamma223332 + dGammadU3338*dMdGamma223333 + dMdPsi22313*dPhidU138 + dMdPsi22323*dPhidU238 + dMdPsi22333*dPhidU338
    dMdU[I3+18,I8] = dGammadU1318*dMdGamma333131 + dGammadU1328*dMdGamma333132 + dGammadU1338*dMdGamma333133 + dGammadU2318*dMdGamma333231 + dGammadU2328*dMdGamma333232 + dGammadU2338*dMdGamma333233 + dGammadU3318*dMdGamma333331 + dGammadU3328*dMdGamma333332 + dGammadU3338*dMdGamma333333 + dMdPsi33313*dPhidU138 + dMdPsi33323*dPhidU238 + dMdPsi33333*dPhidU338
    dMdU[I4+18,I8] = dGammadU1318*dMdGamma233131 + dGammadU1328*dMdGamma233132 + dGammadU1338*dMdGamma233133 + dGammadU2318*dMdGamma233231 + dGammadU2328*dMdGamma233232 + dGammadU2338*dMdGamma233233 + dGammadU3318*dMdGamma233331 + dGammadU3328*dMdGamma233332 + dGammadU3338*dMdGamma233333 + dMdPsi23313*dPhidU138 + dMdPsi23323*dPhidU238 + dMdPsi23333*dPhidU338
    dMdU[I5+18,I8] = dGammadU1318*dMdGamma133131 + dGammadU1328*dMdGamma133132 + dGammadU1338*dMdGamma133133 + dGammadU2318*dMdGamma133231 + dGammadU2328*dMdGamma133232 + dGammadU2338*dMdGamma133233 + dGammadU3318*dMdGamma133331 + dGammadU3328*dMdGamma133332 + dGammadU3338*dMdGamma133333 + dMdPsi13313*dPhidU138 + dMdPsi13323*dPhidU238 + dMdPsi13333*dPhidU338
    dMdU[I6+18,I8] = dGammadU1318*dMdGamma123131 + dGammadU1328*dMdGamma123132 + dGammadU1338*dMdGamma123133 + dGammadU2318*dMdGamma123231 + dGammadU2328*dMdGamma123232 + dGammadU2338*dMdGamma123233 + dGammadU3318*dMdGamma123331 + dGammadU3328*dMdGamma123332 + dGammadU3338*dMdGamma123333 + dMdPsi12313*dPhidU138 + dMdPsi12323*dPhidU238 + dMdPsi12333*dPhidU338
    dMdU[I7+18,I8] = dGammadU1318*dMdGamma323131 + dGammadU1328*dMdGamma323132 + dGammadU1338*dMdGamma323133 + dGammadU2318*dMdGamma323231 + dGammadU2328*dMdGamma323232 + dGammadU2338*dMdGamma323233 + dGammadU3318*dMdGamma323331 + dGammadU3328*dMdGamma323332 + dGammadU3338*dMdGamma323333 + dMdPsi32313*dPhidU138 + dMdPsi32323*dPhidU238 + dMdPsi32333*dPhidU338
    dMdU[I8+18,I8] = dGammadU1318*dMdGamma313131 + dGammadU1328*dMdGamma313132 + dGammadU1338*dMdGamma313133 + dGammadU2318*dMdGamma313231 + dGammadU2328*dMdGamma313232 + dGammadU2338*dMdGamma313233 + dGammadU3318*dMdGamma313331 + dGammadU3328*dMdGamma313332 + dGammadU3338*dMdGamma313333 + dMdPsi31313*dPhidU138 + dMdPsi31323*dPhidU238 + dMdPsi31333*dPhidU338
    dMdU[I9+18,I8] = dGammadU1318*dMdGamma213131 + dGammadU1328*dMdGamma213132 + dGammadU1338*dMdGamma213133 + dGammadU2318*dMdGamma213231 + dGammadU2328*dMdGamma213232 + dGammadU2338*dMdGamma213233 + dGammadU3318*dMdGamma213331 + dGammadU3328*dMdGamma213332 + dGammadU3338*dMdGamma213333 + dMdPsi21313*dPhidU138 + dMdPsi21323*dPhidU238 + dMdPsi21333*dPhidU338

    #Column 9
    dMdU[I1,I9] = dGammadU1219*dMdGamma111121 + dGammadU1229*dMdGamma111122 + dGammadU1239*dMdGamma111123 + dGammadU2219*dMdGamma111221 + dGammadU2229*dMdGamma111222 + dGammadU2239*dMdGamma111223 + dGammadU3219*dMdGamma111321 + dGammadU3229*dMdGamma111322 + dGammadU3239*dMdGamma111323 + dMdPsi11112*dPhidU129 + dMdPsi11122*dPhidU229 + dMdPsi11132*dPhidU329
    dMdU[I2,I9] = dGammadU1219*dMdGamma221121 + dGammadU1229*dMdGamma221122 + dGammadU1239*dMdGamma221123 + dGammadU2219*dMdGamma221221 + dGammadU2229*dMdGamma221222 + dGammadU2239*dMdGamma221223 + dGammadU3219*dMdGamma221321 + dGammadU3229*dMdGamma221322 + dGammadU3239*dMdGamma221323 + dMdPsi22112*dPhidU129 + dMdPsi22122*dPhidU229 + dMdPsi22132*dPhidU329
    dMdU[I3,I9] = dGammadU1219*dMdGamma331121 + dGammadU1229*dMdGamma331122 + dGammadU1239*dMdGamma331123 + dGammadU2219*dMdGamma331221 + dGammadU2229*dMdGamma331222 + dGammadU2239*dMdGamma331223 + dGammadU3219*dMdGamma331321 + dGammadU3229*dMdGamma331322 + dGammadU3239*dMdGamma331323 + dMdPsi33112*dPhidU129 + dMdPsi33122*dPhidU229 + dMdPsi33132*dPhidU329
    dMdU[I4,I9] = dGammadU1219*dMdGamma231121 + dGammadU1229*dMdGamma231122 + dGammadU1239*dMdGamma231123 + dGammadU2219*dMdGamma231221 + dGammadU2229*dMdGamma231222 + dGammadU2239*dMdGamma231223 + dGammadU3219*dMdGamma231321 + dGammadU3229*dMdGamma231322 + dGammadU3239*dMdGamma231323 + dMdPsi23112*dPhidU129 + dMdPsi23122*dPhidU229 + dMdPsi23132*dPhidU329
    dMdU[I5,I9] = dGammadU1219*dMdGamma131121 + dGammadU1229*dMdGamma131122 + dGammadU1239*dMdGamma131123 + dGammadU2219*dMdGamma131221 + dGammadU2229*dMdGamma131222 + dGammadU2239*dMdGamma131223 + dGammadU3219*dMdGamma131321 + dGammadU3229*dMdGamma131322 + dGammadU3239*dMdGamma131323 + dMdPsi13112*dPhidU129 + dMdPsi13122*dPhidU229 + dMdPsi13132*dPhidU329
    dMdU[I6,I9] = dGammadU1219*dMdGamma121121 + dGammadU1229*dMdGamma121122 + dGammadU1239*dMdGamma121123 + dGammadU2219*dMdGamma121221 + dGammadU2229*dMdGamma121222 + dGammadU2239*dMdGamma121223 + dGammadU3219*dMdGamma121321 + dGammadU3229*dMdGamma121322 + dGammadU3239*dMdGamma121323 + dMdPsi12112*dPhidU129 + dMdPsi12122*dPhidU229 + dMdPsi12132*dPhidU329
    dMdU[I7,I9] = dGammadU1219*dMdGamma321121 + dGammadU1229*dMdGamma321122 + dGammadU1239*dMdGamma321123 + dGammadU2219*dMdGamma321221 + dGammadU2229*dMdGamma321222 + dGammadU2239*dMdGamma321223 + dGammadU3219*dMdGamma321321 + dGammadU3229*dMdGamma321322 + dGammadU3239*dMdGamma321323 + dMdPsi32112*dPhidU129 + dMdPsi32122*dPhidU229 + dMdPsi32132*dPhidU329
    dMdU[I8,I9] = dGammadU1219*dMdGamma311121 + dGammadU1229*dMdGamma311122 + dGammadU1239*dMdGamma311123 + dGammadU2219*dMdGamma311221 + dGammadU2229*dMdGamma311222 + dGammadU2239*dMdGamma311223 + dGammadU3219*dMdGamma311321 + dGammadU3229*dMdGamma311322 + dGammadU3239*dMdGamma311323 + dMdPsi31112*dPhidU129 + dMdPsi31122*dPhidU229 + dMdPsi31132*dPhidU329
    dMdU[I9,I9] = dGammadU1219*dMdGamma211121 + dGammadU1229*dMdGamma211122 + dGammadU1239*dMdGamma211123 + dGammadU2219*dMdGamma211221 + dGammadU2229*dMdGamma211222 + dGammadU2239*dMdGamma211223 + dGammadU3219*dMdGamma211321 + dGammadU3229*dMdGamma211322 + dGammadU3239*dMdGamma211323 + dMdPsi21112*dPhidU129 + dMdPsi21122*dPhidU229 + dMdPsi21132*dPhidU329
    dMdU[I1+9,I9] = dGammadU1219*dMdGamma112121 + dGammadU1229*dMdGamma112122 + dGammadU1239*dMdGamma112123 + dGammadU2219*dMdGamma112221 + dGammadU2229*dMdGamma112222 + dGammadU2239*dMdGamma112223 + dGammadU3219*dMdGamma112321 + dGammadU3229*dMdGamma112322 + dGammadU3239*dMdGamma112323 + dMdPsi11212*dPhidU129 + dMdPsi11222*dPhidU229 + dMdPsi11232*dPhidU329
    dMdU[I2+9,I9] = dGammadU1219*dMdGamma222121 + dGammadU1229*dMdGamma222122 + dGammadU1239*dMdGamma222123 + dGammadU2219*dMdGamma222221 + dGammadU2229*dMdGamma222222 + dGammadU2239*dMdGamma222223 + dGammadU3219*dMdGamma222321 + dGammadU3229*dMdGamma222322 + dGammadU3239*dMdGamma222323 + dMdPsi22212*dPhidU129 + dMdPsi22222*dPhidU229 + dMdPsi22232*dPhidU329
    dMdU[I3+9,I9] = dGammadU1219*dMdGamma332121 + dGammadU1229*dMdGamma332122 + dGammadU1239*dMdGamma332123 + dGammadU2219*dMdGamma332221 + dGammadU2229*dMdGamma332222 + dGammadU2239*dMdGamma332223 + dGammadU3219*dMdGamma332321 + dGammadU3229*dMdGamma332322 + dGammadU3239*dMdGamma332323 + dMdPsi33212*dPhidU129 + dMdPsi33222*dPhidU229 + dMdPsi33232*dPhidU329
    dMdU[I4+9,I9] = dGammadU1219*dMdGamma232121 + dGammadU1229*dMdGamma232122 + dGammadU1239*dMdGamma232123 + dGammadU2219*dMdGamma232221 + dGammadU2229*dMdGamma232222 + dGammadU2239*dMdGamma232223 + dGammadU3219*dMdGamma232321 + dGammadU3229*dMdGamma232322 + dGammadU3239*dMdGamma232323 + dMdPsi23212*dPhidU129 + dMdPsi23222*dPhidU229 + dMdPsi23232*dPhidU329
    dMdU[I5+9,I9] = dGammadU1219*dMdGamma132121 + dGammadU1229*dMdGamma132122 + dGammadU1239*dMdGamma132123 + dGammadU2219*dMdGamma132221 + dGammadU2229*dMdGamma132222 + dGammadU2239*dMdGamma132223 + dGammadU3219*dMdGamma132321 + dGammadU3229*dMdGamma132322 + dGammadU3239*dMdGamma132323 + dMdPsi13212*dPhidU129 + dMdPsi13222*dPhidU229 + dMdPsi13232*dPhidU329
    dMdU[I6+9,I9] = dGammadU1219*dMdGamma122121 + dGammadU1229*dMdGamma122122 + dGammadU1239*dMdGamma122123 + dGammadU2219*dMdGamma122221 + dGammadU2229*dMdGamma122222 + dGammadU2239*dMdGamma122223 + dGammadU3219*dMdGamma122321 + dGammadU3229*dMdGamma122322 + dGammadU3239*dMdGamma122323 + dMdPsi12212*dPhidU129 + dMdPsi12222*dPhidU229 + dMdPsi12232*dPhidU329
    dMdU[I7+9,I9] = dGammadU1219*dMdGamma322121 + dGammadU1229*dMdGamma322122 + dGammadU1239*dMdGamma322123 + dGammadU2219*dMdGamma322221 + dGammadU2229*dMdGamma322222 + dGammadU2239*dMdGamma322223 + dGammadU3219*dMdGamma322321 + dGammadU3229*dMdGamma322322 + dGammadU3239*dMdGamma322323 + dMdPsi32212*dPhidU129 + dMdPsi32222*dPhidU229 + dMdPsi32232*dPhidU329
    dMdU[I8+9,I9] = dGammadU1219*dMdGamma312121 + dGammadU1229*dMdGamma312122 + dGammadU1239*dMdGamma312123 + dGammadU2219*dMdGamma312221 + dGammadU2229*dMdGamma312222 + dGammadU2239*dMdGamma312223 + dGammadU3219*dMdGamma312321 + dGammadU3229*dMdGamma312322 + dGammadU3239*dMdGamma312323 + dMdPsi31212*dPhidU129 + dMdPsi31222*dPhidU229 + dMdPsi31232*dPhidU329
    dMdU[I9+9,I9] = dGammadU1219*dMdGamma212121 + dGammadU1229*dMdGamma212122 + dGammadU1239*dMdGamma212123 + dGammadU2219*dMdGamma212221 + dGammadU2229*dMdGamma212222 + dGammadU2239*dMdGamma212223 + dGammadU3219*dMdGamma212321 + dGammadU3229*dMdGamma212322 + dGammadU3239*dMdGamma212323 + dMdPsi21212*dPhidU129 + dMdPsi21222*dPhidU229 + dMdPsi21232*dPhidU329
    dMdU[I1+18,I9] = dGammadU1219*dMdGamma113121 + dGammadU1229*dMdGamma113122 + dGammadU1239*dMdGamma113123 + dGammadU2219*dMdGamma113221 + dGammadU2229*dMdGamma113222 + dGammadU2239*dMdGamma113223 + dGammadU3219*dMdGamma113321 + dGammadU3229*dMdGamma113322 + dGammadU3239*dMdGamma113323 + dMdPsi11312*dPhidU129 + dMdPsi11322*dPhidU229 + dMdPsi11332*dPhidU329
    dMdU[I2+18,I9] = dGammadU1219*dMdGamma223121 + dGammadU1229*dMdGamma223122 + dGammadU1239*dMdGamma223123 + dGammadU2219*dMdGamma223221 + dGammadU2229*dMdGamma223222 + dGammadU2239*dMdGamma223223 + dGammadU3219*dMdGamma223321 + dGammadU3229*dMdGamma223322 + dGammadU3239*dMdGamma223323 + dMdPsi22312*dPhidU129 + dMdPsi22322*dPhidU229 + dMdPsi22332*dPhidU329
    dMdU[I3+18,I9] = dGammadU1219*dMdGamma333121 + dGammadU1229*dMdGamma333122 + dGammadU1239*dMdGamma333123 + dGammadU2219*dMdGamma333221 + dGammadU2229*dMdGamma333222 + dGammadU2239*dMdGamma333223 + dGammadU3219*dMdGamma333321 + dGammadU3229*dMdGamma333322 + dGammadU3239*dMdGamma333323 + dMdPsi33312*dPhidU129 + dMdPsi33322*dPhidU229 + dMdPsi33332*dPhidU329
    dMdU[I4+18,I9] = dGammadU1219*dMdGamma233121 + dGammadU1229*dMdGamma233122 + dGammadU1239*dMdGamma233123 + dGammadU2219*dMdGamma233221 + dGammadU2229*dMdGamma233222 + dGammadU2239*dMdGamma233223 + dGammadU3219*dMdGamma233321 + dGammadU3229*dMdGamma233322 + dGammadU3239*dMdGamma233323 + dMdPsi23312*dPhidU129 + dMdPsi23322*dPhidU229 + dMdPsi23332*dPhidU329
    dMdU[I5+18,I9] = dGammadU1219*dMdGamma133121 + dGammadU1229*dMdGamma133122 + dGammadU1239*dMdGamma133123 + dGammadU2219*dMdGamma133221 + dGammadU2229*dMdGamma133222 + dGammadU2239*dMdGamma133223 + dGammadU3219*dMdGamma133321 + dGammadU3229*dMdGamma133322 + dGammadU3239*dMdGamma133323 + dMdPsi13312*dPhidU129 + dMdPsi13322*dPhidU229 + dMdPsi13332*dPhidU329
    dMdU[I6+18,I9] = dGammadU1219*dMdGamma123121 + dGammadU1229*dMdGamma123122 + dGammadU1239*dMdGamma123123 + dGammadU2219*dMdGamma123221 + dGammadU2229*dMdGamma123222 + dGammadU2239*dMdGamma123223 + dGammadU3219*dMdGamma123321 + dGammadU3229*dMdGamma123322 + dGammadU3239*dMdGamma123323 + dMdPsi12312*dPhidU129 + dMdPsi12322*dPhidU229 + dMdPsi12332*dPhidU329
    dMdU[I7+18,I9] = dGammadU1219*dMdGamma323121 + dGammadU1229*dMdGamma323122 + dGammadU1239*dMdGamma323123 + dGammadU2219*dMdGamma323221 + dGammadU2229*dMdGamma323222 + dGammadU2239*dMdGamma323223 + dGammadU3219*dMdGamma323321 + dGammadU3229*dMdGamma323322 + dGammadU3239*dMdGamma323323 + dMdPsi32312*dPhidU129 + dMdPsi32322*dPhidU229 + dMdPsi32332*dPhidU329
    dMdU[I8+18,I9] = dGammadU1219*dMdGamma313121 + dGammadU1229*dMdGamma313122 + dGammadU1239*dMdGamma313123 + dGammadU2219*dMdGamma313221 + dGammadU2229*dMdGamma313222 + dGammadU2239*dMdGamma313223 + dGammadU3219*dMdGamma313321 + dGammadU3229*dMdGamma313322 + dGammadU3239*dMdGamma313323 + dMdPsi31312*dPhidU129 + dMdPsi31322*dPhidU229 + dMdPsi31332*dPhidU329
    dMdU[I9+18,I9] = dGammadU1219*dMdGamma213121 + dGammadU1229*dMdGamma213122 + dGammadU1239*dMdGamma213123 + dGammadU2219*dMdGamma213221 + dGammadU2229*dMdGamma213222 + dGammadU2239*dMdGamma213223 + dGammadU3219*dMdGamma213321 + dGammadU3229*dMdGamma213322 + dGammadU3239*dMdGamma213323 + dMdPsi21312*dPhidU129 + dMdPsi21322*dPhidU229 + dMdPsi21332*dPhidU329

    #Column 10
    dMdU[I1,I10] = dGammadU12110*dMdGamma111121 + dGammadU12210*dMdGamma111122 + dGammadU12310*dMdGamma111123 + dGammadU22110*dMdGamma111221 + dGammadU22210*dMdGamma111222 + dGammadU22310*dMdGamma111223 + dGammadU32110*dMdGamma111321 + dGammadU32210*dMdGamma111322 + dGammadU32310*dMdGamma111323 + dMdPsi11112*dPhidU1210 + dMdPsi11122*dPhidU2210 + dMdPsi11132*dPhidU3210
    dMdU[I2,I10] = dGammadU12110*dMdGamma221121 + dGammadU12210*dMdGamma221122 + dGammadU12310*dMdGamma221123 + dGammadU22110*dMdGamma221221 + dGammadU22210*dMdGamma221222 + dGammadU22310*dMdGamma221223 + dGammadU32110*dMdGamma221321 + dGammadU32210*dMdGamma221322 + dGammadU32310*dMdGamma221323 + dMdPsi22112*dPhidU1210 + dMdPsi22122*dPhidU2210 + dMdPsi22132*dPhidU3210
    dMdU[I3,I10] = dGammadU12110*dMdGamma331121 + dGammadU12210*dMdGamma331122 + dGammadU12310*dMdGamma331123 + dGammadU22110*dMdGamma331221 + dGammadU22210*dMdGamma331222 + dGammadU22310*dMdGamma331223 + dGammadU32110*dMdGamma331321 + dGammadU32210*dMdGamma331322 + dGammadU32310*dMdGamma331323 + dMdPsi33112*dPhidU1210 + dMdPsi33122*dPhidU2210 + dMdPsi33132*dPhidU3210
    dMdU[I4,I10] = dGammadU12110*dMdGamma231121 + dGammadU12210*dMdGamma231122 + dGammadU12310*dMdGamma231123 + dGammadU22110*dMdGamma231221 + dGammadU22210*dMdGamma231222 + dGammadU22310*dMdGamma231223 + dGammadU32110*dMdGamma231321 + dGammadU32210*dMdGamma231322 + dGammadU32310*dMdGamma231323 + dMdPsi23112*dPhidU1210 + dMdPsi23122*dPhidU2210 + dMdPsi23132*dPhidU3210
    dMdU[I5,I10] = dGammadU12110*dMdGamma131121 + dGammadU12210*dMdGamma131122 + dGammadU12310*dMdGamma131123 + dGammadU22110*dMdGamma131221 + dGammadU22210*dMdGamma131222 + dGammadU22310*dMdGamma131223 + dGammadU32110*dMdGamma131321 + dGammadU32210*dMdGamma131322 + dGammadU32310*dMdGamma131323 + dMdPsi13112*dPhidU1210 + dMdPsi13122*dPhidU2210 + dMdPsi13132*dPhidU3210
    dMdU[I6,I10] = dGammadU12110*dMdGamma121121 + dGammadU12210*dMdGamma121122 + dGammadU12310*dMdGamma121123 + dGammadU22110*dMdGamma121221 + dGammadU22210*dMdGamma121222 + dGammadU22310*dMdGamma121223 + dGammadU32110*dMdGamma121321 + dGammadU32210*dMdGamma121322 + dGammadU32310*dMdGamma121323 + dMdPsi12112*dPhidU1210 + dMdPsi12122*dPhidU2210 + dMdPsi12132*dPhidU3210
    dMdU[I7,I10] = dGammadU12110*dMdGamma321121 + dGammadU12210*dMdGamma321122 + dGammadU12310*dMdGamma321123 + dGammadU22110*dMdGamma321221 + dGammadU22210*dMdGamma321222 + dGammadU22310*dMdGamma321223 + dGammadU32110*dMdGamma321321 + dGammadU32210*dMdGamma321322 + dGammadU32310*dMdGamma321323 + dMdPsi32112*dPhidU1210 + dMdPsi32122*dPhidU2210 + dMdPsi32132*dPhidU3210
    dMdU[I8,I10] = dGammadU12110*dMdGamma311121 + dGammadU12210*dMdGamma311122 + dGammadU12310*dMdGamma311123 + dGammadU22110*dMdGamma311221 + dGammadU22210*dMdGamma311222 + dGammadU22310*dMdGamma311223 + dGammadU32110*dMdGamma311321 + dGammadU32210*dMdGamma311322 + dGammadU32310*dMdGamma311323 + dMdPsi31112*dPhidU1210 + dMdPsi31122*dPhidU2210 + dMdPsi31132*dPhidU3210
    dMdU[I9,I10] = dGammadU12110*dMdGamma211121 + dGammadU12210*dMdGamma211122 + dGammadU12310*dMdGamma211123 + dGammadU22110*dMdGamma211221 + dGammadU22210*dMdGamma211222 + dGammadU22310*dMdGamma211223 + dGammadU32110*dMdGamma211321 + dGammadU32210*dMdGamma211322 + dGammadU32310*dMdGamma211323 + dMdPsi21112*dPhidU1210 + dMdPsi21122*dPhidU2210 + dMdPsi21132*dPhidU3210
    dMdU[I1+9,I10] = dGammadU12110*dMdGamma112121 + dGammadU12210*dMdGamma112122 + dGammadU12310*dMdGamma112123 + dGammadU22110*dMdGamma112221 + dGammadU22210*dMdGamma112222 + dGammadU22310*dMdGamma112223 + dGammadU32110*dMdGamma112321 + dGammadU32210*dMdGamma112322 + dGammadU32310*dMdGamma112323 + dMdPsi11212*dPhidU1210 + dMdPsi11222*dPhidU2210 + dMdPsi11232*dPhidU3210
    dMdU[I2+9,I10] = dGammadU12110*dMdGamma222121 + dGammadU12210*dMdGamma222122 + dGammadU12310*dMdGamma222123 + dGammadU22110*dMdGamma222221 + dGammadU22210*dMdGamma222222 + dGammadU22310*dMdGamma222223 + dGammadU32110*dMdGamma222321 + dGammadU32210*dMdGamma222322 + dGammadU32310*dMdGamma222323 + dMdPsi22212*dPhidU1210 + dMdPsi22222*dPhidU2210 + dMdPsi22232*dPhidU3210
    dMdU[I3+9,I10] = dGammadU12110*dMdGamma332121 + dGammadU12210*dMdGamma332122 + dGammadU12310*dMdGamma332123 + dGammadU22110*dMdGamma332221 + dGammadU22210*dMdGamma332222 + dGammadU22310*dMdGamma332223 + dGammadU32110*dMdGamma332321 + dGammadU32210*dMdGamma332322 + dGammadU32310*dMdGamma332323 + dMdPsi33212*dPhidU1210 + dMdPsi33222*dPhidU2210 + dMdPsi33232*dPhidU3210
    dMdU[I4+9,I10] = dGammadU12110*dMdGamma232121 + dGammadU12210*dMdGamma232122 + dGammadU12310*dMdGamma232123 + dGammadU22110*dMdGamma232221 + dGammadU22210*dMdGamma232222 + dGammadU22310*dMdGamma232223 + dGammadU32110*dMdGamma232321 + dGammadU32210*dMdGamma232322 + dGammadU32310*dMdGamma232323 + dMdPsi23212*dPhidU1210 + dMdPsi23222*dPhidU2210 + dMdPsi23232*dPhidU3210
    dMdU[I5+9,I10] = dGammadU12110*dMdGamma132121 + dGammadU12210*dMdGamma132122 + dGammadU12310*dMdGamma132123 + dGammadU22110*dMdGamma132221 + dGammadU22210*dMdGamma132222 + dGammadU22310*dMdGamma132223 + dGammadU32110*dMdGamma132321 + dGammadU32210*dMdGamma132322 + dGammadU32310*dMdGamma132323 + dMdPsi13212*dPhidU1210 + dMdPsi13222*dPhidU2210 + dMdPsi13232*dPhidU3210
    dMdU[I6+9,I10] = dGammadU12110*dMdGamma122121 + dGammadU12210*dMdGamma122122 + dGammadU12310*dMdGamma122123 + dGammadU22110*dMdGamma122221 + dGammadU22210*dMdGamma122222 + dGammadU22310*dMdGamma122223 + dGammadU32110*dMdGamma122321 + dGammadU32210*dMdGamma122322 + dGammadU32310*dMdGamma122323 + dMdPsi12212*dPhidU1210 + dMdPsi12222*dPhidU2210 + dMdPsi12232*dPhidU3210
    dMdU[I7+9,I10] = dGammadU12110*dMdGamma322121 + dGammadU12210*dMdGamma322122 + dGammadU12310*dMdGamma322123 + dGammadU22110*dMdGamma322221 + dGammadU22210*dMdGamma322222 + dGammadU22310*dMdGamma322223 + dGammadU32110*dMdGamma322321 + dGammadU32210*dMdGamma322322 + dGammadU32310*dMdGamma322323 + dMdPsi32212*dPhidU1210 + dMdPsi32222*dPhidU2210 + dMdPsi32232*dPhidU3210
    dMdU[I8+9,I10] = dGammadU12110*dMdGamma312121 + dGammadU12210*dMdGamma312122 + dGammadU12310*dMdGamma312123 + dGammadU22110*dMdGamma312221 + dGammadU22210*dMdGamma312222 + dGammadU22310*dMdGamma312223 + dGammadU32110*dMdGamma312321 + dGammadU32210*dMdGamma312322 + dGammadU32310*dMdGamma312323 + dMdPsi31212*dPhidU1210 + dMdPsi31222*dPhidU2210 + dMdPsi31232*dPhidU3210
    dMdU[I9+9,I10] = dGammadU12110*dMdGamma212121 + dGammadU12210*dMdGamma212122 + dGammadU12310*dMdGamma212123 + dGammadU22110*dMdGamma212221 + dGammadU22210*dMdGamma212222 + dGammadU22310*dMdGamma212223 + dGammadU32110*dMdGamma212321 + dGammadU32210*dMdGamma212322 + dGammadU32310*dMdGamma212323 + dMdPsi21212*dPhidU1210 + dMdPsi21222*dPhidU2210 + dMdPsi21232*dPhidU3210
    dMdU[I1+18,I10] = dGammadU12110*dMdGamma113121 + dGammadU12210*dMdGamma113122 + dGammadU12310*dMdGamma113123 + dGammadU22110*dMdGamma113221 + dGammadU22210*dMdGamma113222 + dGammadU22310*dMdGamma113223 + dGammadU32110*dMdGamma113321 + dGammadU32210*dMdGamma113322 + dGammadU32310*dMdGamma113323 + dMdPsi11312*dPhidU1210 + dMdPsi11322*dPhidU2210 + dMdPsi11332*dPhidU3210
    dMdU[I2+18,I10] = dGammadU12110*dMdGamma223121 + dGammadU12210*dMdGamma223122 + dGammadU12310*dMdGamma223123 + dGammadU22110*dMdGamma223221 + dGammadU22210*dMdGamma223222 + dGammadU22310*dMdGamma223223 + dGammadU32110*dMdGamma223321 + dGammadU32210*dMdGamma223322 + dGammadU32310*dMdGamma223323 + dMdPsi22312*dPhidU1210 + dMdPsi22322*dPhidU2210 + dMdPsi22332*dPhidU3210
    dMdU[I3+18,I10] = dGammadU12110*dMdGamma333121 + dGammadU12210*dMdGamma333122 + dGammadU12310*dMdGamma333123 + dGammadU22110*dMdGamma333221 + dGammadU22210*dMdGamma333222 + dGammadU22310*dMdGamma333223 + dGammadU32110*dMdGamma333321 + dGammadU32210*dMdGamma333322 + dGammadU32310*dMdGamma333323 + dMdPsi33312*dPhidU1210 + dMdPsi33322*dPhidU2210 + dMdPsi33332*dPhidU3210
    dMdU[I4+18,I10] = dGammadU12110*dMdGamma233121 + dGammadU12210*dMdGamma233122 + dGammadU12310*dMdGamma233123 + dGammadU22110*dMdGamma233221 + dGammadU22210*dMdGamma233222 + dGammadU22310*dMdGamma233223 + dGammadU32110*dMdGamma233321 + dGammadU32210*dMdGamma233322 + dGammadU32310*dMdGamma233323 + dMdPsi23312*dPhidU1210 + dMdPsi23322*dPhidU2210 + dMdPsi23332*dPhidU3210
    dMdU[I5+18,I10] = dGammadU12110*dMdGamma133121 + dGammadU12210*dMdGamma133122 + dGammadU12310*dMdGamma133123 + dGammadU22110*dMdGamma133221 + dGammadU22210*dMdGamma133222 + dGammadU22310*dMdGamma133223 + dGammadU32110*dMdGamma133321 + dGammadU32210*dMdGamma133322 + dGammadU32310*dMdGamma133323 + dMdPsi13312*dPhidU1210 + dMdPsi13322*dPhidU2210 + dMdPsi13332*dPhidU3210
    dMdU[I6+18,I10] = dGammadU12110*dMdGamma123121 + dGammadU12210*dMdGamma123122 + dGammadU12310*dMdGamma123123 + dGammadU22110*dMdGamma123221 + dGammadU22210*dMdGamma123222 + dGammadU22310*dMdGamma123223 + dGammadU32110*dMdGamma123321 + dGammadU32210*dMdGamma123322 + dGammadU32310*dMdGamma123323 + dMdPsi12312*dPhidU1210 + dMdPsi12322*dPhidU2210 + dMdPsi12332*dPhidU3210
    dMdU[I7+18,I10] = dGammadU12110*dMdGamma323121 + dGammadU12210*dMdGamma323122 + dGammadU12310*dMdGamma323123 + dGammadU22110*dMdGamma323221 + dGammadU22210*dMdGamma323222 + dGammadU22310*dMdGamma323223 + dGammadU32110*dMdGamma323321 + dGammadU32210*dMdGamma323322 + dGammadU32310*dMdGamma323323 + dMdPsi32312*dPhidU1210 + dMdPsi32322*dPhidU2210 + dMdPsi32332*dPhidU3210
    dMdU[I8+18,I10] = dGammadU12110*dMdGamma313121 + dGammadU12210*dMdGamma313122 + dGammadU12310*dMdGamma313123 + dGammadU22110*dMdGamma313221 + dGammadU22210*dMdGamma313222 + dGammadU22310*dMdGamma313223 + dGammadU32110*dMdGamma313321 + dGammadU32210*dMdGamma313322 + dGammadU32310*dMdGamma313323 + dMdPsi31312*dPhidU1210 + dMdPsi31322*dPhidU2210 + dMdPsi31332*dPhidU3210
    dMdU[I9+18,I10] = dGammadU12110*dMdGamma213121 + dGammadU12210*dMdGamma213122 + dGammadU12310*dMdGamma213123 + dGammadU22110*dMdGamma213221 + dGammadU22210*dMdGamma213222 + dGammadU22310*dMdGamma213223 + dGammadU32110*dMdGamma213321 + dGammadU32210*dMdGamma213322 + dGammadU32310*dMdGamma213323 + dMdPsi21312*dPhidU1210 + dMdPsi21322*dPhidU2210 + dMdPsi21332*dPhidU3210

    #Column 11
    dMdU[I1,I11] = dGammadU11111*dMdGamma111111 + dGammadU11211*dMdGamma111112 + dGammadU11311*dMdGamma111113 + dGammadU21111*dMdGamma111211 + dGammadU21211*dMdGamma111212 + dGammadU21311*dMdGamma111213 + dGammadU31111*dMdGamma111311 + dGammadU31211*dMdGamma111312 + dGammadU31311*dMdGamma111313 + dMdPsi11111*dPhidU1111 + dMdPsi11121*dPhidU2111 + dMdPsi11131*dPhidU3111
    dMdU[I2,I11] = dGammadU11111*dMdGamma221111 + dGammadU11211*dMdGamma221112 + dGammadU11311*dMdGamma221113 + dGammadU21111*dMdGamma221211 + dGammadU21211*dMdGamma221212 + dGammadU21311*dMdGamma221213 + dGammadU31111*dMdGamma221311 + dGammadU31211*dMdGamma221312 + dGammadU31311*dMdGamma221313 + dMdPsi22111*dPhidU1111 + dMdPsi22121*dPhidU2111 + dMdPsi22131*dPhidU3111
    dMdU[I3,I11] = dGammadU11111*dMdGamma331111 + dGammadU11211*dMdGamma331112 + dGammadU11311*dMdGamma331113 + dGammadU21111*dMdGamma331211 + dGammadU21211*dMdGamma331212 + dGammadU21311*dMdGamma331213 + dGammadU31111*dMdGamma331311 + dGammadU31211*dMdGamma331312 + dGammadU31311*dMdGamma331313 + dMdPsi33111*dPhidU1111 + dMdPsi33121*dPhidU2111 + dMdPsi33131*dPhidU3111
    dMdU[I4,I11] = dGammadU11111*dMdGamma231111 + dGammadU11211*dMdGamma231112 + dGammadU11311*dMdGamma231113 + dGammadU21111*dMdGamma231211 + dGammadU21211*dMdGamma231212 + dGammadU21311*dMdGamma231213 + dGammadU31111*dMdGamma231311 + dGammadU31211*dMdGamma231312 + dGammadU31311*dMdGamma231313 + dMdPsi23111*dPhidU1111 + dMdPsi23121*dPhidU2111 + dMdPsi23131*dPhidU3111
    dMdU[I5,I11] = dGammadU11111*dMdGamma131111 + dGammadU11211*dMdGamma131112 + dGammadU11311*dMdGamma131113 + dGammadU21111*dMdGamma131211 + dGammadU21211*dMdGamma131212 + dGammadU21311*dMdGamma131213 + dGammadU31111*dMdGamma131311 + dGammadU31211*dMdGamma131312 + dGammadU31311*dMdGamma131313 + dMdPsi13111*dPhidU1111 + dMdPsi13121*dPhidU2111 + dMdPsi13131*dPhidU3111
    dMdU[I6,I11] = dGammadU11111*dMdGamma121111 + dGammadU11211*dMdGamma121112 + dGammadU11311*dMdGamma121113 + dGammadU21111*dMdGamma121211 + dGammadU21211*dMdGamma121212 + dGammadU21311*dMdGamma121213 + dGammadU31111*dMdGamma121311 + dGammadU31211*dMdGamma121312 + dGammadU31311*dMdGamma121313 + dMdPsi12111*dPhidU1111 + dMdPsi12121*dPhidU2111 + dMdPsi12131*dPhidU3111
    dMdU[I7,I11] = dGammadU11111*dMdGamma321111 + dGammadU11211*dMdGamma321112 + dGammadU11311*dMdGamma321113 + dGammadU21111*dMdGamma321211 + dGammadU21211*dMdGamma321212 + dGammadU21311*dMdGamma321213 + dGammadU31111*dMdGamma321311 + dGammadU31211*dMdGamma321312 + dGammadU31311*dMdGamma321313 + dMdPsi32111*dPhidU1111 + dMdPsi32121*dPhidU2111 + dMdPsi32131*dPhidU3111
    dMdU[I8,I11] = dGammadU11111*dMdGamma311111 + dGammadU11211*dMdGamma311112 + dGammadU11311*dMdGamma311113 + dGammadU21111*dMdGamma311211 + dGammadU21211*dMdGamma311212 + dGammadU21311*dMdGamma311213 + dGammadU31111*dMdGamma311311 + dGammadU31211*dMdGamma311312 + dGammadU31311*dMdGamma311313 + dMdPsi31111*dPhidU1111 + dMdPsi31121*dPhidU2111 + dMdPsi31131*dPhidU3111
    dMdU[I9,I11] = dGammadU11111*dMdGamma211111 + dGammadU11211*dMdGamma211112 + dGammadU11311*dMdGamma211113 + dGammadU21111*dMdGamma211211 + dGammadU21211*dMdGamma211212 + dGammadU21311*dMdGamma211213 + dGammadU31111*dMdGamma211311 + dGammadU31211*dMdGamma211312 + dGammadU31311*dMdGamma211313 + dMdPsi21111*dPhidU1111 + dMdPsi21121*dPhidU2111 + dMdPsi21131*dPhidU3111
    dMdU[I1+9,I11] = dGammadU11111*dMdGamma112111 + dGammadU11211*dMdGamma112112 + dGammadU11311*dMdGamma112113 + dGammadU21111*dMdGamma112211 + dGammadU21211*dMdGamma112212 + dGammadU21311*dMdGamma112213 + dGammadU31111*dMdGamma112311 + dGammadU31211*dMdGamma112312 + dGammadU31311*dMdGamma112313 + dMdPsi11211*dPhidU1111 + dMdPsi11221*dPhidU2111 + dMdPsi11231*dPhidU3111
    dMdU[I2+9,I11] = dGammadU11111*dMdGamma222111 + dGammadU11211*dMdGamma222112 + dGammadU11311*dMdGamma222113 + dGammadU21111*dMdGamma222211 + dGammadU21211*dMdGamma222212 + dGammadU21311*dMdGamma222213 + dGammadU31111*dMdGamma222311 + dGammadU31211*dMdGamma222312 + dGammadU31311*dMdGamma222313 + dMdPsi22211*dPhidU1111 + dMdPsi22221*dPhidU2111 + dMdPsi22231*dPhidU3111
    dMdU[I3+9,I11] = dGammadU11111*dMdGamma332111 + dGammadU11211*dMdGamma332112 + dGammadU11311*dMdGamma332113 + dGammadU21111*dMdGamma332211 + dGammadU21211*dMdGamma332212 + dGammadU21311*dMdGamma332213 + dGammadU31111*dMdGamma332311 + dGammadU31211*dMdGamma332312 + dGammadU31311*dMdGamma332313 + dMdPsi33211*dPhidU1111 + dMdPsi33221*dPhidU2111 + dMdPsi33231*dPhidU3111
    dMdU[I4+9,I11] = dGammadU11111*dMdGamma232111 + dGammadU11211*dMdGamma232112 + dGammadU11311*dMdGamma232113 + dGammadU21111*dMdGamma232211 + dGammadU21211*dMdGamma232212 + dGammadU21311*dMdGamma232213 + dGammadU31111*dMdGamma232311 + dGammadU31211*dMdGamma232312 + dGammadU31311*dMdGamma232313 + dMdPsi23211*dPhidU1111 + dMdPsi23221*dPhidU2111 + dMdPsi23231*dPhidU3111
    dMdU[I5+9,I11] = dGammadU11111*dMdGamma132111 + dGammadU11211*dMdGamma132112 + dGammadU11311*dMdGamma132113 + dGammadU21111*dMdGamma132211 + dGammadU21211*dMdGamma132212 + dGammadU21311*dMdGamma132213 + dGammadU31111*dMdGamma132311 + dGammadU31211*dMdGamma132312 + dGammadU31311*dMdGamma132313 + dMdPsi13211*dPhidU1111 + dMdPsi13221*dPhidU2111 + dMdPsi13231*dPhidU3111
    dMdU[I6+9,I11] = dGammadU11111*dMdGamma122111 + dGammadU11211*dMdGamma122112 + dGammadU11311*dMdGamma122113 + dGammadU21111*dMdGamma122211 + dGammadU21211*dMdGamma122212 + dGammadU21311*dMdGamma122213 + dGammadU31111*dMdGamma122311 + dGammadU31211*dMdGamma122312 + dGammadU31311*dMdGamma122313 + dMdPsi12211*dPhidU1111 + dMdPsi12221*dPhidU2111 + dMdPsi12231*dPhidU3111
    dMdU[I7+9,I11] = dGammadU11111*dMdGamma322111 + dGammadU11211*dMdGamma322112 + dGammadU11311*dMdGamma322113 + dGammadU21111*dMdGamma322211 + dGammadU21211*dMdGamma322212 + dGammadU21311*dMdGamma322213 + dGammadU31111*dMdGamma322311 + dGammadU31211*dMdGamma322312 + dGammadU31311*dMdGamma322313 + dMdPsi32211*dPhidU1111 + dMdPsi32221*dPhidU2111 + dMdPsi32231*dPhidU3111
    dMdU[I8+9,I11] = dGammadU11111*dMdGamma312111 + dGammadU11211*dMdGamma312112 + dGammadU11311*dMdGamma312113 + dGammadU21111*dMdGamma312211 + dGammadU21211*dMdGamma312212 + dGammadU21311*dMdGamma312213 + dGammadU31111*dMdGamma312311 + dGammadU31211*dMdGamma312312 + dGammadU31311*dMdGamma312313 + dMdPsi31211*dPhidU1111 + dMdPsi31221*dPhidU2111 + dMdPsi31231*dPhidU3111
    dMdU[I9+9,I11] = dGammadU11111*dMdGamma212111 + dGammadU11211*dMdGamma212112 + dGammadU11311*dMdGamma212113 + dGammadU21111*dMdGamma212211 + dGammadU21211*dMdGamma212212 + dGammadU21311*dMdGamma212213 + dGammadU31111*dMdGamma212311 + dGammadU31211*dMdGamma212312 + dGammadU31311*dMdGamma212313 + dMdPsi21211*dPhidU1111 + dMdPsi21221*dPhidU2111 + dMdPsi21231*dPhidU3111
    dMdU[I1+18,I11] = dGammadU11111*dMdGamma113111 + dGammadU11211*dMdGamma113112 + dGammadU11311*dMdGamma113113 + dGammadU21111*dMdGamma113211 + dGammadU21211*dMdGamma113212 + dGammadU21311*dMdGamma113213 + dGammadU31111*dMdGamma113311 + dGammadU31211*dMdGamma113312 + dGammadU31311*dMdGamma113313 + dMdPsi11311*dPhidU1111 + dMdPsi11321*dPhidU2111 + dMdPsi11331*dPhidU3111
    dMdU[I2+18,I11] = dGammadU11111*dMdGamma223111 + dGammadU11211*dMdGamma223112 + dGammadU11311*dMdGamma223113 + dGammadU21111*dMdGamma223211 + dGammadU21211*dMdGamma223212 + dGammadU21311*dMdGamma223213 + dGammadU31111*dMdGamma223311 + dGammadU31211*dMdGamma223312 + dGammadU31311*dMdGamma223313 + dMdPsi22311*dPhidU1111 + dMdPsi22321*dPhidU2111 + dMdPsi22331*dPhidU3111
    dMdU[I3+18,I11] = dGammadU11111*dMdGamma333111 + dGammadU11211*dMdGamma333112 + dGammadU11311*dMdGamma333113 + dGammadU21111*dMdGamma333211 + dGammadU21211*dMdGamma333212 + dGammadU21311*dMdGamma333213 + dGammadU31111*dMdGamma333311 + dGammadU31211*dMdGamma333312 + dGammadU31311*dMdGamma333313 + dMdPsi33311*dPhidU1111 + dMdPsi33321*dPhidU2111 + dMdPsi33331*dPhidU3111
    dMdU[I4+18,I11] = dGammadU11111*dMdGamma233111 + dGammadU11211*dMdGamma233112 + dGammadU11311*dMdGamma233113 + dGammadU21111*dMdGamma233211 + dGammadU21211*dMdGamma233212 + dGammadU21311*dMdGamma233213 + dGammadU31111*dMdGamma233311 + dGammadU31211*dMdGamma233312 + dGammadU31311*dMdGamma233313 + dMdPsi23311*dPhidU1111 + dMdPsi23321*dPhidU2111 + dMdPsi23331*dPhidU3111
    dMdU[I5+18,I11] = dGammadU11111*dMdGamma133111 + dGammadU11211*dMdGamma133112 + dGammadU11311*dMdGamma133113 + dGammadU21111*dMdGamma133211 + dGammadU21211*dMdGamma133212 + dGammadU21311*dMdGamma133213 + dGammadU31111*dMdGamma133311 + dGammadU31211*dMdGamma133312 + dGammadU31311*dMdGamma133313 + dMdPsi13311*dPhidU1111 + dMdPsi13321*dPhidU2111 + dMdPsi13331*dPhidU3111
    dMdU[I6+18,I11] = dGammadU11111*dMdGamma123111 + dGammadU11211*dMdGamma123112 + dGammadU11311*dMdGamma123113 + dGammadU21111*dMdGamma123211 + dGammadU21211*dMdGamma123212 + dGammadU21311*dMdGamma123213 + dGammadU31111*dMdGamma123311 + dGammadU31211*dMdGamma123312 + dGammadU31311*dMdGamma123313 + dMdPsi12311*dPhidU1111 + dMdPsi12321*dPhidU2111 + dMdPsi12331*dPhidU3111
    dMdU[I7+18,I11] = dGammadU11111*dMdGamma323111 + dGammadU11211*dMdGamma323112 + dGammadU11311*dMdGamma323113 + dGammadU21111*dMdGamma323211 + dGammadU21211*dMdGamma323212 + dGammadU21311*dMdGamma323213 + dGammadU31111*dMdGamma323311 + dGammadU31211*dMdGamma323312 + dGammadU31311*dMdGamma323313 + dMdPsi32311*dPhidU1111 + dMdPsi32321*dPhidU2111 + dMdPsi32331*dPhidU3111
    dMdU[I8+18,I11] = dGammadU11111*dMdGamma313111 + dGammadU11211*dMdGamma313112 + dGammadU11311*dMdGamma313113 + dGammadU21111*dMdGamma313211 + dGammadU21211*dMdGamma313212 + dGammadU21311*dMdGamma313213 + dGammadU31111*dMdGamma313311 + dGammadU31211*dMdGamma313312 + dGammadU31311*dMdGamma313313 + dMdPsi31311*dPhidU1111 + dMdPsi31321*dPhidU2111 + dMdPsi31331*dPhidU3111
    dMdU[I9+18,I11] = dGammadU11111*dMdGamma213111 + dGammadU11211*dMdGamma213112 + dGammadU11311*dMdGamma213113 + dGammadU21111*dMdGamma213211 + dGammadU21211*dMdGamma213212 + dGammadU21311*dMdGamma213213 + dGammadU31111*dMdGamma213311 + dGammadU31211*dMdGamma213312 + dGammadU31311*dMdGamma213313 + dMdPsi21311*dPhidU1111 + dMdPsi21321*dPhidU2111 + dMdPsi21331*dPhidU3111

    #Column 12
    dMdU[I1,I12] = dGammadU11112*dMdGamma111111 + dGammadU11212*dMdGamma111112 + dGammadU11312*dMdGamma111113 + dGammadU21112*dMdGamma111211 + dGammadU21212*dMdGamma111212 + dGammadU21312*dMdGamma111213 + dGammadU31112*dMdGamma111311 + dGammadU31212*dMdGamma111312 + dGammadU31312*dMdGamma111313 + dMdPsi11111*dPhidU1112 + dMdPsi11121*dPhidU2112 + dMdPsi11131*dPhidU3112
    dMdU[I2,I12] = dGammadU11112*dMdGamma221111 + dGammadU11212*dMdGamma221112 + dGammadU11312*dMdGamma221113 + dGammadU21112*dMdGamma221211 + dGammadU21212*dMdGamma221212 + dGammadU21312*dMdGamma221213 + dGammadU31112*dMdGamma221311 + dGammadU31212*dMdGamma221312 + dGammadU31312*dMdGamma221313 + dMdPsi22111*dPhidU1112 + dMdPsi22121*dPhidU2112 + dMdPsi22131*dPhidU3112
    dMdU[I3,I12] = dGammadU11112*dMdGamma331111 + dGammadU11212*dMdGamma331112 + dGammadU11312*dMdGamma331113 + dGammadU21112*dMdGamma331211 + dGammadU21212*dMdGamma331212 + dGammadU21312*dMdGamma331213 + dGammadU31112*dMdGamma331311 + dGammadU31212*dMdGamma331312 + dGammadU31312*dMdGamma331313 + dMdPsi33111*dPhidU1112 + dMdPsi33121*dPhidU2112 + dMdPsi33131*dPhidU3112
    dMdU[I4,I12] = dGammadU11112*dMdGamma231111 + dGammadU11212*dMdGamma231112 + dGammadU11312*dMdGamma231113 + dGammadU21112*dMdGamma231211 + dGammadU21212*dMdGamma231212 + dGammadU21312*dMdGamma231213 + dGammadU31112*dMdGamma231311 + dGammadU31212*dMdGamma231312 + dGammadU31312*dMdGamma231313 + dMdPsi23111*dPhidU1112 + dMdPsi23121*dPhidU2112 + dMdPsi23131*dPhidU3112
    dMdU[I5,I12] = dGammadU11112*dMdGamma131111 + dGammadU11212*dMdGamma131112 + dGammadU11312*dMdGamma131113 + dGammadU21112*dMdGamma131211 + dGammadU21212*dMdGamma131212 + dGammadU21312*dMdGamma131213 + dGammadU31112*dMdGamma131311 + dGammadU31212*dMdGamma131312 + dGammadU31312*dMdGamma131313 + dMdPsi13111*dPhidU1112 + dMdPsi13121*dPhidU2112 + dMdPsi13131*dPhidU3112
    dMdU[I6,I12] = dGammadU11112*dMdGamma121111 + dGammadU11212*dMdGamma121112 + dGammadU11312*dMdGamma121113 + dGammadU21112*dMdGamma121211 + dGammadU21212*dMdGamma121212 + dGammadU21312*dMdGamma121213 + dGammadU31112*dMdGamma121311 + dGammadU31212*dMdGamma121312 + dGammadU31312*dMdGamma121313 + dMdPsi12111*dPhidU1112 + dMdPsi12121*dPhidU2112 + dMdPsi12131*dPhidU3112
    dMdU[I7,I12] = dGammadU11112*dMdGamma321111 + dGammadU11212*dMdGamma321112 + dGammadU11312*dMdGamma321113 + dGammadU21112*dMdGamma321211 + dGammadU21212*dMdGamma321212 + dGammadU21312*dMdGamma321213 + dGammadU31112*dMdGamma321311 + dGammadU31212*dMdGamma321312 + dGammadU31312*dMdGamma321313 + dMdPsi32111*dPhidU1112 + dMdPsi32121*dPhidU2112 + dMdPsi32131*dPhidU3112
    dMdU[I8,I12] = dGammadU11112*dMdGamma311111 + dGammadU11212*dMdGamma311112 + dGammadU11312*dMdGamma311113 + dGammadU21112*dMdGamma311211 + dGammadU21212*dMdGamma311212 + dGammadU21312*dMdGamma311213 + dGammadU31112*dMdGamma311311 + dGammadU31212*dMdGamma311312 + dGammadU31312*dMdGamma311313 + dMdPsi31111*dPhidU1112 + dMdPsi31121*dPhidU2112 + dMdPsi31131*dPhidU3112
    dMdU[I9,I12] = dGammadU11112*dMdGamma211111 + dGammadU11212*dMdGamma211112 + dGammadU11312*dMdGamma211113 + dGammadU21112*dMdGamma211211 + dGammadU21212*dMdGamma211212 + dGammadU21312*dMdGamma211213 + dGammadU31112*dMdGamma211311 + dGammadU31212*dMdGamma211312 + dGammadU31312*dMdGamma211313 + dMdPsi21111*dPhidU1112 + dMdPsi21121*dPhidU2112 + dMdPsi21131*dPhidU3112
    dMdU[I1+9,I12] = dGammadU11112*dMdGamma112111 + dGammadU11212*dMdGamma112112 + dGammadU11312*dMdGamma112113 + dGammadU21112*dMdGamma112211 + dGammadU21212*dMdGamma112212 + dGammadU21312*dMdGamma112213 + dGammadU31112*dMdGamma112311 + dGammadU31212*dMdGamma112312 + dGammadU31312*dMdGamma112313 + dMdPsi11211*dPhidU1112 + dMdPsi11221*dPhidU2112 + dMdPsi11231*dPhidU3112
    dMdU[I2+9,I12] = dGammadU11112*dMdGamma222111 + dGammadU11212*dMdGamma222112 + dGammadU11312*dMdGamma222113 + dGammadU21112*dMdGamma222211 + dGammadU21212*dMdGamma222212 + dGammadU21312*dMdGamma222213 + dGammadU31112*dMdGamma222311 + dGammadU31212*dMdGamma222312 + dGammadU31312*dMdGamma222313 + dMdPsi22211*dPhidU1112 + dMdPsi22221*dPhidU2112 + dMdPsi22231*dPhidU3112
    dMdU[I3+9,I12] = dGammadU11112*dMdGamma332111 + dGammadU11212*dMdGamma332112 + dGammadU11312*dMdGamma332113 + dGammadU21112*dMdGamma332211 + dGammadU21212*dMdGamma332212 + dGammadU21312*dMdGamma332213 + dGammadU31112*dMdGamma332311 + dGammadU31212*dMdGamma332312 + dGammadU31312*dMdGamma332313 + dMdPsi33211*dPhidU1112 + dMdPsi33221*dPhidU2112 + dMdPsi33231*dPhidU3112
    dMdU[I4+9,I12] = dGammadU11112*dMdGamma232111 + dGammadU11212*dMdGamma232112 + dGammadU11312*dMdGamma232113 + dGammadU21112*dMdGamma232211 + dGammadU21212*dMdGamma232212 + dGammadU21312*dMdGamma232213 + dGammadU31112*dMdGamma232311 + dGammadU31212*dMdGamma232312 + dGammadU31312*dMdGamma232313 + dMdPsi23211*dPhidU1112 + dMdPsi23221*dPhidU2112 + dMdPsi23231*dPhidU3112
    dMdU[I5+9,I12] = dGammadU11112*dMdGamma132111 + dGammadU11212*dMdGamma132112 + dGammadU11312*dMdGamma132113 + dGammadU21112*dMdGamma132211 + dGammadU21212*dMdGamma132212 + dGammadU21312*dMdGamma132213 + dGammadU31112*dMdGamma132311 + dGammadU31212*dMdGamma132312 + dGammadU31312*dMdGamma132313 + dMdPsi13211*dPhidU1112 + dMdPsi13221*dPhidU2112 + dMdPsi13231*dPhidU3112
    dMdU[I6+9,I12] = dGammadU11112*dMdGamma122111 + dGammadU11212*dMdGamma122112 + dGammadU11312*dMdGamma122113 + dGammadU21112*dMdGamma122211 + dGammadU21212*dMdGamma122212 + dGammadU21312*dMdGamma122213 + dGammadU31112*dMdGamma122311 + dGammadU31212*dMdGamma122312 + dGammadU31312*dMdGamma122313 + dMdPsi12211*dPhidU1112 + dMdPsi12221*dPhidU2112 + dMdPsi12231*dPhidU3112
    dMdU[I7+9,I12] = dGammadU11112*dMdGamma322111 + dGammadU11212*dMdGamma322112 + dGammadU11312*dMdGamma322113 + dGammadU21112*dMdGamma322211 + dGammadU21212*dMdGamma322212 + dGammadU21312*dMdGamma322213 + dGammadU31112*dMdGamma322311 + dGammadU31212*dMdGamma322312 + dGammadU31312*dMdGamma322313 + dMdPsi32211*dPhidU1112 + dMdPsi32221*dPhidU2112 + dMdPsi32231*dPhidU3112
    dMdU[I8+9,I12] = dGammadU11112*dMdGamma312111 + dGammadU11212*dMdGamma312112 + dGammadU11312*dMdGamma312113 + dGammadU21112*dMdGamma312211 + dGammadU21212*dMdGamma312212 + dGammadU21312*dMdGamma312213 + dGammadU31112*dMdGamma312311 + dGammadU31212*dMdGamma312312 + dGammadU31312*dMdGamma312313 + dMdPsi31211*dPhidU1112 + dMdPsi31221*dPhidU2112 + dMdPsi31231*dPhidU3112
    dMdU[I9+9,I12] = dGammadU11112*dMdGamma212111 + dGammadU11212*dMdGamma212112 + dGammadU11312*dMdGamma212113 + dGammadU21112*dMdGamma212211 + dGammadU21212*dMdGamma212212 + dGammadU21312*dMdGamma212213 + dGammadU31112*dMdGamma212311 + dGammadU31212*dMdGamma212312 + dGammadU31312*dMdGamma212313 + dMdPsi21211*dPhidU1112 + dMdPsi21221*dPhidU2112 + dMdPsi21231*dPhidU3112
    dMdU[I1+18,I12] = dGammadU11112*dMdGamma113111 + dGammadU11212*dMdGamma113112 + dGammadU11312*dMdGamma113113 + dGammadU21112*dMdGamma113211 + dGammadU21212*dMdGamma113212 + dGammadU21312*dMdGamma113213 + dGammadU31112*dMdGamma113311 + dGammadU31212*dMdGamma113312 + dGammadU31312*dMdGamma113313 + dMdPsi11311*dPhidU1112 + dMdPsi11321*dPhidU2112 + dMdPsi11331*dPhidU3112
    dMdU[I2+18,I12] = dGammadU11112*dMdGamma223111 + dGammadU11212*dMdGamma223112 + dGammadU11312*dMdGamma223113 + dGammadU21112*dMdGamma223211 + dGammadU21212*dMdGamma223212 + dGammadU21312*dMdGamma223213 + dGammadU31112*dMdGamma223311 + dGammadU31212*dMdGamma223312 + dGammadU31312*dMdGamma223313 + dMdPsi22311*dPhidU1112 + dMdPsi22321*dPhidU2112 + dMdPsi22331*dPhidU3112
    dMdU[I3+18,I12] = dGammadU11112*dMdGamma333111 + dGammadU11212*dMdGamma333112 + dGammadU11312*dMdGamma333113 + dGammadU21112*dMdGamma333211 + dGammadU21212*dMdGamma333212 + dGammadU21312*dMdGamma333213 + dGammadU31112*dMdGamma333311 + dGammadU31212*dMdGamma333312 + dGammadU31312*dMdGamma333313 + dMdPsi33311*dPhidU1112 + dMdPsi33321*dPhidU2112 + dMdPsi33331*dPhidU3112
    dMdU[I4+18,I12] = dGammadU11112*dMdGamma233111 + dGammadU11212*dMdGamma233112 + dGammadU11312*dMdGamma233113 + dGammadU21112*dMdGamma233211 + dGammadU21212*dMdGamma233212 + dGammadU21312*dMdGamma233213 + dGammadU31112*dMdGamma233311 + dGammadU31212*dMdGamma233312 + dGammadU31312*dMdGamma233313 + dMdPsi23311*dPhidU1112 + dMdPsi23321*dPhidU2112 + dMdPsi23331*dPhidU3112
    dMdU[I5+18,I12] = dGammadU11112*dMdGamma133111 + dGammadU11212*dMdGamma133112 + dGammadU11312*dMdGamma133113 + dGammadU21112*dMdGamma133211 + dGammadU21212*dMdGamma133212 + dGammadU21312*dMdGamma133213 + dGammadU31112*dMdGamma133311 + dGammadU31212*dMdGamma133312 + dGammadU31312*dMdGamma133313 + dMdPsi13311*dPhidU1112 + dMdPsi13321*dPhidU2112 + dMdPsi13331*dPhidU3112
    dMdU[I6+18,I12] = dGammadU11112*dMdGamma123111 + dGammadU11212*dMdGamma123112 + dGammadU11312*dMdGamma123113 + dGammadU21112*dMdGamma123211 + dGammadU21212*dMdGamma123212 + dGammadU21312*dMdGamma123213 + dGammadU31112*dMdGamma123311 + dGammadU31212*dMdGamma123312 + dGammadU31312*dMdGamma123313 + dMdPsi12311*dPhidU1112 + dMdPsi12321*dPhidU2112 + dMdPsi12331*dPhidU3112
    dMdU[I7+18,I12] = dGammadU11112*dMdGamma323111 + dGammadU11212*dMdGamma323112 + dGammadU11312*dMdGamma323113 + dGammadU21112*dMdGamma323211 + dGammadU21212*dMdGamma323212 + dGammadU21312*dMdGamma323213 + dGammadU31112*dMdGamma323311 + dGammadU31212*dMdGamma323312 + dGammadU31312*dMdGamma323313 + dMdPsi32311*dPhidU1112 + dMdPsi32321*dPhidU2112 + dMdPsi32331*dPhidU3112
    dMdU[I8+18,I12] = dGammadU11112*dMdGamma313111 + dGammadU11212*dMdGamma313112 + dGammadU11312*dMdGamma313113 + dGammadU21112*dMdGamma313211 + dGammadU21212*dMdGamma313212 + dGammadU21312*dMdGamma313213 + dGammadU31112*dMdGamma313311 + dGammadU31212*dMdGamma313312 + dGammadU31312*dMdGamma313313 + dMdPsi31311*dPhidU1112 + dMdPsi31321*dPhidU2112 + dMdPsi31331*dPhidU3112
    dMdU[I9+18,I12] = dGammadU11112*dMdGamma213111 + dGammadU11212*dMdGamma213112 + dGammadU11312*dMdGamma213113 + dGammadU21112*dMdGamma213211 + dGammadU21212*dMdGamma213212 + dGammadU21312*dMdGamma213213 + dGammadU31112*dMdGamma213311 + dGammadU31212*dMdGamma213312 + dGammadU31312*dMdGamma213313 + dMdPsi21311*dPhidU1112 + dMdPsi21321*dPhidU2112 + dMdPsi21331*dPhidU3112

    return dMdU
    
class TestMicroElement(unittest.TestCase):

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
        
        self.assertEqual(np.allclose(F,Fanalytic),True)
        self.assertEqual(np.allclose(chi,chia),True)
        self.assertEqual(np.allclose(grad_chi,grad_chia),True)
        
        
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
            
        self.assertEqual(np.allclose(F,Fanalytic),True)
        
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
        
        self.assertEqual(np.allclose(chi,chia),True)
        
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
        
        self.assertEqual(np.allclose(grad_chi,grad_chia),True)
        
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
            return self.SOTtensor_to_vector(F)
            
        #Compute the gradients
        dFdUn = fd.numeric_gradient(F_parser,U,1e-6)
        dFdU  = compute_dFdU(xi_vec,rcoords)
        
        #tmp = 4
        #print "Numeric"
        #print dFdUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dFdU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dFdUn.T[:,(12*tmp):12*(tmp+1)] - dFdU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dFdUn.T,dFdU),True)
        
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
            C = hex8.matrix_Tdot(F,F)
            return self.SOTtensor_to_vector(C)
            
        #Compute the gradients
        dCdUn = fd.numeric_gradient(C_parser,U,1e-6)
        F     = compute_F(xi_vec,ccoords,rcoords)
        dFdU  = compute_dFdU(xi_vec,rcoords)
        dCdU  = compute_dCdU(F,dFdU)
        
        #tmp = 0
        #print "Numeric"
        #print dCdUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dCdU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dCdUn.T[:,(12*tmp):12*(tmp+1)] - dCdU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dCdUn.T,dCdU),True)
        
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
            Cinv = hex8.invert_3x3_matrix(hex8.matrix_Tdot(F,F))[1]
            return self.SOTtensor_to_vector(Cinv)
            
        #Compute the gradients
        dCinvdUn = fd.numeric_gradient(Cinv_parser,U,1e-6)
        F        = compute_F(xi_vec,ccoords,rcoords)
        Cinv     = hex8.invert_3x3_matrix(hex8.matrix_Tdot(F,F))[1]
        dFdU     = compute_dFdU(xi_vec,rcoords)
        dCdU     = compute_dCdU(F,dFdU)
        dCinvdU  = compute_dCinvdU(Cinv,dCdU)
        
        #tmp = 0
        #print "Numeric"
        #print dCinvdUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dCinvdU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dCinvdUn.T[:,(12*tmp):12*(tmp+1)] - dCinvdU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dCinvdUn.T,dCinvdU),True)
        
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
            return self.SOTtensor_to_vector(chi)
            
        #Compute the gradients
        dchidUn = fd.numeric_gradient(chi_parser,U,1e-6)
        dchidU  = compute_dchidU(xi_vec)
        
        #tmp = 4
        #print "Numeric"
        #print dchidUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dchidU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dchidUn.T[:,(12*tmp):12*(tmp+1)] - dchidU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dchidUn.T,dchidU),True)
        
    def test_compute_dPhidU(self):
        """Test the computation of the matrix form of the derivative of the
        deformation measure Phi with respect to the degree of freedom vector"""
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
            Phi     = hex8.matrix_Tdot(F,chi)
            return self.SOTtensor_to_vector(Phi)
            
        #Compute required measures
        F      = compute_F(xi_vec,ccoords,rcoords)
        dFdU   = compute_dFdU(xi_vec,rcoords)
        chi    = compute_chi(xi_vec,phi_vectors)
        dchidU = compute_dchidU(xi_vec)
            
        #Compute the gradients
        dPhidUn = fd.numeric_gradient(Phi_parser,U,1e-6)
        dPhidU  = compute_dPhidU(F,chi,dFdU,dchidU)
        
        #tmp = 0
        #print "Numeric"
        #print dPhidUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dPhidU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dPhidUn.T[:,(12*tmp):12*(tmp+1)] - dPhidU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dPhidUn.T,dPhidU),True)
        
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
            return self.TOTtensor_to_vector(grad_chi)
            
        #Compute the gradients
        dgrad_chidUn = fd.numeric_gradient(grad_chi_parser,U,1e-6)
        dgrad_chidU  = compute_dgrad_chidU(xi_vec,rcoords)
        
        #tmp = 4
        #print "Numeric"
        #print dgrad_chidUn.T[:,(12*tmp):12*(tmp+1)]
        #print "Code"
        #print dgrad_chidU[:,(12*tmp):12*(tmp+1)]
        #print "Difference"
        #print dgrad_chidUn.T[:,(12*tmp):12*(tmp+1)] - dgrad_chidU[:,(12*tmp):12*(tmp+1)]
        
        self.assertEqual(np.allclose(dgrad_chidUn.T,dgrad_chidU),True)
        
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
            Gamma    = hex8.matrix_Tdot_TOT(F,grad_chi)
            return self.TOTtensor_to_vector(Gamma)
        
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
        
        self.assertEqual(np.allclose(dGammadUn.T,dGammadU),True)
        
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
    