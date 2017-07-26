import sympy
from sympy.tensor.array import MutableDenseNDimArray as Array
import numpy as np
import operator
import functools

def form_nd_array(symb,dims):
    """Form a symbolic nd array with the given symbol over the given dimensions"""
    indices = create_indices(dims)
    components = get_components(symb,indices,dims)
    return Array(sympy.symbols(components),dims,order='C')
    
def multiply_list(l):
    """Multiply a list of elements together"""
    return functools.reduce(operator.mul, l, 1)
    
def create_indices(dims):
    """Create lists of indices"""
    return [range(1,dim+1) for dim in dims]
        
def dyadic_string_lists(l1,l2):
    """Compute the dyadic product of two lists of strings"""
    if(len(l1)==1):
        l1 = ([l1[0],[]],)
    
    l = [[[str(l1v[0])+str(l2v),[t for t in l1v[1]]+[l2v]] for l1v in l1] for l2v in l2]
    
    return [item for sublist in l for item in sublist]
    
def get_components(symb,indices,dims):
    """Get the components of the array"""
    list = [symb]
    for index in indices:
        list = dyadic_string_lists(list,index)
    
    return order_list(list)
    
def order_list(list):
    """Sort the list into the correct order"""
    names,indices = zip(*list)
    index_vectors = zip(*indices)
    
    #print names
    #print index_vectors
    
    #print names
    #print index_vectors[1]
    #
    #names = sort_list_by_list(names,index_vectors[0])
    #
    #print index_vectors[0]
    #
    #index_vectors[-1] = sort_list_by_list(index_vectors[-1],index_vectors[0])
    #
    #print names
    #print index_vectors[-1]
    #
    #names = sort_list_by_list(names,index_vectors[-1])
    #print names
    #return names
    
    for i in reversed(range(len(index_vectors)-1)):
        names = sort_list_by_list(names,index_vectors[i])
        for j in range(len(index_vectors)):
            if(j!=i):
                index_vectors[j] = sort_list_by_list(index_vectors[j],index_vectors[i])
        
        index_vectors[i] = sort_list_by_list(index_vectors[i],index_vectors[i])
        #print names
        #print index_vectors
    
    #print names
    #raise
    
    return names
    
def get_dFdU():
    """Get dFdU"""
    dFdU = form_nd_array("dFdU",[3,3,12*8])
    zero = sympy.S(0)
    
    #Remove terms associated with chi
    for k in range(3,12*8,12):
        for i in range(3):
            for j in range(3):
                for t in range(k,k+9):
                    dFdU[i,j,t] = zero
    
    #Remove additional zero terms
    for k in range(0,12*8,12):
        for i in range(3):
            #Remove derivatives w.r.t. 2 and 3 from dof 1
            dFdU[1,i,k]   = zero
            dFdU[2,i,k]   = zero
            
            #Remove derivatives w.r.t. 1 and 3 from dof 2
            dFdU[0,i,k+1] = zero
            dFdU[2,i,k+1] = zero
            
            #Remove derivatives w.r.t. 1 and 2 from dof 3
            dFdU[0,i,k+2] = zero
            dFdU[1,i,k+2] = zero
    return dFdU
    
def get_dchidU():
    """Get dPhidU"""
    dchidU = form_nd_array("dchidU",[3,3,12*8])
    zero = sympy.S(0)
    
    #Define the dof to chi mapping
    dof_map = [(-1,-1),(-1,-1),(-1,-1),(0,0),(1,1),(2,2),(1,2),(0,2),(0,1),(2,1),(2,0),(1,0)]
    
    #Remove values associated with x
    #for k in range(0,12*8,12):
    #    for i in range(3):
    #        for j in range(3):
    #            dchidU[i,j,  k] = zero
    #            dchidU[i,j,k+1] = zero
    #            dchidU[i,j,k+2] = zero
    
    #Remove values which are zero
    for k in range(0,12*8):
        idof,jdof = dof_map[k%12]
        for i in range(3):
            for j in range(3):
                if((idof!=i) or (jdof!=j)):
                    dchidU[i,j,k] = zero
    return dchidU
    
def form_symb_dCdU():
    """Form a symbolic version of dCdU"""
    dCdU    = form_nd_array("dCdU",[3,3,8*12])
    for I in range(3):
        for J in range(3):
            for K in range(3,8*12):
                dCdU[I,J,K] = 0
    return dCdU
    
def form_symb_dPhidU():
    """Form a symbolic version of dPhidU"""
    dPhidU = get_dPhidU()
    dPhidUs = form_nd_array("dPhidU",[3,3,8*12])
    for I in range(3):
        for J in range(3):
            for K in range(8*12):
                if((dPhidU[I,J,K]==0) or (dPhidU[I,J,K]==sympy.S(0))):
                    dPhidUs[I,J,K] = 0
    return dPhidUs
    
def form_symb_dGammadU():
    """Form a symbolic version of dGammadU"""
    dGammadU = get_dGammadU()
    dGammadUs = form_nd_array("dGammadU",[3,3,3,8*12])
    for I in range(3):
        for J in range(3):
            for K in range(3):
                for L in range(8*12):
                    if((dGammadU[I,J,K,L]==0) or (dGammadU[I,J,K,L]==sympy.S(0))):
                        dGammadUs[I,J,K,L] = 0
    return dGammadUs
def get_dpk2dU():
    """Compute dpk2dU"""
    dpk2dC     = form_nd_array("dpk2dC",    [3,3,3,3])
    dpk2dPsi   = form_nd_array("dpk2dPsi",  [3,3,3,3])
    dpk2dGamma = form_nd_array("dpk2dGamma",[3,3,3,3,3])
    
    dCdU     = form_symb_dCdU()
    dPsidU   = form_symb_dPhidU()
    dGammadU = form_symb_dGammadU()
    
    dpk2dU = form_nd_array("dpk2dU",[3,3,8*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(8*12):
            
                dpk2dU[I,J,K] = 0
            
                for L in range(3):
                    for M in range(3):
                    
                        dpk2dU[I,J,K] += dpk2dC[I,J,L,M]*dCdU[L,M,K] + dpk2dPsi[I,J,L,M]*dPsidU[L,M,K]
                    
                        for N in range(3):
                            dpk2dU[I,J,K] += dpk2dGamma[I,J,L,M,N]*dGammadU[L,M,N,K]
                            
    print "    #Implementation of tangent"
    
    tmp  = [get_matrix_form_TOT(dCdU)[:,:12],     get_matrix_form_TOT(dPsidU)[:,:12],\
            get_matrix_form_FOT(dGammadU)[:,:12], get_matrix_form_FOT(dpk2dC),\
            get_matrix_form_FOT(dpk2dPsi),        get_matrix_form_VOT(dpk2dGamma)]
    symb = ["dCdU","dPsidU","dGammadU","dpk2dC","dpk2dPsi","dpk2dGamma"]
    [implementation_extract_matrix(t,s,"I","I") for t,s in zip(tmp,symb)]
    implementation_print_matrix(get_matrix_form_TOT(dpk2dU)[:,:12],"dpk2dU","I","I")
    return dpk2dU
    
def get_dSigmadU():
    """Compute dSigmadU
    Note: Identical to dpk2dU except for stress tangents used"""
    dSigmadC     = form_nd_array("dSigmadC",    [3,3,3,3])
    dSigmadPsi   = form_nd_array("dSigmadPsi",  [3,3,3,3])
    dSigmadGamma = form_nd_array("dSigmadGamma",[3,3,3,3,3])
    
    dCdU     = form_symb_dCdU()
    dPsidU   = form_symb_dPhidU()
    dGammadU = form_symb_dGammadU()
    
    dSigmadU = form_nd_array("dSigmadU",[3,3,8*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(8*12):
            
                dSigmadU[I,J,K] = 0
            
                for L in range(3):
                    for M in range(3):
                    
                        dSigmadU[I,J,K] += dSigmadC[I,J,L,M]*dCdU[L,M,K] + dSigmadPsi[I,J,L,M]*dPsidU[L,M,K]
                    
                        for N in range(3):
                            dSigmadU[I,J,K] += dSigmadGamma[I,J,L,M,N]*dGammadU[L,M,N,K]
                            
    print "    #Implementation of tangent"
    
    tmp  = [get_matrix_form_TOT(dCdU)[:,:12],     get_matrix_form_TOT(dPsidU)[:,:12],\
            get_matrix_form_FOT(dGammadU)[:,:12], get_matrix_form_FOT(dSigmadC),\
            get_matrix_form_FOT(dSigmadPsi),        get_matrix_form_VOT(dSigmadGamma)]
    symb = ["dCdU","dPsidU","dGammadU","dSigmadC","dSigmadPsi","dSigmadGamma"]
    [implementation_extract_matrix(t,s,"I","I") for t,s in zip(tmp,symb)]
    implementation_print_matrix(get_matrix_form_TOT(dSigmadU)[:,:12],"dSigmadU","I","I")
    return dSigmadU
    
def get_dMdU():
    """Compute dMdU"""
    
    dMdC     = form_nd_array("dMdC",    [3,3,3,3,3])
    dMdPsi   = form_nd_array("dMdPsi",  [3,3,3,3,3])
    dMdGamma = form_nd_array("dMdGamma",[3,3,3,3,3,3])
    
    dCdU     = form_symb_dCdU()
    dPsidU   = form_symb_dPhidU()
    dGammadU = form_symb_dGammadU()
    
    dMdU = form_nd_array("dMdU",[3,3,3,8*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(3):
                for L in range(8*12):
                    dMdU[I,J,K,L] = 0
                    for O in range(3):
                        for P in range(3):
                            dMdU[I,J,K,L] += dMdC[I,J,K,O,P]*dCdU[O,P,L] + dMdPsi[I,J,K,O,P]*dPsidU[O,P,L]
                            
                            for Q in range(3):
                                dMdU[I,J,K,L] += dMdGamma[I,J,K,O,P,Q]*dGammadU[O,P,Q,L]
    
    tmp  = [get_matrix_form_TOT(dCdU)[:,:12],     get_matrix_form_TOT(dPsidU)[:,:12],\
            get_matrix_form_FOT(dGammadU)[:,:12], get_matrix_form_VOT(dMdC),\
            get_matrix_form_VOT(dMdPsi),        get_matrix_form_VIOT(dMdGamma)]
    symb = ["dCdU","dPsidU","dGammadU","dMdC","dMdPsi","dMdGamma"]
    [implementation_extract_matrix(t,s,"I","I") for t,s in zip(tmp,symb)]
    implementation_print_matrix(get_matrix_form_FOT(dMdU)[:,:12],"dMdU","I","I")
    return dMdU
                
def get_dgrad_chidU():
    """Get dgrad_chidU"""
    dgrad_chidU = form_nd_array("dgrad_chidU",[3,3,3,12*8])
    zero = sympy.S(0)
    
    #Define the dof to chi mapping
    dof_map = [(-1,-1),(-1,-1),(-1,-1),(0,0),(1,1),(2,2),(1,2),(0,2),(0,1),(2,1),(2,0),(1,0)]
    
    for k in range(0,12*8):
        idof,jdof = dof_map[k%12]
        for i in range(3):
            for j in range(3):
                if((idof!=i) or (jdof!=j)):
                    for l in range(3):
                        dgrad_chidU[i,j,l,k] = zero
    return dgrad_chidU
    
def get_dCdU():
    """Compute and return dCdU"""
    F    = form_nd_array("F",[3,3])
    dFdU = get_dFdU()
    
    dCdU = form_nd_array("dCdU",[3,3,12*8])
    
    for I in range(3):
        for J in range(3):
            for K in range(12*8):
                dCdU[I,J,K] = 0
                for i in range(3):
                    dCdU[I,J,K] += F[i,J]*dFdU[i,I,K]+F[i,I]*dFdU[i,J,K]
                    
    return dCdU
    
def get_dPhidU():
    """Compute and return dPhidU"""
    F      = form_nd_array("F",[3,3])
    chi    = form_nd_array("chi",[3,3])
    dFdU   = get_dFdU()
    dchidU = get_dchidU()
    
    dPhidU = form_nd_array("dCdU",[3,3,12*8])
    
    for I in range(3):
        for J in range(3):
            for K in range(12*8):
                dPhidU[I,J,K] = 0
                for i in range(3):
                    dPhidU[I,J,K] += chi[i,J]*dFdU[i,I,K]+F[i,I]*dchidU[i,J,K]
                    
    return dPhidU
    
def get_dGammadU():
    """Compute and return dGammadU"""
    F           = form_nd_array("F",[3,3])
    grad_chi    = form_nd_array("grad_chi",[3,3,3])
    dFdU        = get_dFdU()
    dgrad_chidU = get_dgrad_chidU()
    
    dGammadU   = form_nd_array("dGamma_dU",[3,3,3,8*12])
    
    for I in range(3):
        for J in range(3):
            for L in range(3):
                for K in range(8*12):
                    dGammadU[I,J,L,K] = 0
                    for i in range(3):
                        dGammadU[I,J,L,K] += grad_chi[i,J,L]*dFdU[i,I,K] + F[i,I]*dgrad_chidU[i,J,L,K]
    return dGammadU
    
def get_dCinvdU():
    """Compute and return the derivative of the inverse of 
    C with respect to U"""
    dCinvdU = form_nd_array("dCinvdU",[3,3,8*12])
    Cinv    = form_nd_array("Cinv",[3,3])
    dCdU    = form_nd_array("dCdU",[3,3,8*12])
    
    for I in range(3):
        for J in range(3):
            for K in range(3,8*12):
                dCdU[I,J,K] = 0
    
    for I in range(3):
        for J in range(3):
            for K in range(8*12):
                dCinvdU[I,J,K] = 0
                for L in range(3):
                    for M in range(3):
                        dCinvdU[I,J,K] -= Cinv[I,L]*Cinv[M,J]*dCdU[L,M,K]
                        
    return dCinvdU
    
def implementation_extract_matrix(A,Asymb,I1symb,I2symb,fname='./extract_temp.txt'):
    """Print the extraction of a matrix out so that the implementation in code is easier"""
    
    zero = sympy.S(0)
    M,N = A.shape
    
    f1=open(fname, 'a+')
    f1.write("\n    #Extract components of {0}\n".format(Asymb))
    
    for J in range(N):
        for I in range(M):
            if((A[I,J] != 0) or (A[I,J] != zero)):
                vals = [(I/9)*9, I%9+1, (J/12)*12, J%12+1]
                vals = [str(v) if v>0 else "" for v in vals]
                Itmp1,Itmp2,Jtmp1,Jtmp2 = vals
                if(Itmp1!=""):
                    Itmp1 = "+"+Itmp1
                if(Jtmp1!=""):
                    Jtmp1 = "+"+Jtmp1
                f1.write("    {0} = {1}[{2},{3}]\n".format(A[I,J],Asymb,I1symb+str(Itmp2)+str(Itmp1),I2symb+str(Jtmp2)+str(Jtmp1)))
    
    f1.close()
    
def implementation_print_matrix(A,Asymb,I1symb,I2symb,fname='./implement_temp.txt'):
    """Print a matrix out so that the implementation in code is easier"""
    
    zero = sympy.S(0)
    M,N = A.shape
    
    f1=open(fname,'a+')
    
    for J in range(N):
        f1.write("\n    #Column {0}\n".format(J+1))
        for I in range(M):
            if((A[I,J] != 0) or (A[I,J] != zero)):
                vals = [(I/9)*9, I%9+1, (J/12)*12, J%12+1]
                vals = [str(v) if v>0 else "" for v in vals]
                Itmp1,Itmp2,Jtmp1,Jtmp2 = vals
                if(Itmp1!=""):
                    Itmp1 = "+"+Itmp1
                if(Jtmp1!=""):
                    Jtmp1 = "+"+Jtmp1
                f1.write("    {0}[{1},{2}] = {3}\n".format(Asymb,I1symb+str(Itmp2)+str(Itmp1),I2symb+str(Jtmp2)+str(Jtmp1),A[I,J]))
    f1.close()
    
def get_matrix_form_SOT(SOT):
    """Get the matrix form of a second order tensor"""
    A,B = SOT.shape
    MF = form_nd_array("MF",(A*B,))
    MF[0] = SOT[0,0]
    MF[1] = SOT[1,1]
    MF[2] = SOT[2,2]
    MF[3] = SOT[1,2]
    MF[4] = SOT[0,2]
    MF[5] = SOT[0,1]
    MF[6] = SOT[2,1]
    MF[7] = SOT[2,0]
    MF[8] = SOT[1,0]
    return MF
    
def get_matrix_form_TOT(TOT):
    """Get the matrix form of a third order tensor"""
    A,B,C = TOT.shape
    MF = form_nd_array("MF",(A*B,C))

    for j in range(C):
        VTMP = get_matrix_form_SOT(TOT[:,:,j])
        for i in range(A*B):
            MF[i,j] = VTMP[i]
    return MF
    
def get_matrix_form_FOT(FOT):
    """Get the matrix form of a fourth order tensor"""
    A,B,C,D = FOT.shape
    MF = form_nd_array("MF",(A*B*C,D))
    
    for j in range(D):
        VTMP = form_nd_array("MF",(A*B*C,))
        for k in range(3):
            tmp = get_matrix_form_SOT(FOT[:,:,k,j])
            for i in range(A*B):
                MF[i+k*A*B,j] = tmp[i]
    return MF
    
def get_matrix_form_VOT(VOT):
    """Get the matrix form of a fifth order tensor"""
    A,B,C,D,E = VOT.shape
    MF = form_nd_array("MF",(A*B*C*D,E))
    
    for j in range(E):
        VTMP = form_nd_array("MF",(A*B*C*D,))
        for k in range(C):
            for l in range(D):
                tmp = get_matrix_form_SOT(VOT[:,:,k,l,j])
                
                for i in range(A*B):
                    MF[i+(D*k+l)*A*B,j] = tmp[i]
    return MF
    
def get_matrix_form_VIOT(VIOT):
    """Get the matrix form of a sixth order tensor"""
    A,B,C,D,E,F = VIOT.shape
    MF = form_nd_array("MF",(A*B*C*D*E,F))
    
    for j in range(E):
        VTMP = form_nd_array("MF",(A*B*C*D,))
        for k in range(C):
            for l in range(D):
                for m in range(E):
                    tmp = get_matrix_form_SOT(VIOT[:,:,k,l,m,j])
                    
                    for i in range(A*B):
                        MF[i+A*B*(m+l*E+k*D*E),j] = tmp[i]
    return MF
    
def check_dCdU_vs_matrix_form():
    """Check if the matrix xform of dCdU is correct"""
    dCdU  = get_dCdU()
    F     = form_nd_array("F",[3,3])
    dFdU  = get_dFdU("dFdU",[3,3,12*8])
    VF    = get_vector_form_SOT(F)
    MdFdU = get_matrix_form_TOT(dFdU)
    
def sort_list_by_list(L1,L2):
    """Sort a list by another list"""
    return [x for (y,x) in sorted(zip(L2,L1), key=lambda pair: pair[0])]


    
    