import numpy as np
import hex8
import finite_difference as fd
import unittest
import os

from hex8 import T_to_V_mapping as T2V
from hex8 import V_to_T_mapping as V2T

"""
==============================================
|                                            |
|       Micromorphic Linear Elasticity       |
|                                            |
==============================================
|                                            |
| A constitutive model for micromorphic      |
| linear elasticity. Written in a style      |
| which will be simple to port to Abaqus     |
| for use in finite element simulations.     |
|                                            |
| Model derived from a quadratic form of     |
| the Helmholtz free energy detailed in      |
| "Nonlinear Micromorphic Continuum          |
| Mechanics and Finite Strain                |
| Elastoplasticity" by Dr. Richard Regueiro  |
| with expansions by Nathan Miller           |
==============================================
"""
#Compute deformation measures

def get_deformation_measures(F,chi,grad_chi): #Test function written
    """Compute and return the deformation measures"""
    C   = hex8.matrix_Tdot_V(F,F)
    Psi = hex8.matrix_Tdot_V(F,chi)
    
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
    return C,Psi,Gamma
    
#Compute strain measures

def compute_strain_measures(C,Psi): #Test function written
    """Compute the strain measures"""
    Iten = np.eye(3)
    E_macro = np.zeros([9,])
    E_micro = np.zeros([9,])
    for index in range(len(E_macro)):
        I,J = V2T(index,[3,3])
        E_macro[index] = 0.5*(  C[index]-Iten[I,J])
        E_micro[index] =      Psi[index]-Iten[I,J]
    return E_macro,E_micro
    
def compute_dCinvdC(Cinv): #Test function written
    """Compute the derivative of Cinv w.r.t. C"""
    dCinvdC = np.zeros([3*3*3*3])
    for I in range(3):
        for J in range(3):
            for K in range(3):
                for L in range(3):
                    IJKL = T2V([I,J,K,L],[3,3,3,3])
                    IK   = T2V([I,K],[3,3])
                    LJ   = T2V([L,J],[3,3])
                    dCinvdC[IJKL] = -Cinv[IK]*Cinv[LJ]
    return dCinvdC

###### Form Stiffness Tensors ######

def form_stiffness_tensors(PROPS): #Test function written
    """Form the stiffness tensors"""
    RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11 = PROPS
    AFOT = form_A(LAMBDA,MU)
    BFOT = form_B(ETA,TAU,KAPPA,NU,SIGMA)
    CFOT = form_C(TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11)
    DFOT = form_D(TAU,SIGMA)
    return AFOT,BFOT,CFOT,DFOT
    
def form_A(LAMBDA,MU): #Test function written
    """Form the A stiffness tensor"""
    A = np.zeros([3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    KLMN = T2V([k,l,m,n],[3,3,3,3])
                    A[KLMN] = LAMBDA*I[k,l]*I[m,n] + MU*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return A
    
def form_B(ETA,TAU,KAPPA,NU,SIGMA): #Test function written
    """Form the B stiffness tensor"""
    B = np.zeros([3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    KLMN = T2V([k,l,m,n],[3,3,3,3])
                    B[KLMN] = (ETA-TAU)*I[k,l]*I[m,n]+KAPPA*I[k,m]*I[l,n]+NU*I[k,n]*I[l,m]\
                              -SIGMA*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return B
                                 
def form_C(TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11): #Test function written
    """Form the C stiffness tensor"""
    C = np.zeros([3*3*3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    for p in range(3):
                        for q in range(3):
                            KLMNPQ = T2V([k,l,m,n,p,q],[3,3,3,3,3,3])
                            C[KLMNPQ] = TAU1*(I[k,l]*I[m,n]*I[p,q]+I[k,q]*I[l,m]*I[n,p])\
                                       +TAU2*(I[k,l]*I[m,p]*I[n,q]+I[k,m]*I[l,q]*I[n,p])\
                                       +TAU3*I[k,l]*I[m,q]*I[n,p]+TAU4*I[k,n]*I[l,m]*I[p,q]\
                                       +TAU5*(I[k,m]*I[l,n]*I[p,q]+I[k,p]*I[l,m]*I[n,q])\
                                       +TAU6*I[k,m]*I[l,p]*I[n,q]+TAU7*I[k,n]*I[l,p]*I[m,q]\
                                       +TAU8*(I[k,p]*I[l,q]*I[m,n]+I[k,q]*I[l,n]*I[m,p])+TAU9*I[k,n]*I[l,q]*I[m,p]\
                                       +TAU10*I[k,p]*I[l,n]*I[m,q]+TAU11*I[k,q]*I[l,p]*I[m,n]
    return C

def form_D(TAU,SIGMA): #Test function written
    """Form the D stiffness Matrix"""
    D = np.zeros([3*3*3*3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    KLMN = T2V([k,l,m,n],[3,3,3,3])
                    D[KLMN] = TAU*I[k,l]*I[m,n]+SIGMA*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return D
    
###### Compute Stresses ######

def compute_stresses(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma): #Test function written
    """Compute all of the stress values"""
    CAUCHY,TERM1,TERM2  = compute_pk2_stress(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    SIGMA               = compute_symmetric_stress(TERM1,TERM2)
    M                   = compute_ho_stress(CFOT,Gamma)
    return CAUCHY,SIGMA,M
    
def compute_pk2_stress(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma): #Test function written
    """Compute the Second Piola Kirchhoff Stress returns PK2 and other useful terms"""
    PK2    = np.zeros([9,])
    I      = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
    Cinv   = hex8.invert_3x3_matrix_V(C)[1]
    STRAIN_TERM = E_micro+I#hex8.matrix_dot_V(Cinv,E_micro)+I
    
    TERM1  = np.zeros([9,])
    TERM2  = np.zeros([9,])
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            for k in range(3):
                for l in range(3):
                    IJKL = T2V([i,j,k,l],[3,3,3,3])
                    KL   = T2V([k,l],[3,3])
                    TERM1[IJ] += AFOT[IJKL]*E_macro[KL] + DFOT[IJKL]*E_micro[KL]
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            for k in range(3):
                for l in range(3):
                    KL   = T2V([k,l],[3,3])
                    for q in range(3):
                        for r in range(3):
                            IQKL = T2V([i,q,k,l],[3,3,3,3])
                            KL   = T2V([k,l],[3,3])
                            RQ   = T2V([r,q],[3,3])
                            JR   = T2V([j,r],[3,3])
                            TERM2[IJ] += (BFOT[IQKL]*E_micro[KL]+DFOT[IQKL]*E_macro[KL])*STRAIN_TERM[RQ]*Cinv[JR]
    
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            for l in range(3):
                for m in range(3):
                    for n in range(3):
                        LMN = T2V([l,m,n],[3,3,3])
                        for q in range(3):
                            for r in range(3):
                                for s in range(3):
                                    IQRLMN = T2V([i,q,r,l,m,n],[3,3,3,3,3,3])
                                    JS     = T2V([j,s],[3,3])
                                    SQR    = T2V([s,q,r],[3,3,3])
                                    TERM2[IJ] += CFOT[IQRLMN]*Gamma[LMN]*Cinv[JS]*Gamma[SQR]
    
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            PK2[IJ] = TERM1[IJ]+TERM2[IJ]
    return PK2,TERM1,TERM2
    
def compute_symmetric_stress(TERM1,TERM2): #Test function written
    """Compute the symmetric stress from the given terms"""
    SIGMA = np.zeros([3*3,])
    TERM2_SYMM  = hex8.get_symm_matrix_V(TERM2)
    for i in range(3):
        for j in range(3):
            IJ = T2V([i,j],[3,3])
            SIGMA[IJ] = TERM1[IJ]+2*TERM2_SYMM[IJ]
    return SIGMA
    
def compute_ho_stress(CFOT,Gamma): #Test function written
    """Get the higher order stress"""
    M = np.zeros([27,])
    for k in range(3):
        for l in range(3):
            for m in range(3):
                KLM = T2V([k,l,m],[3,3,3])
                for n in range(3):
                    for p in range(3):
                        for q in range(3):
                            KLMNPQ = T2V([k,l,m,n,p,q],[3,3,3,3,3,3])
                            NPQ    = T2V([n,p,q],[3,3,3])
                            M[KLM] += CFOT[KLMNPQ]*Gamma[NPQ]
    return M
    
###### Compute the tangents ######

def compute_tangents(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma): #Test function written
    """Compute the tangents for a micromorphic linear elastic
    constitutive model"""
    
    #Common Terms
    I    = np.array([1,0,0,0,1,0,0,0,1]).astype(float)#np.eye(3)
    Cinv = hex8.invert_3x3_matrix_V(C)[1]
    
    #Compute dCinvdC
    dCinvdC = compute_dCinvdC(Cinv)
    
    #Compute tangents
    dSdC,    dSigmadC,    dMdC     = compute_stress_derivatives_wrt_C(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,Gamma,I,dCinvdC)
    dSdPsi,  dSigmadPsi,  dMdPsi   = compute_stress_derivatives_wrt_Psi(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I)
    dSdGamma,dSigmadGamma,dMdGamma = compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,Gamma)
    
    return dSdC,dSdPsi,dSdGamma,dSigmadC,dSigmadPsi,dSigmadGamma,dMdC,dMdPsi,dMdGamma

def compute_stress_derivatives_wrt_C(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,Gamma,Iten,dCinvdC): #Test function written
    """Compute the stress derivatives w.r.t. C"""
    #Initialize tangents with respect to C
    dSdC     = np.zeros([3*3*3*3])
    dSigmadC = np.zeros([3*3*3*3])
    dMdC     = np.zeros([3*3*3*3*3])
    
    #Useful terms for computation to prevent repetition
    TERM1 = np.zeros([3*3*3*3])
    TERM2 = np.zeros([3*3*3*3])
    TERM3 = np.zeros([3*3*3*3])
    
    #Compute dCinvdC
    dCinvdC = compute_dCinvdC(Cinv)
    
    #Compute dSdC, dSigmadC, dMdC
    for I in range(3):
        for J in range(3):
            for O in range(3):
                for P in range(3):
                    IJOP = T2V([I,J,O,P],[3,3,3,3])
                    JIOP = T2V([J,I,O,P],[3,3,3,3])
                    TERM1[IJOP] = 0.5*AFOT[IJOP]
                    
    # for I in range(3):
        # for J in range(3):
            # for O in range(3):
                # for P in range(3):
                    # IJOP = T2V([I,J,O,P],[3,3,3,3])
                    # JIOP = T2V([J,I,O,P],[3,3,3,3])
                    for Q in range(3):
                        IQOP = T2V([I,Q,O,P],[3,3,3,3])
                        JQOP = T2V([J,Q,O,P],[3,3,3,3])
                        for R in range(3):
                            JROP = T2V([J,R,O,P],[3,3,3,3])
                            RQ   = T2V([R,Q],[3,3])
                            JR   = T2V([J,R],[3,3])
                            IR   = T2V([I,R],[3,3])
                            TERM2[IJOP] += 0.5*DFOT[IQOP]*(E_micro[RQ]+Iten[RQ])*Cinv[JR]
                            TERM3[IJOP] += 0.5*DFOT[JQOP]*(E_micro[RQ]+Iten[RQ])*Cinv[IR]
                    for Q in range(3):
                        for R in range(3):
                            for K in range(3):
                                for L in range(3):
                                    RQ = T2V([R,Q],[3,3])
                                    IROP = T2V([I,R,O,P],[3,3,3,3])
                                    JROP = T2V([J,R,O,P],[3,3,3,3])
                                    IQKL = T2V([I,Q,K,L],[3,3,3,3])
                                    JQKL = T2V([J,Q,K,L],[3,3,3,3])
                                    KL = T2V([K,L],[3,3])
                                    TERM2[IJOP] += (BFOT[IQKL]*E_micro[KL]+DFOT[IQKL]*E_macro[KL])*(E_micro[RQ]+Iten[RQ])*dCinvdC[JROP]
                                    TERM3[IJOP] += (BFOT[JQKL]*E_micro[KL]+DFOT[JQKL]*E_macro[KL])*(E_micro[RQ]+Iten[RQ])*dCinvdC[IROP]
                                    
                            for L in range(3):
                                for M in range(3):
                                    for N in range(3):
                                        IQRLMN = T2V([I,Q,R,L,M,N],[3,3,3,3,3,3])
                                        JQRLMN = T2V([J,Q,R,L,M,N],[3,3,3,3,3,3])
                                        LMN = T2V([L,M,N],[3,3,3])
                                        for S in range(3):
                                            SQR  = T2V([S,Q,R],[3,3,3])
                                            ISOP = T2V([I,S,O,P],[3,3,3,3])
                                            JSOP = T2V([J,S,O,P],[3,3,3,3])
                                            TERM2[IJOP] += CFOT[IQRLMN]*Gamma[LMN]*Gamma[SQR]*dCinvdC[JSOP]
                                            TERM3[IJOP] += CFOT[JQRLMN]*Gamma[LMN]*Gamma[SQR]*dCinvdC[ISOP]
                    dSdC[IJOP]     = TERM1[IJOP] + TERM2[IJOP]
                    dSigmadC[IJOP] = TERM1[IJOP] + (TERM2[IJOP]+TERM3[IJOP])#(TERM2[IJOP]+TERM2[JIOP])
    #dSdC = hex8.reduce_tensor_to_vector_form(dSdC)
    return dSdC,dSigmadC,dMdC
    
def compute_stress_derivatives_wrt_Psi(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,Iten): #Test function written
    """Compute the stress derivatives w.r.t. Psi"""
    
    #Initialize tangents with respect to Psi
    dSdPsi     = np.zeros([3*3*3*3])
    dSigmadPsi = np.zeros([3*3*3*3])
    dMdPsi     = np.zeros([3*3*3*3*3])
    
    TERM1 = np.zeros([3*3*3*3])
    TERM2 = np.zeros([3*3*3*3])
    TERM3 = np.zeros([3*3*3*3])
    
    #Compute dSdPsi, dSigmadPsi, and dMdPsi
    for I in range(3):
        for J in range(3):
            for O in range(3):
                for P in range(3):
                    IJOP = T2V([I,J,O,P],[3,3,3,3])
                    dSdPsi[IJOP]     = DFOT[IJOP]
                    dSigmadPsi[IJOP] = DFOT[IJOP]
                    
                    for Q in range(3):
                        for R in range(3):
                            IQOP = T2V([I,Q,O,P],[3,3,3,3])
                            JQOP = T2V([J,Q,O,P],[3,3,3,3])
                            RQ   = T2V([R,Q],[3,3])
                            JR   = T2V([J,R],[3,3])
                            IR   = T2V([I,R],[3,3])
                            TERM2[IJOP] += BFOT[IQOP]*(E_micro[RQ]+Iten[RQ])*Cinv[JR]
                            TERM3[IJOP] += BFOT[JQOP]*(E_micro[RQ]+Iten[RQ])*Cinv[IR]
                            
                    for K in range(3):
                        for L in range(3):
                            IPKL = T2V([I,P,K,L],[3,3,3,3])
                            JPKL = T2V([J,P,K,L],[3,3,3,3])
                            KL   = T2V([K,L],[3,3])
                            IO   = T2V([I,O],[3,3])
                            JO   = T2V([J,O],[3,3])
                                    
                            TERM2[IJOP] += (BFOT[IPKL]*E_micro[KL]+DFOT[IPKL]*E_macro[KL])*Cinv[JO]
                            TERM3[IJOP] += (BFOT[JPKL]*E_micro[KL]+DFOT[JPKL]*E_macro[KL])*Cinv[IO]
                                    
                    dSdPsi[IJOP]     += TERM2[IJOP]
                    dSigmadPsi[IJOP] += TERM2[IJOP] + TERM3[IJOP]
                    
    return dSdPsi,dSigmadPsi,dMdPsi
    
def compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,Gamma): #Test function written
    """Compute the stress derivatives w.r.t. Gamma"""
    #Initialize tangents with respect to Gamma
    dSdGamma     = np.zeros([3*3*3*3*3])
    dSigmadGamma = np.zeros([3*3*3*3*3])
    
    #Compute dSdGamma, dSigmadGamma, and dMdGamma
    TERM = np.zeros([3*3*3*3*3])
    for I in range(3):
        for J in range(3):
            for T in range(3):
                IT = T2V([I,T],[3,3])
                JT = T2V([J,T],[3,3])
                for U in range(3):
                    for V in range(3):
                        IJTUV = T2V([I,J,T,U,V],[3,3,3,3,3])
                        for S in range(3):
                            IS = T2V([I,S],[3,3])
                            JS = T2V([J,S],[3,3])
                            for Q in range(3):
                                for R in range(3):
                                    IQRTUV = T2V([I,Q,R,T,U,V],[3,3,3,3,3,3])
                                    JQRTUV = T2V([J,Q,R,T,U,V],[3,3,3,3,3,3])
                                    
                                    IUVSQR = T2V([I,U,V,S,Q,R],[3,3,3,3,3,3])
                                    JUVSQR = T2V([J,U,V,S,Q,R],[3,3,3,3,3,3])
                                    SQR = T2V([S,Q,R],[3,3,3])
                                    
                                    dSdGamma[IJTUV] += CFOT[IQRTUV]*Cinv[JS]*Gamma[SQR] + CFOT[IUVSQR]*Gamma[SQR]*Cinv[JT]
                                    dSigmadGamma[IJTUV] += CFOT[IQRTUV]*Cinv[JS]*Gamma[SQR] + CFOT[IUVSQR]*Gamma[SQR]*Cinv[JT]+CFOT[JQRTUV]*Cinv[IS]*Gamma[SQR] + CFOT[JUVSQR]*Gamma[SQR]*Cinv[IT]
    
    return dSdGamma,dSigmadGamma,CFOT
    
###### Call the model ######

def micromorphic_linear_elasticity(F,chi,grad_chi,params): #Test function written
    """A constitutive model for micromorphic linear elasticity"""
    C,Psi,Gamma         = get_deformation_measures(F,chi,grad_chi)
    E_macro,E_micro     = compute_strain_measures(C,Psi)
    AFOT,BFOT,CFOT,DFOT = form_stiffness_tensors(params)
    PK2,SIGMA,M         = compute_stresses(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    dSdC,dSdPsi,dSdGamma,dSigmadC,dSigmadPsi,dSigmadGamma,dMdC,dMdPsi,dMdGamma = compute_tangents(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    return PK2,SIGMA,M,dSdC,dSdPsi,dSdGamma,dSigmadC,dSigmadPsi,dSigmadGamma,dMdC,dMdPsi,dMdGamma
    
class TestMicro_LE(unittest.TestCase):

    f                    = None
    original_directory   = ""
    module_name           = "micromorphic_linear_elasticity"
    output_file_name      = r"results.txt".format(module_name)
    output_file_location  = r".\tests\unittests\{0}".format(module_name)
    currentResult         = None
    @classmethod
    def setUpClass(self):
        """Setup method"""
        output_file = os.path.join(self.output_file_location,self.output_file_name)
        
        if(not os.path.isdir(self.output_file_location)):
            os.makedirs(self.output_file_location)
        
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
        if(str(ok)):
            str_out = r"\cellcolor{green!25} PASS"
        else:
            str_out = r"\cellcolor{red!25} FAIL"
        
        self.f.write(tname+"\t&\t"+str_out+r"\\"+"\n")
        
    def run(self, result=None):
        """Redefine run to keep track of results"""
        self.currentResult = result
        unittest.TestCase.run(self,result)
        
    def test_form_A(self):
        """Test forming the A stiffness tensor"""
        LAMBDA = 2.4
        MU     = 6.7
        
        A = form_A(LAMBDA,MU)
        
        I = np.eye(3)
        At = np.zeros([3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        At[K,L,M,N] = LAMBDA*I[K,L]*I[M,N] + MU*(I[K,M]*I[L,N] + I[K,N]*I[L,M])
        At = hex8.reduce_tensor_to_vector_form(At)
        
        self.assertEqual(np.allclose(A,At),True)
        
    def test_form_B(self):
        """Test forming the B stiffness tensor"""
        ETA,TAU,KAPPA,NU,SIGMA = 2.4,5.1,5.6,8.2,2.
        
        B  = form_B(ETA,TAU,KAPPA,NU,SIGMA)
        
        I  = np.eye(3)
        Bt = np.zeros([3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        Bt[K,L,M,N] = (ETA-TAU)*I[K,L]*I[M,N] + KAPPA*I[K,M]*I[L,N] + NU*I[K,N]*I[L,M]\
                                      -SIGMA*(I[K,M]*I[L,N]+I[K,N]*I[L,M])
        Bt = hex8.reduce_tensor_to_vector_form(Bt)                        
        self.assertEqual(np.allclose(B,Bt),True)
        
    def test_form_C(self):
        """Test forming the C stiffness tensor"""
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11 =\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
         
        C = form_C(TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11)
        I = np.eye(3)
        
        Ct = np.zeros([3,3,3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        for P in range(3):
                            for Q in range(3):
                                Ct[K,L,M,N,P,Q] = TAU1*(I[K,L]*I[M,N]*I[P,Q] + I[K,Q]*I[L,M]*I[N,P])\
                                                 +TAU2*(I[K,L]*I[M,P]*I[N,Q] + I[K,M]*I[L,Q]*I[N,P])\
                                                 +TAU3*I[K,L]*I[M,Q]*I[N,P]  + TAU4*I[K,N]*I[L,M]*I[P,Q]\
                                                 +TAU5*(I[K,M]*I[L,N]*I[P,Q] + I[K,P]*I[L,M]*I[N,Q])\
                                                 +TAU6*I[K,M]*I[L,P]*I[N,Q]  + TAU7*I[K,N]*I[L,P]*I[M,Q]\
                                                 +TAU8*(I[K,P]*I[L,Q]*I[M,N] + I[K,Q]*I[L,N]*I[M,P])\
                                                 +TAU9*I[K,N]*I[L,Q]*I[M,P]  + TAU10*I[K,P]*I[L,N]*I[M,Q]\
                                                 +TAU11*I[K,Q]*I[L,P]*I[M,N]
        Ct = hex8.reduce_tensor_to_vector_form(Ct)
        
        self.assertEqual(np.allclose(C,Ct),True)

    def test_form_D(self):
        """Test forming the D stiffness tensor"""
        TAU   = 5.1
        SIGMA = 2.
        
        I = np.eye(3)
        
        D = form_D(TAU,SIGMA)
        
        Dt = np.zeros([3,3,3,3])
        
        for K in range(3):
            for L in range(3):
                for M in range(3):
                    for N in range(3):
                        Dt[K,L,M,N] = TAU*I[K,L]*I[M,N]+SIGMA*(I[K,M]*I[L,N]+I[K,N]*I[L,M])
        Dt = hex8.reduce_tensor_to_vector_form(Dt)
        
        self.assertEqual(np.allclose(D,Dt),True)
        
    def test_form_stiffness_tensors(self):
        """Test forming the stiffness tensors"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        A,B,C,D = form_stiffness_tensors(PROPS)
        
        At = form_A(LAMBDA,MU)
        Bt = form_B(ETA,TAU,KAPPA,NU,SIGMA)
        Ct = form_C(TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11)
        Dt = form_D(TAU,SIGMA)
        
        self.assertTrue(np.allclose(A,At))
        self.assertTrue(np.allclose(B,Bt))
        self.assertTrue(np.allclose(C,Ct))
        self.assertTrue(np.allclose(D,Dt))
        
    def test_compute_deformation_measures(self):
        """Test to compute the deformation measures"""
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
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
        
    def test_compute_strain_measures(self):
        """Test the computation of the strain measures"""
        F   = np.array(range(9))
        chi = np.array(range(9,18))
        
        C   = hex8.matrix_Tdot_V(F,F)
        Psi = hex8.matrix_Tdot_V(F,chi)
        I   = hex8.reduce_tensor_to_vector_form(np.eye(3))
        
        E_macro,E_micro = compute_strain_measures(C,Psi)
        
        E_macrot = 0.5*(C-I)
        E_microt = Psi-I
        
        self.assertTrue(np.allclose(E_macro,E_macrot))
        self.assertTrue(np.allclose(E_micro,E_microt))
        
    def test_compute_dCinvdC(self):
        """Test computinig the derivative of Cinv w.r.t. C"""
        Cinv  = np.array(range(9))
        Cinvt = hex8.convert_V_to_T(Cinv,[3,3])
        
        dCinvdC = compute_dCinvdC(Cinv)
        
        dCinvdCt = np.zeros([3,3,3,3])
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    for L in range(3):
                        dCinvdCt[I,J,K,L] = -Cinvt[I,K]*Cinvt[L,J]
        dCinvdCt = hex8.reduce_tensor_to_vector_form(dCinvdCt)
        
        self.assertTrue(np.allclose(dCinvdC,dCinvdCt))
        
    def test_compute_pk2_stress(self):
        """Test the computation of the pk2 stress"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        
        E_micro,E_macro = compute_strain_measures(C,Psi)
        
        pk2   = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)[0]
        
        St  = np.zeros([9,])
        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for K in range(3):
                    for L in range(3):
                        IJKL = T2V([I,J,K,L],[3,3,3,3])
                        KL   = T2V([K,L],        [3,3])
                        St[IJ] +=  AS[IJKL]*E_macro[KL] + DS[IJKL]*E_micro[KL]
                        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for K in range(3):
                    for L in range(3):
                        KL   = T2V([K,L],[3,3])
                        for Q in range(3):
                            for R in range(3):
                                IQKL = T2V([I,Q,K,L],[3,3,3,3])
                                RQ   = T2V([R,Q],[3,3])
                                JR   = T2V([J,R],[3,3])
                                St[IJ] += (BS[IQKL]*E_micro[KL] + DS[IQKL]*E_macro[KL])*(E_micro[RQ] + Iten[RQ])*Cinv[JR]
        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for L in range(3):
                    for M in range(3):
                        for N in range(3):
                            LMN = T2V([L,M,N],[3,3,3])
                            for Q in range(3):
                                for R in range(3):
                                    for S in range(3):
                                        IQRLMN = T2V([I,Q,R,L,M,N],[3,3,3,3,3,3])
                                        JS     = T2V([J,S],[3,3])
                                        SQR    = T2V([S,Q,R],[3,3,3])
                                        St[IJ] += CS[IQRLMN]*Gamma[LMN]*Cinv[JS]*Gamma[SQR]
        
        self.assertTrue(np.allclose(pk2,St))

    def test_compute_symmetric_stress(self):
        """Test the computation of the symmetric stress"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        
        E_micro,E_macro = compute_strain_measures(C,Psi)
        
        _,TEMP1,TEMP2   = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
        symm_stress     = compute_symmetric_stress(TEMP1,TEMP2)
        
        symm_stresst = np.zeros([9,])
        
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for K in range(3):
                    for L in range(3):
                        IJKL = T2V([I,J,K,L],[3,3,3,3])
                        KL   = T2V([K,L],        [3,3])
                        symm_stresst[IJ] += AS[IJKL]*E_macro[KL]+DS[IJKL]*E_micro[KL]
                        for Q in range(3):
                            for R in range(3):
                                IQKL = T2V([I,Q,K,L],[3,3,3,3])
                                JQKL = T2V([J,Q,K,L],[3,3,3,3])
                                RQ = T2V([R,Q],[3,3])
                                JR = T2V([J,R],[3,3])
                                IR = T2V([I,R],[3,3])
                                tmp1 = (BS[IQKL]*E_micro[KL] + DS[IQKL]*E_macro[KL])*(E_micro[RQ]+Iten[RQ])*Cinv[JR]
                                tmp2 = (BS[JQKL]*E_micro[KL] + DS[JQKL]*E_macro[KL])*(E_micro[RQ]+Iten[RQ])*Cinv[IR]
                                symm_stresst[IJ] += tmp1+tmp2
        for I in range(3):
            for J in range(3):
                IJ = T2V([I,J],[3,3])
                for L in range(3):
                    for M in range(3):
                        for N in range(3):
                            LMN = T2V([L,M,N],[3,3,3])
                            for Q in range(3):
                                for R in range(3):
                                    for S in range(3):
                                        IQRLMN = T2V([I,Q,R,L,M,N],[3,3,3,3,3,3])
                                        JQRLMN = T2V([J,Q,R,L,M,N],[3,3,3,3,3,3])
                                        SQR    = T2V([S,Q,R],[3,3,3])
                                        JS     = T2V([J,S],[3,3])
                                        IS     = T2V([I,S],[3,3])
                                        tmp1   = CS[IQRLMN]*Gamma[LMN]*Gamma[SQR]*Cinv[JS]
                                        tmp2   = CS[JQRLMN]*Gamma[LMN]*Gamma[SQR]*Cinv[IS]
                                        symm_stresst[IJ] += tmp1 + tmp2
        self.assertTrue(np.allclose(symm_stress,symm_stresst))
    
    def test_compute_ho_stress(self):
        """Test the computation of the higher order stress"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        Mten  = compute_ho_stress(CS,Gamma)
        Mt = np.zeros([27,])
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for L in range(3):
                        for M in range(3):
                            for N in range(3):
                                IJKLMN = T2V([I,J,K,L,M,N],[3,3,3,3,3,3])
                                LMN    = T2V([L,M,N],[3,3,3])
                                Mt[IJK] += CS[IJKLMN]*Gamma[LMN]
        self.assertTrue(np.allclose(Mten,Mt))

    def test_compute_stresses(self):
        """Test the computation of all of the stress measures"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        E_micro,E_macro = compute_strain_measures(C,Psi)
        
        PK2T,TERM1,TERM2 = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
        SIGMAT           = compute_symmetric_stress(TERM1,TERM2)
        MT               = compute_ho_stress(CS,Gamma)
        PK2,SIGMA,M      = compute_stresses(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
        
        self.assertTrue(np.allclose(PK2T,    PK2))
        self.assertTrue(np.allclose(SIGMAT,SIGMA))
        self.assertTrue(np.allclose(MT,        M))
        
    def test_compute_stress_derivatives_wrt_C(self):
        """Test the computation of the stress derivatives wrt C"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        Iten     = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        #Ften = hex8.convert_V_to_T(F,[3,3])
        #Cten = np.einsum('iI,iJ',Ften,Ften)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        E_macro,E_micro = compute_strain_measures(C,Psi)
        
        def parse_pk2(Ci):
            E_macroi,E_microi = compute_strain_measures(Ci,Psi)
            pk2   = compute_pk2_stress(AS,BS,CS,DS,E_macroi,E_microi,Ci,Gamma)[0]
            return pk2
            
        def parse_symmetric_stress(Ci):
            E_macroi,E_microi = compute_strain_measures(Ci,Psi)
            _,TERM1,TERM2     = compute_pk2_stress(AS,BS,CS,DS,E_macroi,E_microi,Ci,Gamma)
            return compute_symmetric_stress(TERM1,TERM2)
        
        #print "\nNumeric gradient"
        #Test the gradient of PK2
        dpk2dCT = fd.numeric_gradient(parse_pk2,C,1e-6)
        dpk2dCT = self._convert_difference(dpk2dCT)
        dpk2dCT = hex8.reduce_tensor_to_vector_form(dpk2dCT)
        #print hex8.convert_V_to_M(dpk2dCT,[3,3,3,3])       
        
        dCinvdC = compute_dCinvdC(Cinv)
        
        dpk2dC,dSigmadC,dMdC  = compute_stress_derivatives_wrt_C(E_macro,E_micro,AS,BS,CS,DS,Cinv,Gamma,Iten,dCinvdC)
        #print hex8.convert_V_to_M(dpk2dC,[3,3,3,3])
        
        #print "\nSecond Piola Kirchhoff"
        #print "Numeric:\n",dpk2dCT
        #print "Analytic:\n",dpk2dC
        #print "Difference:\n",dpk2dCT-dpk2dC
        
        self.assertTrue(np.allclose(dpk2dCT,dpk2dC))
        
        #Test the gradient of the symmetric stress
        dSigmadCT = fd.numeric_gradient(parse_symmetric_stress,C,1e-6)
        dSigmadCT = self._convert_difference(dSigmadCT)
        dSigmadCT = hex8.reduce_tensor_to_vector_form(dSigmadCT)
        
        # print "\nSymmetric stress"
        # print "Numeric:\n",dSigmadCT
        # print "Analytic:\n",dSigmadC
        # print "Difference:\n",dSigmadCT-dSigmadC
        
        self.assertTrue(np.allclose(dSigmadCT,dSigmadC))
        
        #Test the gradient of the higher order stress
        self.assertTrue(np.allclose(dMdC,np.zeros([3*3*3*3*3])))
    
    def test_compute_stress_derivatives_wrt_Psi(self):
        """Test the computation of the stress derivatives wrt Psi
        Note: finite difference values delta values reduced to allow 
        convergence
        """
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        Iten     = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        #Ften = hex8.convert_V_to_T(F,[3,3])
        #Cten = np.einsum('iI,iJ',Ften,Ften)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        E_macro,E_micro = compute_strain_measures(C,Psi)
        
        def parse_pk2(Psii):
            E_macroi,E_microi = compute_strain_measures(C,Psii)
            pk2   = compute_pk2_stress(AS,BS,CS,DS,E_macroi,E_microi,C,Gamma)[0]
            return pk2
            
        def parse_symmetric_stress(Psii):
            E_macroi,E_microi = compute_strain_measures(C,Psii)
            _,TERM1,TERM2     = compute_pk2_stress(AS,BS,CS,DS,E_macroi,E_microi,C,Gamma)
            return compute_symmetric_stress(TERM1,TERM2)
        
        #print "\nNumeric gradient"
        #Test the gradient of PK2
        dpk2dPsiT = fd.numeric_gradient(parse_pk2,Psi,1e-4)
        dpk2dPsiT = self._convert_difference(dpk2dPsiT)
        dpk2dPsiT = hex8.reduce_tensor_to_vector_form(dpk2dPsiT)
        #print hex8.convert_V_to_M(dpk2dCT,[3,3,3,3])
        
        dpk2dPsi,dSigmadPsi,dMdPsi  = compute_stress_derivatives_wrt_Psi(E_macro,E_micro,AS,BS,CS,DS,Cinv,Iten)
        #print hex8.convert_V_to_M(dpk2dPsi,[3,3,3,3])
        
        #print "\nSecond Piola Kirchhoff"
        #print "Numeric:\n",dpk2dPsiT
        #print "Analytic:\n",dpk2dPsi
        #print "Difference:\n",dpk2dPsiT-dpk2dPsi
        #print "Max difference: {0}".format(max(dpk2dPsiT-dpk2dPsi,key=abs))
        
        self.assertTrue(np.allclose(dpk2dPsiT,dpk2dPsi))
        
        #Test the gradient of the symmetric stress
        dSigmadPsiT = fd.numeric_gradient(parse_symmetric_stress,Psi,1e-4)
        dSigmadPsiT = self._convert_difference(dSigmadPsiT)
        dSigmadPsiT = hex8.reduce_tensor_to_vector_form(dSigmadPsiT)
        
        # print "\nSymmetric stress"
        # print "Numeric:\n",dSigmadPsiT
        # print "Analytic:\n",dSigmadPsi
        # print "Difference:\n",dSigmadPsiT-dSigmadPsi
        
        self.assertTrue(np.allclose(dSigmadPsiT,dSigmadPsi))
        
        #Test the gradient of the higher order stress
        self.assertTrue(np.allclose(dMdPsi,np.zeros([3*3*3*3*3])))
        
    def test_compute_stress_derivatives_wrt_Gamma(self):
        """Test the computation of the stress derivatives wrt Gamma"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        Iten     = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        #Ften = hex8.convert_V_to_T(F,[3,3])
        #Cten = np.einsum('iI,iJ',Ften,Ften)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        E_macro,E_micro = compute_strain_measures(C,Psi)
        
        def parse_pk2(Gammai):
            E_macro,E_micro = compute_strain_measures(C,Psi)
            pk2   = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gammai)[0]
            return pk2
            
        def parse_symmetric_stress(Gammai):
            E_macro,E_micro = compute_strain_measures(C,Psi)
            _,TERM1,TERM2     = compute_pk2_stress(AS,BS,CS,DS,E_macro,E_micro,C,Gammai)
            return compute_symmetric_stress(TERM1,TERM2)
        
        #print "\nNumeric gradient"
        #Test the gradient of PK2
        dpk2dGammaT = fd.numeric_gradient(parse_pk2,Gamma,1e-4)
        dpk2dGammaT = self._convert_difference_TOT(dpk2dGammaT)
        dpk2dGammaT = hex8.reduce_tensor_to_vector_form(dpk2dGammaT)
        #print hex8.convert_V_to_M(dpk2dCT,[3,3,3,3])
        
        dpk2dGamma,dSigmadGamma,dMdGamma  = compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AS,BS,CS,DS,Cinv,Gamma)
        #print hex8.convert_V_to_M(dpk2dPsi,[3,3,3,3])
        
        #print "\nSecond Piola Kirchhoff"
        #print "Numeric:\n",dpk2dGammaT
        #print "Analytic:\n",dpk2dGamma
        #print "Difference:\n",dpk2dGammaT-dpk2dGamma
        #print "Max difference: {0}".format(max(dpk2dGammaT-dpk2dGamma,key=abs))
        
        self.assertTrue(np.allclose(dpk2dGammaT,dpk2dGamma))
        
        #Test the gradient of the symmetric stress
        dSigmadGammaT = fd.numeric_gradient(parse_symmetric_stress,Gamma,1e-4)
        dSigmadGammaT = self._convert_difference_TOT(dSigmadGammaT)
        dSigmadGammaT = hex8.reduce_tensor_to_vector_form(dSigmadGammaT)
        
        #print "\nSymmetric stress"
        #print "Numeric:\n",dSigmadGammaT
        #print "Analytic:\n",dSigmadGamma
        #print "Difference:\n",dSigmadGammaT-dSigmadGamma
        #print "Max difference: {0}".format(max(dSigmadGammaT-dSigmadGamma,key=abs))
        
        self.assertTrue(np.allclose(dSigmadGammaT,dSigmadGamma))
        
        #Test the gradient of the higher order stress
        self.assertTrue(np.allclose(dMdGamma,CS))
    
    def test_compute_stress_derivatives(self):
        """Test the computation of the stress derivatives"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        Iten     = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        #Ften = hex8.convert_V_to_T(F,[3,3])
        #Cten = np.einsum('iI,iJ',Ften,Ften)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        E_macro,E_micro = compute_strain_measures(C,Psi)
        
        dCinvdC = compute_dCinvdC(Cinv)
        
        dpk2dCT,dSigmadCT,dMdCT              = compute_stress_derivatives_wrt_C(E_macro,E_micro,AS,BS,CS,DS,Cinv,Gamma,Iten,dCinvdC)
        dpk2dPsiT,dSigmadPsiT,dMdPsiT        = compute_stress_derivatives_wrt_Psi(E_macro,E_micro,AS,BS,CS,DS,Cinv,Iten)
        dpk2dGammaT,dSigmadGammaT,dMdGammaT  = compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AS,BS,CS,DS,Cinv,Gamma)
        
        dpk2dC,dpk2dPsi,dpk2dGamma,dSigmadC,dSigmadPsi,dSigmadGamma,dMdC,dMdPsi,dMdGamma = compute_tangents(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
        
        self.assertTrue(np.allclose(dpk2dCT,     dpk2dC))
        self.assertTrue(np.allclose(dpk2dPsiT,   dpk2dPsi))
        self.assertTrue(np.allclose(dpk2dGammaT, dpk2dGamma))
        
        self.assertTrue(np.allclose(dSigmadCT,     dSigmadC))
        self.assertTrue(np.allclose(dSigmadPsiT,   dSigmadPsi))
        self.assertTrue(np.allclose(dSigmadGammaT, dSigmadGamma))
        
        self.assertTrue(np.allclose(dMdCT,     dMdC))
        self.assertTrue(np.allclose(dMdPsiT,   dMdPsi))
        self.assertTrue(np.allclose(dMdGammaT, dMdGamma))
        
    def test_micromorphic_linear_elasticity(self):
        """Test function for the micromorphic linear elasticity model"""
        RHO0   = 2.9
        LAMBDA = 2.4
        MU     = 6.7
        ETA    = 2.4
        TAU    = 5.1
        KAPPA  = 5.6
        NU     = 8.2
        SIGMA  = 2.
        TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11=\
        [4.5,1.3,9.2,1.1,6.4,2.4,7.11,5.5,1.5,3.8,2.7]
        PROPS = RHO0,LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11
        
        AS,BS,CS,DS = form_stiffness_tensors(PROPS)
        
        F        = np.array([1,4,2,2,2,3,3,2,1]).astype(float)
        chi      = np.array(range(9,18)).astype(float)
        grad_chi = np.array(range(18,45)).astype(float)
        Iten     = np.array([1,0,0,0,1,0,0,0,1]).astype(float)
        
        C   = hex8.matrix_Tdot_V(F,F)
        Cinv = hex8.invert_3x3_matrix_V(C)[1]
        Psi = hex8.matrix_Tdot_V(F,chi)
        
        #Ften = hex8.convert_V_to_T(F,[3,3])
        #Cten = np.einsum('iI,iJ',Ften,Ften)
        
        Gamma = np.zeros([27,])
        
        for I in range(3):
            for J in range(3):
                for K in range(3):
                    IJK = T2V([I,J,K],[3,3,3])
                    for i in range(3):
                        iI  = T2V([i,I],[3,3])
                        iJK = T2V([i,J,K],[3,3,3])
                        Gamma[IJK] += F[iI]*grad_chi[iJK]
        
        Iten = hex8.reduce_tensor_to_vector_form(np.eye(3))
        E_macro,E_micro = compute_strain_measures(C,Psi)
        
        dCinvdC = compute_dCinvdC(Cinv)
        
        PK2,Sigma,M,dpk2dC,dpk2dPsi,dpk2dGamma,dSigmadC,dSigmadPsi,dSigmadGamma,dMdC,dMdPsi,dMdGamma = micromorphic_linear_elasticity(F,chi,grad_chi,PROPS)
        
        PK2T,SigmaT,MT      = compute_stresses(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
    
        dpk2dCT,dpk2dPsiT,dpk2dGammaT,dSigmadCT,dSigmadPsiT,dSigmadGammaT,dMdCT,dMdPsiT,dMdGammaT = compute_tangents(AS,BS,CS,DS,E_macro,E_micro,C,Gamma)
    
        self.assertTrue(np.allclose(PK2,PK2T))
        self.assertTrue(np.allclose(Sigma,SigmaT))
        self.assertTrue(np.allclose(M,MT))
        
        self.assertTrue(np.allclose(dpk2dCT,     dpk2dC))
        self.assertTrue(np.allclose(dpk2dPsiT,   dpk2dPsi))
        self.assertTrue(np.allclose(dpk2dGammaT, dpk2dGamma))
        
        self.assertTrue(np.allclose(dSigmadCT,     dSigmadC))
        self.assertTrue(np.allclose(dSigmadPsiT,   dSigmadPsi))
        self.assertTrue(np.allclose(dSigmadGammaT, dSigmadGamma))
        
        self.assertTrue(np.allclose(dMdCT,     dMdC))
        self.assertTrue(np.allclose(dMdPsiT,   dMdPsi))
        self.assertTrue(np.allclose(dMdGammaT, dMdGamma))
    
    def _convert_difference(self,A):
        """Convert the finite difference representation to tensor form"""
        #print A.shape
        Aten = np.zeros([3,3,3,3])
        for i in range(9):
            it,jt = V2T(i,[3,3])
            #print "Row: ",mapping[j]
            for j in range(9):
                kt,lt = V2T(j,[3,3])
                #print "Col: ",mapping[j]
                Aten[kt,lt,it,jt] = A[i,j]
        return Aten
        
    def _convert_difference_TOT(self,A):
        """Convert the finite difference representation to tensor form"""
        Aten = np.zeros([3,3,3,3,3])
        for i in range(27):
            it,jt,kt = V2T(i,[3,3,3])
            for j in range(9):
                lt,mt = V2T(j,[3,3])
                Aten[lt,mt,it,jt,kt] = A[i,j]
        return Aten
        
        
if __name__ == '__main__':
    unittest.main()