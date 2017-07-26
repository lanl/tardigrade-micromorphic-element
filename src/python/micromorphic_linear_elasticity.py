import numpy as np
import hex8
import finite_difference as fd
import unittest
import os
"""
==================================
| Micromorphic Linear Elasticity |
==================================

A constitutive model for micromorphic
linear elasticity. Written in a style
which will be simple to port to Abaqus
for use in finite element simulations.

Model derived from a quadratic form of
the Helmholtz free energy detailed in
"Nonlinear Micromorphic Continuum 
Mechanics and Finite Strain 
Elastoplasticity" by Dr. Richard Regueiro
with expansions by Nathan Miller
"""

#Compute strain measures

def compute_strain_measures(C,Psi):
    """Compute the strain measures"""
    I = np.eye(3)
    E_macro = 0.5*(C - I)
    E_micro = Psi - I #TODO I think this has to be Psi-C
    return E_macro,E_micro
    
def compute_dCinvdC(C):
    """Compute the derivative of Cinv w.r.t. C"""
    np.zeros([3,3])
    for I in range(3):
        for J in range(3):
            for K in range(3):
                for L in range(3):
                    dCinvdC[I,J,K,L] = -Cinv[I,K]*Cinv[L,J]
    return dCinvdC

###### Form Stiffness Tensors ######

def form_stiffness_tensors(params):
    """Form the stiffness tensors"""
    LAMBDA,MU,ETA,TAU,KAPPA,NU,SIGMA,TAUS = PARAMS
    AFOT = form_A(LAMBDA,MU)
    BFOT = form_B(ETA,TAU,KAPPA,NU,SIGMA)
    CFOT = form_C(TAUS)
    DFOT = form_D(TAU,SIGMA)
    return AFOT,BFOT,CFOT,DFOT
    
def form_A(LAMBDA,MU):
    """Form the A stiffness tensor"""
    A = np.zeros([3,3,3,3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    A[k,l,m,n] = LAMBDA*I[k,l]*I[m,n] + MU*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
    return A
    
def form_B(ETA,TAU,KAPPA,NU,SIGMA):
    """Form the B stiffness tensor"""
    B = np.zeros([3,3,3,3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    B[k,l,m,n] = (ETA-TAU)*I[k,l]*I[m,n]+KAPPA*I[k,m]*I[l,n]+NU*I[k,n]*I[l,m]\
                                 -SIGMA*(I[k,m]*I[l,n]+I[k,n]*I[l,m])
                                 
def form_C(TAUS):
    """Form the C stiffness tensor"""
    TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,TAU9,TAU10,TAU11 = TAUS
    C = np.zeros([3,3,3,3,3,3])
    I = np.eye(3)
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    for p in range(3):
                        for q in range(3):
                            C[k,l,m,n,p,q] = TAU1*(I[k,l]*I[m,n]*I[p,q]+I[k,q]*I[l,m]*I[n,p])\
                                            +TAU2*(I[k,l]*I[m,p]*I[n,q]+I[k,m]*I[l,q]*I[n,p])\
                                            +TAU3*I[k,l]*I[m,q]*I[n,p]+TAU4*I[k,n]*I[l,m]*I[p,q]\
                                            +TAU5*(I[k,m]*I[l,n]*I[p,q]+I[k,p]*I[l,m]*I[n,q])\
                                            +TAU6*I[k,m]*I[l,p]*I[n,q]+TAU7*I[k,n]*I[l,p]*I[m,q]\
                                            +TAU8*(I[k,p]*I[l,q]*I[m,n]+I[k,q]*I[l,n]*I[m,p])+TAU9*I[k,n]*I[l,q]*I[m,p]\
                                            +TAU10*I[k,p]*I[l,n]*I[m,q]+TAU11*I[k,q]*I[l,p]*I[m,n]

def form_D(TAU,SIGMA):
    """Form the D stiffness Matrix"""
    D = np.zeros([3,3])
    I = np.eye(3)
    for k in range(3):
        for j in range(3):
            for m in range(3):
                for n in range(3):
                    D[k,l,m,n] = TAU*I[k,l]*I[m,n]+SIGMA*(I[k,m]*I[l,n]+I[k,n]*I[l,m])       
    
###### Compute Stresses ######
    
def compute_linear_elastic_pk2_stress(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma):
    """Compute the Second Piola Kirchhoff Stress returns PK2 and other useful terms"""
    PK2    = np.zeros([3,3])
    I      = np.eye(3)
    Cinv   = hex8.invert_3x3_matrix(C)
    STRAIN_TERM = hex8.matrix_dot(Cinv,E_micro)+I
    
    TERM1  = np.zeros([3,3])
    TERM2  = np.zeros([3,3])
    
    for k in range(3):
        for l in range(3):
            
            for m in range(3):
                for n in range(3):
                    TERM1[k,l] += A[k,l,m,n]*E_macro[m,n] + B[k,l,m,n]*E_micro[m,n]
    
    for k in range(3):
        for l in range(3):
            
            for m in range(3):
                for n in range(3):
                    
                    for b in range(3):
                        TERM2[k,l] += (D[k,b,m,n]*E_macro[m,n]+B[k,b,m,n]*E_micro[m,n])*STRAIN_TERM[l,b]
    
    for k in range(3):
        for l in range(3):
            for b in range(3):
                for c in range(3):
                    for n in range(3):
                        for p in range(3):
                            for q in range(3):
                                TERM2[k,l] += C[k,b,c,n,p,q]*Gamma[n,p,q]*Cinv[l,q]*Gamma[q,b,c]
    for i in range(3):
        for j in range(3):
            PK2[i,j] = TERM1[i,j]+TERM2[i,j]
    return PK2,TERM1,TERM2
    
def compute_symmetric_stress(TERM1,TERM2):
    """Compute the symmetric stress from the given terms"""
    SIGMA = np.zeros([3,3])
    TERM2_SYMM  = hex8.get_symm_matrix(TERM2)
    for i in range(3):
        for j in range(3):
            SIGMA[i,j] = TERM1[i,j]+2*TERM2_SYMM[i,j]
    return SIGMA
    
def compute_ho_stress(CFOT,Gamma):
    """Get the higher order stress"""
    M = np.zeros([3,3,3])
    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    for p in range(3):
                        for q in range(3):
                            M[k,l,m] = C[k,l,m,n,p,q]*Gamma[n,p,q]
    return M
    
###### Compute the tangents ######

def compute_linear_elastic_tangents(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma):
    """Compute the tangents for a micromorphic linear elastic
    constitutive model"""
    
    #Common Terms
    I    = np.eye(3)
    Cinv = hex8.invert_3x3_matrix(C)
    
    #Compute dCinvdC
    dCinvdC = compute_dCinvdC(Cinv)
    
    #Compute tangents
    dSdC,    dSigmadC,    dMdC     = compute_stress_derivatives_wrt_C(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I,dCinvdC)
    dSdPhi,  dSigmadPhi,  dMdPhi   = compute_stress_derivatives_wrt_Phi(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I)
    dSdGamma,dSigmadGamma,dMdGamma = compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I)
    
    return dSdC,dSdPhi,dSdGamma,dSigmadC,dSigmadPhi,dSigmadGamma,dMdC,dMdPhi,dMdGamma

def compute_stress_derivatives_wrt_C(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I,dCinvdC):
    """Compute the stress derivatives w.r.t. C"""
    
    #Initialize tangents with respect to C
    dSdC     = np.zeros([3,3,3,3])
    dSigmadC = np.zeros([3,3,3,3])
    dMdC     = np.zeros([3,3,3,3,3])
    
    #Useful terms for computation to prevent repetition
    TERM1 = np.zeros([3,3,3,3])
    TERM2 = np.zeros([3,3,3,3])
    
    #Compute dCinvdC
    dCinvdC = compute_dCinvdC(Cinv)
    
    #Compute dSdC, dSigmadC, dMdC
    for I in range(3):
        for J in range(3):
            for O in range(3):
                for P in range(3):
                    TERM1[I,J,O,P] = AFOT[I,J,O,P]
                    
                    for Q in range(3):
                        for R in range(3):
                            TERM1[I,J,O,P] += DFOT[I,Q,O,P]*(E_micro[R,Q]+I[R,Q])*Cinv[J,R]
                    
                            for K in range(3):
                                for L in range(3):
                                    TERM2[I,J,O,P] += (BFOT[I,Q,K,L]*E_micro[K,L]+DFOT[I,Q,K,L]*E_macro[K,L])*(E_micro[R,Q]+I[R,Q])*dCinvdC[J,R,O,P]
                                    
                            for L in range(3):
                                for M in range(3):
                                    for N in range(3):
                                        TERM2[I,J,O,P] += CFOT[I,Q,R,L,M,N]*Gamma[L,M,N]*Gamma[S,Q,R]*dCinvdC[J,S,O,P]
                                        
                    dSdC[I,J,O,P]     = TERM1[I,J,O,P] + TERM2[I,J,O,P]
                    dSigmadC[I,J,O,P] = TERM1[I,J,O,P] + (TERM2[I,J,O,P]+TERM2[J,I,O,P])
    
    return dSdC,dSigmadC,dMdC
    
def compute_stress_derivatives_wrt_Phi(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I):
    """Compute the stress derivatives w.r.t. Phi"""
    
    #Initialize tangents with respect to Phi
    dSdPhi     = np.zeros([3,3,3,3])
    dSigmadPhi = np.zeros([3,3,3,3])
    dMdPhi     = np.zeros([3,3,3,3,3])
    
    TERM1 = np.zeros([4,4,4,4])
    TERM2 = np.zeros([4,4,4,4])
    
    #Compute dSdPhi, dSigmadPhi, and dMdPhi
    for I in range(3):
        for J in range(3):
            for O in range(3):
                for P in range(3):
                    dSdPhi[I,J,O,P]     = DFOT[I,J,O,P]
                    dSigmadPhi[I,J,O,P] = DFOT[I,J,O,P]
                    
                    for Q in range(3):
                        for R in range(3):
                            TERM2[I,J,O,P] += BFOT[I,Q,O,P]*(E_micro[R,Q]+I[R,Q])*Cinv[J,R]
                            
                            for K in range(3):
                                for L in range(3):
                                    TERM2[I,J,O,P] += (BFOT[I,P,K,L]*E_micro[K,L]+DFOT[I,P,K,L]*E_macro[K,L])*Cinv[J,O]
                                    
                    dSdPhi[I,J,O,P]     += TERM2[I,J,O,P]
                    dSigmadPhi[I,J,O,P] += TERM2[I,J,O,P] + TERM2[J,I,O,P]
                    
    return dSdPhi,dSigmadPhi,dMdPhi
    
def compute_stress_derivatives_wrt_Gamma(E_macro,E_micro,AFOT,BFOT,CFOT,DFOT,Cinv,I):
    """Compute the stress derivatives w.r.t. Gamma"""
    #Initialize tangents with respect to Gamma
    dSdGamma     = np.zeros([3,3,3,3,3])
    dSigmadGamma = np.zeros([3,3,3,3,3])
    
    #Compute dSdGamma, dSigmadGamma, and dMdGamma
    TERM = np.zeros([3,3,3,3,3])
    for I in range(3):
        for J in range(3):
            for T in range(3):
                for U in range(3):
                    for V in range(3):
                        for S in range(3):
                            for Q in range(3):
                                for R in range(3):
                                    dSdGamma[I,J,T,U,V] += CFOT[I,Q,R,T,U,]*Cinv[J,S]*Gamma[S,Q,R] + CFOT[I,U,V,S,Q,R]*Gamma[S,Q,R]*Cinv[J,T]
                        dSigmadGamma[I,J,T,U,V] = dSdGamma[I,J,T,U,V] + dSdGamma[J,I,T,U,V]
    
    return dSdGamma,dSigmadGamma,CFOT
    
    
###### Call the model ######

def micromorphic_linear_elasticity(F,chi,grad_chi,params):
    """A constitutive model for micromorphic linear elasticity"""
    C,Phi,Gamma         = get_deformation_measures(F,chi,grad_chi)
    AFOT,BFOT,CFOT,DFOT = form_stiffness_tensors(params)
    CAUCHY,TERM1,TERM2  = compute_linear_elastic_pk2_stress(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    SIGMA               = compute_linear_elastic_symmetric_stress(TERM1,TERM2)
    M                   = compute_linear_elastic_higher_order_stress(CFOT,Gamma)
    dSdC,dSdPhi,dSdGamma,dSigmadC,dSigmadPhi,dSigmadGamma,dMdC,dMdPhi,dMdGamma = compute_linear_elastic_tangents(AFOT,BFOT,CFOT,DFOT,E_macro,E_micro,C,Gamma)
    return CAUCHY,SIGMA,M,dSdC,dSdPhi,dSdGamma,dSigmadC,dSigmadPhi,dSigmadGamma,dMdC,dMdPhi,dMdGamma
    
class TestMicro_LE(unittest.TestCase):

    f                    = None
    original_directory   = ""
    output_file_name     = r"micromorphic_linear_elasticity.txt"
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