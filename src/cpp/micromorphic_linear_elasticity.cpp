/*!=======================================================
  |                                                     |
  |         micromorphic_linear_elasticity.cpp          |
  |                                                     |
  -------------------------------------------------------
  | The header file for the definition of a             |
  | micromorphic linear elasticity.                     |
  -------------------------------------------------------
  | Notes: Micromorphic constitutive models should be   |
  |        developed in the namespace micro_material    |
  |        and have the function get_stress. This       |
  |        function should read in the float parameter  |
  |        vector, the integer parmeter vector, the     |
  |        right Cauchy-Green deformation tensor, Psi,  |
  |        Gamma, and write the PK2 stress, the         |
  |        symmetric stress in the reference            |
  |        configuration, and the higher order couple   |
  |        stress in the reference configuration.       |
  |        (ADDITIONAL VALUES WILL BE ADDED OVER TIME). |                           |
  =======================================================
  | Dependencies:                                       |
  | tensor:      The class which defines tensor access  |
  |              to an underlying Eigen matrix. This    |
  |              may result in a somewhat slower result |
  |              however it should allow for a faster   |
  |              implementation.                        |
  =======================================================*/
  
#include <vector>
#include <ctime>
#include <tensor.h>
#include <micromorphic_linear_elasticity.h>
  
namespace micro_material{
    
    void get_stress(const std::vector< double > &fparams,  const std::vector< int > &iparams,                                         //Material parameters
                    const tensor::Tensor23&   C_in,        const tensor::Tensor23&   Psi_in,   const tensor::Tensor33& Gamma_in,      //Deformation measures
                            tensor::Tensor23& PK2_stress,  tensor::Tensor23& SIGMA_stress,     tensor::Tensor33& M_stress){           //Stress measures
                            
        /*!========================================
        |              get_stress              |
        ========================================
        
        Compute and update the values of the second Piola 
        Kirchhoff, the symmetric, and higher order stresses.
        
        The incoming values are the floating point parameters 
        for the model, the integer parameters for the model, 
        the deformation metrics C, the right Cauchy-Green 
        deformation tensor; Psi, the microdeformation tensor; 
        and Gamma, the higher order deformation tensor.
        
        Order of fparams:
            fparams = {lambda, mu, eta, tau, kappa, nu, sigma, tau1, ... tau11}
        
        */
        
        //!Common tensors
        tensor::Tensor23 ITEN = tensor::eye(); //!The second order identity tensor
        
        //!Copy over the deformation measures
        tensor::Tensor23 C     = C_in;
        tensor::Tensor23 Cinv  = C.inverse();
        tensor::Tensor23 Psi   = Psi_in;
        tensor::Tensor33 Gamma = Gamma_in;
        
        //!Compute the strain measures
        tensor::Tensor23 macro_E = 0.5*(C - ITEN); //The macro Green-Lagrange strain
        tensor::Tensor23 micro_E = Psi - ITEN;     //The micro equivalent of the Green-Lagrange strain
        
        //!Initialize the stiffness tensors
        tensor::Tensor43 A_stiffness({3,3,3,3});
        tensor::Tensor43 B_stiffness({3,3,3,3});
        tensor::Tensor63 C_stiffness({3,3,3,3,3,3});
        tensor::Tensor43 D_stiffness({3,3,3,3});
        
        //!Compute the stiffness tensors
        generate_A_stiffness(fparams,A_stiffness);
        generate_B_stiffness(fparams,B_stiffness);
        generate_C_stiffness(fparams,C_stiffness);
        generate_D_stiffness(fparams,D_stiffness);
        
        //!Zero the stress measures
        PK2_stress.data.setZero();   //Zero out the Second Piola-Kirchhoff stress
        SIGMA_stress.data.setZero(); //Zero out the symmetric micro-stress tensor
        M_stress.data.setZero();     //Zero out the higher order couple stress
        
        //!Compute the stresses
        compute_stresses(A_stiffness, B_stiffness,  C_stiffness, D_stiffness,
                         C,           Cinv,         Gamma,
                         macro_E,     micro_E,      ITEN,
                         PK2_stress,  SIGMA_stress, M_stress);
    }
    
    void get_stress(const std::vector< double > &fparams,  const std::vector< int > &iparams,                                         //Material parameters
                    const tensor::Tensor23& C_in,          const tensor::Tensor23& Psi_in,      const tensor::Tensor33& Gamma_in,     //Deformation measures
                          tensor::Tensor23& PK2_stress,          tensor::Tensor23& SIGMA_stress,      tensor::Tensor33& M_stress,     //Stress measures
                          tensor::Tensor43& dPK2dC,              tensor::Tensor43& dPK2dPsi,          tensor::Tensor53& dPK2dGamma,   //Tangents of PK2 stress
                          tensor::Tensor43& dSIGMAdC,            tensor::Tensor43& dSIGMAdPsi,        tensor::Tensor53& dSIGMAdGamma, //Tangents of symmetric stress
                          tensor::Tensor53& dMdC,                tensor::Tensor53& dMdPsi,            tensor::Tensor63& dMdGamma      //Tangents of couple stress
                    ){
        /*!========================================
        |              get_stress              |
        ========================================
        
        Compute and update the values of the second Piola 
        Kirchhoff, the symmetric, and higher order stresses.
        
        The incoming values are the floating point parameters 
        for the model, the integer parameters for the model, 
        the deformation metrics C, the right Cauchy-Green 
        deformation tensor; Psi, the microdeformation tensor; 
        and Gamma, the higher order deformation tensor.
        
        This function computes and returns the tangent as well as 
        the stress.
        
        Order of fparams:
            fparams = {lambda, mu, eta, tau, kappa, nu, sigma, tau1, ... tau11}
        
        */
        
        //!Common tensors
        tensor::Tensor23 ITEN = tensor::eye(); //!The second order identity tensor
        
        //!Copy over the deformation measures
        tensor::Tensor23 C     = C_in;
        tensor::Tensor23 Cinv  = C.inverse();
        tensor::Tensor23 Psi   = Psi_in;
        tensor::Tensor33 Gamma = Gamma_in;
        
        //!Compute the strain measures
        tensor::Tensor23 macro_E = 0.5*(C - ITEN); //The macro Green-Lagrange strain
        tensor::Tensor23 micro_E = Psi - ITEN;     //The micro equivalent of the Green-Lagrange strain
        
        //!Initialize the stiffness tensors
        tensor::Tensor43 A_stiffness({3,3,3,3});
        tensor::Tensor43 B_stiffness({3,3,3,3});
        tensor::Tensor63 C_stiffness({3,3,3,3,3,3});
        tensor::Tensor43 D_stiffness({3,3,3,3});
        
        //!Compute the stiffness tensors
        generate_A_stiffness(fparams,A_stiffness);
        generate_B_stiffness(fparams,B_stiffness);
        generate_C_stiffness(fparams,C_stiffness);
        generate_D_stiffness(fparams,D_stiffness);
        
        //!Zero the stress measures
        PK2_stress.data.setZero();   //Zero out the Second Piola-Kirchhoff stress
        SIGMA_stress.data.setZero(); //Zero out the symmetric micro-stress tensor
        M_stress.data.setZero();     //Zero out the higher order couple stress
        
        //!Compute the stresses
        compute_stresses(A_stiffness, B_stiffness,  C_stiffness, D_stiffness,
                         C,           Cinv,         Gamma,
                         macro_E,     micro_E,      ITEN,
                         PK2_stress,  SIGMA_stress, M_stress);
        
        //!=
        //!| Compute the tangent
        //!=
        
        //!Common terms
        tensor::Tensor43 term1T({3,3,3,3});
        tensor::Tensor43 term2T({3,3,3,3});
        tensor::Tensor43 term3T({3,3,3,3});
        
        tensor::Tensor43 dCinvdC({3,3,3,3});
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    for(int L=0; L<3; L++){
                        dCinvdC(I,J,K,L) = -Cinv(I,K)*Cinv(L,J);
                    }
                }
            }
        }
        
        //!Compute tangents w.r.t. C
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int O=0; O<3; O++){
                    for(int P=0; P<3; P++){
                        term1T(I,J,O,P) = 0.5*A_stiffness(I,J,O,P);
                        
                        for(int Q=0; Q<3; Q++){
                            for(int R=0; R<3; R++){
                                
                                term2T(I,J,O,P) += 0.5*D_stiffness(I,Q,O,P)*(micro_E(R,Q)+ITEN(R,Q))*Cinv(J,R);
                                term3T(I,J,O,P) += 0.5*D_stiffness(J,Q,O,P)*(micro_E(R,Q)+ITEN(R,Q))*Cinv(I,R);
                            }
                        }
                            
                        for(int Q=0; Q<3; Q++){
                            for(int R=0; R<3; R++){
                                for(int K=0; K<3; K++){
                                    for(int L=0; L<3; L++){
                                        term2T(I,J,O,P) += (B_stiffness(I,Q,K,L)*micro_E(K,L)+D_stiffness(I,Q,K,L)*macro_E(K,L))*(micro_E(R,Q)+ITEN(R,Q))*dCinvdC(J,R,O,P);
                                        term3T(I,J,O,P) += (B_stiffness(J,Q,K,L)*micro_E(K,L)+D_stiffness(J,Q,K,L)*macro_E(K,L))*(micro_E(R,Q)+ITEN(R,Q))*dCinvdC(I,R,O,P);
                                    }
                                }

                                for(int L=0; L<3; L++){
                                    for(int M=0; M<3; M++){
                                        for(int N=0; N<3; N++){
                                            for(int S=0; S<3; S++){
                                                term2T(I,J,O,P) += C_stiffness(I,Q,R,L,M,N)*Gamma(L,M,N)*Gamma(S,Q,R)*dCinvdC(J,S,O,P);
                                                term3T(I,J,O,P) += C_stiffness(J,Q,R,L,M,N)*Gamma(L,M,N)*Gamma(S,Q,R)*dCinvdC(I,S,O,P);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        dPK2dC(I,J,O,P)   = term1T(I,J,O,P) + term2T(I,J,O,P);
                        dSIGMAdC(I,J,O,P) = term1T(I,J,O,P) + (term2T(I,J,O,P)+term3T(I,J,O,P));
                    }
                }
            }
        }
            
        //!Compute tangents w.r.t. Psi
        term1T.data.setZero();
        term2T.data.setZero();
        term3T.data.setZero();
            
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int O=0; O<3; O++){
                    for(int P=0; P<3; P++){
                        dPK2dPsi(I,J,O,P)   = D_stiffness(I,J,O,P);
                        dSIGMAdPsi(I,J,O,P) = D_stiffness(I,J,O,P);
                    
                        for(int Q=3; Q<3; Q++){
                            for(int R=0; R<3; R++){
                                term2T(I,J,O,P) += B_stiffness(I,Q,O,P)*(micro_E(R,Q)+ITEN(R,Q))*Cinv(J,R);
                                term3T(I,J,O,P) += B_stiffness(J,Q,O,P)*(micro_E(R,Q)+ITEN(R,Q))*Cinv(I,R);
                            }
                        }
                            
                        for(int K=0; K<3; K++){
                            for(int L=0; L<3; L++){
                                term2T(I,J,O,P) += (B_stiffness(I,P,K,L)*micro_E(K,L)+D_stiffness(I,P,K,L)*macro_E(K,L))*Cinv(J,O);
                                term3T(I,J,O,P) += (B_stiffness(J,P,K,L)*micro_E(K,L)+D_stiffness(J,P,K,L)*macro_E(K,L))*Cinv(I,O);
                            }
                        }
                                    
                        dPK2dPsi(I,J,O,P)   = term2T(I,J,O,P);
                        dSIGMAdPsi(I,J,O,P) = term2T(I,J,O,P) + term3T(I,J,O,P);
                    }
                }
            }
        }
        
        //!Compute tangents w.r.t. Gamma
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int T=0; T<3; T++){
                    for(int U=0; U<3; U++){
                        for(int V=0; V<3; V++){
                            for(int S=0; S<3; S++){
                                for(int Q=0; Q<3; Q++){
                                    for(int R=0; R<3; R++){
                                    
                                        dPK2dGamma(I,J,T,U,V)   += C_stiffness(I,Q,R,T,U,V)*Cinv(J,S)*Gamma(S,Q,R) + C_stiffness(I,U,V,S,Q,R)*Gamma(S,Q,R)*Cinv(J,T);
                                        dSIGMAdGamma(I,J,T,U,V) += C_stiffness(I,Q,R,T,U,V)*Cinv(J,S)*Gamma(S,Q,R) + C_stiffness(I,U,V,S,Q,R)*Gamma(S,Q,R)*Cinv(J,T)
                                                                  +C_stiffness(J,Q,R,T,U,V)*Cinv(I,S)*Gamma(S,Q,R) + C_stiffness(J,U,V,S,Q,R)*Gamma(S,Q,R)*Cinv(I,T);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    dMdGamma.data = C_stiffness.data;
    return;
}
    
    void compute_stresses(const tensor::Tensor43& A_stiffness, const tensor::Tensor43& B_stiffness,  const tensor::Tensor63& C_stiffness, const tensor::Tensor43& D_stiffness,
                          const tensor::Tensor23& C,           const tensor::Tensor23& Cinv,         const tensor::Tensor33& Gamma,
                          const tensor::Tensor23& macro_E,     const tensor::Tensor23& micro_E,      const tensor::Tensor23& ITEN,
                                tensor::Tensor23& PK2_stress,        tensor::Tensor23& SIGMA_stress,       tensor::Tensor33& M_stress){
        /*!==========================
        |    compute_stresses    |
        ==========================
        
        Compute the stresses from the stiffness and strain tensors
        
        */
        
        //!Create temporary storage tensors for the common terms
        tensor::Tensor23 term1({3,3});
        tensor::Tensor23 term2({3,3});
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    for(int L=0; L<3; L++){
                        term1(I,J) += A_stiffness(I,J,K,L)*macro_E(K,L) + D_stiffness(I,J,K,L)*micro_E(K,L);
                    }
                }
                for(int K=0; K<3; K++){
                    for(int L=0; L<3; L++){
                        for(int Q=0; Q<3; Q++){
                            for(int R=0; R<3; R++){
                                term2(I,J) +=  (B_stiffness(I,Q,K,L)*micro_E(K,L)+D_stiffness(I,Q,K,L)*macro_E(K,L))*(micro_E(R,Q)+ITEN(R,Q))*Cinv(J,R);
                            }
                        }
                    }
                }
                
                for(int Q=0; Q<3; Q++){
                    for(int R=0; R<3; R++){
                        for(int M=0; M<3; M++){
                            for(int N=0; N<3; N++){
                                for(int L=0; L<3; L++){
                                    for(int S=0; S<3; S++){
                                        term2(I,J) += C_stiffness(I,Q,R,L,M,N)*Gamma(L,M,N)*Cinv(J,S)*Gamma(S,Q,R);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //!Compute the stress measures
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                PK2_stress(I,J)   = term1(I,J)+term2(I,J);
                SIGMA_stress(I,J) = term1(I,J)+(term2(I,J)+term2(J,I));
            }
        }
        
        for(int I=0; I<3; I++){
            for(int J=0; J<3; J++){
                for(int K=0; K<3; K++){
                    for(int M=0; M<3; M++){
                        for(int N=0; N<3; N++){
                            for(int O=0; O<3; O++){
                                M_stress(I,J,K) += C_stiffness(I,J,K,M,N,O)*Gamma(M,N,O);
                            }
                        }
                    }
                }
            }
        }
        
        return;
    }
    
    void generate_A_stiffness(const std::vector< double > &fparams, tensor::Tensor43& A_stiffness){
        /*!==============================
        |    generate_A_stiffness    |
        ==============================
        
        Generate the fourth order stiffness tensor 
        A from Regueiro 2010 Eqn 2.106
        
        */
        
        tensor::Tensor23 I = tensor::eye();      //!Second order identity tensor
        double lambda = fparams[0];              //!lambda micromorphic material parameter
        double mu     = fparams[1];              //!mu micromorphic material parameter
        
        for(int K=0; K<3; K++){
            for(int L=0; L<3; L++){
                for(int M=0; M<3; M++){
                    for(int N=0; N<3; N++){
                        A_stiffness(K,L,M,N) = lambda*I(K,L)*I(M,N) + mu*(I(K,M)*I(L,N) + I(K,N)*I(L,M));
                    }
                }
            }
        }
        
        return;
    }
    
    void generate_B_stiffness(const std::vector< double > &fparams, tensor::Tensor43& B_stiffness){
        /*!==============================
        |    generate_B_stiffness    |
        ==============================
        
        Generate the fourth order stiffness tensor 
        B from Regueiro 2010 Eqn 2.107
        
        */
        
        tensor::Tensor23 I = tensor::eye();      //!Second order identity tensor
        double eta   = fparams[2];               //!eta micromorphic material parameter
        double tau   = fparams[3];               //!tau micromorphic material parameter
        double kappa = fparams[4];               //!kappa micromorphic material parameter
        double nu    = fparams[5];               //!nu micromorphic material parameter
        double sigma = fparams[6];               //!sigma micromorphic material parameter
        
        for(int K=0; K<3; K++){
            for(int L=0; L<3; L++){
                for(int M=0; M<3; M++){
                    for(int N=0; N<3; N++){
                        B_stiffness(K,L,M,N) = (eta-tau)*I(K,L)*I(M,N) + kappa*I(K,M)*I(L,N) + nu*I(K,N)*I(L,M)
                                               -sigma*(I(K,M)*I(L,N) + I(K,N)*I(L,M));
                    }
                }
            }
        }
        
        return;
    }
    
    void generate_C_stiffness(const std::vector< double > &fparams, tensor::Tensor63& C_stiffness){
        /*!==============================
        |    generate_C_stiffness    |
        ==============================
        
        Generate the sixth order stiffness tensor 
        C from Regueiro 2010 Eqn 2.108
        
        */
        
        tensor::Tensor23 I = tensor::eye();          //!Second order identity tensor
        double tau1  = fparams[ 7];                  //!tau1  micromorphic material parameter
        double tau2  = fparams[ 8];                  //!tau2  micromorphic material parameter
        double tau3  = fparams[ 9];                  //!tau3  micromorphic material parameter
        double tau4  = fparams[10];                  //!tau4  micromorphic material parameter
        double tau5  = fparams[11];                  //!tau5  micromorphic material parameter
        double tau6  = fparams[12];                  //!tau6  micromorphic material parameter
        double tau7  = fparams[13];                  //!tau7  micromorphic material parameter
        double tau8  = fparams[14];                  //!tau8  micromorphic material parameter
        double tau9  = fparams[15];                  //!tau9  micromorphic material parameter
        double tau10 = fparams[16];                  //!tau10 micromorphic material parameter
        double tau11 = fparams[17];                  //!tau11 micromorphic material parameter
        
        for(int K=0; K<3; K++){
            for(int L=0; L<3; L++){
                for(int M=0; M<3; M++){
                    for(int N=0; N<3; N++){
                        for(int P=0; P<3; P++){
                            for(int Q=0; Q<3; Q++){
                                C_stiffness(K,L,M,N,P,Q) =  tau1*(I(K,L)*I(M,N)*I(P,Q) + I(K,Q)*I(L,M)*I(N,P))
                                                           +tau2*(I(K,L)*I(M,P)*I(N,Q) + I(K,M)*I(L,Q)*I(N,P))
                                                           +tau3*I(K,L)*I(M,Q)*I(N,P)  + tau4*I(K,N)*I(L,M)*I(P,Q)
                                                           +tau5*(I(K,M)*I(L,N)*I(P,Q) + I(K,P)*I(L,M)*I(N,Q))
                                                           +tau6*I(K,M)*I(L,P)*I(N,Q)  + tau7*I(K,N)*I(L,P)*I(M,Q)
                                                           +tau8*(I(K,P)*I(L,Q)*I(M,N) + I(K,Q)*I(L,N)*I(M,P)) + tau9*I(K,N)*I(L,Q)*I(M,P)
                                                           +tau10*I(K,P)*I(L,N)*I(M,Q) + tau11*I(K,Q)*I(L,P)*I(M,N);
                            }
                        }
                    }
                }
            }
        }
        
        return;
    }
    
    void generate_D_stiffness(const std::vector< double > &fparams, tensor::Tensor43& D_stiffness){
        /*!==============================
        |    generate_D_stiffness    |
        ==============================
        
        Generate the fourth order stiffness tensor 
        D from Regueiro 2010 Eqn 2.109
        
        */
        
        tensor::Tensor23 I = tensor::eye();      //!Second order identity tensor
        double tau   = fparams[3];               //!tau micromorphic material parameter
        double sigma = fparams[6];               //!sigma micromorphic material parameter
        
        for(int K=0; K<3; K++){
            for(int L=0; L<3; L++){
                for(int M=0; M<3; M++){
                    for(int N=0; N<3; N++){
                        D_stiffness(K,L,M,N) = tau*I(K,L)*I(M,N) + sigma*(I(K,M)*I(L,N) + I(K,N)*I(L,M));
                    }
                }
            }
        }
        
        return;
    }
}