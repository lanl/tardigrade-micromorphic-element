/*
 * micromorphic_linear_elasticity.h
 *
 * An implimentation of linear elasticity in the micromorphic context.
 *
 * Based around a quadratic form of the Helmholtz free energy:
 * \rho \psi = \frac{1}{2} E_{IJ} A_{IJKL} E_{KL} + \frac{1}{2} \mathcal{E}_{IJ} B_{IJKL} \mathcal{E}_{KL} 
 *           + \frac{1}{2} \Gamma_{IJK} C_{IJKLMN} \Gamma_{LMN} + E_{IJ} D_{IJKL} \mathcal{E}_{KL}
 */

#ifndef MICROMORPHIC_LINEAR_ELASTICITY_H
#define MICROMORPHIC_LINEAR_ELASTICITY_H

#include<error_tools.h>
#define USE_EIGEN
#include<vector_tools.h>
#include<constitutive_tools.h>
#include<micromorphic_tools.h>
#include<micromorphic_material_library.h>

namespace micromorphicLinearElasticity{

    typedef micromorphicTools::variableType variableType;
    typedef micromorphicTools::variableVector variableVector;
    typedef micromorphicTools::variableMatrix variableMatrix;

    typedef micromorphicTools::parameterType parameterType;
    typedef micromorphicTools::parameterVector parameterVector;
    typedef micromorphicTools::parameterMatrix parameterMatrix;

    typedef micromorphicTools::constantType constantType;
    typedef micromorphicTools::constantVector constantVector;
    typedef micromorphicTools::constantMatrix constantMatrix;

    typedef errorTools::Node errorNode;
    typedef errorNode* errorOut;

    struct cout_redirect{
        cout_redirect( std::streambuf * new_buffer)
            : old( std::cout.rdbuf( new_buffer ) )
        { }

        ~cout_redirect( ) {
            std::cout.rdbuf( old );
        }

        private:
            std::streambuf * old;
    };

    struct cerr_redirect{
        cerr_redirect( std::streambuf * new_buffer)
            : old( std::cerr.rdbuf( new_buffer ) )
        { }

        ~cerr_redirect( ) {
            std::cerr.rdbuf( old );
        }

        private:
            std::streambuf * old;
    };

    errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress );

    errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress,
                               variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdChi, variableMatrix &dCauchyStressdGradChi,
                               variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdChi, variableMatrix &dMicroStressdGradChi,
                               variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                               variableMatrix &dHigherOrderStressdGradChi );

    errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress );

    errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdChi, variableMatrix &dPK2StressdGradChi,
                                        variableMatrix &dReferenceMicroStressdF, variableMatrix &dReferenceMicroStressdChi,
                                        variableMatrix &dReferenceMicroStressdGradChi, variableMatrix &dMdF, variableMatrix &dMdGradChi );

    errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                       const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress );

    errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                       const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress,
                                                       variableMatrix &dPK2StressdRCG, variableMatrix &dPK2StressdPsi,
                                                       variableMatrix &dPK2StressdGamma,
                                                       variableMatrix &dReferenceMicroStressdRCG,
                                                       variableMatrix &dReferenceMicroStressdPsi,
                                                       variableMatrix &dReferenceMicroStressdGamma,
                                                       variableMatrix &dMdGamma );

    errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress );

    errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress,
                                         variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdPK2Stress,
                                         variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdReferenceMicroStress,
                                         variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                         variableMatrix &dHigherOrderStressdReferenceHigherOrderStress );

    errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma );

    errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                         variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidChi, 
                                         variableMatrix &dGammadF, variableMatrix &dGammadGradChi );

    errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const parameterVector &A, const parameterVector &D, variableVector &term1 );

    errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const parameterVector &A, const parameterVector &D, variableVector &term1,
                                        variableMatrix &dTerm1dGreenLagrangeStrain, variableMatrix &dTerm1dMicroStrain );

    errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &incCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2 );

    errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2, variableMatrix &dTerm2dGreenLagrangeStrain,
                                        variableMatrix &dTerm2dMicroStrain, variableMatrix &dTerm2dInvCPsi );

    errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                variableVector &referenceHigherOrderStress );

    errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C, 
                                                variableVector &referenceHigherOrderStress, 
                                                variableMatrix &dHigherOrderStressdGamma );

    errorOut computeLinearElasticTerm3( const variableVector &invCGamma,    
                                        const variableVector &referenceHigherOrderStress, variableVector &term3 );

    errorOut computeLinearElasticTerm3( const variableVector &invCGamma,
                                        const variableVector &referenceHigherOrderStress, variableVector &term3,
                                        variableMatrix &dTerm3dInvCGamma, variableMatrix &dTerm3dReferenceHigherOrderStress );

    errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi );

    errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                               variableMatrix &dInvRCGPsidRGG, variableMatrix &dInvRCGPsidPsi );

    errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma );

    errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                 variableMatrix &dInvRCGGammadRCG, variableMatrix &dInvRCGGammadGamma );

    errorOut formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A );

    errorOut formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                             const parameterType &nu,  const parameterType &sigma, parameterVector &B );

    errorOut formIsotropicC( const parameterVector &taus, parameterVector &C );

    errorOut formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D );

    errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation );

    errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation,
                                                     variableMatrix &dFdGradU, variableMatrix &dChidPhi,
                                                     variableMatrix &dGradChidGradPhi );

    errorOut extractMaterialParameters( const std::vector< double > &fparams,
                                        parameterVector &Amatrix, parameterVector &Bmatrix,
                                        parameterVector &Cmatrix, parameterVector &Dmatrix );

    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::string &output_message
                      );

    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &DPK2Dgrad_u,   std::vector< std::vector< double > > &DPK2Dphi,
                        std::vector< std::vector< double > > &DPK2Dgrad_phi,
                        std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi,
                        std::vector< std::vector< double > > &DSIGMADgrad_phi,
                        std::vector< std::vector< double > > &DMDgrad_u,     std::vector< std::vector< double > > &DMDphi,
                        std::vector< std::vector< double > > &DMDgrad_phi,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                        std::string &output_message
                      );

    class LinearElasticity: public micromorphic_material_library::IMaterial{
        /*!
         * The class which is called when evaluating a
         * linear elastic micromorphic constitutive model.
         *
         * This class registers linear elasticity into a library of models which can then
         * be called by name.
         */

        public:
            int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                std::vector< double > &SDVS,
                                const std::vector< double > &current_ADD_DOF,
                                const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                const std::vector< double > &previous_ADD_DOF,
                                const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                                std::vector< std::vector< double > > &ADD_TERMS,
                                std::string &output_message
                              )
            {
#ifdef DEBUG_MODE
                std::map< std::string, std::vector< double > > DEBUG;
#endif
                return micromorphicLinearElasticity::evaluate_model(
                                         time, fparams,
                                         current_grad_u, current_phi, current_grad_phi,
                                         previous_grad_u, previous_phi, previous_grad_phi,
                                         SDVS,
                                         current_ADD_DOF, current_ADD_grad_DOF,
                                         previous_ADD_DOF, previous_ADD_grad_DOF,
                                         PK2, SIGMA, M,
                                         ADD_TERMS,
                                         output_message
#ifdef DEBUG_MODE
                                         , DEBUG
#endif
                                         );
            }
        
            int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                std::vector< double > &SDVS,
                                const std::vector< double > &current_ADD_DOF,
                                const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                const std::vector< double > &previous_ADD_DOF,
                                const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                                std::vector< std::vector< double > > &DPK2Dgrad_u,
                                std::vector< std::vector< double > > &DPK2Dphi,
                                std::vector< std::vector< double > > &DPK2Dgrad_phi,
                                std::vector< std::vector< double > > &DSIGMADgrad_u,
                                std::vector< std::vector< double > > &DSIGMADphi,
                                std::vector< std::vector< double > > &DSIGMADgrad_phi,
                                std::vector< std::vector< double > > &DMDgrad_u,
                                std::vector< std::vector< double > > &DMDphi,
                                std::vector< std::vector< double > > &DMDgrad_phi,
                                std::vector< std::vector< double > > &ADD_TERMS,
                                std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                                std::string &output_message
                              )
            {
#ifdef DEBUG_MODE
                std::map< std::string, std::vector< double > > DEBUG;
#endif
                return micromorphicLinearElasticity::evaluate_model(
                                         time, fparams,
                                         current_grad_u, current_phi, current_grad_phi,
                                         previous_grad_u, previous_phi, previous_grad_phi,
                                         SDVS,
                                         current_ADD_DOF, current_ADD_grad_DOF,
                                         previous_ADD_DOF, previous_ADD_grad_DOF,
                                         PK2, SIGMA, M,
                                         DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                         DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                         DMDgrad_u, DMDphi, DMDgrad_phi,
                                         ADD_TERMS, ADD_JACOBIANS,
                                         output_message
#ifdef DEBUG_MODE
                                         , DEBUG
#endif
                                         );
            }
    };

    REGISTER_MATERIAL(LinearElasticity)
}

#endif
