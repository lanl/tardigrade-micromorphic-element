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

    errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &microDeformationGradient,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress );

    errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &microDeformationGradient,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdXi, variableMatrix &dPK2StressdGradXi,
                                        variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdXi,
                                        variableMatrix &dMicroStressdGradXi, variableMatrix &dMdF, variableMatrix &dMdGradXi );

    errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma );

    errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                         variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidXi, 
                                         variableMatrix &dGammadF, variableMatrix &dGammadGradXi );

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

}

#endif
