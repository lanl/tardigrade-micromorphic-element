#include<material_python_interface.h>

namespace materialPythonInterface{

    int evaluate_model( const std::string &model_name,
                        const std::vector< double > &time, const std::vector< double > &fparams,
                        const double ( &current_grad_u )[ 3 ][ 3 ], const double ( &current_phi )[ 9 ], const double ( &current_grad_phi )[ 9 ][ 3 ],
                        const double ( &previous_grad_u )[ 3 ][ 3 ], const double ( &previous_phi )[ 9 ], const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &DPK2Dgrad_u, std::vector< std::vector< double > > &DPK2Dphi, std::vector< std::vector< double > > &DPK2Dgrad_phi,
                        std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi, std::vector< std::vector< double > > &DSIGMADgrad_phi,
                        std::vector< std::vector< double > > &DMDgrad_u, std::vector< std::vector< double > > &DMDphi, std::vector< std::vector< double > > &DMDgrad_phi,
                        std::vector< std::vector< double > > &ADD_TERMS, std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                        std::string &output_message ){
        /*!
         * A wrapper for the evaluate model method of micromorphic_material_library::IMaterial which also includes the call to the material name
         * 
         * Evaluate the jacobian of the material model using a numeric gradient
         *
	 * :param const std::string &model_name: The name of the model to be evaluated
         * :param const std::vector< double > &time: The current time and the timestep
         *     [ current_t, dt ]
         * :param const std::vector< double > ( &fparams ): The parameters for the constitutive model
         * :param const double ( &current_grad_u )[ 3 ][ 3 ]: The current displacement gradient
         *     Assumed to be of the form [ [ u_{1,1}, u_{1,2}, u_{1,3} ],
         *                                 [ u_{2,1}, u_{2,2}, u_{2,3} ],
         *                                 [ u_{3,1}, u_{3,2}, u_{3,3} ] ]
         * :param const double ( &current_phi )[ 9 ]: The current micro displacment values.
         *     Assumed to be of the form [ \phi_{11}, \phi_{12}, \phi_{13}, \phi_{21}, \phi_{22}, \phi_{23}, \phi_{31}, \phi_{32}, \phi_{33} ]
         * :param const double ( &current_grad_phi )[ 9 ][ 3 ]: The current micro displacement gradient
         *     Assumed to be of the form [ [ \phi_{11,1}, \phi_{11,2}, \phi_{11,3} ],
         *                                 [ \phi_{12,1}, \phi_{12,2}, \phi_{12,3} ],
         *                                 [ \phi_{13,1}, \phi_{13,2}, \phi_{13,3} ],
         *                                 [ \phi_{21,1}, \phi_{21,2}, \phi_{21,3} ],
         *                                 [ \phi_{22,1}, \phi_{22,2}, \phi_{22,3} ],
         *                                 [ \phi_{23,1}, \phi_{23,2}, \phi_{23,3} ],
         *                                 [ \phi_{31,1}, \phi_{31,2}, \phi_{31,3} ],
         *                                 [ \phi_{32,1}, \phi_{32,2}, \phi_{32,3} ],
         *                                 [ \phi_{33,1}, \phi_{33,2}, \phi_{33,3} ] ]
         * :param const double ( &previous_grad_u )[ 3 ][ 3 ]: The previous displacement gradient.
         * :param const double ( &previous_phi )[ 9 ]: The previous micro displacement.
         * :param const double ( &previous_grad_phi )[ 9 ][ 3 ]: The previous micro displacement gradient.
         * :param std::vector< double > &SDVS: The previously converged values of the state variables
         * :param std::vector< double > &current_ADD_DOF: The current values of the additional degrees of freedom
         * :param std::vector< std::vector< double > > &current_ADD_grad_DOF: The current values of the gradients of the
         *     additional degrees of freedom
         * :param std::vector< double > &current_PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ]
         * :param std::vector< double > &current_SIGMA: The current value of the reference micro stress. The format is
         *     [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ]
         * :param std::vector< double > &current_M: The current value of the reference higher order stress. The format is
         *     [ M_{111}, M_{112}, M_{113}, M_{121}, M_{122}, M_{123}, M_{131}, M_{132}, M_{133},
         *       M_{211}, M_{212}, M_{213}, M_{221}, M_{222}, M_{223}, M_{231}, M_{232}, M_{233},
         *       M_{311}, M_{312}, M_{313}, M_{321}, M_{322}, M_{323}, M_{331}, M_{332}, M_{333} ]
         * :param std::vector< std::vector< double > > &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the
         *     gradient of macro displacement.
         * :param std::vector< std::vector< double > > &DPK2Dphi: The Jacobian of the PK2 stress w.r.t. the
         *     micro displacement.
         * :param std::vector< std::vector< double > > &DPK2Dgrad_phi: The Jacobian of the PK2 stress w.r.t.
         *     the gradient of the micro displacement.
         * :param std::vector< std::vector< double > > &DSIGMAdgrad_u: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the macro displacement.
         * :param std::vector< std::vector< double > > &DSIGMAdphi: The Jacobian of the reference symmetric micro
         *     stress w.r.t. the micro displacement.
         * :param std::vector< std::vector< double > > &DSIGMAdgrad_phi: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the micro displacement.
         * :param std::vector< std::vector< double > > &DMDgrad_u: The Jacobian of the reference higher order
         *     stress w.r.t. the gradient of the macro displacement.
         * :param std::vector< std::vector< double > > &DMDphi: The Jacobian of the reference higher order stress
         *     w.r.t. the micro displacement.
         * :param std::vector< std::vector< double > > &DMDgrad_phi: The Jacobian of the reference higher order stress
         *     w.r.t. the gradient of the micro displacement.
         * :param std::vector< std::vector< double > > &ADD_TERMS: Additional terms ( unused )
         * :param std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS: The jacobians of the additional
         *     terms w.r.t. the deformation ( unused )
         * :param std::string &output_message: The output message string.
         * :param double delta = 1e-6: The perturbation to be applied to the incoming degrees of freedom.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */
    
        auto &factory = micromorphic_material_library::MaterialFactory::Instance( );
        auto material = factory.GetMaterial( model_name );
    
        int errorCode = material->evaluate_model( time, fparams,
                                                  current_grad_u, current_phi, current_grad_phi,
                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                  SDVS,
                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                  PK2,           SIGMA,         M,
                                                  DPK2Dgrad_u,   DPK2Dphi,      DPK2Dgrad_phi,
                                                  DSIGMADgrad_u, DSIGMADphi,    DSIGMADgrad_phi,
                                                  DMDgrad_u,     DMDphi,        DMDgrad_phi,
                                                  ADD_TERMS,     ADD_JACOBIANS, output_message );
    
        return errorCode;
    
    }
}
