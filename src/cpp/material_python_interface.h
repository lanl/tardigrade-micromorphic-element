#include<micromorphic_material_library.h>

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
                        std::string &output_message );
}
