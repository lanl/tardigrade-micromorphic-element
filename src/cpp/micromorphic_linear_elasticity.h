/*!=======================================================
  |                                                     |
  |         micromorphic_linear_elasticity.h            |
  |                                                     |
  -------------------------------------------------------
  | The header file for the definition of a             |
  | micromorphic linear elasticity.                     |
  -------------------------------------------------------
  | Notes: Micromorphic constitutive models should be   |
  |        developed in the namespace micro_material    |
  |        and have the function get_stress. This       |
  |        function should read in the right Cauchy     |
  |        green deformation tensor, Psi, Gamma, and    |
  |        write the PK2 stress, the symmetric stress   |
  |        in the reference configuration, and the      |
  |        higher order couple stress in the reference  |
  |        configuration. (ADDITIONAL VALUES WILL BE    |
  |        ADDED OVER TIME).                            |
  =======================================================
  | Dependencies:                                       |
  | tensor:      The class which defines tensor access  |
  |              to an underlying Eigen matrix. This    |
  |              may result in a somewhat slower result |
  |              however it should allow for a faster   |
  |              implementation.                        |
  =======================================================*/
  
#include <vector>
  
namespace micro_material{
    
    void get_stress( const std::vector< double >&, const std::vector< int >&,
                     const tensor::Tensor23&,      const tensor::Tensor23&,   const tensor::Tensor33&,
                           tensor::Tensor23&,            tensor::Tensor23&,         tensor::Tensor33&
                   );
    
    void get_stress( const std::vector< double >&, const std::vector< int >&,
                     const tensor::Tensor23&,      const tensor::Tensor23&,   const tensor::Tensor33&,
                           tensor::Tensor23&,            tensor::Tensor23&,         tensor::Tensor33&,
                           tensor::Tensor43&,            tensor::Tensor43&,         tensor::Tensor53&,
                           tensor::Tensor43&,            tensor::Tensor43&,         tensor::Tensor53&,
                           tensor::Tensor53&,            tensor::Tensor53&,         tensor::Tensor63&
                   );
    
    void compute_stresses(const tensor::Tensor43& A_stiffness, const tensor::Tensor43& B_stiffness,  const tensor::Tensor63& C_stiffness, const tensor::Tensor43& D_stiffness,
                          const tensor::Tensor23& C,           const tensor::Tensor23& Cinv,         const tensor::Tensor33& Gamma,
                          const tensor::Tensor23& macro_E,     const tensor::Tensor23& micro_E,      const tensor::Tensor23& ITEN,
                                tensor::Tensor23& PK2_stress,        tensor::Tensor23& SIGMA_stress,       tensor::Tensor33& M_stress);
    
    void generate_A_stiffness(const std::vector< double > &, tensor::Tensor43&);
    void generate_B_stiffness(const std::vector< double > &, tensor::Tensor43&);
    void generate_C_stiffness(const std::vector< double > &, tensor::Tensor63&);
    void generate_D_stiffness(const std::vector< double > &, tensor::Tensor43&);
}