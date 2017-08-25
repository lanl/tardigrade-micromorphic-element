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
    
    void get_stress(const std::vector< double >&, const std::vector< int >&,
                          const tensor::Tensor&,  const tensor::Tensor&, const tensor::Tensor&,
                          tensor::Tensor&,        tensor::Tensor&,       tensor::Tensor&);
    
    tensor::Tensor generate_A_stiffness(const std::vector< double > &);
    tensor::Tensor generate_B_stiffness(const std::vector< double > &);
    tensor::Tensor generate_C_stiffness(const std::vector< double > &);
    tensor::Tensor generate_D_stiffness(const std::vector< double > &);
}