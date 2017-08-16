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
#include <tensor.h>
  
namespace micro_material{
    
    void get_stress(const tensor::Tensor&   C, const tensor::Tensor&   Psi, const tensor::Tensor& Gamma,
                          tensor::Tensor& PK2,       tensor::Tensor& SIGMA,       tensor::Tensor& M){
        /*!========================================
        |              get_stress              |
        ========================================
        
        Compute and update the values of the second Piola 
        Kirchhoff, the symmetric, and higher order stresses.
        
        The incoming values are the deformation metrics C, the 
        right Cauchy-Green deformation tensor; Psi, the micro 
        deformation tensor; and Gamma, the higher order deformation 
        tensor.
        
        */
    }
}