/*!===========================================================================
   |                                                                         |
   |                         micromorphic_measures.h                         |
   |                                                                         |
   ===========================================================================
   | The header file for a wrapper that converts variables and their         |
   | gradients generated by MOOSE into deformation measures which can be     |
   | used to compute the micromorphic stress measures and their tangents.    |
   | This is done to avoid any possible assumptions of symmetry which could  |
   | be present in the Tensor Mechanics physics module.                      |
   ===========================================================================
   | Dependencies:                                                           |
   |     Eigen: A matrix library available at eigen.tuxfamily.org            |
   ===========================================================================
   */
   
   #include <Eigen/Dense>
   
   