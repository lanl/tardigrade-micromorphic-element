/*!=======================================================
  |                                                     |
  |                  micro_element.h                    |
  |                                                     |
  -------------------------------------------------------
  | The header file for the definition of a             |
  | micromorphic continuum element.                     |
  =======================================================
  | Dependencies:                                       |
  | tensor:      The class which defines tensor access  |
  |              to an underlying Eigen matrix. This    |
  |              may result in a somewhat slower result |
  |              however it should allow for a faster   |
  |              implementation.                        |
  =======================================================*/
  
#include <iostream>
#include <vector>
#include <tensor.h>
  
namespace micro_element
{
    class Hex8{
        /*!===
         |
         | H E X 8
         |
        ===
        
        An eight noded hexehedral element.
        
        The nodes have all degrees of freedom i.e. the 
        deformation u and the micro-deformation phi are 
        defined at all nodes.
        
        The element is defined such that different constitutive 
        equations can be implemented without extreme difficulty.
        
        */
        
        public:
        
            //!==
            //!|
            //!| Attribute Definitions
            //!|
            //!==
            
            std::vector< std::vector< double > > reference_coords; //!The reference coordinates of the element's nodes
            std::vector< std::vector< double > > current_coords;   //!The current coordinates of the element's nodes
            
            //!==
            //!|
            //!| Constructors
            //!|
            //!==
        
            //Constructors
            Hex8();
            Hex8(std::vector< std::vector< double > >);
            
            //!==
            //!|
            //!| Operators
            //!|
            //!==
    };
}