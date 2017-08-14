/*!=======================================================
  |                                                     |
  |                 micro_element.cpp                   |
  |                                                     |
  -------------------------------------------------------
  | The source file for the definition of a             |
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
    
    /*!===
     |
     | Constructors
     |
    ===*/
    
    Hex8::Hex8(){
        
    }
    
    Hex8::Hex8(std::vector< std::vector< double > > rcs){
        /*!====================
        |        Hex8       |
        =====================
        
        The constructor for a hexehedral element when 
        given the nodal reference coordinates.
        
        Input:
           rcs: A vector of vectors of doubles which are the
                coordinates of the nodes.
                
                The nodes are ordered in a counter clockwise 
                manner i.e.
                
               4,8      3,7
                o--------o
                |        |
                |        |
                |        |
                o--------o
               1,5      2,6
               
               where the comma indicates the ``upper layer'' of the
               hexehedral element.
        */
    }
    class Hex8{
        /*!==
        |
        | H E X 8
        |
        ==
        
        An eight noded hexehedral element.
        
        The nodes have all degrees of freedom i.e. the 
        deformation u and the micro-deformation phi are 
        defined at all nodes.
        
        The element is defined such that different constitutive 
        equations can be implemented without extreme difficulty.
        
        */
          
        //Constructors
        Hex8();
        Hex8(std::vector< std::vector< double > >);
    };
}