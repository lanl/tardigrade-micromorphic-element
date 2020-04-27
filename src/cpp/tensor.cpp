/*!=======================================================
  |                                                     |
  |                     tensor.cpp                      |
  |                                                     |
  -------------------------------------------------------
  | The source file for a class definition that allows  |
  | n-dimensional tensors. The data for the tensor is   |
  | stored in a 2 dimensional array and there is logic  |
  | which allows referencing to the correct indices to  |
  | occur.                                              |
  =======================================================
  | Dependencies:                                       |
  | Eigen: An implementation of various matrix          |
  |        commands. The implementation of the data     |
  |        matrix uses such a matrix.                   |
  =======================================================*/
  
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>
#include <Eigen/Dense>
#include <tensor.h>

namespace tensor{
    //!==
    //!|
    //!| Functions
    //!|
    //!==
        
    Tensor23 eye(){
        /*!Return the second order identity tensor
        
        Populate the second order identity tensor.
        This function is useful so that a consistent identity tensor can be used 
        with different storage schemes.*/
        
        //Initialize the matrix
        
        Tensor23 I({3,3});
        
        for(int i=0; i<3; i++){
            I(i,i) = 1.;
        }
        return I;
    }
    Tensor43 FOT_eye(){
        /*!Return the fourth order identity tensor
        
        Populate the fourth order identity tensor.
        This function is useful so that a consistent identity tensor can be used
        with different storage schemes.*/
        
        //Initialize the tensor
        Tensor43 FOTI({3,3,3,3});
        Tensor23 I    = eye(); //Get the second order identity tensor
        
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int k=0; k<3; k++){
                    for(int l=0; l<3; l++){
                        FOTI(i,j,k,l) = I(i,k)*I(j,l);
                    }
                }
            }
        }
        
        return FOTI;
    }
}