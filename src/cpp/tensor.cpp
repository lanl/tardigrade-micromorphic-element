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
#include <cstdarg>
#include <Eigen/Dense>
#include <tensor.h>

namespace tensor{
    
    void Tensor::Tensor(){
        /*!The default constructor for the tensor. 
        Zeros the sizes of the shape and data attributes
        to save space*/
        
        shape.resize(0);
        data.resize(0,0);
        data_dimensions[0] = 0;
        data_dimensions[1] = 0;
    }
    
    void Tensor::Tensor(std::vector< unsigned int > dimensions){
        /*!A constructor for a tensor where the dimensions of the tensor are read in
        The tensor will be initialized to 0.*/
        
        //Assign the incoming dimensions to shape
        shape = dimensions;
        
        //Set the data matrix dimensions
        set_dimensions()
        
    }
    
    void Tensor::set_dimensions(){
        /*!Set the dimensions of the data storage of the tensor.
        
        This function initializes the contents of the tensor to zero.
        
        This is used in initialization to size the self.data matrix.
        It is not particularly useful after the initialization task.
        
        This is (one location) where the self.format string is important 
        as different formats could produce drastically different matrices.
        
        */
        
        /*Check if the format is the default format*/
        if(format.compare("default")){
            /*!Default format takes the first shape.size()/2 indices
            and places them vertically (rows) and then the remaining 
            shape.size()/2+shape.size()%2 indices and places them 
            horizontally (columns)
            
            EX: A 9x1 vector
            
            ->data_dimensions[0] = 9
            ->data_dimensions[1] = 1
            
            
            EX: A 10x3x3 tensor
            
            ->data_dimensions[0] = 10
            ->data_dimensions[1] = 3*3 = 9
            */
            
            //Check the length of the shape
            if(shape.size()==0){
                //TODO: Should raise an error
            }
            else if(shape.size()==1){//If the shape is a vector
                
                data_dimensions[0] = shape[0];
                data_dimensions[1] = 1;
            }
            else{//If the shape is a tensor of rank >=2
                
                index_split    = shape.size()/2;           //Locate the split index
                iterator_split = shape.begin()+index_split //Construct the iterator
                
                //Multiply each element in the subvectors to obtain the matrix dimensions
                data_dimensions[0] = std::accumulate(shape.begin(),        iterator_split, 1, std::multiplies<double>());
                data_dimensions[1] = std::accumulate(shape.iterator_split, shape.end(),    1, std::multiplies<double>());
                
            }
            
        }
        
        //Set the dimensions of the data matrix and initialize to zero
        data = Eigen::MatrixXd::Zeros(data_dimensions[0],data_dimensions[1]);
    }
    
    
    
}