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
        /*!===========================
        |     set_dimensions       |
        ============================
        
        Description:
        Set the dimensions of the data storage of the tensor.
        
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
    
    
    
    double& Tensor::operator()(const unsigned int &num, ...){
        /*!================================================
        |               Tensor::operator()              |
        =================================================
        
        Description:
        Operator which returns the *address* of the underlying data matrix 
        using index notation. Allows the user to set the value of the 
        data matrix.
        
        The definition of the indexing operator that allows index notation
        access to the underlying data matrix. This allows the user to use a 
        more natural index notation. It is hoped that the somewhat higher 
        computational cost will be offset by the ease of understanding and 
        developing the code.
        
        Input:
            num: The number of arguments provided
        */
        
        
        //Variable definitions
        std::va_list indices;                     //! The incoming indices
        std::va_start (indices, num);             //! Store all of the arguments, 
                                                  //! other than num, into indices
        
        unsigned int[2]   data_indices  = {0,0};  //! The row and column numbers
        
        //Error Handling
        if(num != shape.length){
            //TODO: Should raise an error!
        }
        
        data_indices = map_index(indices);
        
        return data(data_indices[0],data_indices[1])
        
    }
    
    double Tensor::operator()(const unsigned int &num, ...) const{
        /*!================================================
        |               Tensor::operator()              |
        =================================================
        
        Description:
        Operator which returns the *value* of the underlying data matrix 
        using index notation. Allows the user to access the value of the
        data matrix.
        
        The definition of the indexing operator that allows index notation
        access to the underlying data matrix. This allows the user to use a 
        more natural index notation. It is hoped that the somewhat higher 
        computational cost will be offset by the ease of understanding and 
        developing the code.
        
        Input:
            num: The number of arguments provided
        */
        
        
        //Variable definitions
        std::va_list indices;                     //! The incoming indices
        std::va_start (indices, num);             //! Store all of the arguments, 
                                                  //! other than num, into indices
        
        unsigned int[2]   data_indices  = {0,0};  //! The row and column numbers
        
        //Error Handling
        if(num != shape.length){
            //TODO: Should raise an error!
        }
        
        data_indices = map_index(indices);
        
        return data(data_indices[0],data_indices[1])
        
    }
    
    unsigned int[2] map_index(std::va_list indices){
        /*!===================================
        |            map_index             |
        ====================================
        
        Map the indices of the tensor to the 
        row and column indices of the matrix
        used for data storage
        
        */
        
        unsigned int[2] data_indices = {0,0}; //!The mapped row and column indices corresponding
                                              //!to the tensor indices
        unsigned int val; //!Temporary variable
        
        //Get the row index
        for(int i=0; i<index_split; i++){//Increment through the indices before the split
            val = 1; //Set the initial value of val
            for(int j=(i+1); j<index_split; j++){//Increment through the remaining indices before the split
                val *= shape[j];
            }
            data_indices[0] += i*val;//Add the remaining values to row
        }
        
        //Get the column index
        for(int i=index_split; i<shape.size(); i++){//Increment through the indices after the split
            val = 1; //Set the initial value of val
            for(int j=(i+1); j<shape.size(); j++){//Increment through the remaining indices before the split
                val *= shape[j];
            }
            data_indices[1] += i*val;//Add the remaining values to row
        }
        
        return data_indices;
    }
    
}