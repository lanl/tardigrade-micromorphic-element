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
//#include <cstdarg>
#include <Eigen/Dense>
#include <tensor.h>

namespace tensor{
    
    Tensor::Tensor(){
        /*!The default constructor for the tensor. 
        Zeros the sizes of the shape and data attributes
        to save space*/
        
        shape.resize(0);
        data.resize(0,0);
        data_dimensions[0] = 0;
        data_dimensions[1] = 0;
    }
    
    Tensor::Tensor(std::vector< int > dimensions){
        /*!A constructor for a tensor where the dimensions of the tensor are read in
        The tensor will be initialized to 0.*/
        
        //Assign the incoming dimensions to shape
        shape = dimensions;
        
        //Set the data matrix dimensions
        set_dimensions();
        
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
        if(!format.compare("default")){
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
                assert(1==0);
            }
            else if(shape.size()==1){//If the shape is a vector
                
                std::cout << "In vector definition\n";
                data_dimensions[0] = shape[0];
                data_dimensions[1] = 1;
            }
            else{//If the shape is a tensor of rank >=2
                
                index_split    = shape.size()/2;            //Locate the split index
                iterator_split = shape.begin()+index_split; //Construct the iterator
                
                //Multiply each element in the subvectors to obtain the matrix dimensions
                data_dimensions[0] = 1;
                for(int i=0; i<index_split; i++){
                    data_dimensions[0] *= shape[i];
                }
                
                data_dimensions[1] = 1;
                for(int i=index_split; i<shape.size(); i++){
                    data_dimensions[1] *= shape[i];
                }
            }
            
        }
        else{
            std::cout << "Error: format not recognized\n";
            std::cout << "       format: " << format << "\n";
            assert(1==0);
        }
        
        //Set the dimensions of the data matrix and initialize to zero
        std::cout << "\ndata_dimensions[0]: " << data_dimensions[0] << "\ndata_dimensions[1]: " << data_dimensions[1] << "\n";
        data.resize(data_dimensions[0],data_dimensions[1]);
        data = Eigen::MatrixXd::Zero(data_dimensions[0],data_dimensions[1]);
    }
    
    Tensor& Tensor::operator=(const Tensor& T){
        /*!================================================
        |               Tensor::operator=               |
        =================================================
        
        Redefine the copy operator
        
        */
        
        shape           = T.shape;
        data_dimensions = T.data_dimensions;
        data            = T.data;
        format          = T.format;
        index_split     = T.index_split;
        iterator_split  = T.iterator_split;
    }
    
    std::array< int, 2 > Tensor::map_index(std::initializer_list< int > indices) const{
        /*!===================================
        |            map_index             |
        ====================================
        
        Map the indices of the tensor to the 
        row and column indices of the matrix
        used for data storage
        
        */
        
        std::array< int ,2> data_indices = {0,0}; //!The mapped row and column indices corresponding
                                                  //!to the tensor indices
        int val;                                  //!Temporary variable
        std::initializer_list<int>::iterator it;  //!Iterator variable
        int inc;                                  //!The incrementation variable
        
        std::cout << "indices.size(): " << indices.size() << "\n" << "shape.size(): " << shape.size() << "\n";
        
        //Error Handling
        if(indices.size() != shape.size()){
            //TODO: Should raise an error!
            assert(1==0);
        }
        
        //Get the row index
        inc = 0;
        for(it=indices.begin(); inc<=index_split; it++){//Iterate through the indices before the split
            val = 1;
            
            for(int j=(it-indices.begin()+1); j<index_split; j++){//Increment through the remaining indices before the split
                val *= shape[j];
            }
            
            data_indices[0] += *it*val; //Add the remaining values to row
            
            inc++;                      //Increment the increment variable
        }
        //Get the column index
        inc = index_split;
        for(it=indices.begin(); inc<shape.size(); it++){//Iterate through the indices before the split
            val = 1;
            
            for(int j=(it-indices.begin()+1); j<index_split; j++){//Increment through the remaining indices before the split
                val *= shape[j];
            }
            
            data_indices[1] += *it*val; //Add the remaining values to row
            
            inc++;                      //Increment the increment variable
        }
        /*for(int i=0; i<index_split; i++){//Increment through the indices before the split
            val = 1; //Set the initial value of val
            for(int j=(i+1); j<index_split; j++){//Increment through the remaining indices before the split
                val *= shape[j];
            }
            data_indices[0] += indices[i]*val;//Add the remaining values to row
        }
        
        //Get the column index
        for(int i=index_split; i<shape.size(); i++){//Increment through the indices after the split
            val = 1; //Set the initial value of val
            for(int j=(i+1); j<shape.size(); j++){//Increment through the remaining indices before the split
                val *= shape[j];
            }
            data_indices[1] += indices[i]*val;//Add the remaining values to row
        }
        
        */
        
        std::cout << "data_indices[0]: " << data_indices[0] << "\n";
        std::cout << "data_indices[1]: " << data_indices[1] << "\n";
        
        return data_indices;
    }
    
}