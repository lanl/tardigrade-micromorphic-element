/*!=======================================================
  |                                                     |
  |                     tensor.h                        |
  |                                                     |
  -------------------------------------------------------
  | The header file for a class definition that allows  |
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
  
#include <iostream>
#include <cstdarg>
#include <Eigen/Dense>
  
namespace tensor
{
    class Tensor{
        /*!===
           |
           | T e n s o r
           |
          ===
        
        The base class for the tensor implementation
        
        Defines and allows access to a tensor using indicial 
        notation but the data is stored as an array.
        
        The storage format of Tensor.data is controlled by 
        the method Tensor.format(string) which takes in a 
        string that indicates how the information should be 
        stored
        
        */
        
        public:
        
            //!==
            //!|
            //!| Attribute Definitions
            //!|
            //!==
        
            std::vector< unsigned int >     shape;              //!The shape of the tensor a vector of integers indicating each dimensions length
            unsigned int[2]                 data_dimensions;    //!The dimensions of the data matrix
            Eigen::MatrixXd                 data;               //!The data object for the tensor a Eigen::MatrixXd
            std::string                     format = "default"; //!The format of the data object. This modifies the way the tensor is stored
                                                                //!options: (default)
            unsigned int                    index_split;        //!Where, in the indices for the tensor, the split between rows and
                                                                //!Columns occurs for the storage array
            std::vector<T>::const_iterator  iterator_split;     //!The location of the iterator which marks the split in the array
        
            //Constructors
            void Tensor();
            void Tensor(std::vector< double >);

        private:
            void map_index(int i, ...); //Map the tensor index to the data matrix
            void set_dimensions();      //Set the dimensions of the data matrix
    }
}