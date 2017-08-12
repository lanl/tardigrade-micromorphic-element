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
//#include <cstdarg>
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
        
            std::vector< int >                  shape;              //!The shape of the tensor a vector of integers indicating each dimensions length
            std::array< int,2 >                 data_dimensions;    //!The dimensions of the data matrix
            Eigen::MatrixXd                     data;               //!The data object for the tensor a Eigen::MatrixXd
            std::string                         format = "default"; //!The format of the data object. This modifies the way the tensor is stored
                                                                    //!options: (default)
            int                                 index_split;        //!Where, in the indices for the tensor, the split between rows and
                                                                    //!Columns occurs for the storage array
            std::vector< int >::const_iterator  iterator_split;     //!The location of the iterator which marks the split in the array
        
            //Constructors
            Tensor();
            Tensor(std::vector< int >);
            
            //Operators
            /*Operators*/
            Tensor& operator=(const Tensor& T);
            template <typename ...ArgsT>
            double& operator()(ArgsT ...indices){
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
                    indices: The indices provided by the user. Should be of type int
                */
        
                std::array< int , 2>   data_indices  = map_index({indices...});  //! The row and column numbers
        
                return data(data_indices[0],data_indices[1]);
        
            }
    
            template <typename ...ArgsT>
            double operator()(ArgsT ...indices) const{
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
                    indices: The indices provided by the user. Should be of type int
                */
        
                std::array< int, 2>  data_indices  = map_index({indices...});  //! The row and column numbers
        
                return data(data_indices[0],data_indices[1]);
        
            }  

        private:
            std::array<int, 2> map_index(std::initializer_list< int > indices) const; //Map the tensor index to the data matrix
            void set_dimensions();                 //Set the dimensions of the data matrix
    };
}