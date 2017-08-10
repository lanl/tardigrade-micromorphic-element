/*=======================================================
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
  
namespace tensor
{
    class Tensor{
        /*===
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
        
            //==
            //|
            //| Parameter defintions
            //|
            //==
        
            std::vector< int > shape;
        
            //Constructors
            Tensor(){
                
            }
    }
}