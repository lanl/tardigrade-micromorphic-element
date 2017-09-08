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
#include <Eigen/Dense>
  
namespace tensor
{
    template<int m_b, int n_b> class BaseTensor{
        /*!===
           |
           | B a s e T e n s o r
           |
          ===
        
        The base class for the tensor implementation
        
        Defines and allows access to a tensor using indicial 
        notation but the data is stored as an array.
        
        The storage format of Tensor.data is controlled by 
        the method Tensor.format(string) which takes in a 
        string that indicates how the information should be 
        stored
        
        The class is a template that requires the number of rows and 
        columns to be input. The default values are for a dynamic 
        matrix size which is the most general case. For better speed 
        on smaller tensors (>16 terms) the number of terms in the 
        matrix should be specified unless the size of the tensor will 
        be changing over time.
        
        */
        
        public:
        
            //!==
            //!|
            //!| Attribute Definitions
            //!|
            //!==
        
            std::vector< int >                  shape;                          //!The shape of the tensor a vector of integers indicating each dimensions length
            std::array< int,2 >                 data_dimensions = {m_b,n_b};    //!The dimensions of the data matrix
            Eigen::Matrix<double, m_b, n_b>     data;                           //!The data object for the tensor a Eigen::Matrix
            std::string                         format = "default";             //!The format of the data object. This modifies the way the tensor is stored
                                                                                //!options: (default)
            int                                 index_split;                    //!Where, in the indices for the tensor, the split between rows and
                                                                                //!Columns occurs for the storage array. This is the first index 
                                                                                //!associated with the columns.
            std::vector< int >::const_iterator  iterator_split;                 //!The location of the iterator which marks the split in the array
        
            //!==
            //!|
            //!| Constructors
            //!|
            //!==
            
            BaseTensor(){
                /*!The default constructor for the tensor.*/
            }
            
            
            BaseTensor(std::vector< int > dimensions){
                /*!A constructor for a tensor where the dimensions of the tensor are read in
                The tensor will be initialized to 0.*/
        
                //Assign the incoming dimensions to shape
                shape = dimensions;
        
                //Set the multiplication factors
                set_factors();
        
            }
            
            BaseTensor(std::vector< int > dimensions, Eigen::MatrixXd data_in){
                /*!A constructor for a tensor where the dimensions of the tensor are read in
                along with the data values.
        
                dimensions and data must be consistent!*/
        
                //Assign the incoming dimensions to shape
                shape = dimensions;
        
                //Set the multiplication factors
                set_factors();
        
                if((data.rows()==data_in.rows()) && (data.cols()==data_in.cols())){
                    data = data_in;
                }
                else{
                    std::cout << "\nError: Tensor dimensions and the provided data matrix are\n";
                    std::cout << "         not consistent.";
                    assert(1==0);
                }
            }

            //!==
            //!|
            //!| Operators
            //!|
            //!==
            template< int m, int n> class BaseTensor& operator=(const BaseTensor<m,n>& T);
            template< int m, int n> class BaseTensor operator+(const BaseTensor<m,n>& T1);
            void operator+=(const BaseTensor& T2);
            
            template< int m, int n> BaseTensor operator-(const BaseTensor<m,n>& T1);
            template< int m, int n> BaseTensor operator-();
            void operator-=(const BaseTensor &T1);
            
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
            
            //!==
            //!|
            //!| Methods
            //!|
            //!==
            
            template< int m, int n> BaseTensor inverse();
            double det();
            void zero();

        private:
            std::vector< int > index_factors;                                         //!The factors by which each index is multiplied  
                                                                                      //!by in the mapping from the tensor indices to the 
                                                                                      //!matrix form.
            void set_factors();                                                       //!Set the multiplication factors for the tensor 
                                                                                      //!index to matrix form
            std::array<int, 2> map_index(std::initializer_list< int > indices) const; //Map the tensor index to the data matrix
            void set_dimensions();                                                    //Set the dimensions of the data matrix
    };
    
    //!==
    //!|
    //!| Operators
    //!|
    //!==
    
    template< int m, int n> BaseTensor<m,n> operator*(const double&,const BaseTensor<m,n>&);
    template< int m, int n> BaseTensor<m,n> operator*(const int&,   const BaseTensor<m,n>&);
    
    template< int m, int n> BaseTensor<m,n> operator*(const BaseTensor<m,n>&,const double&);
    template< int m, int n> BaseTensor<m,n> operator*(const BaseTensor<m,n>&,const int&);
    
    template< int m, int n> BaseTensor<m,n> operator/(const BaseTensor<m,n>&,const double&);
    template< int m, int n> BaseTensor<m,n> operator/(const BaseTensor<m,n>&,const int&);
    
        
    //!==
    //!|
    //!| Type definitions
    //!|
    //!==
    
    typedef BaseTensor<3,3> Tensor23;
    typedef BaseTensor<9,9> Tensor43;
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
        
    Tensor23 eye();     //!Return the second order identity tensor in three dimensions
    Tensor43 FOT_eye(); //!Return the fourth order identity tensor in four dimensions
}