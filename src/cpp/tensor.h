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
    template<int m_b = Eigen::Dynamic, int n_b = Eigen::Dynamic> class BaseTensor{
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
        
            std::vector< int >                  shape;                           //!The shape of the tensor a vector of integers indicating each dimensions length
            std::array< int,2 >                 data_dimensions = { {m_b,n_b} }; //!The dimensions of the data matrix (note, double braces for gcc compiler)
            Eigen::Matrix<double, m_b, n_b>     data;                            //!The data object for the tensor a Eigen::Matrix
            std::string                         format = "default";              //!The format of the data object. This modifies the way the tensor is stored
                                                                                 //!options: (default)
            int                                 index_split;                     //!Where, in the indices for the tensor, the split between rows and
                                                                                 //!Columns occurs for the storage array. This is the first index 
                                                                                 //!associated with the columns.
            std::vector< int >::const_iterator  iterator_split;                  //!The location of the iterator which marks the split in the array
        
            //!==
            //!|
            //!| Constructors
            //!|
            //!==
            
            BaseTensor(){
                /*!The default constructor for the tensor.*/
                data.setZero();
            }
            
            
            BaseTensor(std::vector< int > _shape){
                /*!A constructor for a tensor where the dimensions of the tensor are read in
                The tensor will be initialized to 0.*/
                
                int temp_rows; //!A temporary variable indicating the expected number of rows 
                               //!Of the matrix given the current value of index_split
                
                data.setZero(); //!Set the data matrix to zero
                //Assign the incoming dimensions to shape
                shape = _shape;
                
                if((data_dimensions[0]==Eigen::Dynamic) || (data_dimensions[1]==Eigen::Dynamic)){
                    set_dimensions();
                }
                else{
                    //Assign the value to index_split which 
                    //Corresponds with the assigned data 
                    //matrix dimensions.
                    index_split = 1;
                    temp_rows   = shape[0];
                    while((temp_rows<data.rows()) && (index_split<shape.size())){
                        temp_rows *= shape[index_split];
                        index_split++;
                    }
                    if(index_split>=shape.size()){
                        std::cout << "Error: value of index_split is outside of allowable range.\n";
                        assert(1==0);
                    }
                }
        
                //Set the multiplication factors
                set_factors();
        
            }
            
            BaseTensor(std::vector< int > _shape, Eigen::MatrixXd data_in){
                /*!A constructor for a tensor where the dimensions of the tensor are read in
                along with the data values.
        
                dimensions and data must be consistent!*/
                
                int temp_rows; //!A temporary variable indicating the expected number of rows 
                               //!Of the matrix given the current value of index_split
                
                data.setZero(); //!Set the data matrix to zero
                //Assign the incoming dimensions to shape
                shape = _shape;
                
                if((data_dimensions[0]==Eigen::Dynamic) || (data_dimensions[1]==Eigen::Dynamic)){
                    set_dimensions();
                }
                else{
                    //Assign the value to index_split which 
                    //Corresponds with the assigned data 
                    //matrix dimensions.
                    index_split = 1;
                    temp_rows   = shape[0];
                    while((temp_rows<data.rows()) && (index_split<shape.size())){
                        temp_rows *= shape[index_split];
                        index_split++;
                    }
                    if(index_split>=shape.size()){
                        std::cout << "Error: value of index_split is outside of allowable range.\n";
                        assert(1==0);
                    }
                }
        
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
            BaseTensor<m_b,n_b>& operator=(const BaseTensor<m_b,n_b>& T){
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
                index_factors   = T.index_factors;
            }
    
            BaseTensor<m_b,n_b> operator+(const BaseTensor<m_b,n_b>& T1){
                /*!================================================
                |               Tensor::operator+               |
                =================================================
        
                Redefine the addition operator
        
                */
        
                if(T1.shape != shape){
                    std::cout << "\nError: tensors must have the same shape\n";
                    std::cout << "first tensor shape: {";
                    for(int i=0; i<shape.size(); i++){std::cout << " " << shape[i];}
                    std::cout << "}\n";
                    std::cout << "second matrix shape: {";
                    for(int i=0; i<T1.shape.size(); i++){std::cout << " " << T1.shape[i];}
                    std::cout << " }\n";
                    assert(1==0); //TODO: allow to raise error
                }
                
                if(T1.format.compare(format)){//Note, this is not strictly accurate but it makes sense to only use tensors with the same storage format
                    std::cout << "\nError: tensors must have the same storage format\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                Eigen::Matrix<double,m_b,n_b> new_data = T1.data + data;
        
                BaseTensor<m_b,n_b> T(T1.shape,new_data);
                return T;
        
            }
            
            //template<int m, int n>
            void operator+=(const BaseTensor<m_b,n_b>& T1){
                /*!================================================
                |               Tensor::operator+=              |
                =================================================
        
                Redefine the addition equals operator
        
                */
        
                if(T1.shape != shape){
                    std::cout << "\nError: tensors must have the same shape\n";
                    std::cout << "base tensor shape: {";
                    for(int i=0; i<shape.size(); i++){std::cout << " " << shape[i];}
                    std::cout << "}\n";
                    std::cout << "added tensor shape: {";
                    for(int i=0; i<T1.shape.size(); i++){std::cout << " " << T1.shape[i];}
                    std::cout << " }\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                if(T1.format.compare(format)){//Note, this is not strictly accurate but it makes sense to only use tensors with the same storage format
                    std::cout << "\nError: tensors must have the same storage format\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                data += T1.data;
        
                return;
        
            }
            
            BaseTensor<m_b,n_b> operator-(const BaseTensor<m_b,n_b>& T1){
                /*!================================================
                |               Tensor::operator-               |
                =================================================
        
                Redefine the subtraction operator
        
                */
        
                if(T1.shape != shape){
                    std::cout << "\nError: tensors must have the same shape\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                if(T1.format.compare(format)){//Note, this is not strictly accurate but it makes sense to only use tensors with the same storage format
                    std::cout << "\nError: tensors must have the same storage format\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                Eigen::Matrix<double,m_b,n_b> new_data = data - T1.data;
        
                BaseTensor<m_b,n_b> T(T1.shape,new_data);
                return T;
        
            }
    
            BaseTensor<m_b,n_b> operator-(){
                /*!================================================
                |               Tensor::operator-               |
                =================================================
        
                Redefine the negative operator
        
                */
        
                Eigen::Matrix<double, m_b, n_b> new_data = -data;
        
                BaseTensor<m_b,n_b> T(shape,new_data);
                return T;
        
            }
    
            void operator-=(const BaseTensor<m_b,n_b> &T1){
                /*!================================================
                |               Tensor::operator-=              |
                =================================================
        
                Redefine the subtraction equals operator
        
                */
        
                if(T1.shape != shape){
                    std::cout << "\nError: tensors must have the same shape\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                if(T1.format.compare(format)){//Note, this is not strictly accurate but it makes sense to only use tensors with the same storage format
                    std::cout << "\nError: tensors must have the same storage format\n";
                    assert(1==0); //TODO: allow to raise error
                }
        
                data -= T1.data;
                return;
        
            }
            
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
            
            BaseTensor<m_b,n_b> inverse(){
                /*!Invert a tensor of even order. 
        
                This method directly inverts the data matrix which returns an inverse
                of the tensor. This method will only work with tensors of even order 
                such that the storage matrix is square.*/
        
                Eigen::Matrix<double,m_b,n_b> inv_data = data.inverse();
                BaseTensor<m_b,n_b> T_out(shape,inv_data);
                return T_out;
            }
            
            double det(){
                /*!Get the determinant of a tensor of even order
        
                This method gets the determinant of the data matrix and returns this 
                as the determinant of the tensor. This method will only work with 
                tensors of even order such that the storage matrix is square.*/
        
                return data.determinant();
            }
            
            void zero(){
                /*! Set all values of the tensor to zero */
        
                for(int i=0; i<data.rows(); i++){
                    for(int j=0; j<data.rows(); j++){
                        data(i,j) = 0.;
                    }
                }
            }

        private:
            std::vector< int > index_factors;                                         //!The factors by which each index is multiplied  
                                                                                      //!by in the mapping from the tensor indices to the 
                                                                                      //!matrix form.
            
            void set_factors(){
                /*!=====================
                |    set_factors    |
                =====================
        
                Set the multiplication factors for each of the indices
        
                */
        
                index_factors    = std::vector< int >(shape.size(),1);
        
                //Get the row index
                for(int i=0; i<index_split; i++){//Iterate through the pre-split indices
        
                    for(int j=i+1; j<index_split; j++){//Assemble the index multiplier
                        index_factors[i] *= shape[j];
                    }
                }
                //Get the column index
                for(int i=index_split; i<shape.size(); i++){//Iterate through the post-split indices
            
                    for(int j=i+1; j<shape.size(); j++){//Assemble the index multipler
                        index_factors[i] *= shape[j];
                    }
                }
            }
            
            std::array< int, 2 > map_index(std::initializer_list< int > indices) const{
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
        
                //Error Handling
                if(indices.size() != shape.size()){//Make sure the number of indices is equal to the tensor order
                    //TODO: Should raise an error!
                    std::cout << "Error: indices size does not equal the order of the tensor.\n";
                    assert(1==0);
                }
        
                inc = 0;
                for(it=indices.begin(); inc<shape.size(); it++){//Make sure the dimension of each index is consistent with the shape
            
                    if((*it>=shape[inc]) || (*it<0)){
                        //TODO: Should raise an error!
                        std::cout << "Error: index " << inc << " with value " << *it << "\noutside allowed bounds.\n";
                        assert(1==0);
                    }
                    inc++;
                }
        
                //Get the row index
                inc = 0;
                for(it=indices.begin(); inc<index_split; it++){//Iterate through the pre-split indices
            
                    data_indices[0] += *it*index_factors[inc];
            
                    inc++;
            
                }
                //Get the column index
                inc = index_split;
                for(it=(indices.begin()+inc); inc<shape.size(); it++){//Iterate through the post-split indices
            
                    data_indices[1] += *it*index_factors[inc];
            
                    inc++;
                }
                
                return data_indices;
            }
    
            void set_dimensions(){
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
                
                        index_split        = 1;
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
                data.resize(data_dimensions[0],data_dimensions[1]);
                data = Eigen::MatrixXd::Zero(data_dimensions[0],data_dimensions[1]);
            }
            
            
    };
    
    //!==
    //!|
    //!| Operators
    //!|
    //!==
    
    template< int m, int n> BaseTensor<m,n> operator*(const double& a, const BaseTensor<m,n>& T){
        /*!================================================
        |               Tensor::operator*               |
        =================================================
        
        Redefine the multiplication operator for a scalar 
        double.
        
        */
        
        Eigen::Matrix<double,m,n> new_data = a*T.data;
        BaseTensor<m,n> Tout(T.shape,new_data);
        
        return Tout;
    }
    
    template< int m, int n> BaseTensor<m,n> operator*(const int& a, const BaseTensor<m,n>& T){
        /*!================================================
        |               Tensor::operator*               |
        =================================================
        
        Redefine the multiplication operator for a scalar 
        integer.
        
        */
        
        Eigen::Matrix<double,m,n> new_data = a*T.data;
        BaseTensor<m,n> Tout(T.shape,new_data);
        
        return Tout;
    }
    
    template< int m, int n> BaseTensor<m,n> operator*(const BaseTensor<m,n>& T, const double& a){
        /*!================================================
        |               Tensor::operator*               |
        =================================================
        
        Redefine the multiplication operator for a scalar 
        double.
        
        */
        return a*T;
    }
    
    template< int m, int n> BaseTensor<m,n> operator*(const BaseTensor<m,n>& T, const int& a){
        /*!================================================
        |               Tensor::operator*               |
        =================================================
        
        Redefine the multiplication operator for a scalar 
        integer.
        
        */
        return a*T;
    }
    
    template< int m, int n> BaseTensor<m,n> operator/(const BaseTensor<m,n>& T,const double& a){
        /*!================================================
        |               Tensor::operator/               |
        =================================================
        
        Redefine the division operator for a scalar 
        integer.
        
        */
        
        return T*(1./a);
    }
    
    template< int m, int n> BaseTensor<m,n> operator/(const BaseTensor<m,n>& T,const int& a){
        /*!================================================
        |               Tensor::operator/               |
        =================================================
        
        Redefine the division operator for a scalar 
        integer.
        
        */
        
        Eigen::Matrix<double,m,n> new_data = T.data/a;
        BaseTensor<m,n> Tout(T.shape,new_data);
        
        return Tout;
    }
    
        
    //!==
    //!|
    //!| Type definitions
    //!|
    //!==
    
    typedef BaseTensor<3,3>   Tensor23;
    typedef BaseTensor<9,3>   Tensor33;
    typedef BaseTensor<9,9>   Tensor43;
    typedef BaseTensor<27,27> Tensor63;
    typedef BaseTensor<Eigen::Dynamic, Eigen::Dynamic> Tensor;
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
        
    Tensor23 eye();
    Tensor43 FOT_eye();
}