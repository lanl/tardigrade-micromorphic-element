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
    //!| Constructors
    //!|
    //!==
    
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
    
    Tensor::Tensor(std::vector< int > dimensions, Eigen::MatrixXd data_in){
        /*!A constructor for a tensor where the dimensions of the tensor are read in
        along with the data values.
        
        dimensions and data must be consistent!*/
        
        //Assign the incoming dimensions to shape
        shape = dimensions;
        
        //Set the data matrix dimensions
        set_dimensions();
        
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
    
    Tensor Tensor::operator+(const Tensor& T1){
        /*!================================================
        |               Tensor::operator+               |
        =================================================
        
        Redefine the addition operator
        
        */
        
        if(T1.shape != shape){
            std::cout << "\nError: tensors must have the same shape\n";
            assert(1==0); //TODO: allow to raise error
        }
                
        if(T1.format.compare(format)){//Note, this is not strictly accurate but it makes sense to only use tensors with the same storage format
            std::cout << "\nError: tensors must have the same storage format\n";
            assert(1==0); //TODO: allow to raise error
        }
        
        Eigen::MatrixXd new_data = T1.data + data;
        
        tensor::Tensor T = tensor::Tensor(T1.shape,new_data);
        return T;
        
    }
    
    void Tensor::operator+=(const Tensor& T1){
        /*!================================================
        |               Tensor::operator+=              |
        =================================================
        
        Redefine the addition equals operator
        
        */
        
        if(T1.shape != shape){
            std::cout << "\nError: tensors must have the same shape\n";
            assert(1==0); //TODO: allow to raise error
        }
        
        if(T1.format.compare(format)){//Note, this is not strictly accurate but it makes sense to only use tensors with the same storage format
            std::cout << "\nError: tensors must have the same storage format\n";
            assert(1==0); //TODO: allow to raise error
        }
        
        data += T1.data;
        
        return;
        
    }
    
    Tensor Tensor::operator-(const Tensor& T1){
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
        
        Eigen::MatrixXd new_data = data - T1.data;
        
        tensor::Tensor T = tensor::Tensor(T1.shape,new_data);
        return T;
        
    }
    
    Tensor Tensor::operator-(){
        /*!================================================
        |               Tensor::operator-               |
        =================================================
        
        Redefine the negative operator
        
        */
        
        Eigen::MatrixXd new_data = -data;
        
        tensor::Tensor T = tensor::Tensor(shape,new_data);
        return T;
        
    }
    
    void Tensor::operator-=(const Tensor& T1){
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
    
    //!==
    //!|
    //!| Methods
    //!|
    //!==
    
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
        
        //Error Handling
        if(indices.size() != shape.size()){//Make sure the number of indices is equal to the tensor order
            //TODO: Should raise an error!
            std::cout << "Error: indices size does not equal the order of the tensor.";
            assert(1==0);
        }
        
        inc = 0;
        for(it=indices.begin(); inc<shape.size(); it++){//Make sure the dimension of each index is consistent with the shape
            
            if((*it>=shape[inc]) || (*it<0)){
                //TODO: Should raise an error!
                std::cout << "Error: index " << inc << "with value " << *it << "\noutside allowed bounds.\n";
                assert(1==0);
            }
            inc++;
        }
        
        //Get the row index
        inc = 0;
        for(it=indices.begin(); inc<index_split; it++){//Iterate through the pre-split indices
            val = 1;
            
            for(int j=inc+1; j<index_split; j++){//Assemble the index multiplier
                val *= shape[j];
            }
            
            data_indices[0] += *it*val;
            
            inc++;
            
        }
        //Get the column index
        inc = index_split;
        for(it=(indices.begin()+inc); inc<shape.size(); it++){//Iterate through the post-split indices
            val = 1;
            
            for(int j=inc+1; j<shape.size(); j++){//Assemble the index multipler
                val *= shape[j];
            }
            
            data_indices[1] += *it*val;
            
            inc++;
        }
                
        return data_indices;
    }
    
    Tensor Tensor::inverse(){
        /*!Invert a tensor of even order. 
        
        This method directly inverts the data matrix which returns an inverse
        of the tensor. This method will only work with tensors of even order 
        such that the storage matrix is square.*/
        
        Eigen::MatrixXd inv_data = data.inverse();
        Tensor T_out = Tensor(shape,inv_data);
        return T_out;
    }
    
    double Tensor::det(){
        /*!Get the determinant of a tensor of even order
        
        This method gets the determinant of the data matrix and returns this 
        as the determinant of the tensor. This method will only work with 
        tensors of even order such that the storage matrix is square.*/
        
        return data.determinant();
    }
    
    void Tensor::zero(){
        /*! Set all values of the tensor to zero */
        
        for(int i=0; i<data.rows(); i++){
            for(int j=0; j<data.rows(); j++){
                data(i,j) = 0.;
            }
        }
    }
    
    //!==
    //!|
    //!| Functions
    //!|
    //!==
    
    Tensor eye(){
        /*!Return the second order identity tensor
        
        Populate the second order identity tensor.
        This function is useful so that a consistent identity tensor can be used 
        with different storage schemes.*/
        
        //Initialize the matrix shape
        std::vector< int > shape;
        shape.resize(2);
        shape[0] = 3;
        shape[1] = 3;
        
        //Initialize the matrix
        Tensor I = Tensor(shape);
        
        for(int i=0; i<3; i++){
            I(i,i) = 1.;
        }
        return I;
    }
    
    Tensor FOT_eye(){
        /*!Return the fourth order identity tensor
        
        Populate the fourth order identity tensor.
        This function is useful so that a consistent identity tensor can be used
        with different storage schemes.*/
        
        //Initialize the matrix shape
        std::vector< int > shape;
        shape.resize(4);
        shape[0] = 3;
        shape[1] = 3;
        shape[2] = 3;
        shape[3] = 3;
        
        //Initialize the matrix
        Tensor FOTI = Tensor(shape);
        Tensor I    = eye(); //Get the second order identity tensor
        
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