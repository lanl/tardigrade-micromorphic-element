/*!=======================================================
  |                                                     |
  |                 micro_element.cpp                   |
  |                                                     |
  -------------------------------------------------------
  | The source file for the definition of a             |
  | micromorphic continuum element.                     |
  =======================================================
  | Dependencies:                                       |
  | tensor:      The class which defines tensor access  |
  |              to an underlying Eigen matrix. This    |
  |              may result in a somewhat slower result |
  |              however it should allow for a faster   |
  |              implementation.                        |
  =======================================================*/
  
#include <iostream>
#include <vector>
#include <tensor.h>
  
namespace micro_element
{
    
    /*!==
    |
    | Constructors
    |
    ==*/
    
    Hex8::Hex8(){
        
    }
    
    Hex8::Hex8(std::vector< std::vector< double > > rcs){
        /*!====================
        |        Hex8       |
        =====================
        
        The constructor for a hexehedral element when 
        given the nodal reference coordinates.
        
        Input:
           rcs: A vector of vectors of doubles which are the
                coordinates of the nodes.
                
                The nodes are ordered in a counter clockwise 
                manner i.e.
                
               4,8      3,7
                o--------o
                |        |
                |        |
                |        |
                o--------o
               1,5      2,6
               
               where the comma indicates the ``upper layer'' of the
               hexehedral element.
        */
        
        reference_coords = rcs;
        current_coords   = rcs;
    }
    
    //!==
    //!|
    //!| Operators
    //!|
    //!==
    
    Hex8::operator=(const Hex8& hex8_in){
        /*========================
        |      operator=       |
        ========================
        
        Copy operator to allow for copying 
        hexahedral elements*/
        
        reference_coords = hex8_in.reference_coords;
        current_coords   = hex8_in.current_coords;
    }
    
    //!==
    //!|
    //!| Methods
    //!|
    //!==
    
    //!=
    //!| Shape Functions
    //!=
    
    Hex8::shape_function(int n, const std::vector< double > &xi){
        /*==========================
        |      shape_function    |
        ==========================
        
        Compute the value of the shape function
        at a given node and xi (local coordinate) 
        locations
        
        Input:
            n  : Node number
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        return 0.125*(1+xi[0]*local_coords[n][0])*(1+xi[1]*local_coords[n][1])*(1+xi[2]*local_coords[n][2])
    }
    
    std::vector< double > Hex8::local_gradient_shape_function(int n, const std::vector< double > &xi){
        /*========================================
        |    local_gradient_shape_function    |
        =======================================
        
        Compute the value of the gradient of the 
        shape function at a given node and xi 
        (local coordinate) locations locally.
        
        Input:
            n  : Node number
            xi : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        //Initialize the shape of the output tensor
        std::vector< double > dNdxi;
        dNdxi.resize(3);
        
        dNdxi[0] = 0.125*local_coords[n][0]*(1+xi[1]*local_coords[n][1])*(1+xi[2]*local_coords[n][2]);
        dNdxi[1] = 0.125*(1+xi[0]*local_coords[n][0])*local_coords[n][1]*(1+xi[2]*local_coords[n][2]);
        dNdxi[2] = 0.125*(1+xi[0]*local_coords[n][0])*(1+xi[1]*local_coords[n][1])*local_coords[n][2];
        
        return dNdxi
    }
    
    std::vector< double > Hex8::global_gradient_shape_function(bool mode, int n, const std::vector< double > &xi){
        /*========================================
        |    global_gradient_shape_function    |
        ========================================
        
        Compute the value of the gradient of the 
        shape function at a given node and xi 
        (local coordinate) locations globally.
        
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
            n   : Node number
            xi  : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        //Initialize the output
        std::vector< double > dNdx;
        dNdx.resize(3);
        
        tensor::Tensor J     = compute_jacobian(mode,xi);              //Compute the jacobian
        tensor::Tensor Jinv  = J.inverse();                            //Invert the jacobian
        double         Jdet  = J.det();                                //Get the determinant of the jacobian
        tensor::Tensor dNdxi = local_gradient_of_shape_function(n,xi); //Compute the local gradient of the shape function
        
        //Compute the global gradient w.r.t. either the reference or global x
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                dNdx[i] += dNdxi[j]*Jinv(j,i);
            }
        }                       
        return dNdx;
    }
    
    tensor::Tensor& compute_jacobian(bool mode, const std::vector< double > xi){
        /*==========================
        |    compute_jacobian    |
        ==========================
        
        Compute the value of the jacobian of 
        transformation between the given 
        coordinates and the local coordinates.
        Input:
            mode: Selection between the gradient w.r.t.
                  the reference coordinates (0) or the 
                  current coordinates (1)
            xi  : xi vector xi = [xi_1, xi_2, xi_3]
        
        */
        
        //Initialize the output
        std::vector< double > shape = {3,3};
        tensor::Tensor J = tensor::Tensor.zeros(shape);
        
        //Initialize local gradient value
        std::vector< double > local_gradient;
        local_gradient.resize(3);

        //Set up a temporary renaming variable depending on the mode
        std::vector< std::vector< double > > coordinates;
        
        if(mode){
            coordinates = current_coords;
        }
        else{
            coordinates = reference_coords;
        }
        
        //Add the contributions of each node in the element
        for(int n = 0; n<reference_coords.size; n++){
            local_gradient = local_gradient_shape_function(n,xi); //Get the local gradient of the given nodes shape function
            J += vector_dyadic_product(local_gradient,coordinates[n][j]); //Compute the vector dyadic product and add it to the jacobian
        }
        
        return J
    }
    
    vector_dyadic_product(const std::vector< double > & V1, const std::vector< double >& V2){
        /*================================
        |    vector_dyadic_product    |
        ===============================
        
        Compute the dyadic product of two vectors.
        
        i.e. A_{ij} = V1_i V2_j
        
        Input:
            V1: Vector 1
            V2: Vector 2
        
        */
        
        //Compute the dimension of the resulting tensor
        int rows = V1.rows();
        int cols = V2.cols();
        
        //Initialize the tensor data and shape variables
        Eigen::MatrixXd data = Eigen::MatrixXd::Zero(rows,cols);
        std::vector< int > shape = {rows, cols};
        
        //Compute the dyadic product of the two vectors
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                data(i,j) = V1[i]*V2[j];
            }
        }
        
        tensor::Tensor T = tensor::Tensor(shape,data);
        return T
    }
}