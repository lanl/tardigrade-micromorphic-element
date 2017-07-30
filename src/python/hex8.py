import numpy as np
import unittest
import os

"""Definition of a Hex8 element"""

def Hex8_shape_function_loc(xi_vec,node_local_coords): #Test function written
    """Compute the value of the shape function at xi_vec
    for a node with local coordinates node_local_coords"""
    xi1,xi2,xi3 = xi_vec
    nc1,nc2,nc3 = node_local_coords
    return 0.125*(1+xi1*nc1)*(1+xi2*nc2)*(1+xi3*nc3)
    
def Hex8_shape_function(n,xi_vec): #Test function written
    """Compute the value of the shape function from node n at local position xi_vec"""
    xi1,xi2,xi3 = xi_vec
    nc1,nc2,nc3 = Hex8_node_coords(n)
    return 0.125*(1+xi1*nc1)*(1+xi2*nc2)*(1+xi3*nc3)
    
def Hex8_local_grad_shape_function(n,xi_vec): #Test function written
    """Compute the gradient of the shape function w.r.t. the local coordinates"""
    xi1,xi2,xi3 = xi_vec
    nc1,nc2,nc3 = Hex8_node_coords(n)
    return 0.125*nc1*(1+xi2*nc2)*(1+xi3*nc3),\
           0.125*nc2*(1+xi1*nc1)*(1+xi3*nc3),\
           0.125*nc3*(1+xi1*nc1)*(1+xi2*nc2)
    
def Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords): #Test function written
    """Compute the gradient of the shape function w.r.t. the global coordinates
    along with detJ
    n:                   node number
    nodal_global_coords: global coordinates of the nodes
    xi_vec:              local coordinates of interest
    """
    J                 = get_jacobian(nodal_global_coords,xi_vec)
    detJ,Jinv         = invert_3x3_matrix(J)
    local_grad        = Hex8_local_grad_shape_function(n,xi_vec)
    return np.dot(Jinv,local_grad),detJ
    
def Hex8_get_shape_function_info(n,xi_vec,nodal_global_coords):
    """Compute and return the shape function value, the global reference gradient of the shape function"""
    N           = Hex8_shape_function(n,xi_vec)
    grad_N,detJ = Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords)
    return N,grad_N,detJ
    
def Hex8_node_coords(n): #Test function written
    """Get the local coordinates for node n (0-7)"""
    coords = [[-1,-1,-1],[ 1,-1,-1],[ 1, 1,-1],[-1, 1,-1],\
              [-1,-1, 1],[ 1,-1, 1],[ 1, 1, 1],[-1, 1, 1]]
    return coords[n]
    
def Hex8_interpolation_matrix(nodal_global_coords,nodal_dof_values): #Test function written
    """Form the interpolation matrix for the Hex8 element
    nodal_global_coords: A list (8 values) of coordinates (x,y,z)
    nodal_dof_values: A list (n values) of degrees of freedom (u1,u2,u3,...)"""
    
    column_vectors = [np.concatenate(([1],nodal_global_coords[n],nodal_dof_values[n])) for n in range(8)]
    column_vectors = [np.reshape(vector,[len(vector),1]) for vector in column_vectors]
    return np.hstack(column_vectors)
    
def get_jacobian(nodal_global_coords,xi_vec): #Test function written
    """Compute the jacobian of a given node at a location
    nodal_global_coords: Global coordinates of the nodes
    xi_vec :             local coordinates of point of interest"""
    
    J = np.zeros([3,3])
    for n in range(8):
        node_local_grad  = Hex8_local_grad_shape_function(n,xi_vec)
        J += vector_dyadic_product(node_local_grad,nodal_global_coords[n])
    return J
    
def invert_3x3_matrix_V(V):
    """Invert a 3x3 matrix where the matrix is
    given in vector form"""
    A = convert_V_to_T(V,[3,3])
    deta,Ainv = invert_3x3_matrix(A)
    return deta,reduce_tensor_to_vector_form(Ainv)
    
def invert_3x3_matrix(a): #Test function written
    """Invert a 3x3 matrix"""
    A      = np.empty([3,3])    
    A[0,0] = a[1,1]*a[2,2] - a[1,2]*a[2,1]
    A[1,1] = a[0,0]*a[2,2] - a[0,2]*a[2,0]
    A[2,2] = a[0,0]*a[1,1] - a[0,1]*a[1,0]
    A[1,0] = a[1,2]*a[2,0] - a[1,0]*a[2,2]
    A[0,1] = a[0,2]*a[2,1] - a[0,1]*a[2,2]
    A[2,0] = a[1,0]*a[2,1] - a[1,1]*a[2,0]
    A[0,2] = a[0,1]*a[1,2] - a[0,2]*a[1,1]
    A[2,1] = a[0,1]*a[2,0] - a[0,0]*a[2,1]
    A[1,2] = a[0,2]*a[1,0] - a[0,0]*a[1,2]
    
    deta = a[0,0]*A[0,0]+a[0,1]*A[1,0]+a[0,2]*A[2,0]
    return deta,A/deta
    
def compute_BVoigt(xi_vec,nodal_global_coords): #Test function written
    """Compute the strain displacement matrix
    xi_vec:              local location
    nodal_global_coords: The global coordinates of the element nodes"""
    
    Bsubs = [BVoigt_submatrix(n,xi_vec,nodal_global_coords) for n in range(8)]
    return np.hstack(Bsubs)
    
def BVoigt_submatrix(n,xi_vec,nodal_global_coords): #Test function written
    """Compute a submatrix associated with node n for the strain displacement matrix
    n:                   node number
    nodal_global_coords: Global coordinates of the nodes
    xi_vec:              local coordinates of interest"""
    
    dNdx,detJ = Hex8_global_grad_shape_function(n,xi_vec,nodal_global_coords)
    
    Bsub = np.zeros([6,3])
    Bsub[0,0] = dNdx[0]
    Bsub[1,1] = dNdx[1]
    Bsub[2,2] = dNdx[2]
    Bsub[3,0] = dNdx[1]
    Bsub[3,1] = dNdx[0]
    Bsub[4,1] = dNdx[2]
    Bsub[4,2] = dNdx[1]
    Bsub[5,0] = dNdx[2]
    Bsub[5,2] = dNdx[0]
    return Bsub
    
def get_gpw(order): #Test function written
    """Get the gauss points and weights"""
    p,w = get_1D_gpw(order)
    
    P = []
    W = []
    
    gpw = []
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                P.append([p[i],p[j],p[k]])
                W.append(w[i]*w[j]*w[k])
    return P,W
    
def get_1D_gpw(order): #Test function written
    """Get the gauss points and weights in 1D"""
    if(order==1):
        return [0.],[2.]
    if(order==2):
        return [-np.sqrt(1./3), np.sqrt(1./3)],[1.,1.]
    if(order==3):
        return [-np.sqrt(3./5), 0., np.sqrt(3./5)],[5./9, 8./9, 5./9]
    
def vector_dyadic_product(a,b): #Test function written
    """Compute the vector dyadic product of two vectors"""
    A = np.zeros([len(a),len(b)])
    
    for i in range(len(a)):
        for j in range(len(b)):
            A[i,j] = a[i]*b[j]
    return A
    
def matrix_dot(A,B): #Test function written
    """Compute the dot product of two matrices"""
    C = np.zeros([3,3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i,j] += A[i,k]*B[k,j]
    return C
    
def matrix_dot_V(A,B): #Test function written
    """Compute the dot product of two matrices in vector form"""
    C = np.zeros([9,])
    for i in range(3):
        for j in range(3):
            Cindx = T_to_V_mapping([i,j],[3,3])
            
            for k in range(3):
                Aindx = T_to_V_mapping([i,k],[3,3])
                Bindx = T_to_V_mapping([k,j],[3,3])
                C[Cindx] += A[Aindx]*B[Bindx]
    return C
    
def matrix_Tdot(A,B): #Test function written
    """Compute the dot product of a transposed matrix and a matrix"""
    C = np.zeros([3,3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i,j] += A[k,i]*B[k,j]
    return C
    
def matrix_Tdot_V(A,B): #Test function written
    """Compute the dot product of a transposed matrix and a matrix in vector form"""
    C = np.zeros([9,])
    for i in range(3):
        for j in range(3):
            Cindx = T_to_V_mapping([i,j],[3,3])
            
            for k in range(3):
                Aindx = T_to_V_mapping([k,i],[3,3])
                Bindx = T_to_V_mapping([k,j],[3,3])
                C[Cindx] += A[Aindx]*B[Bindx]
    return C
    
def matrix_Tdot_TOT(A,B): #Test function written
    """Compute the dot product of a transposed matrix and a third order matrix"""
    C = np.zeros([3,3,3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    C[i,j,k] += A[l,i]*B[l,j,k]
    return C
    
def get_symm_matrix(A): #Test function written
    """Compute the symmetric part of A"""
    Asymm = np.zeros([3,3])
    for i in range(3):
        for j in range(3):
            Asymm[i,j] = 0.5*(A[i,j]+A[j,i])
    return Asymm
    
def get_symm_matrix_V(V):
    """Get the symmetric part of a matrix in vector form"""
    A = convert_V_to_T(V,[3,3])
    Asymm = get_symm_matrix(A)
    return reduce_tensor_to_vector_form(Asymm)
    
def vector_dot_matrix(v,A): #Test function written
    """Compute the vector resulting from a vector dotted with a matrix"""
    w = np.zeros([3,])
    
    for i in range(3):
        for j in range(3):
            w[j] += v[i]*A[i,j]
    return w
    
def reduce_tensor_to_matrix_form(A): #Test function written
    """Reduce a tensor to matrix form"""
    
    #Get the shape of the tensor
    shape = A.shape
    
    #Put the tensor in vector form
    V = reduce_tensor_to_vector_form(A)
    
    #Convert vector to matrix
    return convert_V_to_M(V,shape)
    
def reduce_tensor_to_vector_form(A): #Test function written (included in reduce_tensor_to_matrix_form)
    """Reduce a tensor to its vector form"""
    return np.reshape(A,[A.size])
    
def T_to_V_mapping(indices,shape): #Test function written
    """Map the indices to the vector index"""
    
    if(len(indices)!=len(shape)):
        print "\nError: The number of indices should equal the\n"+\
              "       dimension of the tensor.\n"
        raise ValueError()
        
    #Specific Cases
    if(len(shape)==1):
        return indices[0]
    if(len(shape)==2):
        return   indices[0]*shape[1]\
               + indices[1]
    if(len(shape)==3):
        return   indices[0]*shape[1]*shape[2]\
               + indices[1]*shape[2]\
               + indices[2]
    if(len(shape)==4):
        return   indices[0]*shape[1]*shape[2]*shape[3]\
               + indices[1]*shape[2]*shape[3]\
               + indices[2]*shape[3]\
               + indices[3]
    if(len(shape)==5):
        return   indices[0]*shape[1]*shape[2]*shape[3]*shape[4]\
               + indices[1]*shape[2]*shape[3]*shape[4]\
               + indices[2]*shape[3]*shape[4]\
               + indices[3]*shape[4]\
               + indices[4]
    if(len(shape)==6):
        return   indices[0]*shape[1]*shape[2]*shape[3]*shape[4]*shape[5]\
               + indices[1]*shape[2]*shape[3]*shape[4]*shape[5]\
               + indices[2]*shape[3]*shape[4]*shape[5]\
               + indices[3]*shape[4]*shape[5]\
               + indices[4]*shape[5]\
               + indices[5]
    
    #General case
    index = 0
    for n,i in enumerate(indices):
        index += i*reduce(lambda x, y: x*y, shape[(n+1):],1)
        
    return int(index)
    
def V_to_T_mapping(index,shape): #Test function written
    """Map a vector index to the tensor indices"""
    
    indices = np.zeros([len(shape),]).astype(int)
    
    for n,i in enumerate(range(len(shape))):
        factor     = int(reduce(lambda x, y: x*y, shape[(n+1):],1))
        indices[n] = int(index)/int(reduce(lambda x, y: x*y, shape[(n+1):],1))
        index  %= int(factor)
    return indices
    
def V_to_M_mapping(index,shape): #Test function written (contained in convert_V_to_M)
    """Map a vector index to the matrix indices"""
    #Define the mapping of the vector
    mapping = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1),(2,1),(2,0),(1,0)]
    #Get tensor indices
    indices = V_to_T_mapping(index,shape)
    #Get mapping index
    map_index = [i for i in range(9) if ((mapping[i][0]==indices[0]) and (mapping[i][1]==indices[1]))][0]
    #Get mapping factors
    factors = [int(reduce(lambda x, y: x*y, shape[(n+1):-1],1)) for n in range(len(indices))[2:-1]]
    #Get matrix indices
    i = int(map_index+len(mapping)*sum([f*i for f,i in zip(factors,indices[2:-1])]))
    j = int(indices[-1] if len(indices)>2 else 0)
    return i,j
    
    
def convert_V_to_M(V,shape): #Test function written
    """Convert a vector to its matrix form"""
    
    #Compute the size of the matrix
    Ndim = int(shape[-1] if len(shape)>2 else 1)
    Mdim = len(V)/Ndim
    
    #Initialize the empty matrix
    M = np.zeros([Mdim,Ndim])
    
    #Iterate through all members of V
    for index in range(len(V)):
        i,j = V_to_M_mapping(index,shape)
        M[i,j] = V[index]
    return M
    
def convert_M_to_V(M,shape):
    """Convert a matrix to its vector form"""
    
    #Compute the size of the vector
    Vdim = M.shape[0]*M.shape[1]
    
    #Initialize the empty vector
    V = np.zeros([Vdim,])
    
    #Iterate through all members of V
    for index in range(Vdim):
        #Get tensor indices
        i,j = V_to_M_mapping(index,shape)
        V[index] = M[i,j]
    return V
    
def convert_V_to_T(V,shape):
    """Convert a vector to its tensor form"""
    return np.reshape(V,shape)
    
#Testing suite
    
class TestHex8(unittest.TestCase):

    f                    = None
    original_directory   = ""
    output_file_name     = r"hex8_unittests.txt"
    output_file_location = r".\tests\unittests"
    currentResult        = None
    @classmethod
    def setUpClass(self):
        """Setup method"""
        output_file = os.path.join(self.output_file_location,self.output_file_name)
        if(os.path.isfile(output_file)):
            os.remove(output_file)
        self.f = open(output_file,"w+")
    @classmethod
    def tearDownClass(self):
        """Teardown method"""
        self.f.close()
        
    def setUp(self):
        pass
        
    def tearDown(self):
        ok = self.currentResult.wasSuccessful()
        tname = self.id().split(".")[-1]
        self.f.write(tname+"\t&\t"+str(ok)+"\n")
        
    def run(self, result=None):
        """Redefine run to keep track of results"""
        self.currentResult = result
        unittest.TestCase.run(self,result)
            
        
    def test_Hex8_shape_function(self):
        """Test the shape function for a Hex8 element"""
        #Test 1
        xi_vec = [0.,0.,0.]
        result = [Hex8_shape_function(n,xi_vec) for n in range(8)]
        self.assertEqual(np.allclose(result,0.125),True)
        
        #Test 2
        coords = [Hex8_node_coords(n) for n in range(8)]
        temp = [[Hex8_shape_function(n,nlc) for nlc in coords ] for n in range(8)]
        result = [ len([val for val in t if abs(val)>1e-9]) for t in temp]
        self.assertEqual(np.allclose(result,1),True)
        
    def test_Hex8_shape_function_loc(self):
        """Test the shape function for a Hex8 element"""
        #Test 1
        xi_vec = [0.,0.,0.]
        coords = [Hex8_node_coords(n) for n in range(8)]
        result = [Hex8_shape_function_loc(xi_vec,nlc) for nlc in coords]
        self.assertEqual(np.allclose(result,0.125),True)
        
        #Test 2
        temp = [[Hex8_shape_function_loc(coords[n],nlc) for nlc in coords ] for n in range(8)]
        result = [ len([val for val in t if abs(val)>1e-9]) for t in temp]
        self.assertEqual(np.allclose(result,1),True)
        
    def test_Hex8_local_grad_shape_function(self):
        """Test the local gradient of the shape function for a Hex8 element"""
        #Test 1
        xi_vec = [-.65,.23,.9]
        gradients = [Hex8_local_grad_shape_function(n,xi_vec) for n in range(8)]
        self.assertEqual(np.allclose(gradients,self._numeric_sflocal_gradient(xi_vec)),True)        
        
    def test_Hex8_global_grad_shape_function(self):
        """Test the global gradient of the shape function for a Hex8 element"""
        #Test 1 (note, had to drop tolerance on np.allclose)
        xi_vec = [-.65,.23,.9]
        gcoords = [[0,0,0],[1.5,0,0],[1.0,1,0],[0,1,0],\
                   [0,0,1],[1.5,0,1],[1.0,1,1],[0,1,1]]
        gradients = [Hex8_global_grad_shape_function(n,xi_vec,gcoords)[0] for n in range(len(gcoords))]
        
        self.assertEqual(np.allclose(gradients,self._numeric_sfglobal_gradient(xi_vec,gcoords),atol=1e-4,rtol=1e-4),True)
        
    def test_Hex8_get_shape_function_info(self):
        """Test the Hex8_get_shape_function_info function"""
        xi_vec = [-.65,.23,.9]
        gcoords = [[0,0,0],[1.5,0,0],[1.0,1,0],[0,1,0],\
                   [0,0,1],[1.5,0,1],[1.0,1,1],[0,1,1]]
        vals = [Hex8_get_shape_function_info(n,xi_vec,gcoords) for n in range(len(gcoords))]
        
        Ns,gradients,detJs = zip(*vals)
        
        self.assertEqual(np.allclose(Ns,[Hex8_shape_function(n,xi_vec) for n in range(8)]),True)
        self.assertEqual(np.allclose(gradients,self._numeric_sfglobal_gradient(xi_vec,gcoords),atol=1e-4,rtol=1e-4),True)
        
    def test_invert_3x3_matrix(self):
        """Test the matrix inversion function"""
        A = np.reshape(np.array([1,2,3,4,5,6,7,8,6]).astype(float),[3,3])
        self.assertEqual(np.isclose(np.linalg.det(A),invert_3x3_matrix(A)[0]),True)
        self.assertEqual(np.allclose(np.linalg.inv(A),invert_3x3_matrix(A)[1]),True)
        
    def test_invert_3x3_matrix_V(self):
        """Test the matrix in vector form inversion function"""
        A = np.reshape(np.array([1,2,3,4,5,6,7,8,6]).astype(float),[3,3])
        V = reduce_tensor_to_vector_form(A)
        self.assertEqual(np.isclose(np.linalg.det(A),invert_3x3_matrix_V(V)[0]),True)
        self.assertEqual(np.allclose(reduce_tensor_to_vector_form(np.linalg.inv(A)),invert_3x3_matrix_V(V)[1]),True)
        
    def test_get_jacobian(self):
        """Test the get_jacobian function"""
        xi_vec  = [-.65,.23,.9]
        gcoords = [[0,0,0],[1.5,0,0],[1.0,1,0],[0,1,0],\
                   [0,0,1],[1.5,0,1],[1.5,1,1],[0,1,1]]
        dNdxis  = [Hex8_local_grad_shape_function(n,xi_vec) for n in range(8)]
        dNdxi   = np.zeros([3,8])
        xmat    = np.zeros([8,3])
        
        for i in range(8):
            dNdxi[0,i] =  dNdxis[i][0]
            dNdxi[1,i] =  dNdxis[i][1]
            dNdxi[2,i] =  dNdxis[i][2]
            xmat[i,0]  = gcoords[i][0]
            xmat[i,1]  = gcoords[i][1]
            xmat[i,2]  = gcoords[i][2]

        self.assertEqual(np.allclose(np.dot(dNdxi,xmat),get_jacobian(gcoords,xi_vec)),True)
        
        
    def test_vector_dyadic_product(self):
        """Test the vector dyadic product function"""
        a = [ 1, 2, 3]
        b = [-2, 7, 8]
        A = np.array([[a[0]*b[0], a[0]*b[1], a[0]*b[2]],\
                      [a[1]*b[0], a[1]*b[1], a[1]*b[2]],\
                      [a[2]*b[0], a[2]*b[1], a[2]*b[2]]])
        self.assertEqual(np.allclose(vector_dyadic_product(a,b),A),True)
        
    def test_Hex8_interpolation_matrix(self):
        """Test the generation of the interpolation matrix"""
        
        nodal_global_coords = np.reshape(range(24),[8,3])
        
        nodal_dof_values    = np.reshape(range(24,24+8*3),[8,3])
        
        M = Hex8_interpolation_matrix(nodal_global_coords,nodal_dof_values)
        
        Msoln = np.reshape([ 1, 1, 1, 1, 1, 1, 1, 1,\
                             0, 3, 6, 9,12,15,18,21,\
                             1, 4, 7,10,13,16,19,22,\
                             2, 5, 8,11,14,17,20,23,\
                            24,27,30,33,36,39,42,45,\
                            25,28,31,34,37,40,43,46,\
                            26,29,32,35,38,41,44,47],[7,8])
        self.assertEqual(np.allclose(M,Msoln),True)
        
    def test_get_1D_gpw(self):
        """Test the returned gauss points from get_1D_gpw"""
        f1 = lambda x:         .7*x     +    1.7
        a1 = lambda x:    0.5*0.7*x**2. +    1.7*x
        f2 = lambda x:       0.12*x**2. +    3.4*x     - 6
        a2 = lambda x:    0.12/3.*x**3. +0.5*3.4*x**2. - 6*x
        f3 = lambda x:      -2.13*x**3. +    2.1*x**2. - 6.3*x         + 2.1
        a3 = lambda x: -0.25*2.13*x**4. + 2.1/3.*x**3. - 0.5*6.3*x**2. + 2.1*x
        
        gpws = [get_1D_gpw(n) for n in range(1,4)]
        
        res1 = sum([f1(x)*w for x,w in zip(gpws[0][0],gpws[0][1])])
        res2 = sum([f2(x)*w for x,w in zip(gpws[1][0],gpws[1][1])])
        res3 = sum([f3(x)*w for x,w in zip(gpws[2][0],gpws[2][1])])
        
        int1 = a1(1)-a1(-1)
        int2 = a2(1)-a2(-1)
        int3 = a3(1)-a3(-1)
        
        #print "\nres1: {0}\tres2: {1}\tres3: {2}\nint1: {3}\tint2: {4}\tint3: {5}".format(res1,res2,res3,int1,int2,int3)
        
        self.assertEqual(np.isclose(res1,int1),True)
        self.assertEqual(np.isclose(res2,int2),True)
        self.assertEqual(np.isclose(res3,int3),True)
        
    def test_get_gpq(self):
        """Test the returned gauss points from get_gpw"""
        f = lambda x,y,z: x**2 - y*z**2. + z
        
        ps,ws = get_gpw(3)
        
        res1 = sum([f(p[0],p[1],p[2])*w for p,w in zip(ps,ws)])
        int1 = 8/3.
        
        #print "res: {0}\nint: {1}".format(res1,int1)
        
        self.assertEqual(np.isclose(res1,int1),True)
        
    def test_BVoigt_submatrix(self):
        """Test the BVoigt submatrix function"""
        
        n = 0
        xi_vec  = [-.65,.23,.9]
        gcoords = [[0,0,0],[1.5,0,0],[1.0,1,0],[0,1,0],\
                   [0,0,1],[1.5,0,1],[1.5,1,1],[0,1,1]]
        
        BVcalc = BVoigt_submatrix(n,xi_vec,gcoords)
        
        dNdx,detJ = Hex8_global_grad_shape_function(n,xi_vec,gcoords)
        
        B = np.zeros([6,3])
        B[0,0] = dNdx[0]
        B[1,1] = dNdx[1]
        B[2,2] = dNdx[2]
        B[3,0] = dNdx[1]
        B[3,1] = dNdx[0]
        B[4,1] = dNdx[2]
        B[4,2] = dNdx[1]
        B[5,0] = dNdx[2]
        B[5,2] = dNdx[0]
        
        self.assertEqual(np.allclose(BVcalc,B),True)
        
    def test_compute_BVoigt(self):
        """Test for computing the full strain displacement matrix using compute_BVoigt"""
        xi_vec  = [-.65,.23,.9]
        gcoords = [[0,0,0],[1.5,0,0],[1.0,1,0],[0,1,0],\
                   [0,0,1],[1.5,0,1],[1.5,1,1],[0,1,1]]
        B = compute_BVoigt(xi_vec,gcoords)
        
        Bsubs = [BVoigt_submatrix(n,xi_vec,gcoords) for n in range(8)]
        
        tests = [np.allclose(B[:,3*n:(3*n+3)],Bsubs[n]) for n in range(8)]
        
        self.assertEqual(sum(tests),8)
        
    def test_matrix_dot(self):
        """Test the matrix dot product"""
        A = np.reshape(np.random.rand(9),[3,3])
        B = np.reshape(np.random.rand(9),[3,3])
        C = A.dot(B)
        Cdot = matrix_dot(A,B)
        self.assertEqual(np.allclose(C,Cdot),True)
        
    def test_matrix_dot_V(self):
        """Test the matrix dot product in vector form"""
        A = np.reshape(np.random.rand(9),[3,3])
        B = np.reshape(np.random.rand(9),[3,3])
        C = np.reshape(A.dot(B),[9,])
        Cdot = matrix_dot_V(np.reshape(A,[9,]),np.reshape(B,[9,]))
        self.assertEqual(np.allclose(C,Cdot),True)
        
    def test_matrix_Tdot(self):
        """Test the matrix transpose dot product"""
        A = np.reshape(np.random.rand(9),[3,3])
        B = np.reshape(np.random.rand(9),[3,3])
        C = (A.T).dot(B)
        Cdot = matrix_Tdot(A,B)
        self.assertEqual(np.allclose(C,Cdot),True)
        
    def test_matrix_Tdot_V(self):
        """Test the matrix transpose dot product in vector form"""
        A = np.reshape(np.random.rand(9),[3,3])
        B = np.reshape(np.random.rand(9),[3,3])
        C = np.reshape((A.T).dot(B),[9,])
        Cdot = matrix_dot_V(np.reshape(A.T,[9,]),np.reshape(B,[9,]))
        self.assertEqual(np.allclose(C,Cdot),True)
        
    def test_matrix_Tdot_TOT(self):
        """Test the matrix transpose dot product with a third order tensor"""
        A = np.reshape(np.random.rand(9),[3,3])
        B = np.reshape(np.random.rand(27),[3,3,3])
        C = np.einsum('li,ljk',A,B)
        Cdot = matrix_Tdot_TOT(A,B)
        self.assertEqual(np.allclose(C,Cdot),True)
        
    def test_vector_dot_matrix(self):
        """Test the vector dot product"""
        v = np.reshape(np.random.rand(3),[3,])
        A = np.reshape(np.random.rand(9),[3,3])
        w = v.dot(A)
        wdot = vector_dot_matrix(v,A)
        self.assertEqual(np.allclose(w,wdot),True)
        
    def test_get_symm_matrix(self):
        """Test getting the symmetric part of a matrix A"""
        A = np.reshape(np.random.rand(9),[3,3])
        Asymm = 0.5*(A+A.T)
        AsymmA = get_symm_matrix(A)
        self.assertEqual(np.allclose(Asymm,AsymmA),True)
        
    def test_get_symm_matrix_V(self):
        """Test getting the symmetric part of a matrix in vector form"""
        V      = np.array(range(9))
        Asymm  = get_symm_matrix_V(V)
        At     = convert_V_to_T(V,[3,3])
        AsymmT = 0.5*(At+At.T)
        AsymmT = reduce_tensor_to_vector_form(AsymmT)
        self.assertTrue(np.allclose(Asymm,AsymmT))
        
    def test_T_to_V_mapping(self):
        """Test the mapping from the tensor indices to the vector form index"""
        A = np.reshape(range(2*7*3*5),[2,7,3,5])
        indices = [(0,0,0,2),(0,0,2,1),(1,0,0,2),(1,1,1,1),(1,6,0,4)]
        results = [T_to_V_mapping(i,A.shape) for i in indices]
        answers = [2,11,107,126,199]
        
        #print results
        #print answers
        
        self.assertEqual(np.allclose(results,answers),True)
        
        #Test 2
        a = range(9)
        A = np.reshape(a,[3,3])
        indices = [(1,2),(1,1),(1,0)]
        mapping = [T_to_V_mapping(i,A.shape) for i in indices]
        results = [A[i,j] for i,j in indices]
        answers = [a[m] for m in mapping]
        self.assertTrue(np.allclose(results,answers))
        
        #Test 3
        a = range(3)
        A = np.reshape(a,[3])
        indices = [(1,),(0,)]
        mapping = [T_to_V_mapping(i,A.shape) for i in indices]
        results = [A[i] for i in indices]
        answers = [a[m] for m in mapping]
        self.assertTrue(np.allclose(results,answers))
        
        #Test 4
        a = range(2*6*3)
        A = np.reshape(a,[2,6,3])
        indices = [(1,4,0),(0,0,2),(1,5,0)]
        mapping = [T_to_V_mapping(i,A.shape) for i in indices]
        results = [A[i,j,k] for i,j,k in indices]
        answers = [a[m] for m in mapping]
        self.assertTrue(np.allclose(results,answers))
        
        #Test 5
        a = range(5*2*6*9*3)
        A = np.reshape(a,[5,2,6,9,3])
        indices = [(4,1,5,8,0),(0,0,3,3,2),(1,1,0,2,2)]
        mapping = [T_to_V_mapping(i,A.shape) for i in indices]
        results = [A[i,j,k,l,m] for i,j,k,l,m in indices]
        answers = [a[m] for m in mapping]
        self.assertTrue(np.allclose(results,answers))
        
        #Test 6
        a = range(4*6*2*1*7*9)
        A = np.reshape(a,[4,6,2,1,7,9])
        indices = [(2,1,1,0,3,5),(0,0,1,0,6,6),(1,1,1,0,5,7)]
        mapping = [T_to_V_mapping(i,A.shape) for i in indices]
        results = [A[i,j,k,l,m,n] for i,j,k,l,m,n in indices]
        answers = [a[m] for m in mapping]
        self.assertTrue(np.allclose(results,answers))
        
        #Test 7
        a = range(6*4*5*6*3*9*3)
        A = np.reshape(a,[6,4,5,6,3,9,3])
        indices = [(2,1,1,0,2,5,1),(0,0,1,0,0,6,1),(1,1,1,0,1,7,0)]
        mapping = [T_to_V_mapping(i,A.shape) for i in indices]
        results = [A[i,j,k,l,m,n,o] for i,j,k,l,m,n,o in indices]
        answers = [a[m] for m in mapping]
        self.assertTrue(np.allclose(results,answers))
        
        #Test 8 (should fail)
        print "\nExecuting Test 8: Should observe Error!\n"
        try:
            T_to_V_mapping([2,1,3],[1,0])
            self.assertTrue(False)
        except:
            self.assertTrue(True)
            print "\nTest 8 passed!\n"
        
    def test_V_to_T_mapping(self):
        """Test the mapping from the vector index to the tensor indices"""
        A = np.reshape(range(2*7*3*5),[2,7,3,5])
        indices = [2,11,107,126,199]
        results = [V_to_T_mapping(i,A.shape) for i in indices]
        answers = [(0,0,0,2),(0,0,2,1),(1,0,0,2),(1,1,1,1),(1,6,0,4)]
        
        #print results
        #print answers
        
        self.assertEqual(np.allclose(results,answers),True)
        
    def test_convert_V_to_M(self):
        """Test the conversion of a vector to a matrix"""
        
        #Test 1
        V        = range(3*3*4*6*5)
        shape    = [3,3,4,6,5]
        A        = np.reshape(V,shape)
        results  = convert_V_to_M(V,shape)
        
        answers = np.zeros([3*3*4*6,5])
        
        mapping = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1),(2,1),(2,0),(1,0)]
        
        for j in range(5):
            for i in range(9):
                mi,mj = mapping[i]
                for k in range(4):
                    for l in range(6):
                        answers[i+l*9+k*6*9,j] = A[mi,mj,k,l,j]
                        
        self.assertEqual(np.allclose(results,answers),True)
        
        #Test 2
        V = range(3*3)
        shape = [3,3]
        A = np.reshape(V,shape)
        results = convert_V_to_M(V,shape)
        
        answers = np.zeros([9,1])
        
        for j in range(1):
            for i in range(9):
                mi,mj = mapping[i]
                answers[i,j] = A[mi,mj]
        self.assertEqual(np.allclose(results,answers),True)
        
    def test_convert_M_to_V(self):
        """Test the conversion of a vector to a matrix"""
        
        #Test 1
        V        = range(3*3*4*6*5)
        shape    = [3,3,4,6,5]
        A        = np.reshape(V,shape)
        results  = convert_V_to_M(V,shape)
        
        answers  = reduce_tensor_to_matrix_form(A)
                        
        self.assertEqual(np.allclose(results,answers),True)
        
        #Test 2
        V = range(3*3)
        shape = [3,3]
        A = np.reshape(V,shape)
        results = convert_V_to_M(V,shape)
        answers = reduce_tensor_to_matrix_form(A)
        self.assertEqual(np.allclose(results,answers),True)
        
    def test_convert_V_to_T(self):
        #Test 1
        V        = range(3*3*4*6*5)
        shape    = [3,3,4,6,5]
        A        = np.reshape(V,shape)
        results  = convert_V_to_T(V,shape)
        
        answers  = A
                        
        self.assertEqual(np.allclose(results,answers),True)
        
        #Test 2
        V = range(3*3)
        shape = [3,3]
        A = np.reshape(V,shape)
        results = convert_V_to_T(V,shape)
        answers = A
        self.assertEqual(np.allclose(results,answers),True)
        
    def test_reduce_tensor_to_matrix_form(self):
        """Test the conversion of a tensor to a matrix"""
        #Test 1
        V        = range(3*3*4*6*5)
        shape    = [3,3,4,6,5]
        A        = np.reshape(V,shape)
        results  = reduce_tensor_to_matrix_form(A)
        
        answers = np.zeros([3*3*4*6,5])
        
        mapping = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1),(2,1),(2,0),(1,0)]
        
        for j in range(5):
            for i in range(9):
                mi,mj = mapping[i]
                for k in range(4):
                    for l in range(6):
                        answers[i+l*9+k*6*9,j] = A[mi,mj,k,l,j]
                        
        self.assertEqual(np.allclose(results,answers),True)
        
        #Test 2
        V = range(3*3)
        shape = [3,3]
        A = np.reshape(V,shape)
        results = reduce_tensor_to_matrix_form(A)
        
        answers = np.zeros([9,1])
        
        for j in range(1):
            for i in range(9):
                mi,mj = mapping[i]
                answers[i,j] = A[mi,mj]
        self.assertEqual(np.allclose(results,answers),True)
        
    def _numeric_sflocal_gradient(self,v):
        """Compute a numeric gradient of the shape function in local coordinates"""
        h = 1e-7
        varray = [[[v[0]+h,  v[1],  v[2]],[v[0]-h,  v[1],  v[2]]],\
                  [[  v[0],v[1]+h,  v[2]],[  v[0],v[1]-h,  v[2]]],\
                  [[  v[0],  v[1],v[2]+h],[  v[0],  v[1],v[2]-h]]]
                  
        coords = [Hex8_node_coords(i) for i in range(8)]
        
        results = [[(Hex8_shape_function_loc(v2,nlc)-Hex8_shape_function_loc(v1,nlc))/(2*h) for v2,v1 in varray] for nlc in coords]
        return results
        
    def _numeric_sfglobal_gradient(self,v,gcoords):
        """Compute a numeric gradient of the shape function in global coordinates"""
        h = 1e-5
        
        J = get_jacobian(gcoords,v)
        detJ,invJ = invert_3x3_matrix(J)
        
        #Retrieving the global coordinates of the point
        M = Hex8_interpolation_matrix(gcoords,[[],[],[],[],[],[],[],[]])
        N = np.reshape(np.array([Hex8_shape_function(n,v) for n in range(8)]),[8,1])
        vg = np.dot(M,N)[1:4]
        
        #Global perturbation
        vgarray = [np.array([ vg[0]+h,   vg[1],   vg[2]]), np.array([ vg[0]-h,     vg[1],     vg[2]]),\
                   np.array([   vg[0], vg[1]+h,   vg[2]]), np.array([   vg[0],   vg[1]-h,     vg[2]]),\
                   np.array([   vg[0],   vg[1], vg[2]+h]), np.array([   vg[0],     vg[1],   vg[2]-h])]
                  
        varray  = []
                  
        #Solve for local perturbed node coordinates (Newton-Raphson)
                  
        for vp in vgarray:
            vc = np.reshape(np.copy(v),[3,1])
            
            N  = np.reshape(np.array([Hex8_shape_function(n,vc) for n in range(8)]),[8,1])
            R  = vp - np.dot(M,N)[1:4]
            
            r  = np.linalg.norm(R)
            r0 = r
            niter = 0
            maxiter = 20
            tol = 1e-9
            
            while(((r/r0)>tol) and (r>tol) and niter<maxiter):
                J = get_jacobian(gcoords,vc)
                detJ,invJ = invert_3x3_matrix(J)
                vc += np.dot(invJ,R)
                N  = np.reshape(np.array([Hex8_shape_function(n,vc) for n in range(8)]),[8,1])
                R  = vp - np.dot(M,N)[1:4]
                r  = np.linalg.norm(R)
                niter += 1
                
            if((r/r0)>tol and (r>tol) and niter>=maxiter):
                print "\nError: Newton Raphson failed to converge\n"
                raise
                
            else:
                varray.append(np.copy(vc))
                  
        Narray = [np.reshape(np.array([Hex8_shape_function(n,vt) for n in range(8)]),[8,1]) for vt in varray]
        Iarray = [np.dot(M,N)[1:4] for N in Narray]
        varr = [[varray[i],varray[i+1]] for i in range(0,6,2)]
        
        coords = [Hex8_node_coords(n) for n in range(8)]
        return [np.concatenate([(Hex8_shape_function_loc(v2,nlc)-Hex8_shape_function_loc(v1,nlc))/(2*h) for v2,v1 in varr]) for n,nlc in enumerate(coords)]
        
    def _flatten_list(self,l):
        """Flatten list of lists"""
        return [item for sublist in l for item in sublist]
        
if __name__ == '__main__':
    unittest.main()