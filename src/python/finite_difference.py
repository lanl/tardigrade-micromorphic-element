import numpy as np
import unittest
import os

def finite_difference(fxn,x0,h,accuracy_order=2): #Test Function Written
    """Compute the finite difference derivative of a function in the direction of h"""
    
    #x0 = np.reshape(x0,[len(x0),1])
    deltah_multipliers = get_deltah_multipliers(accuracy_order)
    coefficients  = get_coefficients(accuracy_order)
    
    if(type(h)!=float):
        hden = np.linalg.norm(h)
    else:
        hden = h
    
    return sum([C*fxn(x0+h*dhm) for C,dhm in zip(coefficients,deltah_multipliers) if abs(C)>0])/hden
    
def numeric_gradient(fxn,x0,h,accuracy_order=2):
    """Compute the numeric gradient of a function of many variables.
    fxn:            The function in question
    x0:             The point to compute the gradient at
    h:              A scalar value for the perturbation to x0
    accuracy_order: The accuracy of the solution
    """
    
    #Form the perturbation vectors
    hvecs = np.zeros([len(x0),len(x0)])
    for i in range(len(x0)):
        hvecs[i,i] = h
        
    #Compute the derivatives
    return np.array([finite_difference(fxn,x0,hvecs[:,i],accuracy_order=accuracy_order) for i in range(len(x0))])
    
    
    
def get_deltah_multipliers(accuracy_order): #Test Function Written
    """Get the deltah multipliers"""
    deltah_multipliers = [[-1.,0.,1.],[-2.,-1.,0.,1.,2.],\
                          [-3.,-2.,-1.,0.,1.,2.,3.],[-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.]]
    return deltah_multipliers[accuracy_order/2-1]

def get_coefficients(accuracy_order):
    """Get the coefficients of the finite difference"""
    first_derivative_coeffs = [[-0.5,0,0.5],[1./12,-2./3,0.,2./3,-1./12],\
                               [-1./60,3./20,-3./4,0.,3./4,-3./20,1./60],\
                               [-1./280,-4./105,-1./5,-4./5,0,4./5,-1./5,4./105,-1./280]]
                               
    return first_derivative_coeffs[accuracy_order/2-1]
    
class TestMicroElement(unittest.TestCase):

    f                    = None
    original_directory   = ""
    output_file_name     = r"finite_difference_unittests.txt"
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

    def test_finite_difference(self):
        """Test the finite difference function"""
        #Test 1
        f1 = lambda x: x**2.
        a1 = lambda x: 2.*x
        
        #Test 2
        f2 = lambda x: np.array([x[0]**2.,x[1]-1])
        a2 = lambda x: np.array([2*x[0], 1])
        
        x01 = 2.6
        x02 = np.array([2.3,-5.7])
        
        result1 = finite_difference(f1,x01,1e-6)
        answer1 = a1(x01)
        
        result2 = finite_difference(f2,x02,np.array([1e-6,0]))[0]
        answer2 = a2(x02)[0]
        
        self.assertEqual(np.allclose(result1,answer1),True)
        self.assertEqual(np.allclose(result2,answer2),True)
        
    def test_numeric_gradient(self):
        """Test the numeric gradient function"""
        #Test 1
        f1 = lambda x: np.array([x[0]**2.,x[1]-1])
        a1 = lambda x: np.array([[2*x[0], 0,],[0, 1]])
        
        #Test 2
        f2 = lambda x: np.array([x[1]*x[0]**2,x[1]+x[2]*x[0], x[0]])
        a2 = lambda x: np.array([[2*x[0]*x[1], x[2], 1.],\
                                 [   x[0]**2.,   1., 0.],\
                                 [         0., x[0], 0.]])
        
        x01 = np.array([2.3,-5.7])
        result1 = numeric_gradient(f1,x01,1e-6)
        answer1 = a1(x01)
        
        x02 = np.array([1.4,-2.,3.4])
        result2 = numeric_gradient(f2,x02,1e-6)
        answer2 = a2(x02)
        
        self.assertEqual(np.allclose(result1,answer1),True)
        self.assertEqual(np.allclose(result2,answer2),True)

    def test_get_deltah_multipliers(self):
        """Test the get_deltah_multipliers subroutine"""
        answers = [[-1,0,1],[-2,-1,0,1,2],[-3,-2,-1,0,1,2,3],[-4,-3,-2,-1,0,1,2,3,4]]
        results = [get_deltah_multipliers(i+1) for i in range(2,10,2)]
        self.assertEqual(np.allclose(self._flatten_list(answers),self._flatten_list(results)),True)
        
    def test_get_coefficients(self):
        """Test the get_coefficients subroutine"""
        answers = [[-0.5,0,0.5],[1./12,-2./3,0.,2./3,-1./12],\
                               [-1./60,3./20,-3./4,0.,3./4,-3./20,1./60],\
                               [-1./280,-4./105,-1./5,-4./5,0,4./5,-1./5,4./105,-1./280]]
        results = [get_coefficients(i+1) for i in range(2,10,2)]
        self.assertEqual(np.allclose(self._flatten_list(answers),self._flatten_list(results)),True)
        
    def _flatten_list(self,l):
        """Flatten list of lists"""
        return [item for sublist in l for item in sublist]
        
    
if __name__ == '__main__':
    unittest.main()