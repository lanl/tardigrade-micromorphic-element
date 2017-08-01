import numpy as np
import unittest
import os

import micro_element as mel

class InputParser:
    """The input parser class which reads in data from a text file"""
    
    filename      = "" #String giving the filename
    latex_string  = "" #The description of the input deck
    nodes_coords  = [] #List of node numbers and their reference coordinates
    dirichlet_bcs = [] #List of the dirichlet boundary conditions
    ndof          = 0. #The number of degrees of freedom at each node
    nodes_dof     = [] #List of lists of degrees of freedom assocated with each node
    element_nodes = [] #List of nodes associated with each element
    properties    = [] #List of the material properties
    
    def __init__(self,filename_):
        """Initialize the input parser object"""
        self.filename = filename_
        
    def __repr__(self):
        """Define the repr string"""
        return "InputParser({0})".format(self.filename)
    
    def read_input(self):
    
        print "\n=================================================\n"+\
                "|                                               |\n"+\
                "|                  INPUT PARSER                 |\n"+\
                "|                                               |\n"+\
                "=================================================\n"
            
        #Read in the input deck 
        print "Reading data from: {0}".format(self.filename)
        
        with open(self.filename) as f:
            for line in f.readlines():
                line = line.strip()
                if((len(line)>0) and (line[0]!="#")):
                    bool,fxn,line = self.check_keywords(line)
                    if(bool):
                        parse_fxn = fxn
                    parse_fxn(line)
        
        print "\n=================================================\n"+\
                "|                                               |\n"+\
                "|             INPUT PARSER COMPLETED            |\n"+\
                "|                                               |\n"+\
                "=================================================\n"
        
    def check_keywords(self,line):
        """Check a line for keywords"""
        keywords   = ["*NODES","*DIRICHLET_BCS","*ELEMENTS","*PROPERTIES","*LATEX"]
        parse_fxns = [self.parse_nodes, self.parse_dirichlet_bcs, self.parse_elements, self.parse_props, self.parse_latex]
        
        for key,fxn in zip(keywords,parse_fxns):
            if key in line:
                print "Keyword \"{0}\" found".format(key)
                line = line.replace(key,'')
                return True,fxn,line
        return False,None,line
        
    def parse_latex(self,line):
        """Parse the file when triggered by a latex keyword"""
        self.latex_string += line
        
    def parse_nodes(self,line):
        """Parse the file when triggered by a node keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)==2):
                self.ndof = int(sline[1])
                self.nodes_dof = [range(n*self.ndof,(n+1)*self.ndof) for n in range(len(self.nodes_coords))]
            elif(len(sline)==4):
                self.nodes_coords.append([int(sline[0]),float(sline[1]),float(sline[2]),float(sline[3])])
            else:
                print "Error: A node must be defined by its number,"+\
                      "       followed by its x,y,z coordinates."
                raise ValueError()
                
    def parse_dirichlet_bcs(self,line):
        """Parse the dirichlet boundary conditions applied at the nodes"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)==3):
                self.dirichlet_bcs.append([int(sline[0]),int(sline[1]),float(sline[2])])
            else:
                print "Error: A boundary condition is defined by,"+\
                      "       the node it is applied to, the dof"+\
                      "       on the node, and the value of the"+\
                      "       displacement."
 
    def parse_elements(self,line):
        """Parse the file when triggered by a element keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)<8):
                print "Error: An element must be defined by 8 nodes."
                raise ValueError()
            else:
                self.element_nodes.append([int(sline[0]), int(sline[1]), int(sline[2]), int(sline[3]), int(sline[4]),\
                                                          int(sline[5]), int(sline[6]), int(sline[7]), int(sline[8])])
        
    def parse_props(self,line):
        """Parse the file when triggered by a properties keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(self.properties)==0):
                self.properties = [float(i) for i in sline]
            else:
                print "Error: Properties can only be defined once."
                raise ValueError()
        
class TestFEA(unittest.TestCase):

    f                    = None
    original_directory   = ""
    output_file_name     = r"fea_driver_unittests.txt"
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
        
    def test_read_input(self):
        """Test the read input command"""
        
        IP = InputParser("unittest.inp")
        
        IP.read_input()
      
class FEAModel():
    """The class which defines the FEA model"""
    input      = None #Input of type InputParser
    t0         = 0.   #Previous value of time
    ti         = 1.   #Current value of time
    x0         = None #Previous value of x in Ax = b
    xi         = None #Current increment value of x in solver
    total_ndof = None #Total number of degrees of freedom
    RHS        = None #The right hand side vector
    AMATRX     = None #The A matrix in Ax = b
    RHS_bc     = None #The right hand side vector with dirichlet bcs applied
    AMATRX_bc  = NONE #The A matrix with dirichlet bcs applied
    gdof       = None
    def __init__(sefl,IP_):
        """Initialize the finite element model"""
        self.input = IP_
        total_ndof = len(self.input.nodes_dof)*self.input.ndof
        self.x0    = np.zeros([self.total_ndof])
        self.xi    = np.zeros([self.total_ndof])
    
    def __repr__(self):    
        """Define the repr string"""
        return "FEAModel({0})".format(self.input)
        
    def assemble_RHS_and_jacobian_matrix(self):
        """Assemble the global right hand side vector and jacobian matrix"""
        
        self.RHS    = np.zeros([self.total_ndof]) #Initialize the right hand side
        self.AMATRX = np.zeros([self.total_ndof,self.total_ndof]) #Initialize the A matrix
        
        for e in self.input.element_nodes:
            #Get the reference coordinates of the nodes
            coords = [self.input.nodes_coords[i][1:] for ei in e if self.input.nodes_coords[i][0]==ei]
            #Get the dof associated with each node in the element
            edof = [self.input.nodes_dof[ei] for ei in e]
            flat_edof = flatten_list(edof)
            #Get the previous values of the dof
            dof0  = [[self.x0[d] for d in eidof] for eidof in edof]
            #Get the current values of the dof
            dofi  = [[self.xi[d] for d in eidof] for eidof in edof]
            #Create the U and DU vectors
            U  = flatten_list(dofi)
            U0 = flatten_list(dof0)
            DU = [Ui-U0i for Ui,U0i in zip(U,U0)]
            #Call the element to provide the right hand side and jacobian matrix
            RHSe,AMATRXe = UEL(self.input.properties,len(self.input.properties),\
                               coords,None,U,DU,self.t0,self.ti-self.t0)
            #Assemble the global matrix
            for i in range(self.input.ndof*8):
                Iglob = flat_edof[i] #Get the mapping to the global degree of freedom index for i
                self.RHS[Iglob] += RHSe[i]
                for j in range(self.input.ndof*8):
                    Jglob = flat_edof[j] #Get the mapping to the global degree of freedom index for j
                    self.AMATRX[Iglob,Jglob] += AMATRXe[i,j]
        
    def apply_dirichlet_bcs(self):
        """Apply the dirichlet boundary conditions to the system of equations"""
        if(gdof==None):
            gdof = np.zeros([len(self.input.dirichlet_bcs),]).astype(int)
        for i,dbc in enumerate(self.input.dirichlet_bcs):
            #Get the global dof number and value
            gdof[i] = self.input.nodes_dof[dbc[0]][dbc[1]]
            val     = self.input.nodes_dof[dlb[0]][dbc[2]]
            
            self.RHS -= self.AMATRX[:,gdof[i]]*val #Subtract this force from the right hand side
            
        for dof in reversed(gdof):
            self.RHS_bc    = np.delete(self.RHS,dof,0)
            self.AMATRX_bc = np.delete(self.AMATRX,0)
            self.AMATRX_bc = np.delete(self.AMATRX,1)
            
    def solve_increment(self):
        """Solve the currently defined increment and update xi"""
        self.assemble_RHS_and_jacobian_matrix()
        self.apply_dirichlet_bcs()
        self.xi = np.linalg.solve(self.AMATRX_bc,self.RHS,bc)
        
def run_finite_elmement_model(input_filename):
    """Run the finite element model identified by the input filename"""
    IP = InputParser(input_filename)
    IP.read_input()
    
def flatten_list(l):
    """Flatten list of lists"""
    return [item for sublist in l for item in sublist]
    
        
if __name__ == '__main__':
    unittest.main()
        