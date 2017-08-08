import numpy as np
import unittest
import os

import micro_element as mel

import cProfile

class InputParser(object):
    """The input parser class which reads in data from a text file"""

    def __init__(self,filename_):
        """Initialize the input parser object"""
        self.filename = filename_
        temp = os.path.split(self.filename)
        self.filename = temp[1]              #String giving the filename
        self.path_to_file = temp[0]          #String giving the path to the file
        
        self.latex_string  = "" #The description of the input deck
        self.nodes_coords  = [] #List of node numbers and their reference coordinates
        self.dirichlet_bcs = [] #List of the dirichlet boundary conditions
        self.ndof          = 0. #The number of degrees of freedom at each node
        self.element_nodes = [] #List of nodes associated with each element
        self.properties    = [] #List of the material properties
        self.svars         = [] #List of the state variables
        self.nodesets      = [] #List of the nodesets [name, node1, node2, ...]
    
        self.mms_name      = None #The name of the manufactured solution being tested
        self.mms_fxn       = None #The function to compute U values for the method
                                  #of manufactured solutions
        self.mms_dir_set   = None #The nodeset to be used for the method of manufactured
                                  #solutions
        
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
        
        with open(os.path.join(self.path_to_file,self.filename)) as f:
            for linenum,line in enumerate(f.readlines()):
                line = line.strip()
                if((len(line)>0) and (line[0]!="#")):
                    bool,fxn,line = self.check_keywords(linenum,line)
                    if(bool):
                        parse_fxn = fxn
                    parse_fxn(linenum+1,line)
        
        f.close()
        print "\n=================================================\n"+\
                "|                                               |\n"+\
                "|             INPUT PARSER COMPLETED            |\n"+\
                "|                                               |\n"+\
                "=================================================\n"
        
    def check_keywords(self,linenum,line):
        """Check a line for keywords"""
        keywords   = ["*NODES","*DIRICHLET_BCS","*ELEMENTS","*PROPERTIES","*LATEX","*NSET","*MMS"]
        parse_fxns = [self.parse_nodes, self.parse_dirichlet_bcs, self.parse_elements, self.parse_props, self.parse_latex, self.parse_nodeset, self.parse_mms]
        for key,fxn in zip(keywords,parse_fxns):
            if key in line:
                print "Keyword \"{0}\" found".format(key)
                line = line.replace(key,'')
                return True,fxn,line
        return False,None,line
        
    def parse_latex(self,lineum,line):
        """Parse the file when triggered by a latex keyword"""
        self.latex_string += line + "\n"
        
    def parse_nodes(self,linenum,line):
        """Parse the file when triggered by a node keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)==2):
                self.ndof = int(sline[1])
            elif(len(sline)==4):
                self.nodes_coords.append([int(sline[0]),float(sline[1]),float(sline[2]),float(sline[3])])
            else:
                print "Error: On line {0}, a node must be defined by its number,"+\
                      "       followed by its x,y,z coordinates.".format(linenum)
                raise ValueError()
                
    def parse_dirichlet_bcs(self,linenum,line):
        """Parse the dirichlet boundary conditions applied at the nodes"""
        if(len(line)>0):
            sline = line.split(',')
            #Parse a definition of a boundary condition that specifies a single node
            if(len(sline)==3):
                try:
                    self.dirichlet_bcs.append([int(sline[0]),int(sline[1]),float(sline[2])])
                except:
                    print "Error: On line {0}, problem parsing line.".format(linenum)
            #Parse a definition of a boundary condition that specifies a nodeset
            if(len(sline)==4):
                #Remove whitespace from the nodeset name
                sline[1] = sline[1].strip()
                #Find the index of the node in self.nodeset
                index = [i for i in range(len(self.nodesets)) if self.nodesets[i][0]==sline[1]]
                #Error handing of the index to prevent multiple nodesets with the same name, undefined nodesets, etc.
                if (len(index)==1):
                    for node in self.nodesets[index[0]][1:]:
                        self.dirichlet_bcs.append([node, int(sline[2]), float(sline[3])])
                elif (len(index)>1):
                    print "Error: On line {0}, multiple sidesets with the same name.".format(linenum)
                    raise ValueError()
                else:
                    print "Error: On line{0}, nodeset {1} not found".format(linenum,sline[1])
                    raise ValueError()
                
            else:
                print "Error: On line {0}, a boundary condition is defined by,\n".format(linenum)+\
                      "       the node it is applied to, the dof\n"+\
                      "       on the node, and the value of the\n"+\
                      "       displacement or the command:\n"+\
                      "       sideset,name,local_dof,value.\n"+\
                      "       note, bc's must come AFTER the nodesets.\n"
                raise
 
    def parse_elements(self,linenum,line):
        """Parse the file when triggered by a element keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)<8):
                print "Error: On line {0}, an element must be defined by 8 nodes.".format(linenum)
                raise ValueError()
            else:
                self.element_nodes.append([int(sline[0]), int(sline[1]), int(sline[2]), int(sline[3]), int(sline[4]),\
                                                          int(sline[5]), int(sline[6]), int(sline[7]), int(sline[8])])
        
    def parse_props(self,linenum,line):
        """Parse the file when triggered by a properties keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(self.properties)==0):
                self.properties = [float(i) for i in sline]
            else:
                print "Error: On line {0}, properties can only be defined once.".format(linenum)
                raise ValueError()
                
    def parse_nodeset(self,linenum,line):
        """Parse the incoming nodesets when triggered by the keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)<2):
                print "Error: On line {0}, a nodeset must have at least one node.".format(linenum)
                raise ValueError()
            else:
                self.nodesets.append([sline[0]]+[int(i) for i in sline[1:]])
    
    def parse_mms(self,linenum,line):
        """Parse the method of manufactured solutions keyword"""
        mms_fxn_names = ["const_u","linear_u"]
        mms_fxns      = [self.mms_const_u,self.mms_linear_u]
        if(len(line)>0):
            sline =  line.split(',')
            #Get the method of manufactured solutions u function
            if((len(sline)==2) and (sline[0].strip()=="")):
                try:
                    i = mms_fxn_names.index(sline[1].strip())
                    self.mms_fxn = mms_fxns[i]
                    self.mms_fxn_name = sline[1].strip()
                except:
                    print "Error: On line {0}, Method of Manufactured Solutions function name not found".format(linenum)
                    raise ValueError()
            #Read in the method of manufactured solutions nodeset
            elif(len(line)<2):
                print "Error: on line {0}, the nodeset for the method of manufactured\n".format(linenum)+\
                      "       solutions mush have at least one node."
                raise ValueError()
            else:
                self.mms_dir_set = [sline[0]]+[int(i) for i in sline[1:]]
                
    #MMS functions
    def mms_const_u(self,coords):
        """Compute the method of manufactured solutions for a constant displacement of u
        with the phi values fixed at zero"""
        return .1,.2,.3,0.,0.,0.,0.,0.,0.,0.,0.,0.
        
    def mms_linear_u(self,coords):
        """Compute the method of manufactured solutions for a linear displacement of u 
        with the phi values fixed at zero"""
        a = .021
        b = .013
        c = .034
        return .1+a*coords[0], .2+b*coords[1],.3+c*coords[2],0.,0.,0.,0.,0.,0.,0.,0.,0.
                
class FEAModel(object):
    """The class which defines the FEA model"""
    
                      
    def __init__(self,IP_):
        """Initialize the finite element model"""
        self.input      = IP_                                          #Input of type InputParser
        self.total_ndof = len(self.input.nodes_coords)*self.input.ndof #Total number of degrees of freedom
        self.x0         = np.zeros([self.total_ndof])                  #Previous value of x in Ax = b
        self.xi         = np.zeros([self.total_ndof])                  #Current increment value of x in solver
        self.dxp        = np.zeros([self.total_ndof])                  #The previous change in x
        
        #Time bounds (TODO: should be in input file)
        self.ttot       = 0.3  #Total time
        self.t0         = 0.   #Previous value of time
        self.ti         = 0.1  #Current value of time
        self.dt         = 0.3  #The timestep
        
        self.nodes_dof  = []   #List of lists of degrees of freedom assocated with each node
        self.RHS        = None #The right hand side vector
        self.AMATRX     = None #The A matrix in Ax = b
        self.RHS_bc     = None #The right hand side vector with dirichlet bcs applied
        self.AMATRX_bc  = None #The A matrix with dirichlet bcs applied
        self.dbcdof     = None #The global degrees of freedom that have dirichlet boundary
                          #conditions applied and the value of the boundary condition
        self.delta_dbcdof = None #The change in dbcdof from this timestep to the previous
    
        self.maxiter    = 20   #The maximum number of iterations allowed in the Newton-Raphson solver
        self.atol       = 1e-8 #The absolute tolerance on the Newton-raphson iterations
        self.rtol       = 1e-5 #The relative tolerance on the Newton-Raphson iterations
        self.inc_number = 0    #The current increment number
    
        self.F          = None #The forcing function for the method of manufactured solutions
        self.alpha      = 1.0  #The relaxation parameter
    
    def __repr__(self):    
        """Define the repr string"""
        return "FEAModel({0})".format(self.input)
        
    def assemble_RHS_and_jacobian_matrix(self):
        """Assemble the global right hand side vector and jacobian matrix"""
        
        print self.total_ndof
        
        self.RHS    = np.zeros([self.total_ndof]) #Initialize the right hand side
        self.AMATRX = np.zeros([self.total_ndof,self.total_ndof]) #Initialize the A matrix
        
        print "=\n"+\
              "| Computing global stiffness matrix\n"+\
              "=";
        
        for e in self.input.element_nodes:
            print "=> Processing element: {0}".format(e[0])
            #Get the reference coordinates of the nodes
            coords = [[nc[1:] for nc in self.input.nodes_coords if nc[0]==e[i]][0] for i in range(1,9)]
            #Get the dof numbers associated with each node in the element
            edof = [[ndof[1:] for ndof in self.nodes_dof if ndof[0]==e[i]][0] for i in range(1,9)]
            flat_edof = flatten_list(edof)
            #Get the previous values of the dof
            dof0  = np.array([self.x0[d] for d in flat_edof])
            #Get the current values of the dof (U)
            U  = np.array([self.xi[d] for d in flat_edof])
            #print U
            #Create the DU vector
            DU = U-dof0
            #Call the element to provide the right hand side and jacobian matrix
            RHSe,AMATRXe = mel.UEL(self.input.properties,len(self.input.properties),\
                                   self.input.svars,len(self.input.svars),coords,None,U,DU,self.t0,self.ti-self.t0)
            
            #print "RHS:\n",RHSe
            #print "max RHS: {0}\nmin RHS: {1}".format(max(RHSe),min(RHSe))
            #ftot = np.zeros([3,])
            #mtot = np.zeros([9,])
            #for n in range(8):
            #    print "Node: {0}".format(n+1)
            #    print "Forces: {0}\nMoments: {1}".format(RHSe[n*12:(n*12+3)],RHSe[(n*12+3):((n+1)*12)])
            #    for i in range(3):
            #        ftot[i] += RHSe[n*12:(n*12+3)][i]
            #    for i in range(9):
            #        mtot[i] += RHSe[(n*12+3):((n+1)*12)][i]
            #print "Total Forces: {0}\nTotal Moments: {1}".format(ftot,mtot)
            
            #Assemble the global matrix
            for i in range(self.input.ndof*8):
                Iglob = flat_edof[i] #Get the mapping to the global degree of freedom index for i
                self.RHS[Iglob] += RHSe[i]
                for j in range(self.input.ndof*8):
                    Jglob = flat_edof[j] #Get the mapping to the global degree of freedom index for j
                    self.AMATRX[Iglob,Jglob] += AMATRXe[i,j]
        
        print "=\n"+\
              "| Completed global stiffness matrix formation\n"+\
              "=";
        #Include the method of manufactured solutions part if required  
        if(self.F!=None):
            print "=\n|=> Including method of manufactured solutions forcing function."
            for i in range(len(self.RHS)):
                self.RHS[i] += self.F[i]
        
    def form_increment_dof_vector(self):
        """Form the incremented dof vector at a given increment"""
        #Compute the scaling factor of the dirichlet boundary conditions
        factor = (self.ti-self.t0)/self.ttot
        #Copy the previous step to the current
        self.xi = np.copy(self.x0)
        for i in range(len(self.dbcdof)):
            #Set the change in the dirichlet boundary condition
            self.delta_dbcdof[i][1] = factor*self.dbcdof[i][1]
            #Add the change in dirichlet boundary conditions to the dof vector
            self.xi[self.dbcdof[i][0]] += factor*self.dbcdof[i][1]
        
    def initialize_dof(self):
        """Initialize the degrees of freedom for the FEA problem"""
        
        print "=\n"+\
              "|=> Initializing the degrees of freedom\n"+\
              "=";
        
        #Assign the dof
        self.assign_dof()
        #Compute the required values for the method of maufactured solutions if required
        if(self.input.mms_fxn!=None):
            self.apply_manufactured_solution()
        self.convert_local_dbc_to_global()
        
        print "=\n"+\
              "| Degrees of freedom intialized\n"+\
              "=";
              
    def initialize_timestep(self):
        """Initialize the timestep for the FEA problem"""
        #Create self.xi
        self.form_increment_dof_vector()
        self.inc_number += 1
        
    def assign_dof(self):
        "Assign dof to the nodes"
        for n in range(len(self.input.nodes_coords)):
            #Assign the number of local degrees of freedom required to each of the nodes
            self.nodes_dof = [[self.input.nodes_coords[n][0]]+range(n*self.input.ndof,(n+1)*self.input.ndof)\
                              for n in range(len(self.input.nodes_coords))]
        
    def convert_local_dbc_to_global(self):
        """Convert the locally defined dirichlet boundary conditions to their global values"""
        self.dbcdof = []
        for index,[node,ldof,val] in enumerate(self.input.dirichlet_bcs):
            #Find the global dof number from the node and local dof numbers
            fixed_dof = [l[ldof] for i,l in enumerate(self.nodes_dof) if l[0]==node]
            #Return an error if two nodes have the same number
            if(len(fixed_dof)==0):
                print "Error: Node does not exist for dirichlet boundary condition!"
                raise ValueError()
            elif(len(fixed_dof)>1):
                print "Error: Two nodes cannot have the same number!"
                raise ValueError()
            #Append the fixed dof to the list
            self.dbcdof.append([fixed_dof[0],val])
        #Remove any duplicates which may exist
        self.dbcdof = list(set([tuple(k) for k in self.dbcdof]))
        #Sort the boundary conditions by dof
        self.dbcdof.sort(key=lambda x: x[0])
        self.delta_dbcdof = [[dof[0],0.] for dof in self.dbcdof]
        
    def apply_dirichlet_bcs(self):
        """Apply the dirichlet boundary conditions to the system of equations"""
        
        #print self.dbcdof
        #print self.delta_dbcdof
        #raise
        
        for i,dbc in enumerate(self.delta_dbcdof):
            #Subtract the resulting force from the applied dof from the RHS
            self.RHS -= self.AMATRX[:,dbc[0]]*dbc[1]
            
        #Include the relaxation modification to the RHS
        #self.RHS  = (self.RHS - (1.-self.alpha)*np.dot(self.AMATRX,self.dxp))/self.alpha
        
        #Copy over the RHS and AMATRX terms to have the dirichlet bc terms removed
        self.RHS_bc    = np.copy(self.RHS)
        self.AMATRX_bc = np.copy(self.AMATRX)
        #Create reduced matrices for inversion
        for dbc in reversed(self.dbcdof):
            self.RHS_bc    = np.delete(self.RHS_bc,dbc[0])
            self.AMATRX_bc = np.delete(self.AMATRX_bc,dbc[0],axis=0)
            self.AMATRX_bc = np.delete(self.AMATRX_bc,dbc[0],axis=1)
            
    def update_increment(self):
        """update the currently defined increment and update xi"""
        self.assemble_RHS_and_jacobian_matrix()
        self.apply_dirichlet_bcs()
        #Use the numpy solver for speed
        #print self.RHS_bc
        dxi = np.linalg.solve(self.AMATRX_bc,self.RHS_bc)
        #print self.dbcdof
        for ddof in self.delta_dbcdof:
            dxi = np.insert(dxi,ddof[0],0.)#ddof[1])
        #print dxi
        self.xi = self.xi+dxi#self.x0+dxi #Add the change to the dof vector
        self.dxp = np.copy(dxi)
        #for n in range(len(self.nodes_dof)):
        #    temp = self.xi[(n*12):((n+1)*12)]
        #    print "Node: {0}\nu: {1}\nphi: {2}\n".format(self.nodes_dof[n],temp[:3],temp[3:])
        #raise
        
    def increment_solution(self):
        """Perform a Newton-Raphson iteration to find the solution"""
        self.initialize_timestep()
        self.update_increment()
        
        R  = np.linalg.norm(self.RHS_bc)
        R0 = R
        niter = 1
        
        print "===========================\n"+\
              " At increment:      {0}    \n".format(self.inc_number)+\
              " Iteration:         {0}    \n".format(niter)+\
              " Residual:          {0}    \n".format(R)+\
              " Relative Residual: {0}    \n".format(R/R0)+\
              "==========================="
        
        while ((R/R0)>self.rtol and R>self.atol and niter<self.maxiter):
            self.update_increment()
            R = np.linalg.norm(self.RHS_bc)
            
            print "===========================\n"+\
                  " At increment:      {0}    \n".format(self.inc_number)+\
                  " Iteration:         {0}    \n".format(niter)+\
                  " Residual:          {0}    \n".format(R)+\
                  " Relative Residual: {0}    \n".format(R/R0)+\
                  "==========================="
            niter +=1
            
        if((niter>=self.maxiter) and (R/R0>self.rtol) and R>self.atol):
            print "Error: Newton-Raphson solution did not converge for increment"
            return False
            
        print "=========================================\n"+\
              " At increment:      {0}                  \n".format(self.inc_number)+\
              " Solution Converged after {0} iterations \n".format(niter)+\
              "=========================================\n"
            
        return True
        
    def solve(self):
        """Solve the defined finite element problem"""
        #Initialize the degrees of freedom
        
        print "\n=================================================\n"+\
                "|                                               |\n"+\
                "|                BEGINNING SOLVER               |\n"+\
                "|                                               |\n"+\
                "=================================================\n"
        
        self.initialize_dof()
        self.ti = self.t0+self.dt
        
        while (self.t0<self.ttot):
            #Increment the solution at the timestep
            result = self.increment_solution()
            #Increment the timestep if solution converged, break otherwise
            if(result):
                self.t0 = self.ti
                self.ti = self.ti+self.dt
                if(self.ti>self.ttot):
                    self.ti = self.ttot
                self.x0 = np.copy(self.xi)
            else:
                if(self.input.mms_fxn!=None):
                    self.compare_manufactured_solution()
                return False
                
        print "\n=================================================\n"+\
                "|                                               |\n"+\
                "|                SOLVER COMPLETED               |\n"+\
                "|                                               |\n"+\
                "=================================================\n"
        
        if(self.input.mms_fxn!=None):
            self.compare_manufactured_solution()
        
        return True
        
    def compare_manufactured_solution(self):
        """Compare the manufactured solution to the computed solution"""
        mms = self.get_mms_dof_vector()
        
        if(np.allclose(mms,self.xi)):
            result = "Manufactured solution passed"
        else:
            result = "Error: Manufactured solution did not pass"
            
        #Write the test description file
        fname = os.path.join(self.input.path_to_file,"description.tex")
        if(os.path.isfile(fname)):
            os.remove(fname)
        fout = open(fname,"w+")
        fout.write(self.input.latex_string)
        fout.close()
        
        #Write the results of the test
        fname = os.path.join(self.input.path_to_file,"results.tex")
        if(os.path.isfile(fname)):
            os.remove(fname)
        fout = open(fname,"w+")
        fout.write(result)
        fout.write("\nManufactured Solution:\n"+str(mms)+"\n")
        fout.write("\n\FEA Solution:\n"+str(self.xi)+"\n")
        fout.write("\nDifference:\n"+str(mms-self.xi)+"\n")
        fout.close()
        
        print result
        
    def apply_manufactured_solution(self):
        """Apply the manufactured solution indicated"""
        print "=\n|=> Computing forcing function for method of manufactured solutions\n="
        self.compute_mms_forcing_function()
        self.compute_mms_bc_values()
        
    def get_mms_dof_vector(self):
        """Compute the displacement vector for the method of manufactured solutions"""
        U = np.zeros([12*len(self.input.nodes_coords)])
        for i,node in enumerate(self.input.nodes_coords):
            Us = self.input.mms_fxn(node[1:])
            for d in range(12):
                U[d+i*12] = Us[d]
        return U
        
    def compute_mms_forcing_function(self):
        """Compute the forcing function for the method of manufactured solutions"""
        
        #Get the manufactured solutions dof vector
        U = self.get_mms_dof_vector()
        #Set the degree of freedom to the computed U
        xi_old = np.copy(self.xi) #Copy initial value
        self.xi = U #Set the DOF vector to the computed U
        self.assemble_RHS_and_jacobian_matrix() #Compute the forcing function
        self.mms_F = np.copy(self.RHS) #Copy the right hand side to the forcing function value
        
        
    def compute_mms_bc_values(self):
        """Compute the boundary condition values for the method of manufactured solutions.
        Currently only works for dirichlet boundary conditions."""
        
        #Clear any dirichlet boundary conditions defined
        self.input.dirichlet_bcs = []
        for node in self.input.mms_dir_set[1:]:
            
            coords = [i[1:] for i in self.input.nodes_coords if i[0]==node]
            if(len(coords)==1):
                #Compute the applied displacements
                [ux,uy,uz,phi11,phi22,phi33,phi23,phi13,phi12,phi32,phi31,phi21] = self.input.mms_fxn(coords[0])
                #Apply the boundary condition values
                self.input.dirichlet_bcs.append([node,  1,    ux])
                self.input.dirichlet_bcs.append([node,  2,    uy])
                self.input.dirichlet_bcs.append([node,  3,    uz])
                self.input.dirichlet_bcs.append([node,  4, phi11])
                self.input.dirichlet_bcs.append([node,  5, phi22])
                self.input.dirichlet_bcs.append([node,  6, phi33])
                self.input.dirichlet_bcs.append([node,  7, phi23])
                self.input.dirichlet_bcs.append([node,  8, phi13])
                self.input.dirichlet_bcs.append([node,  9, phi12])
                self.input.dirichlet_bcs.append([node, 10, phi32])
                self.input.dirichlet_bcs.append([node, 11, phi31])
                self.input.dirichlet_bcs.append([node, 12, phi21])
            else:
                print "Error: Two nodes have the same number"
                raise ValueError()
        
        
        
def run_finite_element_model(input_filename):
    """Run the finite element model identified by the input filename"""
    IP = InputParser(input_filename)
    IP.read_input()
    FM = FEAModel(IP)
    FM.solve()
    
def profile_solve():
    """Solver run to profile"""
    IP = InputParser("unittest.inp")
    IP.read_input()
    FE = FEAModel(IP)
    FE.solve()
    
def flatten_list(l):
    """Flatten list of lists"""
    return [item for sublist in l for item in sublist]
    
class TestFEA(unittest.TestCase):

    f                    = None
    original_directory   = ""
    module_name           = "fea_driver"
    output_file_name      = r"results.tex".format(module_name)
    output_file_location  = r".\tests\unittests\{0}".format(module_name)
    currentResult         = None
    @classmethod
    def setUpClass(self):
        """Setup method"""
        #Define the results output format
        output_file = os.path.join(self.output_file_location,self.output_file_name)
        
        if(not os.path.isdir(self.output_file_location)):
            os.makedirs(self.output_file_location)
        
        #Write the description output
        description_file = os.path.join(self.output_file_location,r"description.tex")
        
        if(os.path.isfile(description_file)):
            os.remove(description_file)
        
        df = open(description_file,'w+')
        description_string = r"Unit tests of the \verb|{0}.py| module.".format(self.module_name)+"\n"
        df.write(description_string)
        df.close()
        
        #Write the results output
        if(os.path.isfile(output_file)):
            os.remove(output_file)
        self.f = open(output_file,"w+")
        table_open = r"\begin{table}[htb!]"+"\n"+\
                     r"\centering" +"\n"+\
                     r"\begin{tabular}{|l|c|}"+"\n"+\
                     r"\hline"+"\n"+\
                     r"module name & status\\"+"\n"+\
                     r"\hline"+"\n"+\
                     r"\hline"+"\n"
        print table_open
        self.f.write(table_open)
    @classmethod
    def tearDownClass(self):
        """Teardown method"""
        table_close = r"\hline"+"\n"+r"\end{tabular}"+"\n"+\
                      r"\end{table}" +"\n"+\
                      r"\FloatBarrier" + "\n"
        self.f.write(table_close)
        self.f.close()
        
    def setUp(self):
        pass
        
    def tearDown(self):
        ok = self.currentResult.wasSuccessful()
        tname = self.id().split(".")[-1].replace("_","\_")
        if(ok):
            str_out = r"\cellcolor{green!25} PASS"
        else:
            str_out = r"\cellcolor{red!25} FAIL"
        
        self.f.write(tname+"\t&\t"+str_out+r"\\"+"\n")
        
    def run(self, result=None):
        """Redefine run to keep track of results"""
        self.currentResult = result
        unittest.TestCase.run(self,result)
        
    def _test_read_input(self):
        """Test the read input command"""
        IP = InputParser("unittest.inp")
        IP.read_input()
        
    def test_fea_solver(self):
        """Perform test of the FEA solver"""
        IP = InputParser("unittest.inp")
        IP.read_input()
        FE = FEAModel(IP)
        self.assertTrue(FE.solve())
if __name__ == '__main__':
    unittest.main()
        