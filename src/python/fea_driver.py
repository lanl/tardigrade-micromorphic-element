import numpy as np
import unittest
import os

class InputParser:
    """The input parser class which reads in data from a text file"""
    
    filename      = "" #String giving the filename
    latex_string  = "" #The description of the input deck
    nodes_coords  = [] #List of node numbers and their coordinates
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
        keywords   = ["*NODES","*DOF","*ELEMENTS","*PROPERTIES","*LATEX"]
        parse_fxns = [self.parse_nodes, self.parse_dof, self.parse_elements, self.parse_props, self.parse_latex]
        
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
            if(len(sline)<4):
                print "Error: A node must be defined by its number,"+\
                      "       followed by its x,y,z coordinates."
                raise ValueError()
            else:
                self.nodes_coords.append([int(sline[0]),float(sline[1]),float(sline[2]),float(sline[3])])
        
    def parse_dof(self,line):
        """Parse the file when triggered by the dof keyword"""
        if(len(line)>0):
            sline = line.split(',')
            if(len(sline)<13):
                print "Error: A node must be defined by its number,"+\
                      "       followed by its 12 dof."
                raise ValueError()
            else:
                self.nodes_dof.append([int(sline[0]), float(sline[1]), float(sline[2]), float(sline[3]),\
                                                      float(sline[4]), float(sline[5]), float(sline[6]),\
                                                      float(sline[7]), float(sline[8]), float(sline[9]),\
                                                     float(sline[10]),float(sline[11]),float(sline[12])])
        
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
        
        
if __name__ == '__main__':
    unittest.main()
        