import os
import sys
import glob

class TestHarness():
    """Class which defines the test harness"""
    
    def __init__(self):
        """Initialize the test harness"""
        self.orig_stdout = sys.stdout
        
        self.filename = "results.txt"
        if(os.path.isfile(self.filename)):
            os.remove(self.filename)
        self.f = open(self.filename,"w")
        sys.stdout = self.f
        self.owd   = os.getcwd()
        
        #Identify module locations
        self.module_location = r".\src\python"
        self.fformat = "./*.py"
        
        #Add the module locations to the system path for now
        sys.path.append(os.path.abspath(self.module_location))
        
        #Identify regression tests
        self.regression_tests = [r"tests/regression_tests/const_u/const_u.inp",\
                                 r"tests/regression_tests/linear_u/linear_u.inp"]
        
    def __del__(self):
        print "======================================\n"+\
              "|                                    |\n"+\
              "|          Tests Completed!          |\n"+\
              "|                                    |\n"+\
              "======================================\n"
        self.f.close()
        sys.stdout = self.orig_stdout
    
    def run_function_tests(self):
        """Run the tests associated with verifying the lower-level functions"""
        
        print "======================================\n"+\
              "|                                    |\n"+\
              "|         Run function tests         |\n"+\
              "|                                    |\n"+\
              "======================================\n"+\
              "| These tests are of lower level     |\n"+\
              "| functions which can be verified    |\n"+\
              "| best individually and depend only  |\n"+\
              "| a limit and well known set of      |\n"+\
              "| functions.                         |\n"+\
              "======================================\n"
        
        #Run the unit tests for modules located in the module location
        os.chdir(os.path.join(self.owd,self.module_location))
        files = glob.glob(self.fformat)
        
        #Run unit tests on modules
        for f in files:
            print "\nTesting: {0}\n".format(f)
            os.system("python {0} -v".format(f))

        print "======================================\n"+\
              "|                                    |\n"+\
              "|        Run regression tests        |\n"+\
              "|                                    |\n"+\
              "======================================\n"+\
              "| These tests are of higher level    |\n"+\
              "| functions which must be verified   |\n"+\
              "| against the code as a whole.       |\n"+\
              "======================================\n"
        
        import fea_driver
        for f in self.regression_tests:
            try:
                sys.stdout = self.orig_stdout                           #Change the output to the screen instead of the results file
                fea_driver.run_finite_element_model(os.path.abspath(f)) #Run the finite element model indicated
            except:
                print "!!! Error in execution !!!\nError: {0}".format(sys.exc_info()[0]) #Raise exception
                sys.stdout = self.f                                                      #Set the standard out back to the results file
                print "!!! Error in execution !!!\nError: {0}".format(sys.exc_info()[0]) #Repeat the error message into the file
        os.chdir(self.owd) #Reset the working directory
        
        
        
    def _flatten_list(self,l):
        """Flatten list of lists"""
        return [item for sublist in l for item in sublist]
        
if __name__ == '__main__':
    TestHarness().run_function_tests()