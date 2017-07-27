import os
import sys
import glob

class TestHarness():
    """Class which defines the test harness"""
    
    def __init__(self):
        """Initialize the test harness"""
        self.orig_stdout = sys.stdout
        self.filename = "test_results.txt"
        os.remove(self.filename)
        self.f = open(self.filename,"w")
        sys.stdout = self.f
        self.owd   = os.getcwd()
        self.file_locations = [r".\src\python",r".\src\fortran"]
        self.fformat = "./*.py"
        
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
        
        os.chdir(os.path.join(self.owd,self.file_locations[0]))
        files = glob.glob(self.fformat)
        
        #Parse through files
        for f in files:
            print "\nTesting: {0}\n".format(f)
            os.system("python {0} -v".format(f))
            
        os.chdir(self.owd)
        
        print "======================================\n"+\
              "|                                    |\n"+\
              "|        Run regression tests        |\n"+\
              "|                                    |\n"+\
              "======================================\n"+\
              "| These tests are of higher level    |\n"+\
              "| functions which must be verified   |\n"+\
              "| against the code as a whole.       |\n"+\
              "======================================\n"
        
        os.chdir(os.path.join(self.owd,self.file_locations[1]))
        files = glob.glob(self.fformat)
        
        for f in files:
            print "\nTesting: {0}\n".format(f)
            os.system("python {0} -v".format(f))
            
        os.chdir(self.owd)
        
        
        
    def _flatten_list(self,l):
        """Flatten list of lists"""
        return [item for sublist in l for item in sublist]
        
if __name__ == '__main__':
    TestHarness().run_function_tests()