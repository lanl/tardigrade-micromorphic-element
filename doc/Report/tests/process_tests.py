import os
import sys
import shutil

"""===================================================
| process_test.py                                   |
| A collection of routines which take results from  |
| a test and then generate the required LaTeX       |
| documentation.                                    |
|                                                   |
| The expected files are a main LaTeX file which is |
| generated from the input deck:                    |
| description.tex : describes the test              |
| results.tex     : contains the results of the     |
|                   test in LaTeX format as they    |
|                   should be put into the report.  |
| Any other files associated with the test should   |
| be stored in the test directory.                  |
|                                                   |
| All unit tests should be put into the report with |
| a table which lists the function being tested and |
| whether it passed or failed. The cells should be  |
| colored green for pass and red for fail using the |
| command \cellcolor{color!25} where color is       |
| red or green                                      |
====================================================="""

class LatexParser(object):
    """A class which parses the indicated LaTeX files from 
    the tests and adds them to the Report compile chain"""
    
    report_location = r"../.."                     #The location of the report
    unittest_dir    = r"unittests"                 #The location of the unittest LaTeX dir
    regression_dir  = r"regression_tests"          #The location of the regression LaTeX dir
    test_dir        = r"../../../src/python/tests" #The location of the test files
    
    unittest_file        = open("../tex/unittest_results.tex","w+")        #The LaTeX results of the unit tests
    regression_test_file = open("../tex/regression_test_results.tex","w+") #The LaTeX results of the regression tests
    
    def __init__(self):
        """Initialization of the LatexParser class"""
        
        self.get_test_names()
        self.copy_tex_files()
        self.form_test_section()
        
        
    def get_test_names(self):
        """Get the test names required for the LaTeX report"""
        self.unittests        = [t[1] for t in os.walk(self.unittest_dir)][0]
        self.regression_tests = [t[1] for t in os.walk(self.regression_dir)][0]
        
    def copy_tex_files(self):
        """Copy the required LaTeX files from the test directory to the LaTeX directory"""
        
        #copy the unit tests
        unittest_path  = os.path.join(self.test_dir,self.unittest_dir)
                
        for test_name in self.unittests:
            #Identify the file sources
            test_path = os.path.join(unittest_path,test_name)
            results_source_fname           = os.path.join(test_path,'results.tex')
            description_source_fname       = os.path.join(test_path,'description.tex')
            
            #Identify the file destinations
            latex_path = os.path.join(self.unittest_dir,test_name)
            results_destination_fname      = os.path.join(latex_path,'results.tex')
            description_destination_fname  = os.path.join(latex_path,'description.tex')
            
            #Copy the files
            shutil.copyfile(results_source_fname    , results_destination_fname)
            shutil.copyfile(description_source_fname, description_destination_fname)
            
        #Copy the regression tests
        regression_test_path  = os.path.join(self.test_dir,self.regression_dir)
        
        for test_name in self.regression_tests:
            #Identify the file sources
            test_path = os.path.join(regression_test_path,test_name)
            results_source_fname           = os.path.join(test_path,'results.tex')
            description_source_fname       = os.path.join(test_path,'description.tex')
            
            #Identify the file destinations
            latex_path = os.path.join(self.regression_dir,test_name)
            results_destination_fname      = os.path.join(latex_path,'results.tex')
            description_destination_fname  = os.path.join(latex_path,'description.tex')
            
            #Copy the files
            shutil.copyfile(results_source_fname    , results_destination_fname)
            shutil.copyfile(description_source_fname, description_destination_fname)
        
    def form_test_section(self):
        """Form the section of the LaTeX document that details the tests"""
        
        unittest_header        = "\section{Unit Test Results}\n\nThis section details the results of the function evaluations. Details on each simulation being evaluated are included in each subsection."
        regression_test_header = "\section{Regression Test Results}\n\n This section details the results of the regression tests. Details on each simulation being evaluated are included in each subsection."
        
        self.unittest_file.write(unittest_header)
        self.regression_test_file.write(regression_test_header)
        
        for module in self.unittests:
            path = os.path.join("./tests",self.unittest_dir,module)
            
            self.unittest_file.write("\n\n\\subsection{{{0}}}\n\\input{{{1}}}".format(module.replace("_","\_"),os.path.join(path,"description.tex").replace("\\","/")))
            self.unittest_file.write("\n\\input{{{1}}}\n".format(module,os.path.join(path,"results.tex").replace("\\","/")))
        
        for test in self.regression_tests:
            path = os.path.join("./tests",self.regression_dir,test)
        
            self.regression_test_file.write("\n\n\\subsection{{{0}}}\n\\input{{{1}}}".format(test.replace("_","\_"),os.path.join(path,"description.tex").replace("\\","/")))
            self.regression_test_file.write("\n\\input{{{1}}}\n\\clearpage".format(test,os.path.join(path,"results.tex").replace("\\","/")))
        
        self.unittest_file.close()
        self.regression_test_file.close()
if __name__ == '__main__':
    LP = LatexParser()