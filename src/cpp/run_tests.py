###========================================
#|                                      |
#|           run_tests.py               |
#|                                      |
#========================================
#| Compiles and runs the tests of the   |
#| individual functions.                |
#|                                      |
#| These tests should output a LaTeX    |
#| file "results.tex" which contains    |
#| function/method name & True/False    |
#| where True indicates that the method |
#| or function passed the test and      |
#| False indicates the converse.        |
#========================================"""

import os
import sys
import subprocess

make_command     = 'mingw32-make'                         #The system command to run make
results_filename = "results.tex"                          #The results textfile

def print_char(num,char="#"):
    ###Print out the number of char symbols indicated by num
    str = ""
    for i in range(num):
        str += char
    
    print str
    
def print_bound_whitespace(num,bound_symbol = "#"):
    ###Print out the number of " " indicated by num bounded by bound_symbol
    str = bound_symbol
    for i in range(num):
        str += " "
        
    print str+bound_symbol

def print_test_header(dir):
    ###Print the header associated with the directory dir
    str = "#    Processing test in: " + os.path.join('tests',dir) + "    #"
    print ""
    print_char(len(str))
    print str
    print_char(len(str))
    print ""
    
def print_header():
    ###Print the test script header
    
    str = "#    Compiling and Executing tests located in \"tests\" subdirectory    #"
    
    print_char(len(str))
    print_bound_whitespace(len(str)-2)
    print str
    print_bound_whitespace(len(str)-2)
    print_char(len(str))
    
def process_results_file():
    ###Process the results for a unit test
    if(os.isfile(results_filename)):
        
        f = open(results_filename)
        
        
    else:
        str = "!!!    Error: Results file not found    !!!"
        print_char(len(str),"!")
        print str
        print_char(len(str),"!")
    
    
#Main code execution block
if(__name__=="__main__"):

    print_header()

    initial_directory = os.getcwd()                       #Get the current working directory

    sub_directories = [x for x in os.walk('tests')][0][1] #Get the subdirectories

    str = "#    Executing Unit Tests    #"
    print_char(len(str))
    print str
    print_char(len(str))
    
    for dir in sub_directories:                           #Iterate through the directories
        if(dir!='regression_tests'):
        
            print_test_header(dir)                             #Print the header
            os.chdir(os.path.join('tests',dir))           #Move to the directory
            exe_name = 'test_' + dir + ".exe"             #The resulting executable is assumed 
                                                      #to have the name test_directory.exe
            try:
                subprocess.call(make_command)                 #Call make in the current directory
            except:
                print "Error: Problem encountered during compilation."
        
            try:
                subprocess.call(exe_name)                     #Call the resulting executable
            except:
                print "Error: Problem encountered during execution."
            os.chdir(initial_directory)                   #Change back to the inital directory