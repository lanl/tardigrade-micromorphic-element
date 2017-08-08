import os
import sys

pdflatex_command = "pdflatex"
bibtex_command   = "bibtex"

def generate_manuals():
    """Generate the Theory, Programmers, and Users manuals"""
    #Save the current directory
    original_directory = os.getcwd()
    
    #Generate the manuals
    generate_theory_manual()
    
    #Return to the original directory
    os.chdir(original_directory)
    
def generate_theory_manual():
    #Change directories
    os.chdir("./doc/TheoryManual")
    
    #Generate the theory manual
    os.system(r"{0} TheoryManual.tex".format(pdflatex_command))
    os.system(r"{0} TheoryManual.tex".format(pdflatex_command))
    
    return
    
if __name__ == '__main__':
    generate_manuals()