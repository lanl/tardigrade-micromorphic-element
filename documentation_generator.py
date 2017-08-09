import os
import sys

pdflatex_command = "pdflatex"
bibtex_command   = "bibtex"

def generate_manuals():
    """Generate the Theory, Programmers, and Users manuals"""
    #Save the current directory
    original_directory = os.getcwd()
    
    #Generate the manuals
    generate_theory_manual(original_directory)
    generate_users_manual(original_directory)
    
    #Return to the original directory
    os.chdir(original_directory)
    
def generate_theory_manual(original_directory):
    #Change directories
    os.chdir("./doc/TheoryManual")
    
    #Generate the theory manual
    os.system(r"{0} TheoryManual.tex".format(pdflatex_command))
    os.system(r"{0} TheoryManual.tex".format(pdflatex_command))
    
    os.chdir(original_directory)
    return
    
def generate_users_manual(original_directory):
    #Change directories
    os.chdir("./doc/UsersManual")
    
    #Generate the theory manual
    os.system(r"{0} UsersManual.tex".format(pdflatex_command))
    os.system(r"{0} UsersManual.tex".format(pdflatex_command))
    
    os.chdir(original_directory)
    return
    
if __name__ == '__main__':
    generate_manuals()