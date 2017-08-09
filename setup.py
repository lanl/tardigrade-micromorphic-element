import os
import sys
import shutil

"""A script which runs the test suite and generates the various manuals.
 
See the report (.\doc\Report) for verification that installation was 
completed successfully"""
 
#Copy the bibliography file to the required LaTeX locations
bibliography_file       = "micromorphic.bib"
documentation_directory = "doc"
manual_locations = ["Report","TheoryManual","UsersManual","ProgrammersManual"]

#Form the source path
source_path      = os.path.join(documentation_directory,bibliography_file)

for location in manual_locations:
    #Form the destination path
    destination_path = os.path.join(documentation_directory,location,bibliography_file)
    #Copy the micromorphic bibliography
    shutil.copyfile(source_path,destination_path)
    
os.system("python run_tests.py")
os.system("python documentation_generator.py")