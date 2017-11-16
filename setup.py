import os
import sys
import shutil
import fileinput
import config

"""A script which runs the test suite and generates the various manuals.
 
See the report (.\doc\Report) for verification that installation was 
completed successfully"""
 
#Copy the bibliography file to the required LaTeX locations
bibliography_file       = "micromorphic.bib"
documentation_directory = "doc"
manual_locations = ["Report","TheoryManual","UsersManual","ProgrammersManual"]

#Define the absolute paths to the required libraries
eigen_location   = config.Configuration().eigen

#Define the compiler command
compiler_command = config.Configuration().compiler

#Define the make command
make_command     = config.Configuration().make

#Form the source path
source_path      = os.path.join(documentation_directory,bibliography_file)

for location in manual_locations:
    #Form the destination path
    destination_path = os.path.join(documentation_directory,location,bibliography_file)
    #Copy the micromorphic bibliography
    shutil.copyfile(source_path,destination_path)

#Find all of the makefiles in src/cpp
makefiles = []
for dirpath,_,filenames in os.walk("."):
    for f in filenames:
        if("makefile" in f.lower()):
            makefiles.append( os.path.abspath(os.path.join(dirpath,f)))

#Replace the keystring with the location of Eigen
for file in makefiles:
    output_file = file.replace(".default","")
    f = open(output_file,"w")
    for line in fileinput.input(file, inplace = True):
        line = line.replace("EIGEN_LOCATION",eigen_location)
        line = line.replace("COMPILER_COMMAND",compiler_command)
        #sys.stdout.write(line)
        f.write(line)
    f.close()

output_file = "./src/cpp/run_tests.py"
f = open(output_file,"w")
for line in fileinput.input("./src/cpp/run_tests.py.default", inplace=True):
    line = line.replace("MAKE_COMMAND",make_command)
    #sys.stdout.write(line)
    f.write(line)
f.close()
    
#os.system("python run_tests.py")
#os.system("python documentation_generator.py")
