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

#Define the python command
python_command   = config.Configuration().python

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
        if(("makefile" in f.lower()) and ("default" in f.lower())):
            makefiles.append( os.path.abspath(os.path.join(dirpath,f)))

#Replace the keystring with the location of Eigen
for file in makefiles:
    output_file = file.replace(".default","")
    f = open(output_file,"w+")
    df = open(file,"r+")
    default = df.readlines()
    for line in default:
        line = line.replace("EIGEN_LOCATION",eigen_location)
        line = line.replace("COMPILER_COMMAND",compiler_command)
        f.write(line)
    df.close()
    f.close()

output_file = "./src/cpp/run_tests.py"
f = open(output_file,"w+")
df = open("./src/cpp/run_tests.py.default","r+")
default = df.readlines()
for line in default:
    line = line.replace("MAKE_COMMAND",make_command)
    f.write(line)
f.close()
df.close()

output_file = "run_tests.py"
f = open(output_file,"w+")
df = open("run_tests.py.default","r+")
default = df.readlines()
for line in default:
    line = line.replace("MAKE_COMMAND",make_command)
    line = line.replace("PYTHON_COMMAND",python_command)
    f.write(line)
f.close()
df.close()


#os.system("python run_tests.py")
#os.system("python documentation_generator.py")
