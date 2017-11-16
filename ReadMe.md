####################
# Micromorphic_UEL
####################
Developed by:

Nathan A. Miller (nathanm@lanl.gov)



Under the supervision of:

Dr. Richard Regueiro

Professor of Civil Engineering

The University of   
Colorado at Boulder



** Description **

An implementation of a micromorphic hex8 element. The form is 
intended to be used for easy implementation into an Abaqus 
user element (UEL). The current implementation is in python 
and the intent is to make sure the method is sound in this 
format before transitioning to either a Fortran or C++ code 
which will be used as a user element for Abaqus. At this time 
the C++ code appears the most advantageous due to the object 
oriented nature but the Fortran approach is more well trod.

** Upcoming work **

- Debug the python code finding the source of non-convergence
- Port the python code to C++
- Continue to improve the documentation

** Description of directories **

- .\doc: Where all of the genereated LaTeX documentation will be located. This includes 
       the main report of the testing of the code, the users manual, the programmers 
       manual and the theory manual. The resulting manuals are in PDF form.
       
- .\src: Where all of the source code for the micromorphic element lies. There are three
       subdirectories, python, cpp, and fortran which contain all python, C++, and 
       Fortran files required to compile the code.

** Python Requirements **

Requires default python libraries (python.org) as well as numpy (numpy.org). These packages are 
available for download in a convenient package at [continuum](www.continuum.io/downloads).

** LaTeX Requirements **

The documentation requires an installation of LaTeX and Bibtex. It is assumed that the commands for these 
functions are pdflatex and bibtex respectively

*** LaTeX packages ***

- Report/Users Manual/Programmers Manual:
    - \usepackage{listings, xcolor, subcaption, placeins}
    - \usepackage{undertilde}
    - \usepackage{algorithm,algpseudocode}
    - \usepackage{multicol}
    - \usepackage{makecell}
    - \usepackage[table]{colortbl}

- Theory Manual:
    - beamer theme Pittsburg
    - \usepackage[utf8]{inputenc}
    - \usepackage{amsmath}
    - \usepackage{amsfonts}
    - \usepackage{amssymb}
    - \usepackage{undertilde}
    - \usepackage{bm}
    - \usepackage{subcaption}

** C++ Compiler Requirements **

Requires the library [Eigen](http://eigen.tuxfamily.org) which is a collection of header files and does
not require any compilation. The user must define the path to this library in '''setup.py'''
('''eigen_location = /absolute/path/to/eigen''')

** Code Setup **

The code can be tested and the documentation generated by issuing the command >python setup.py
where "python" is the call to the python 2.7 installation. The results of the tests can be viewed in
the directory '''.\doc\Report\'''

A full test of the code can be executed with the command
> '''python run_tests.py -v'''

where "python" is the call to the python 2.7 installation. This will run the unittest module 
which will execute the test functions in the code. The -v option allows more information to 
be printed to the screen. This is useful to understand which tests are being executed.

** CPP Code **

Requires the GCC compiler (or other) though it defaults to gcc.

A local installation of GCC can be used by:

1. Downloading the gcc compiler and untaring it
2. Changing to the directory and running: '''./contrib/download_prerequisites'''
3. Running the command '''./configure --prefix=/absolute/path/to/install/directory'''
4. Running the command '''make'''
5. Before running a program compiled with this compiler one must set the environment 
   variable '''export LD_LIBRARY_PATH=/absolute/path/to/install/directory'''

Also requires the library [Eigen](http://eigen.tuxfamily.org) which requires that 
the path is defined in '''setup.py''' ('''eigen_location = /absolute/path/to/eigen''')