An implementation of a micromorphic hex8 element. The form is 
intended to be used for easy implementation into an Abaqus 
user element (UEL).


Requires default python libraries (python.org) as well as numpy (numpy.org). These packages are 
available for download in a convenient package at www.continuum.io/downloads.

A full test of the code can be executed with the command >python test.py -v

where "python" is the call to the python 2.7 installation. This will run the unittest module which will execute the test functions in the code.
The -v option allows more information to be printed to the screen. This is useful to understand which tests are being executed.

Also requires an installation of LaTeX and Bibtex. It is assumed that the commands for these functions are pdflatex and bibtex respectively