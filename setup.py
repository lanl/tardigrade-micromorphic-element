import os
import sys

"""A script which runs the test suite and generates the various manuals.
 
See the report (.\doc\Report) for verification that installation was 
completed successfully"""
 
os.system("python test.py")
os.system("python documentation_generator.py")
