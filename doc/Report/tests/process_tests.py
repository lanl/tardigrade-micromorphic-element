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
=====================================================


 