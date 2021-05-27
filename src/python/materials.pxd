from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

import numpy as np
cimport numpy as np

# Function definitions

cdef extern from "material_python_interface.h" namespace "materialPythonInterface":

    int evaluate_model(const string &model_name,
                       const vector[double] &time, const vector[double] &fparams,
                       const double ( &current_grad_u )[ 3 ][ 3 ], const double ( &current_phi )[ 9 ], const double ( &current_grad_phi )[ 9 ][ 3 ],
                       const double ( &previous_grad_u )[ 3 ][ 3 ], const double ( &previous_phi )[ 9 ], const double ( &previous_grad_phi )[ 9 ][ 3 ],
                       vector[double] &SDVS,
                       const vector[double] &current_ADD_DOF,
                       const vector[vector[double]] &current_ADD_grad_DOF,
                       const vector[double] &previous_ADD_DOF,
                       const vector[vector[double]] &previous_ADD_grad_DOF,
                       vector[double] &PK2, vector[double] &SIGMA, vector[double] &M,
                       vector[vector[double]] &DPK2Dgrad_u,   vector[vector[double]] &DPK2Dphi,   vector[vector[double]] &DPK2Dgrad_phi,
                       vector[vector[double]] &DSIGMADgrad_u, vector[vector[double]] &DSIGMADphi, vector[vector[double]] &DSIGMADgrad_phi,
                       vector[vector[double]] &DMDgrad_u,     vector[vector[double]] &DMDphi,     vector[vector[double]] &DMDgrad_phi,
                       vector[vector[double]] &ADD_TERMS,     vector[vector[vector[double]]] &ADD_JACOBIANS,
                       string &output_message);
