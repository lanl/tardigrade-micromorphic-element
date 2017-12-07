/*!==============================================
|                    uel.h                   |
==============================================
| Defintions of functions and other required |
| information for uel.cpp                    |
==============================================
*/

#include <iostream>
#include <Eigen/Dense>

extern "C" void UEL(double *RHS,double *AMATRX,double *SVARS,double *ENERGY,
                    int NDOFEL,int NRHS,int NSVARS,double *PROPS,int NPROPS,
                    double *COORDS,int MCRD,int NNODE,double *U,double *DU,
                    double *V,double *A,int JTYPE,double TIME[2],double DTIME,
                    int KSTEP,int KINC,int JELEM,double *PARAMS,int NDLOAD,
                    int *JDLTYP,double *ADLMAG,double *PREDEF,int NPREDF,
                    int *LFLAGS,int MLVARX,double *DDLMAG,int MDLOAD,
                    double PNEWDT,int *JPROPS,int NJPROP,double PERIOD);
                    
void compute_hex8(double *RHS,          double *AMATRX,     Vector &SVARS,  energy_vector &ENERGY,
                 Vector &PROPS,         Matrix_RM &COORDS,  Vector &U,      Vector &DU,
                 Vector &V,             Vector &A,          double TIME[2], double DTIME, 
                 int KSTEP,             int KINC,           int JELEM,      params_vector &PARAMS,
                 Matrixi_RM &JDLTYP,    Vector &ADLMAG,     double *PREDEF, int NPREDF,
                 lflags_vector &LFLAGS, Matrix_RM &DDLMAG,  double PNEWDT,  Vectori &JPROPS,
                 double PERIOD,         int NDOFEL,         int NRHS,       double *SVARS_ptr,
                 std::string output_fn);//,       std::ofstream&);
