/*!
=====================================================================
|                  micromorphic_material_library.h                  |
=====================================================================
| A header file which defines a class which registers all of the    |
| available micromorphic material models and their interface. The   |
| library is implemented in such a way that someone can register a  |
| new micromorphic constitutive model and have it available for use |
| merely by calling it by name. This allows us to re-use code as    |
| much as possible without needing to do time consuming rewrites of |
| already existing code.                                            |
---------------------------------------------------------------------
| Note: Registration approach taken from stackoverflow question     |
|       compile time plugin system 2                                |
=====================================================================
*/

#include "micromorphic_material_library.h"

namespace micromorphic_material_library {

    //IMaterial::~Imaterial(){}
    
    //IMaterialRegistrar::~IMaterialRegistrar(){}


    void IMaterial::evaluate_model_numeric_gradients(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                   const double (&grad_u)[3][3],           const double (&phi)[9],
                                   const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                   const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                   std::vector<double> &cauchy,                      std::vector<double> &s,                        std::vector<double> &m,
                                   std::vector<std::vector<double>> &DcauchyDgrad_u, std::vector<std::vector<double>> &DcauchyDphi, std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                   std::vector<std::vector<double>> &DsDgrad_u,      std::vector<std::vector<double>> &DsDphi,      std::vector<std::vector<double>> &DsDgrad_phi,
                                   std::vector<std::vector<double>> &DmDgrad_u,      std::vector<std::vector<double>> &DmDphi,      std::vector<std::vector<double>> &DmDgrad_phi,
                                   std::vector<std::vector<double>> &ADD_TERMS,      std::vector<std::vector<std::vector<double>>> &ADD_JACOBIANS, double delta){
        /*!
        ==========================================
        |    evaluate_model_numeric_gradients    |
        ==========================================

        Evaluates the model and computes the numeric gradients of the stress measures w.r.t the degrees of freedom and their gradients.

        TODO: Add gradients of the additional terms.

        */

        //Create temporary matrices and vectors
        Vector_9  _cauchy;
        Vector_9  _s;
        Vector_27 _m;

        Vector_9  _cauchy_p;
        Vector_9  _s_p;
        Vector_27 _m_p;

        Vector_9  _cauchy_n;
        Vector_9  _s_n;
        Vector_27 _m_n;

        std::vector<Eigen::VectorXd> _ADD_TERMS;

        //Evaluate the model at the reference point
        evaluate_model(time,    fparams,  grad_u, phi, grad_phi, SDVS, ADD_DOF, ADD_grad_DOF,
                       _cauchy, _s, _m, _ADD_TERMS);

        //Populate the stress outputs
        map_eigen_to_vector(_cauchy,cauchy);
        map_eigen_to_vector(_s,s);
        map_eigen_to_vector(_m,m);

        //Populate the additional terms
        ADD_TERMS.resize(_ADD_TERMS.size());
        for (unsigned int i=0; i<_ADD_TERMS.size(); i++){map_eigen_to_vector(_ADD_TERMS[i], ADD_TERMS[i]);}

        //Form the total vector
        std::vector<double> U(45);

        //Transfer the reference point for grad_u
        U[0] = grad_u[0][0];
        U[1] = grad_u[1][1];
        U[2] = grad_u[2][2];
        U[3] = grad_u[1][2];
        U[4] = grad_u[0][2];
        U[5] = grad_u[0][1];
        U[6] = grad_u[2][1];
        U[7] = grad_u[2][0];
        U[8] = grad_u[1][0];

        //Transfer the reference point for phi
        for (int i=0; i<9; i++){
            U[i+9] = phi[i];
        }

        //Transfer the reference point for grad_phi
        int sot_to_voigt_map[3][3]    = {{0,5,4},
                                         {8,1,3},
                                         {7,6,2}};
        int tot_to_voigt_map[3][3][3];
        deformation_measures::get_tot_to_voigt_map(tot_to_voigt_map);
        int voigt_to_tot_map[81];
        deformation_measures::get_voigt_to_tot_map(voigt_to_tot_map);

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                for (int k=0; k<3; k++){
                    U[tot_to_voigt_map[i][j][k]+18] = grad_phi[sot_to_voigt_map[i][j]][k];
                }
            }
        }

        double _grad_u[3][3];
        double _phi[9];
        double _grad_phi[9][3];

        std::vector<std::vector<double>> STRESS_MATRIX(45);
        for (int I=0; I<45; I++){STRESS_MATRIX[I].resize(45);}

        for (int J=0; J<45; J++){

            //Perturb the state in the positive direction
            U[J] += delta;

            //Extract the perturbed state
            _grad_u[0][0] = U[0];
            _grad_u[1][1] = U[1];
            _grad_u[2][2] = U[2];
            _grad_u[1][2] = U[3];
            _grad_u[0][2] = U[4];
            _grad_u[0][1] = U[5];
            _grad_u[2][1] = U[6];
            _grad_u[2][0] = U[7];
            _grad_u[1][0] = U[8];

            for (int i=0; i<9; i++){_phi[i] = U[i+9];}

            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    for (int k=0; k<3; k++){
                        _grad_phi[sot_to_voigt_map[i][j]][k] = U[tot_to_voigt_map[i][j][k]+18];
                    }
                }
            }

            //Evaluate the model (positive shift)
            evaluate_model(time,      fparams,  _grad_u, _phi, _grad_phi, SDVS, ADD_DOF, ADD_grad_DOF,
                           _cauchy_p, _s_p,     _m_p,    _ADD_TERMS);

            //Perturb the state in the positive direction
            U[J] -= 2*delta;

            //Extract the perturbed state
            _grad_u[0][0] = U[0];
            _grad_u[1][1] = U[1];
            _grad_u[2][2] = U[2];
            _grad_u[1][2] = U[3];
            _grad_u[0][2] = U[4];
            _grad_u[0][1] = U[5];
            _grad_u[2][1] = U[6];
            _grad_u[2][0] = U[7];
            _grad_u[1][0] = U[8];

            for (int i=0; i<9; i++){_phi[i] = U[i+9];}

            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    for (int k=0; k<3; k++){
                        _grad_phi[sot_to_voigt_map[i][j]][k] = U[tot_to_voigt_map[i][j][k]+18];
                    }
                }
            }

            //Evaluate the model (positive shift)
            evaluate_model(time,      fparams,  _grad_u, _phi, _grad_phi, SDVS, ADD_DOF, ADD_grad_DOF,
                           _cauchy_n, _s_n,     _m_n,    _ADD_TERMS);

            //Extract the response
            for (int I=0; I<9;  I++){STRESS_MATRIX[I   ][J] = 0.5*(_cauchy_p(I) - _cauchy_n[I])/delta;}
            for (int I=0; I<9;  I++){STRESS_MATRIX[I+9 ][J] = 0.5*(_s_p(I)      - _s_n[I])/delta;}
            for (int I=0; I<27; I++){STRESS_MATRIX[I+18][J] = 0.5*(_m_p(I)      - _m_n[I])/delta;}

            U[J] += delta;

        }

        //Output the matrix to the terminal
//        std::cout << "STRESS_MATRIX:\n";
//        for (int I=0; I<45; I++){
//            for (int J=0; J<45; J++){
//                std::cout << STRESS_MATRIX[I][J] << " ";
//            }
//            std::cout << "\n";
//        }

        //Extract the stress jacobians
        //std::cout << "dcauchydgrad_u\n";
        DcauchyDgrad_u.resize(9);
        for (int I=0; I<9; I++){
            DcauchyDgrad_u[I].resize(9);
            for (int J=0; J<9; J++){
                DcauchyDgrad_u[I][J] = STRESS_MATRIX[I][J];
            }
        }
        //std::cout << "dcauchydphi\n";
        DcauchyDphi.resize(9);
        for (int I=0; I<9; I++){
            DcauchyDphi[I].resize(9);
            for (int J=0; J<9; J++){
                DcauchyDphi[I][J] = STRESS_MATRIX[I][J+9];
            }
        }
        //std::cout << "dcauchydgrad_phi\n";
        DcauchyDgrad_phi.resize(9);
        for (int I=0; I<9; I++){
            DcauchyDgrad_phi[I].resize(27);
            for (int J=0; J<27; J++){
                DcauchyDgrad_phi[I][J] = STRESS_MATRIX[I][J+18];
            }
        }
        //std::cout << "dsdgrad_u\n";
        DsDgrad_u.resize(9);
        for (int I=0; I<9; I++){
            DsDgrad_u[I].resize(9);
            for (int J=0; J<9; J++){
                DsDgrad_u[I][J] = STRESS_MATRIX[I+9][J];
            }
        }
        //std::cout << "dsdphi\n";
        DsDphi.resize(9);
        for (int I=0; I<9; I++){
            DsDphi[I].resize(9);
            for (int J=0; J<9; J++){
                DsDphi[I][J] = STRESS_MATRIX[I+9][J+9];
            }
        }
        //std::cout << "dsdgrad_phi\n";
        DsDgrad_phi.resize(9);
        for (int I=0; I<9; I++){
            DsDgrad_phi[I].resize(27);
            for (int J=0; J<27; J++){
                DsDgrad_phi[I][J] = STRESS_MATRIX[I+9][J+18];
            }
        }
        //std::cout << "dmdgrad_u\n";
        DmDgrad_u.resize(27);
        for (int I=0; I<27; I++){
            DmDgrad_u[I].resize(9);
            for (int J=0; J<9; J++){
                DmDgrad_u[I][J] = STRESS_MATRIX[I+18][J];
            }
        }
        //std::cout << "dmdphi\n";
        DmDphi.resize(27);
        for (int I=0; I<27; I++){
            DmDphi[I].resize(9);
            for (int J=0; J<9; J++){
                DmDphi[I][J] = STRESS_MATRIX[I+18][J+9];
            }
        }
        //std::cout << "dmdgrad_phi\n";
        DmDgrad_phi.resize(27);
        for (int I=0; I<27; I++){
            DmDgrad_phi[I].resize(27);
            for (int J=0; J<27; J++){
                DmDgrad_phi[I][J] = STRESS_MATRIX[I+18][J+18];
            }
        }

        return;

    }

    void IMaterial::evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                   const double (&grad_u)[3][3],           const double (&phi)[9],
                                   const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                   const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                   std::vector<double> &cauchy,                      std::vector<double> &s,                        std::vector<double> &m,
                                   std::vector<std::vector<double>> &ADD_TERMS){

        /*!
        ========================
        |    evaluate_model    |
        ========================

        Copies values from vectors to eigen matrices, 
        submits to the Eigen::Matrix based function, 
        and then passes those values back.

        */

        //Create temporary matrices and vectors
        Vector_9  _cauchy;
        Vector_9  _s;
        Vector_27 _m;

        std::vector<Eigen::VectorXd> _ADD_TERMS;

        //Evaluate the model
        evaluate_model(time,    fparams,  grad_u, phi, grad_phi, SDVS, ADD_DOF, ADD_grad_DOF,
                       _cauchy, _s, _m, _ADD_TERMS);

        //Populate the stress outputs
        map_eigen_to_vector(_cauchy,cauchy);
        map_eigen_to_vector(_s,s);
        map_eigen_to_vector(_m,m);

        //Populate the additional terms
        ADD_TERMS.resize(_ADD_TERMS.size());
        for (unsigned int i=0; i<_ADD_TERMS.size(); i++){map_eigen_to_vector(_ADD_TERMS[i], ADD_TERMS[i]);}

    }

    void IMaterial::evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                   const double (&grad_u)[3][3],           const double (&phi)[9],
                                   const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                   const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                   std::vector<double> &cauchy,                      std::vector<double> &s,                        std::vector<double> &m,
                                   std::vector<std::vector<double>> &DcauchyDgrad_u, std::vector<std::vector<double>> &DcauchyDphi, std::vector<std::vector<double>> &DcauchyDgrad_phi,
                                   std::vector<std::vector<double>> &DsDgrad_u,      std::vector<std::vector<double>> &DsDphi,      std::vector<std::vector<double>> &DsDgrad_phi,
                                   std::vector<std::vector<double>> &DmDgrad_u,      std::vector<std::vector<double>> &DmDphi,      std::vector<std::vector<double>> &DmDgrad_phi,
                                   std::vector<std::vector<double>> &ADD_TERMS,      std::vector<std::vector<std::vector<double>>> &ADD_JACOBIANS){

        /*!
        ========================
        |    evaluate_model    |
        ========================

        Copies values from vectors to eigen matrices, 
        submits to the Eigen::Matrix based function, 
        and then passes those values back.

        */

        //std::cout << "In evaluate_model";

        //Create temporary matrices and vectors
        Vector_9  _cauchy;
        Vector_9  _s;
        Vector_27 _m;

        Matrix_9x9  _DcauchyDgrad_u;
        Matrix_9x9  _DcauchyDphi;
        Matrix_9x27 _DcauchyDgrad_phi;

        Matrix_9x9  _DsDgrad_u;
        Matrix_9x9  _DsDphi;
        Matrix_9x27 _DsDgrad_phi;

        Matrix_27x9  _DmDgrad_u;
        Matrix_27x9  _DmDphi;
        Matrix_27x27 _DmDgrad_phi;

        std::vector<Eigen::VectorXd> _ADD_TERMS;
        std::vector<Eigen::MatrixXd> _ADD_JACOBIANS;

        //Evaluate the model
        //std::cout << "temporary matrices created\n";
        //assert(-1==0);
        evaluate_model(time,    fparams,  grad_u, phi, grad_phi, SDVS, ADD_DOF, ADD_grad_DOF,
                       _cauchy,                      _s, _m,
                       _DcauchyDgrad_u, _DcauchyDphi, _DcauchyDgrad_phi,
                       _DsDgrad_u,      _DsDphi,      _DsDgrad_phi,
                       _DmDgrad_u,      _DmDphi,      _DmDgrad_phi,
                       _ADD_TERMS,      _ADD_JACOBIANS);

        //assert(0==1);

        //Populate the stress outputs
        map_eigen_to_vector(_cauchy,cauchy);
        map_eigen_to_vector(_s,s);
        map_eigen_to_vector(_m,m);

        //assert(1==2);

        //Populate the jacobian outputs
        map_eigen_to_vector(_DcauchyDgrad_u,   DcauchyDgrad_u);
        map_eigen_to_vector(_DcauchyDphi,      DcauchyDphi);
        map_eigen_to_vector(_DcauchyDgrad_phi, DcauchyDgrad_phi);

        //assert(2==3);

        map_eigen_to_vector(_DsDgrad_u,   DsDgrad_u);
        map_eigen_to_vector(_DsDphi,      DsDphi);
        map_eigen_to_vector(_DsDgrad_phi, DsDgrad_phi);

        //assert(3==4);

        map_eigen_to_vector(_DmDgrad_u,   DmDgrad_u);
        map_eigen_to_vector(_DmDphi,      DmDphi);
        map_eigen_to_vector(_DmDgrad_phi, DmDgrad_phi);

        //assert(4==5);
        ADD_TERMS.resize(_ADD_TERMS.size());
        for (unsigned int i=0; i<_ADD_TERMS.size(); i++){map_eigen_to_vector(_ADD_TERMS[i], ADD_TERMS[i]);}

        ADD_JACOBIANS.resize(_ADD_JACOBIANS.size());
        for (unsigned int i=0; i<_ADD_JACOBIANS.size(); i++){map_eigen_to_vector(_ADD_JACOBIANS[i], ADD_JACOBIANS[i]);}

        //assert(5==6);

        return;
    }

    void IMaterial::map_eigen_to_vector(const Vector_9 &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = V.size();
        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void IMaterial::map_eigen_to_vector(const Vector_27 &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = V.size();
        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void IMaterial::map_eigen_to_vector(const Eigen::VectorXd &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = V.size();
        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

   void IMaterial::map_eigen_to_vector(const Matrix_9x9 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void IMaterial::map_eigen_to_vector(const Matrix_9x27 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void IMaterial::map_eigen_to_vector(const Matrix_27x9 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void IMaterial::map_eigen_to_vector(const Matrix_27x27 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */
        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    void IMaterial::map_eigen_to_vector(const Eigen::MatrixXd &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        unsigned int A = M.rows();
        unsigned int B = M.cols();

        unsigned int a = v.size();

        if(a<A){v.resize(A);}

        for (unsigned int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (unsigned int j=0; j<B; j++){
                v[i][j] = M(i,j);
            }
        }

        return;
    }

    MaterialFactory& MaterialFactory::Instance() {
        static MaterialFactory instance;
        return instance;
    }

    void MaterialFactory::Register(IMaterialRegistrar* registrar, std::string name) {
        registry_[name] = registrar;
    }

    std::unique_ptr<IMaterial> MaterialFactory::GetMaterial(std::string name) {
        /* throws out_of_range if material unknown */
        IMaterialRegistrar* registrar = registry_.at(name);
        return registrar->GetMaterial();
    }

}
