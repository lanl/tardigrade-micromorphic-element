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
        for (int i=0; i<_ADD_TERMS.size(); i++){map_eigen_to_vector(_ADD_TERMS[i], ADD_TERMS[i]);}

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

        std::cout << "In evaluate_model";

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
        std::cout << "temporary matrices created\n";
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
        for (int i=0; i<_ADD_TERMS.size(); i++){map_eigen_to_vector(_ADD_TERMS[i], ADD_TERMS[i]);}

        ADD_JACOBIANS.resize(_ADD_JACOBIANS.size());
        for (int i=0; i<_ADD_JACOBIANS.size(); i++){map_eigen_to_vector(_ADD_JACOBIANS[i], ADD_JACOBIANS[i]);}

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

        int A = V.size();
        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void IMaterial::map_eigen_to_vector(const Vector_27 &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        int A = V.size();
        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

    void IMaterial::map_eigen_to_vector(const Eigen::VectorXd &V, std::vector<double> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        int A = V.size();
        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){v[i] = V[i];}

        return;
    }

   void IMaterial::map_eigen_to_vector(const Matrix_9x9 &M, std::vector<std::vector<double>> &v){
        /*!
        =============================
        |    map_eigen_to_vector    |
        =============================

        Map an eigen matrix to a standard vector.

        */

        int A = M.rows();
        int B = M.cols();

        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (int j=0; j<B; j++){
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

        int A = M.rows();
        int B = M.cols();

        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (int j=0; j<B; j++){
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

        int A = M.rows();
        int B = M.cols();

        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (int j=0; j<B; j++){
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
        int A = M.rows();
        int B = M.cols();

        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (int j=0; j<B; j++){
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

        int A = M.rows();
        int B = M.cols();

        int a = v.size();

        if(a<A){v.resize(A);}

        for (int i=0; i<A; i++){
            if(v[i].size()<B){v[i].resize(B);}
            for (int j=0; j<B; j++){
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
