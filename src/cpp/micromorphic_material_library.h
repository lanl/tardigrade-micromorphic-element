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

#include<vector>
#include <list>
#include <string>
#include <map>
#include <memory>
#include<deformation_measures.h>

#ifndef MICROMORPHIC_MATERIAL_LIBRARY_H
#define MICROMORPHIC_MATERIAL_LIBRARY_H

namespace micromorphic_material_library{

    /* Base class for materials */
    class IMaterial {
    public:
    
        //virtual ~IMaterial();
        
        virtual void evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                    const double (&grad_u)[3][3],           const double (&phi)[9],
                                    const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                    const std::vector<double> &ADDDOF,      const std::vector<std::vector<double>> &ADD_grad_DOF,
                                    Vector_9 &cauchy, Vector_9 &s, Vector_27 &m, std::vector<Eigen::VectorXd> &ADD_TERMS) = 0;
                                    
        virtual void evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                    const double (&grad_u)[3][3],           const double (&phi)[9],
                                    const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                    const std::vector<double> &ADDDOF,      const std::vector<std::vector<double>> &ADD_grad_DOF,
                                    Vector_9    &cauchy,    Vector_9    &s,           Vector_27    &m,
                                    Matrix_9x9  &DcauchyDgrad_u, Matrix_9x9  &DcauchyDphi, Matrix_9x27  &DcauchyDgrad_phi,
                                    Matrix_9x9  &DsDgrad_u,      Matrix_9x9  &DsDphi,      Matrix_9x27  &DsDgrad_phi,
                                    Matrix_27x9 &DmDgrad_u,      Matrix_27x9 &DmDphi,      Matrix_27x27 &DmDgrad_phi,
                                    std::vector<Eigen::VectorXd> &ADD_TERMS,               std::vector<Eigen::MatrixXd> &ADD_JACOBIANS) = 0;
    };

    /* 
     * Base class for MaterialRegistrar
     * See MaterialRegistrar below for explanations
     */
    class IMaterialRegistrar {
    public:
        //virtual ~IMaterialRegistrar();
    
        virtual std::unique_ptr<IMaterial> GetMaterial() = 0;
    };

    /* 
     * This is the factory, the common interface to "materials".
     * Materials registers themselves here and the factory can serve them on
     * demand.
     * It is a Singleton
     */
    class MaterialFactory {
    public:
        /* Get Singleton instance */
        static MaterialFactory& Instance();
        /* Register a new material */
        void Register(IMaterialRegistrar* registrar, std::string name);
        /* Get an instance of a material based on its name */
        /* throws out_of_range if material not found */
        std::unique_ptr<IMaterial> GetMaterial(std::string name);

    private:
        /* Holds pointers to material registrars */
        std::map<std::string, IMaterialRegistrar*> registry_;
        /* Make constructors private and forbid cloning */
        MaterialFactory(): registry_() {};
        MaterialFactory(MaterialFactory const&) = delete;
        void operator=(MaterialFactory const&) = delete;
    };

    /* 
     * Helper class that registers a material upon construction.
     * Actually, the registrar registers itself, and the proxied material is only
     * created on-demand. This mechanism can be shortened by directly 
     * registering and instance of the material, but the assumption here is that
     * instanciating the material can be heavy and not necessary.
     */
    template<class TMaterial>
    class MaterialRegistrar: public IMaterialRegistrar {
    public:
        MaterialRegistrar(std::string classname);
        std::unique_ptr<IMaterial> GetMaterial();
    private:
        /* That is not really used there, but could be useful */
        std::string classname_;
    };

    /* template functions in header */

    template<class TMaterial>
    MaterialRegistrar<TMaterial>::MaterialRegistrar(std::string classname): classname_(classname) {
        MaterialFactory &factory = MaterialFactory::Instance();
        factory.Register(this, classname);
    }

    template<class TMaterial>
    std::unique_ptr<IMaterial>
    MaterialRegistrar<TMaterial>::GetMaterial() {
        std::unique_ptr<IMaterial> material(new TMaterial());
        return material;
    }
}

/*
 * Here is the trick: upon creation of the global variable, the class created
 * out of the template will get instanciated once, and will register itself.
 * The template contains the information to create a material instance.
 * An unnamed namespace is used to enclose this later unused variable in the
 * compilation unit.
 */
#define REGISTER_MATERIAL(CLASSNAME) \
    namespace { \
        static micromorphic_material_library::MaterialRegistrar<CLASSNAME> \
        _registrar( #CLASSNAME ); \
    };
    
#endif