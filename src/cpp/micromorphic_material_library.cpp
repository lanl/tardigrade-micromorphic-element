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