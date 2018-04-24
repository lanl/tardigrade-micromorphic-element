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
=====================================================================
*/

#include<vector>
#include<deformation_measures.h>

namespace micromorphic_material_library{

    class MicromorphicMaterial{
        /*!
        ==============================
        |    MicromorphicMaterial    |
        ==============================

        The base micromorphic material which defines the 
        fundamental functions. 

        */

        public:
            MicromorphicMaterialLibraryProxyBase();

            //Expose the functions we will require each material model to have
            virtual void evaluate_stress(const double &t, const double &dt, std::vector<double> (&fparams),
                                         )
            virtual void evaluate_stress_jacobian() 

    };

    class MaterialLibrary{
        /*!
        =================
        |    Library    |
        =================

        The class which contains all of the registered material models.

        */

        std::vector<std::string>          model_names;
        std::vector<MicromorphicMaterial> models;

        public:
            
            //Constructors
            MaterialLibrary(int n = 0){
                /*!
                =========================
                |    MaterialLibrary    |
                =========================

                The constructor for MaterialLibrary.

                n: The number of material models.

                */

                model_names.resize(n);
                models.resize(n);
            }

            //Methods
            void add_material(const std::string &name, MaterialLibrary &model){
                /*!
                ======================
                |    add_material    |
                ======================

                Add a material to the library.

                */

                model_names.pushback(name);
                models.pushback(model);
            }

            void get_material(const std::string &name, MaterialLibrary &model){
                /*!
                ======================
                |    get_material    |
                ======================

                Get the material model indicated by the model name

                */
                
                int model_num = -1;

                for (int i=0; i<model_names.size(); i++){
                    if(model_names[i].compare(name)){
                        model_num = 1;
                        model = models[i];
                        break;
                    }
                }
                if (model_num < 0){
                    std::cout << "Error: Material model not found.\n";
                    assert(1==0);
                }

            }

            MaterialLibrary library;

    };
}
