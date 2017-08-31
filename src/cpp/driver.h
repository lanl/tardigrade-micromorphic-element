/*!=======================================================
  |                                                     |
  |                     driver.h                        |
  |                                                     |
  -------------------------------------------------------
  | The driver file to exercise the micromorphic finite |
  | element.                                            |
  =======================================================
  | Dependencies:                                       |
  | tensor:        The class which defines tensor       |
  |                access to an underlying Eigen        |
  |                matrix. This may result in a         |
  |                somewhat slower result however it    |
  |                should allow for a faster            |
  |                implementation time.                 |
  | micro_element: The definition of the micromorphic   |
  |                finite element. Currently, only a    |
  |                hexehedral element is defined.       |
  | micromorphic_linear_elasticity:                     |
  |                An implementation of the linear      |
  |                elastic model developed by Richard   |
  |                Regueiro (2010).                     |
  | newton_krylov: The definition of a newton-krylov    |
  |                solver. This solver is used as a     |
  |                proof of concept prior to the        |
  |                implementation of a more traditional |
  |                Newton-Raphson method.               |
  =======================================================*/
  
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <ctime>

std::string trim(const std::string& str, const std::string& whitespace = " \t");

std::vector< std::string > split(const std::string &str, const std::string& delimiter = ",");

class DirichletBC{
    /*!===
       |
       | D i r i c h l e t B C
       |
      ===
        
        A class for storing a dirichlet boundary condition 
        at a node.
        
    */
    
    public:
        unsigned int node_number; //!The node number at which the BC is applied
        int          dof_number;  //!The local degree of freedom which is constrained
        double       value;       //!The value of the boundary condition
        
    DirichletBC();
    
    DirichletBC(unsigned int _node_number, int _dof_number, double _value);
};

class Node{
    /*!  ===
       |
       | N o d e
       |
      ===
        
        A class which defines a node
        
    */
    
    public:
        unsigned int          number; //!The node number
        std::array<double,3>  coordinates; //!The coordinates of the node
        
    Node();
    
    Node(unsigned int _number, float x, float y, float z);
};

class Element{
    /*!  ===
       |
       | E l e m e n t
       |
      ===
        
        A class which defines an element
        
    */
    
    public:
        unsigned int               number; //!The element number
        std::array<unsigned int,8>  nodes; //!The nodes which make up the element
        
        Element();
        Element(unsigned int _number, unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4,
                                              unsigned int n5, unsigned int n6, unsigned int n7, unsigned int n8);
        
        Element(unsigned int _number, std::vector< unsigned int> _nodes);
};

class NodeSet{
    
    /*!  ===
       |
       | N o d e S e t
       |
      ===
        
        A class which defines a collection of nodes
        
    */
    
    public:
        std::string                 name;  //!The nodeset name
        std::vector< unsigned int > nodes; //!The node numbers
        
    NodeSet();
    
    NodeSet(std::string _name, std::vector< unsigned int> _nodes);
};

std::vector< double > mms_const_u(std::array< double, 3 > coords, double t);

std::vector< double > mms_linear_u(std::array< double, 3 > coords, double t);

double vector_norm(const std::vector< double > &x);

class InputParser{
    /*!===
       |
       | I n p u t P a r s e r
       |
      ===
        
        The base class for the input parser class.
        
        This class allows the user to read in information 
        from an input deck and parse it for use in the 
        finite element model.
        
    */
        
        public:
            std::string filename;                                             //!The location and filename of the input deck
            std::string path_to_file;                                         //!String giving the path to the file
            
            std::string latex_string;                                         //!The description of the input deck
            std::vector< Node > nodes;                                        //!List of the nodes in the finite element model
                                                                              //!and their coordinates
            std::vector< DirichletBC > dirichlet_bcs;                         //!A list of the dirichlet boundary conditions
            
            unsigned int node_dof;                                            //!The number of degrees of freedom at a node
            std::vector< Element > elements;                                  //!A list of the elements as defined by their nodes in the model
            std::vector< double > fprops;                                     //!The floating point properties of the material model
            std::vector< int   > iprops;                                      //!The integer properties of the material model
            std::vector< std::vector< float > > svars;                        //!A list of the state variables for each element
            std::vector< NodeSet > nodesets;                                  //!A list of the nodesets
            
            std::string mms_name = "";                                        //!The name of the manufactured solution being tested
            unsigned int mms_dirichlet_set_number = -1;                       //!The number of the dirichlet nodeset name in the nodesets vector
            std::vector< double > (*mms_fxn)(std::array< double, 3 >, double) = NULL; //!The function to compute U values for the method of 
                                                                                      //!manufactured solutions.
            
            double total_time = 1.0;                                          //!The total time of the simulation
            double tp         = 0.0;                                          //!The previous time
            double t          = 0.0;                                          //!The current value of the time
            double dt         = 0.3;                                          //!The current timestep
            
            bool verbose = true;                                              //!The verbosity of the output
            void (InputParser::* keyword_fxn)(unsigned int, std::string);     //!The keyword processing function
            
            //!=
            //!|=> Constructors
            //!=
            
            InputParser();
            
            InputParser(std::string _filename);
            
            //!=
            //!|=> Methods
            //!=
            
            void read_input();
            
        private:
            //!Private attributes
            char comment = '#'; //!Character which indicates a comment
            char keyword = '*'; //!Character which indicates a keyword
            
            //!Private methods
            void process_keyword(const unsigned int &line_number, std::string &line);
            
            void default_function(unsigned int line_number, std::string line);
            
            void parse_latex(unsigned int line_number, std::string line);
            
            void parse_nodes(unsigned int line_number, std::string line);
            
            void parse_elements(unsigned int line_number, std::string line);
            
            void parse_properties(unsigned int line_number, std::string line);
            
            void parse_nodesets(unsigned int line_number, std::string line);
            
            void parse_manufactured_solution(unsigned int line_number, std::string line);
};

class FEAModel{
    /*!===
       |
       | F E A M o d e l
       |
      ===
        
        The base class for the finite element model class.
        
        This stores the solution and runs the model.
        
    */
    
    public:
        InputParser input;                        //!The input parser class associated with the simulation
        unsigned int total_ndof;                  //!The total number of degrees of freedom
        std::vector< double > up;                 //!The previous value of the degree of freedom vector
        std::vector< double > u;                  //!The current value of the degree of freedom vector
        std::vector< double > du;                 //!The change in the degree of freedom vector in the increment
        
        std::vector< std::vector< unsigned int > > internal_nodes_dof; //!The global degrees of freedom associated with the nodes
        std::vector< DirichletBC > dbcdof;        //!The dirichlet boundary conditions defined w.r.t. the global degrees of freedom
        
        std::vector< Element > mapped_elements;   //!The elements defined with internal node numbering
        std::vector< NodeSet > mapped_nodesets;   //!A list of the nodesets mapped to the internal node numbering
        std::vector< unsigned int > unbound_dof;  //!A list of all of the unbound degrees of freedom
        
        std::vector< double > RHS;                //!The right hand side vector
        
        int maxiter = 20;                         //!The maximum number of iterations allowed at each timestep
        double atol = 1e-9;                       //!The absolute tolerance on the solver
        double rtol = 1e-7;                       //!The relative tolerance on the solver
        double mms_tol = 1e-6;                    //!The tolerance on the method of manufactured solutions comparison
        unsigned int increment_number = 0;        //!The current increment number
        
        std::vector< double > F;                  //!The forcing function for the method of manufactured solutions
        std::vector< double > mms_u;              //!The manufactured solution.
        
        double alpha = 1.0;                       //!The relaxation parameter
        std::string solver = "NewtonKrylov";      //!The solution to use
    
    FEAModel();
    
    FEAModel(InputParser _input);
    
    void map_element_nodes();
    
    void map_nodesets();
    
    /*!=
    |=> Degrees of freedom methods
    =*/
    
    void convert_local_dbc_to_global();
    
    void form_increment_dof_vector();
    
    void initialize_dof();
    
    void assign_dof();
    
    void id_unbound_dof();
    
    /*!=
    |=> Solver methods
    =*/
    
    bool solve();
    
    bool increment_solution();
    
    void initialize_timestep();
    
    void update_increment();
    
    void assemble_RHS_and_jacobian_matrix();
    
    void run_newton_krylov();
    
    std::vector< double > krylov_residual(std::vector<double> _du);
    
    std::vector< double > get_unbound_du();
    
    /*!=
    |=> Manufactured solutions methods
    =*/
    
    void apply_manufactured_solution();
    
    void compute_mms_forcing_function();
    
    void set_mms_dof_vector();
    
    void compute_mms_bc_values();
    
    void compare_manufactured_solution();
};

typedef std::vector<double> (FEAModel::*residual_function)(const std::vector<double>&); //!Type definition for a newton-krylov residual function
    
class KrylovSolver: public krylov::Solver{
    /*Define a Newton-Krylov solver for MPM functions*/
        
    public:
        
        //Overload residual function R_fxn
        FEAModel *model;
        residual_function R_fxn;
        
        KrylovSolver():krylov::Solver(){
            /*Constructor for blank Krylov Solver*/
        }
            
        KrylovSolver(FEAModel *_model, std::vector< double > _u,
            unsigned int _imax, unsigned int _kmax, bool _verbose){
            /*Initializer*/
            model   = _model;
            imax    = _imax;
            kmax    = _kmax;
            verbose = _verbose;
        
            u       = _u;
            norm_u  = vector_norm(u);
        
            if(u.size()<kmax){
                std::cout << "\nWarning: Reducing gkmax to size of u\n";
                kmax = u.size();
            }
        }
};

class FEAKrylovSolver: public KrylovSolver{
    /*!=========================
    |    FEAKrylovSolver    |
    =========================
    
    The newton krylov solver for the 
    finite element model.
    
    */
    public:
        FEAKrylovSolver(FEAModel &model, std::vector< double > _u,
                        unsigned int _imax, unsigned int _kmax, bool _verbose):KrylovSolver(&model, _u, _imax, _kmax, _verbose){};
                    
        std::vector< double > get_residual(std::vector< double > du){
            /*!Redefine the get_residual method to use the desired method of model*/
            return model->krylov_residual(du);
        }
};