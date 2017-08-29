/*!=======================================================
  |                                                     |
  |                     driver.cpp                      |
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
  =======================================================*/
  
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <tensor.h>
#include <micro_element.h>
#include <micromorphic_linear_elasticity.h>
#include <ctime>

std::string trim(const std::string& str, const std::string& whitespace = " \t"){
    /*!==============
    |    trim    |
    ==============
    
    Remove the leading and trailing whitespace on a string
    
    */
    
    const auto strBegin = str.find_first_not_of(whitespace); //Find the first non-whitespace string
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace); //Find the last non-whitespace string
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange); //Return the trimmed string
}

std::vector< std::string > split(const std::string &str, const std::string& delimiter = ","){
    /*!===============
    |    split    |
    ===============
    
    Split the string into substrings
    defined by the delimiter.
    
    */
    
    std::vector< std::string > split_str(0);
    std::size_t found_p = 0;
    std::size_t found   = str.find(delimiter); //Search the string for the delimiter
   
    
    while((found != std::string::npos) && (found != found_p)){
        if(found_p>0){split_str.push_back(str.substr(found_p+1,found-found_p-1));}//Get the substring
        else{split_str.push_back(str.substr(found_p,found-found_p));}
        found_p = found;                                 //Set found_p to found
        found = str.find(delimiter,found_p+1);           //Find the next delimiter
    }
    
    split_str.push_back(str.substr(found_p+1));
    return split_str;
}

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
        
    DirichletBC(){
        /*!Default constructor*/
    }
    
    DirichletBC(unsigned int _node_number, int _dof_number, double _value){
        /*!Full constructor*/
        node_number = _node_number;
        dof_number  = _dof_number;
        value       = _value;
    }
};

class Node{
    /*!===
       |
       | N o d e
       |
      ===
        
        A class which defines a node
        
    */
    
    public:
        unsigned int         node_number; //!The node number
        std::array<float,3>  coordinates; //!The coordinates of the node
        
    Node(){
        /*!Default constructor*/
    }
    
    Node(unsigned int _node_number, float x, float y, float z){
        /*!Full constructor*/
        node_number    = _node_number;
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }
};

class Element{
    /*!===
       |
       | E l e m e n t
       |
      ===
        
        A class which defines an element
        
    */
    
    public:
        unsigned int         element_number; //!The element number
        std::array<unsigned int,8>  nodes; //!The nodes which make up the element
        
        Element(){
            /*!Default constructor*/
        }
        Element(unsigned int _element_number, unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4,
                                              unsigned int n5, unsigned int n6, unsigned int n7, unsigned int n8){
            /*!Full constructor*/
            element_number = _element_number;
            nodes[0] = n1;
            nodes[1] = n2;
            nodes[2] = n3;
            nodes[3] = n4;
            nodes[4] = n5;
            nodes[5] = n6; 
            nodes[6] = n7;
            nodes[7] = n8;
        }
};

class NodeSet{
    
    /*!===
       |
       | N o d e S e t
       |
      ===
        
        A class which defines a collection of nodes
        
    */
    
    public:
        std::string                 name;  //!The nodeset name
        std::vector< unsigned int > nodes; //!The node numbers
        
    NodeSet(){
        /*!Default constructor*/
        name = "";
        nodes.resize(0);
    }
    
    NodeSet(std::string _name, std::vector< unsigned int> _nodes){
        /*!Full constructor*/
        name           = _name;
        nodes          = _nodes;
    }
};

std::vector< double > mms_const_u(std::vector< double > coords){
    /*!=====================
    |    mms_const_u    |
    =====================
    
    Compute the method of manufactured solutions for a 
    constant displacement of u with the phi values 
    fixed at zero.
    
    */
    
    return {0.1, 0.2, 0.3, 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

std::vector< double > mms_linear_u(std::vector< double > coords){
    /*!======================
    |    mms_linear_u    |
    ======================
    
    Compute the method of manufactured solutions for a 
    linear displacement of u with the phi values 
    fixed at zero.
    
    */
    
    double a = 0.021;
    double b = 0.013;
    double c = 0.034;
    
    return {0.1+a*coords[0], 0.2+b*coords[1], 0.3+c*coords[2], 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

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
            std::string filename;                                           //!The location and filename of the input deck
            std::string path_to_file;                                       //!String giving the path to the file
            
            std::string latex_string;                                       //!The description of the input deck
            std::vector< Node > nodes;                                      //!List of the nodes in the finite element model
                                                                            //!and their coordinates
            std::vector< DirichletBC > dirichlet_bcs;                       //!A list of the dirichlet boundary conditions
            
            unsigned int node_dof;                                          //!The number of degrees of freedom at a node
            std::vector< Element > elements;                                //!A list of the elements as defined by their nodes in the model
            std::vector< float > fprops;                                    //!The floating point properties of the material model
            std::vector< std::vector< float > > svars;                      //!A list of the state variables for each element
            std::vector< NodeSet > nodesets;                                //!A list of the nodesets
            
            std::string mms_name = "";                                      //!The name of the manufactured solution being tested
            std::vector< double > (*mms_fxn)(std::vector< double >) = NULL; //!The function to compute U values for the method of 
                                                                            //!manufactured solutions.
            NodeSet mms_dirichlet_set;                                      //!The nodeset for the dirichlet boundary conditions for 
            
            double total_time = 1.0;                                        //!The total time of the simulation
            double t0         = 0.0;                                        //!The initial time
            double t          = 0.0;                                        //!The current value of the time
            double dt         = 0.3;                                        //!The current timestep
            
            bool verbose = 1;                                               //!The verbosity of the output
            void (InputParser::* keyword_fxn)(unsigned int, std::string);   //!The keyword processing function
            
            //!=
            //!|=> Constructors
            //!=
            
            InputParser(std::string _filename){
                filename    = _filename;
                keyword_fxn = &default_function;
                nodes.resize(0);
                elements.resize(0);
                fprops.resize(0);
                svars.resize(0);
                nodesets.resize(0);
            }
            
            //!=
            //!|=> Methods
            //!=
            
            void read_input(){
                /*!====================
                |    read_input    |
                ====================
                
                Read the input file (filename) and parse the 
                information into the class structure.
                
                */
                
                std::cout << "\n=================================================\n"<<
                               "|                                               |\n"<<
                               "|                  INPUT PARSER                 |\n"<<
                               "|                                               |\n"<<
                               "=================================================\n";
                
                //Initialize line string
                std::string line;
                //Initialize the line numbers
                unsigned int line_number=0;
                
                //Read in the input deck
                std::cout << "Reading data from " << filename << "\n";
                
                //Open the file
                std::ifstream f(filename);
                
                //Read the lines from the file
                if (f.is_open()){
                    while (std::getline(f,line)){
                        line = trim(line); //Remove leading and trailing whitespace
                        if(!(line[0]==comment) && (line.length()>0)){//Ignore comment lines
                            if(line[0]==keyword){
                                process_keyword(line_number,line);
                            }
                            (*this.*keyword_fxn)(line_number,line); //Call the current line parser
                        }
                        line_number++;
                    }
                }
                else{//Error if file cannot be opened
                    std::cout << "Error: Could not open file\n";
                }
                
                //Close the file
                f.close();
                
                std::cout << "\n=================================================\n"<<
                               "|                                               |\n"<<
                               "|             INPUT PARSER COMPLETED            |\n"<<
                               "|                                               |\n"<<
                               "=================================================\n";
                
            }
            
            
        private:
            //!Private attributes
            char comment = '#'; //!Character which indicates a comment
            char keyword = '*'; //!Character which indicates a keyword
            
            //!Private methods
            void process_keyword(const unsigned int &line_number, std::string &line){
                /*!=========================
                |    process_keyword    |
                =========================
                
                Process a line which has a keyword 
                defined in it.
                
                */
                
                if (line.find("*LATEX") != std::string::npos){
                    std::cout << "Keyword *LATEX found\n";
                    keyword_fxn = &parse_latex;
                    line.erase(line.begin(), line.begin()+6);
                }
                else if (line.find("*NODES") != std::string::npos){
                    std::cout << "Keyword *NODES found\n";
                    keyword_fxn = &parse_nodes;
                    line.erase(line.begin(), line.begin()+6);
                }
                else if (line.find("*DIRICHLET_BCS") != std::string::npos){
                    std::cout << "Keyword *DIRICHLET_BCS found\n";
                    std::cout << "Error: Not currently implemented.\n";
                    assert(1==0);
                }
                else if (line.find("*ELEMENTS") != std::string::npos){
                    std::cout << "Keyword *ELEMENTS found\n";
                    keyword_fxn = &parse_elements;
                    line.erase(line.begin(), line.begin()+9);
                }
                else if (line.find("*PROPERTIES") != std::string::npos){
                    std::cout << "Keyword *PROPERTIES found\n";
                    keyword_fxn = &parse_properties;
                    line.erase(line.begin(), line.begin()+11);
                }
                else if (line.find("*NSET") != std::string::npos){
                    std::cout << "Keyword *NSET found\n";
                    keyword_fxn = &parse_nodesets;
                    line.erase(line.begin(), line.begin()+5);
                }
                else if (line.find("*MMS") != std::string::npos){
                    std::cout << "Keyword *MMS found\n";
                    keyword_fxn = &parse_manufactured_solution;
                    line.erase(line.begin(), line.begin()+4);
                }
                else{
                    std::cout << "Error: Keyword not recognized\n";
                    assert(1==0);
                }
                
            }
            
            void default_function(unsigned int line_number, std::string line){
                /*!==========================
                |    default_function    |
                ==========================
                
                The default function which executes if a 
                line is read in prior to a keyword.
                
                */
                
                std::cout << "Warning: Line " << line_number << " is not associated\n"
                          << "         with a keyword. It will be ignored.\n";
            }
            
            void parse_latex(unsigned int line_number, std::string line){
                /*!=====================
                |    parse_latex    |
                =====================
                
                Parse the fine when triggered by a latex keyword
                
                Appends the line whole to latex_string
                
                input:
                    line_number: The number of the line (used primarily for error handling)
                    line:        The line read from the file
                
                */
                
                latex_string += line + "\n";
            }
            
            void parse_nodes(unsigned int line_number, std::string line){
                /*!=====================
                |    parse_nodes    |
                =====================
                
                Parse the fine when triggered by a node keyword
                
                Appends a new node to the node vector
                
                input:
                    line_number: The number of the line (used primarily for error handling)
                    line:        The line read from the file
                
                */
                
                //Initialize the split line
                std::vector< std::string> split_line;
                
                //Check if the line has length
                if(line.length()>0){
                    //Split the line at the commas
                    split_line =  split(line);
                    
                    if(split_line.size()==1){node_dof = std::strtoul(split_line[0].c_str(),NULL,10);} //Handle the case when the original line was *NODE,node_dof
                    else if(split_line.size()==4){
                        nodes.push_back(Node(std::strtoul(split_line[0].c_str(),NULL,10),
                                             std::strtod( split_line[1].c_str(),NULL),
                                             std::strtod( split_line[2].c_str(),NULL),
                                             std::strtod( split_line[3].c_str(),NULL))); //Convert to uint and double and create a node
                    }
                    else{
                        std::cout << "Error: On line " << line_number << ", a node must be defined by its number,"<<
                                     "       followed by its x,y,z coordinates.\n";
                        assert(1==0);
                    }
                    
                    if(verbose){
                        if(nodes.size()>0){
                            std::cout << "node number: " << nodes[nodes.size()-1].node_number;
                            std::cout << " coordinates: ";
                            for(int i=0; i<3; i++){std::cout << " " << nodes[nodes.size()-1].coordinates[i];}
                            std::cout << "\n";
                        }
                    }
                    
                }
            }
            
            void parse_elements(unsigned int line_number, std::string line){
                /*!========================
                |    parse_elements    |
                ========================
                
                Parse the fine when triggered by a element keyword
                
                Appends a new element to the element vector
                
                input:
                    line_number: The number of the line (used primarily for error handling)
                    line:        The line read from the file
                
                */
                
                //Initialize the split line
                std::vector< std::string> split_line;
                
                if(line.size()>0){ //Check if the line has useful information
                
                    split_line = split(line); //Split the line
                
                    //Only allow elements defined by the element number and eight nodes
                    if(split_line.size()<9){
                        std::cout << line;
                        std::cout << "Error: On line " << line_number << ", an element must be defined by the number and 8 nodes.";
                        assert(1==0);
                    }
                    else{
                        elements.push_back(Element(std::strtoul(split_line[0].c_str(),NULL,10),
                                                   std::strtoul(split_line[1].c_str(),NULL,10),
                                                   std::strtoul(split_line[2].c_str(),NULL,10),
                                                   std::strtoul(split_line[3].c_str(),NULL,10),
                                                   std::strtoul(split_line[4].c_str(),NULL,10),
                                                   std::strtoul(split_line[5].c_str(),NULL,10),
                                                   std::strtoul(split_line[6].c_str(),NULL,10),
                                                   std::strtoul(split_line[7].c_str(),NULL,10),
                                                   std::strtoul(split_line[8].c_str(),NULL,10)));
                    }
                    if(verbose){
                        if(elements.size()>0){
                            std::cout << "element number: " << elements[elements.size()-1].element_number;
                            std::cout << " nodes: ";
                            for(int i=0; i<8; i++){std::cout << " " << elements[elements.size()-1].nodes[i];}
                            std::cout << "\n";
                        }
                    }
                    
                }
                
            }
            
            void parse_properties(unsigned int line_number, std::string line){
                /*!==========================
                |    parse_properties    |
                ==========================
                
                Parse the properties when triggered by a element keyword
                
                input:
                    line_number: The number of the line (used primarily for error handling)
                    line:        The line read from the file
                
                */
                
                //Initialize the split line
                std::vector< std::string> split_line;
                
                if(line.size()>0){//Check if the line has useful information
                
                    split_line = split(line); //Split the line
                    
                    //Append the properties if they have not been previously defined
                    if(fprops.size()==0){
                        fprops.resize(split_line.size());
                        for(int i=0; i<split_line.size(); i++){fprops[i] = std::strtod( split_line[i].c_str(),NULL);}
                    }
                    else{
                        std::cout << "Error: On line " << line_number << ", properties can only be defined once.";
                        assert(1==0);
                    }
                    
                    if(verbose){
                        std::cout << "fprops: ";
                        for(int i=0; i<fprops.size(); i++){std::cout << " " << fprops[i];}
                        std::cout << "\n";
                    }
                }
            }
            
            void parse_nodesets(unsigned int line_number, std::string line){
                /*!========================
                |    parse_nodesets    |
                ========================
                
                Parse the line when triggered by a nodeset keyword
                
                Appends a new nodeset to the nodeset vector
                
                input:
                    line_number: The number of the line (used primarily for error handling)
                    line:        The line read from the file
                
                */
            
                //Initialize the split line
                std::vector< std::string> split_line;
            
                //Initialize the nodes
                std::vector< unsigned int > nset_nodes;
                
                if(line.size()>0){//Check if the line has useful information
                
                    split_line = split(line); //Split the line
                    
                    //Append the properties if they have not been previously defined
                    if(split_line.size()<2){
                        std::cout << "Error: On line " << line_number << ", a nodeset must have at least one node.";
                        assert(1==0);
                    }
                    else{
                        for(unsigned int i=1; i<split_line.size(); i++){nset_nodes.push_back(std::strtoul(split_line[i].c_str(),NULL,10));}
                            nodesets.push_back(NodeSet(trim(split_line[0]),nset_nodes)); 
                    }
                    
                    
                    if(verbose){
                        std::cout << "nodeset: " << nodesets[nodesets.size()-1].name <<"\n";
                        for(int i=0; i<nodesets[nodesets.size()-1].nodes.size(); i++){std::cout << " " << nodesets[nodesets.size()-1].nodes[i];}
                        std::cout << "\n";
                    }
                }
            }
            
            void parse_manufactured_solution(unsigned int line_number, std::string line){
                /*!=====================================
                |    parse_manufactured_solution    |
                =====================================
                
                Parse the line when triggered by a manufactured solutions keyword
                
                Defines the identified manufactured solution
                
                input:
                    line_number: The number of the line (used primarily for error handling)
                    line:        The line read from the file
                
                */
                
                //Initialize the split line
                std::vector< std::string> split_line;
                std::string fxn_name;
                std::vector< unsigned int > nset_nodes;
                
                if(line.length()>0){
                    
                    split_line = split(line); //Split the line
                    
                    if((split_line.size()==1) && (mms_fxn==NULL)){ //Try to find the manufactured solution function name
                        fxn_name = trim(split_line[0]); //Remove white space
                        if(fxn_name.compare("const_u")){mms_fxn = &mms_const_u;}
                        else if(fxn_name.compare("linear_u")){mms_fxn = &mms_linear_u;}
                        else{
                            std::cout << "Error: On line " << line_number << ", Method of Manufactured Solutions function name not found.";
                            assert(1==0);
                        }
                        if(verbose){
                            std::cout << "Manufactured Solution Function: " << fxn_name << "\n";
                        }
                    }
                    else if(split_line.size()<2){
                        std::cout << "Error: On line " << line_number << ", the nodeset for the method of manufactured\n"<<
                                     "       solutions must have at least one node.\n";
                        assert(1==0);
                    }
                    else{
                        for(unsigned int i=1; i<split_line.size(); i++){nset_nodes.push_back(std::strtoul(split_line[i].c_str(),NULL,10));}
                            mms_dirichlet_set = NodeSet(trim(split_line[0]),nset_nodes);
                            
                            if(verbose){
                                std::cout << "nodeset: " << mms_dirichlet_set.name << "\n";
                                for(int i=0; i<mms_dirichlet_set.nodes.size(); i++){std::cout << " " << mms_dirichlet_set.nodes[i];}
                                std::cout << "\n";
                            }
                            
                    }
                }
            }
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
        InputParser input;                       //!The input parser class associated with the simulation
        unsigned int total_ndof;                 //!The total number of degrees of freedom
        std::vector< double > up;                //!The previous value of the degree of freedom vector
        std::vector< double > u;                 //!The current value of the degree of freedom vector
        std::vector< double > du;                //!The change in the degree of freedom vector in the increment
        
        std::vector< std::vector< unsigned int > > internal_nodes_dof; //!The global degrees of freedom associated with the nodes
        
        std::vector< double > RHS;               //!The right hand side vector
        
        int maxiter = 20;                        //!The maximum number of iterations allowed at each timestep
        double atol = 1e-8;                      //!The absolute tolerance on the solver
        double rtol = 1e-6;                      //!The relative tolerance on the solver
        unsigned int increment_number = 0;       //!The current increment number
        
        std::vector< double > F;                 //!The forcing function for the method of manufactured solutions
        double alpha = 1.0;                      //!The relaxation parameter
    
}

int main( int argc, char *argv[] ){
    /*!===
       |
       | D r i v e r
       |
      ===
        
        The driver code which can solve finite element formulations using 
        micromorphic finite elements written for use in Abaqus UEL 
        subroutines.
        
    */
    
    if (argc != 2){ // We expect two arguments for use of the code
        std::cout << "usage: " << argv[0] << " <filename>\n";
    }
    else{
        // The first argument is assumed to be a filename to open
        InputParser IP(argv[1]);
        IP.read_input();
    }
        
}