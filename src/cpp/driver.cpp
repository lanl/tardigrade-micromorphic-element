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
    
    while(found != std::string::npos){
        split_str.push_back(str.substr(found_p,found)); //Get the substring
        found_p = found;                                //set found_p to found
        found = str.find(delimiter,found_p+1);          //Find the next delimiter
    }
    
    split_str.push_back(str.substr(found_p));
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
    }
    
    NodeSet(std::string _name, std::vector< unsigned int> _nodes){
        /*!Full constructor*/
        name           = _name;
        nodes          = _nodes;
    }
};
  
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
            std::string filename;                                   //!The location and filename of the input deck
            std::string path_to_file;                               //!String giving the path to the file
            
            std::string latex_string;                               //!The description of the input deck
            std::vector< Node > nodes;                              //!List of the nodes in the finite element model
                                                                    //!and their coordinates
            std::vector< DirichletBC > dirichlet_bcs;               //!A list of the dirichlet boundary conditions
            
            unsigned int node_dof;                                  //!The number of degrees of freedom at a node
            std::vector< std::vector< unsigned int > > elements;    //!A list of the elements as defined by their nodes in the model
            std::vector< float > fprops;                            //!The floating point properties of the material model
            std::vector< std::vector< float > > svars;              //!A list of the state variables for each element
            std::vector< NodeSet > nodesets;                        //!A list of the nodesets
            
            std::string mms_name = "";                              //!The name of the manufactured solution being tested
            void (*mms_fxn)(unsigned int, std::string);             //!The function to compute U values for the method of 
                                                                    //!manufactured solutions.
                                                                    
            std::vector< std::string > keywords = {"*NODES", "*DIRICHLET_BCS", "*ELEMENTS", "*PROPERTIES", "*LATEX", "*NSET", "*MMS"}; //!All of the currently defined keywords
            void (InputParser::* keyword_fxn)(unsigned int, std::string);         //!The keyword processing function
            
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
            void process_keyword(unsigned int line_number, std::string line){
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
                }
                else if (line.find("*ELEMENTS") != std::string::npos){
                    std::cout << "Keyword *ELEMENTS found\n";
                }
                else if (line.find("*PROPERTIES") != std::string::npos){
                    std::cout << "Keyword *PROPERTIES found\n";
                }
                else if (line.find("*NSET") != std::string::npos){
                    std::cout << "Keyword *NSET found\n";
                }
                else if (line.find("*MMS") != std::string::npos){
                    std::cout << "Keyword *MMS found\n";
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
                    for(int i=0; i<split_line.size(); i++){std::cout << "split_line["<<i<<"]: " << split_line[i] << "\n";}
                }
            }
            
            
};

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