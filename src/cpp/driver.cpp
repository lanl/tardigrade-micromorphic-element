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
#include <tensor.h>
#include <micro_element.h>
#include <micromorphic_linear_elasticity.h>
#include <newton_krylov.h>
#include <driver.h>
#include <ctime>
#include <algorithm>

std::string trim(const std::string& str, const std::string& whitespace){
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

std::vector< std::string > split(const std::string &str, const std::string& delimiter){
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

void print_vector(std::string name, std::vector< double > V){
    /*!======================
    |    print_vector    |
    ======================
    
    Print a vector to the screen
    
    */
    
    std::cout << name << ": ";
    
    for(int i=0; i<V.size(); i++){
        std::cout << V[i] << " ";
    }
    std::cout << "\n";
}

void print_vector_of_vectors(std::string name, std::vector< std::vector< double > > V){
    /*!=================================
    |    print_vector_of_vectors    |
    =================================
    
    Print a vector of vectors to the screen
    
    */
    
    for(int i=0; i<V.size(); i++){
        std::cout << name << "[" << i << "]: ";
        for(int j=0; j<V[i].size(); j++){
            std::cout << V[i][j] << " ";
        }
        std::cout << "\n";
    }
}

DirichletBC::DirichletBC(){
    /*!Default constructor*/
}
    
DirichletBC::DirichletBC(unsigned int _node_number, int _dof_number, double _value){
    /*!Full constructor*/
    node_number = _node_number;
    dof_number  = _dof_number;
    value       = _value;
}

Node::Node(){
    /*!Default constructor*/
}
    
Node::Node(unsigned int _number, float x, float y, float z){
    /*!Full constructor*/
    number    = _number;
    coordinates[0] = x;
    coordinates[1] = y;
    coordinates[2] = z;
}


Element::Element(){
    /*!Default constructor*/
}
Element::Element(unsigned int _number, unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4,
                                       unsigned int n5, unsigned int n6, unsigned int n7, unsigned int n8){
    /*!Full constructor*/
    number   = _number;
    nodes[0] = n1;
    nodes[1] = n2;
    nodes[2] = n3;
    nodes[3] = n4;
    nodes[4] = n5;
    nodes[5] = n6; 
    nodes[6] = n7;
    nodes[7] = n8;
}
        
Element::Element(unsigned int _number, std::vector< unsigned int> _nodes){
    /*!Constructor for vector form of node numbers*/
    number = _number;
    for(int i=0; i<8; i++){
        nodes[i] = _nodes[i];
    }
}

NodeSet::NodeSet(){
    /*!Default constructor*/
    name = "";
    nodes.resize(0);
}
    
NodeSet::NodeSet(std::string _name, std::vector< unsigned int> _nodes){
    /*!Full constructor*/
    name           = _name;
    nodes          = _nodes;
}

std::vector< double > mms_const_u(std::array< double, 3 > coords, double t){
    /*!=====================
    |    mms_const_u    |
    =====================
    
    Compute the method of manufactured solutions for a 
    constant displacement of u with the phi values 
    fixed at zero.
    
    */
    
    return {0.1*t, 0.2*t, 0.3*t, 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

std::vector< double > mms_linear_u(std::array< double, 3 > coords, double t){
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
    
    return {(0.1+a*coords[0])*t, (0.2+b*coords[1])*t, (0.3+c*coords[2])*t, 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

double vector_norm(const std::vector< double > &x){
    /*!=====================
    |    vector_norm    |
    =====================
    
    Compute the L2 norm of a vector.
    
    */
    
    double norm = 0;
    
    for(int i=0; i<x.size(); i++){
        norm += x[i]*x[i];
    }
    
    norm = sqrt(norm);
    
    return norm;
}

InputParser::InputParser(){
    /*Default constructor*/
    filename    = "";
    keyword_fxn = &InputParser::default_function;
    nodes.resize(0);
    elements.resize(0);
    fprops.resize(0);
    svars.resize(0);
    nodesets.resize(0);
}
            
InputParser::InputParser(std::string _filename){
    filename    = _filename;
    keyword_fxn = &InputParser::default_function;
    nodes.resize(0);
    elements.resize(0);
    fprops.resize(0);
    svars.resize(0);
    nodesets.resize(0);
}

//!=
//!|=> Methods
//!=
            
void InputParser::read_input(){
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
    return;
}

void InputParser::process_keyword(const unsigned int &line_number, std::string &line){
    /*!=========================
    |    process_keyword    |
    =========================
    
    Process a line which has a keyword 
    defined in it.
    
    */
    
    if (line.find("*LATEX") != std::string::npos){
        std::cout << "Keyword *LATEX found\n";
        keyword_fxn = &InputParser::parse_latex;
        line.erase(line.begin(), line.begin()+6);
    }
    else if (line.find("*NODES") != std::string::npos){
        std::cout << "Keyword *NODES found\n";
        keyword_fxn = &InputParser::parse_nodes;
        line.erase(line.begin(), line.begin()+6);
    }
    else if (line.find("*DIRICHLET_BCS") != std::string::npos){
        std::cout << "Keyword *DIRICHLET_BCS found\n";
        std::cout << "Error: Not currently implemented.\n";
        assert(1==0);
    }
    else if (line.find("*ELEMENTS") != std::string::npos){
        std::cout << "Keyword *ELEMENTS found\n";
        keyword_fxn = &InputParser::parse_elements;
        line.erase(line.begin(), line.begin()+9);
    }
    else if (line.find("*PROPERTIES") != std::string::npos){
        std::cout << "Keyword *PROPERTIES found\n";
        keyword_fxn = &InputParser::parse_properties;
        line.erase(line.begin(), line.begin()+11);
    }
    else if (line.find("*NSET") != std::string::npos){
        std::cout << "Keyword *NSET found\n";
        keyword_fxn = &InputParser::parse_nodesets;
        line.erase(line.begin(), line.begin()+5);
    }
    else if (line.find("*MMS") != std::string::npos){
        std::cout << "Keyword *MMS found\n";
        keyword_fxn = &InputParser::parse_manufactured_solution;
        line.erase(line.begin(), line.begin()+4);
    }
    else{
        std::cout << "Error: Keyword not recognized\n";
        assert(1==0);
    }
    return;
}

void InputParser::default_function(unsigned int line_number, std::string line){
    /*!==========================
    |    default_function    |
    ==========================
    
    The default function which executes if a 
    line is read in prior to a keyword.
    
    */
    
    std::cout << "Warning: Line " << line_number << " is not associated\n"
              << "         with a keyword. It will be ignored.\n";
    return;
}

void InputParser::parse_latex(unsigned int line_number, std::string line){
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
    return;
}

void InputParser::parse_nodes(unsigned int line_number, std::string line){
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
                std::cout << "node number: " << nodes[nodes.size()-1].number;
                std::cout << " coordinates: ";
                for(int i=0; i<3; i++){std::cout << " " << nodes[nodes.size()-1].coordinates[i];}
                std::cout << "\n";
            }
        } 
    }
    return;
}

void InputParser::parse_elements(unsigned int line_number, std::string line){
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
                std::cout << "element number: " << elements[elements.size()-1].number;
                std::cout << " nodes: ";
                for(int i=0; i<8; i++){std::cout << " " << elements[elements.size()-1].nodes[i];}
                    std::cout << "\n";
                }
            }
        }
}
            
void InputParser::parse_properties(unsigned int line_number, std::string line){
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
            
void InputParser::parse_nodesets(unsigned int line_number, std::string line){
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
    return;
}
            
void InputParser::parse_manufactured_solution(unsigned int line_number, std::string line){
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
            if(!fxn_name.compare("const_u")){mms_fxn = &mms_const_u;}
            else if(!fxn_name.compare("linear_u")){mms_fxn = &mms_linear_u;}
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
            nodesets.push_back(NodeSet(trim(split_line[0]),nset_nodes));
            mms_dirichlet_set_number = nodesets.size()-1;
                            
            if(verbose){
                std::cout << "nodeset: " << nodesets[mms_dirichlet_set_number].name << "\n";
                for(int i=0; i<nodesets[mms_dirichlet_set_number].nodes.size(); i++){std::cout << " " << nodesets[mms_dirichlet_set_number].nodes[i];}
                std::cout << "\n";
            }  
        }
    }
}

FEAModel::FEAModel(){
    /*!Default constructor*/
}
    
FEAModel::FEAModel(InputParser _input){
    /*!The full constructor*/
    
    std::cout << "=\n|=> Constructing FEA Model\n=\n";
    
    input = _input; //!Copy input
    
    //!Set the vector sizes
    total_ndof = input.nodes.size()*input.node_dof;
    up         = std::vector< double >(total_ndof,0.);
    u          = std::vector< double >(total_ndof,0.);
    du         = std::vector< double >(total_ndof,0.);
    
    //!Initialize the size of the dirichlet boundary condition vector
    dbcdof.resize(0);
    
    //!Define the internal node numbering of the elements
    map_element_nodes();
    
    map_nodesets();
    
    //!Initialize the degree of freedom vector
    initialize_dof();
    
    if(input.verbose){
        for(int i=0; i<mapped_elements.size(); i++){
            std::cout << "Element " << mapped_elements[i].number << " internal node numbers:";
            for(int j=0; j<mapped_elements[i].nodes.size(); j++){std::cout << " " << mapped_elements[i].nodes[j];}
            std::cout << "\n";
        }
        
        for(int i=0; i<mapped_nodesets.size(); i++){
            std::cout << "Nodeset: " << mapped_nodesets[i].name << " internal node numbers:";
            for(int j=0; j<mapped_nodesets[i].nodes.size(); j++){std::cout << " " << mapped_nodesets[i].nodes[j];}
            std::cout << "\n";
        }
    }
}
    
void FEAModel::map_element_nodes(){
    /*!===========================
    |    map_element_nodes    |
    ===========================
    
    Map the nodes that define each element 
    as the user defined them to the internal 
    node numbering.
    
    */
    
    std::cout << "\n|=> Mapping element nodes\n";
    
    bool node_number_check = false;                         //!Boolean which checks that the node numbers are 
                                                            //!set coorectly
    unsigned int current_node_number;                       //!The current node number as defined by the user
    std::vector< unsigned int > internal_node_numbers(8,0); //!The internal node numbers for a given element
    
    for(int e=0; e<input.elements.size(); e++){//Iterate through the elements
        for(int el_n=0; el_n<input.elements[e].nodes.size(); el_n++){//Iterate through the nodes in the element
            current_node_number = input.elements[e].nodes[el_n]; //Set the current user defined node number
            node_number_check   = false;                         //Reset the node number check boolean
            for(int n=0; n<input.nodes.size(); n++){             //Search the list of nodes to find the indicated number
                if(input.nodes[n].number==current_node_number){
                    internal_node_numbers[el_n] = n;
                    node_number_check = true;
                }
            }
            if(!node_number_check){//Check for if the indicated node is not defined
                std::cout << "Error: Element "<<input.elements[e].number<< " calls for node number " << input.elements[e].nodes[el_n] <<
                             "       which is not defined.\n";
                assert(1==0);
            }
        }
        mapped_elements.push_back(Element(input.elements[e].number, internal_node_numbers));
    } 
    std::cout << "\n|=> Mapping complete\n";
}
    
void FEAModel::map_nodesets(){
    /*!======================
    |    map_nodesets    |
    ======================
    
    Map the nodesets from the user 
    defined numbering to the internal 
    numbering.
    
    */
    
    std::cout << "\n|=> Mapping nodesets\n";
    
    unsigned int node_number; //!The current node
    
    mapped_nodesets = input.nodesets; //!Copy over the nodesets
    
    for(int n=0; n<input.nodes.size(); n++){//Iterate through all of the defined nodes
        
        node_number = input.nodes[n].number;
        
        for(int s=0; s<input.nodesets.size(); s++){//Iterate through all of the nodesets
            
            for(int m=0; m<input.nodesets[s].nodes.size(); m++){//Iterate through the nodes in the nodeset
                
                if(input.nodesets[s].nodes[m]==node_number){
                    mapped_nodesets[s].nodes[m] = n; //Map the user defined node numbering to the global node numbering in the 
                                                     //Nodesets
                }
                
            }
            
        }
    
    }
        
    std::cout << "\n|=> Mapping complete\n";
}
    
/*!=
|=> Degrees of freedom methods
=*/
    
void FEAModel::initialize_dof(){
    /*!========================
    |    initialize_dof    |
    ========================
    
    Initialize the degree of freedom vector
        
    */
        
    std::cout << "=\n"<<
                 "|=> Initializing degrees of freedom\n"<<
                 "=\n";
                     
    assign_dof();                      //!Assign the degrees of freedom to the nodes
    
    convert_local_dbc_to_global();     //!Convert the local dirichlet boundary conditions to 
                                       //!global boundary conditions.
    
    id_unbound_dof();                  //!Identify all of the unbound dof
    std::cout << "=\n"<<
                 "|=> Degrees of freedom initialized\n"<<
                 "=\n";
    return;
}
    
void FEAModel::convert_local_dbc_to_global(){
    /*!=====================================
    |    convert_local_dbc_to_global    |
    =====================================
        
    Convert the locally defined dirichlet boundary 
    conditions to the global dof number
        
    */
        
    bool node_number_check = false;   //!Check to make sure the node is defined correctly
    unsigned int current_node_number; //!The current node number
    dbcdof = input.dirichlet_bcs;     //!Copy over the dirichlet boundary condition vector
        
    for(int dbc=0; dbc<input.dirichlet_bcs.size(); dbc++){//Iterate through the dirichlet boundary conditions
        node_number_check = false;
        //std::cout << "node_number: " << input.dirichlet_bcs[dbc].node_number << " dof number: " << input.dirichlet_bcs[dbc].dof_number << " value: " << input.dirichlet_bcs[dbc].value << "\n";
        for(int n=0; n<input.nodes.size(); n++){//Iterate through the nodes
            if(input.dirichlet_bcs[dbc].node_number==input.nodes[n].number){
                dbcdof[dbc] = DirichletBC(n,internal_nodes_dof[n][input.dirichlet_bcs[dbc].dof_number-1],input.dirichlet_bcs[dbc].value); //Append the dirichlet boundary condition to the list
                node_number_check = true;
            }
        }
        if(!node_number_check){
            std::cout << "Error: Node number " << input.dirichlet_bcs[dbc].node_number << " in the dirichlet boundary condition on"<<
                         "       local degree of freedom " << input.dirichlet_bcs[dbc].dof_number << " with a "<<
                         "       value of " << input.dirichlet_bcs[dbc].value << " is not defined.\n";
            assert(1==0);
        }
    }
    return;
}
    
/*!=
|=> Solver methods
=*/
    
bool FEAModel::solve(){
    /*!===============
    |    solve    |
    ===============
        
    Solve the finite element problem
        
    */
        
    std::cout << "\n=================================================\n"<<
                   "|                                               |\n"<<
                   "|                BEGINNING SOLVER               |\n"<<
                   "|                                               |\n"<<
                   "=================================================\n";
    
    input.t += input.tp+input.dt; //!Set initial timestep increment
    bool result;                  //!The result of the simulation
        
    while(input.tp<input.total_time){  //Iterate through the timesteps
    
        if(input.mms_fxn!=NULL){
            apply_manufactured_solution(); //!Apply the manufactured solution forcing function
        }
            
        result = increment_solution(); //Increment the solution at the timestep
            
        if(result){
            input.tp = input.t;          //Increment time
            input.t  = input.t+input.dt;
            if(input.t>input.total_time){
                input.t  = input.total_time;
                input.dt = input.t - input.tp;
            }
            
            for(int i=0; i<u.size(); i++){
                up[i] = u[i]; //Set the previous dof vector to the current
            }
                
        }
        else{
            if(input.mms_fxn!=NULL){
                compare_manufactured_solution();
            }
            return false;
        }
            
    }
        
    std::cout << "\n=================================================\n"<<
                   "|                                               |\n"<<
                   "|                SOLVER COMPLETED               |\n"<<
                   "|                                               |\n"<<
                   "=================================================\n";
        
    if(input.mms_fxn!=NULL){
        compare_manufactured_solution();
    }
        
    return true;
}
    
bool FEAModel::increment_solution(){
    /*!============================
    |    increment_solution    |
    ============================
        
    Perform a timestep of the solution
        
    */
        
    initialize_timestep();
    update_increment();
        
    return true;
}
    
void FEAModel::initialize_timestep(){
    /*!=============================
    |    initialize_timestep    |
    =============================
        
    Initialize the timestep for the 
    current iteration.
        
    */
        
    form_increment_dof_vector();
    increment_number += 1;
        
    return;
}
    
void FEAModel::update_increment(){
    /*!==========================
    |    update_increment    |
    ==========================
        
    Update the increment using the 
    chosen solution technique.
        
    */
        
    std::cout << "=\n|=> Begining solution of increment " << increment_number << "\n=\n";
        
    if(input.verbose){
        std::cout << "Solution prior to iteration\n";
        for(int n=0; n<internal_nodes_dof.size(); n++){
            std::cout << "Internal node " << n << " dof ";
            for(int i=0; i<input.node_dof; i++){
                std::cout << " (" << internal_nodes_dof[n][i] << ", " << u[internal_nodes_dof[n][i]]<<")";
            }
            std::cout << "\n";
        }
        std::cout << "Change in solution prior to iteration\n";
        for(int n=0; n<internal_nodes_dof.size(); n++){
            std::cout << "Internal node " << n << " dof ";
            for(int i=0; i<input.node_dof; i++){
                std::cout << " (" << internal_nodes_dof[n][i] << ", " << du[internal_nodes_dof[n][i]]<<")";
            }
            std::cout << "\n";
        }
    }
        
    if(!solver.compare("NewtonKrylov")){//Solve the equations using a Jacobian free Newton-Krylov method
        run_newton_krylov();
    }
    
    if(input.verbose){
        std::cout << "Solution after iteration\n";
        for(int n=0; n<internal_nodes_dof.size(); n++){
            std::cout << "Internal node " << n << " dof ";
            for(int i=0; i<input.node_dof; i++){
                std::cout << " " << u[internal_nodes_dof[n][i]];
            }
            std::cout << "\n";
        }
    }
    
    for(int i=0; i<u.size(); i++){
        up[i] = u[i];              //Update the previous value of u
    }
    
    return;
}
    
void FEAModel::run_newton_krylov(){
    /*!===========================
    |    run_newton_krylov    |
    ===========================
        
    Run the Newton-Krylov solver.
        
    */
    
    std::vector< double > ub_du = get_unbound_du();
    FEAKrylovSolver krylov_solver = FEAKrylovSolver(*this, ub_du, 10, 40, true);
    krylov_solver.solve();
    return;
}

void FEAModel::id_unbound_dof(){
    /*!========================
    |    id_unbound_dof    |
    ========================
    
    Identify the degrees of freedom 
    which are unbound.
    
    */
    
    std::vector< bool > unbound_dof_bool(total_ndof,true); //!A vector indicating if a 
                                                           //!degree of freedom is unbound or not
    std::vector< unsigned int > dof_tmp(input.node_dof,0); //!A temporary vector of degrees of freedom
    
    if(input.mms_dirichlet_set_number>0){ //!Identify all of the dof associated with the dirichlet set if it exists
        for(int i=0; i<mapped_nodesets[input.mms_dirichlet_set_number].nodes.size(); i++){//Index through the dirichlet bc set
            dof_tmp = internal_nodes_dof[mapped_nodesets[input.mms_dirichlet_set_number].nodes[i]];
            for(int j=0; j<dof_tmp.size(); j++){
                unbound_dof_bool[dof_tmp[j]] = false;
            }
        }
    }
    
    for(int i=0; dbcdof.size(); i++){//Identify all of the dof associated with Dirichlet boundary conditions
        unbound_dof_bool[dbcdof[i].dof_number] = false;
    }
    
    for(int i=0; i<unbound_dof_bool.size(); i++){//Set the unbound degrees of freedom
        if(unbound_dof_bool[i]){
            unbound_dof.push_back(i);
        }
    }
    
    if(input.verbose){
        std::cout << "unbound dof: ";
        for(int i=0; i<unbound_dof.size(); i++){std::cout << " " << unbound_dof[i];}
        std::cout << "\n";
    }
}

std::vector< double > FEAModel::get_unbound_du(){
    /*!========================
    |    get_unbound_du    |
    ========================
    
    Get the change in all of the 
    unbound degrees of freedom
    
    */
    
    std::vector< double > ub_du(0,0);
    
    for(int i=0; i<unbound_dof.size(); i++){
        ub_du.push_back(du[unbound_dof[i]]);
    }
    
    return ub_du;
}
    
std::vector< double > FEAModel::krylov_residual(std::vector<double> ub_du){
    /*!=========================
    |    krylov_residual    |
    =========================
        
    Access the RHS vector in a way useful to the 
    Krylov solver.
        
    */
    
    unsigned int dof_num;
        
    for(int i=0; i<ub_du.size(); i++){
        dof_num = unbound_dof[i];
        du[dof_num]  = ub_du[i];                       //Update du
        u[dof_num]   = up[dof_num] + du[dof_num];      //Update u
    }
    
    assemble_RHS_and_jacobian_matrix(); //Compute the RHS
    
    std::vector< double > sub_RHS(ub_du.size(),0.); //Get the residual values associated 
                                                    //with the degrees of freedom
    
    for(int i=0; i<unbound_dof.size(); i++){
        sub_RHS[i] = RHS[unbound_dof[i]];           //Get the RHS values not on the boundary
    }
    
    return sub_RHS;                         //Return the residual
}
    
void FEAModel::form_increment_dof_vector(){
    /*!===================================
    |    form_increment_dof_vector    |
    ===================================
        
    Form the incremented dof vector at a given increment
        
    */
    
    double factor = input.dt/input.total_time; //!The scaling factor for the dirichlet boundary conditions
    
    for(int dbc=0; dbc<dbcdof.size(); dbc++){//Apply the dirichlet boundary conditions
        du[dbcdof[dbc].dof_number]  = factor*dbcdof[dbc].value;
        u[dbcdof[dbc].dof_number]  += du[dbcdof[dbc].dof_number];
    }
    return;
}
    
void FEAModel::assign_dof(){
    /*!====================
    |    assign_dof    |
    ====================
        
    Assign degrees of freedom to the nodes
        
    */
    internal_nodes_dof.resize(input.nodes.size()); //!Resize the internal nodes to dof vector
    unsigned int dof_val = 0;                      //!The current degree of freedom value
        
    for(unsigned int n=0; n<input.nodes.size(); n++){
        internal_nodes_dof[n].resize(input.node_dof);
        for(int i=0; i<input.node_dof; i++){
            internal_nodes_dof[n][i] = dof_val;
            dof_val++;
        }
    }
        
    if(input.verbose){
        std::cout << "Global degrees of freedom at nodes (internal numbering)\n";
        for(int n=0; n<internal_nodes_dof.size(); n++){
            std::cout << "Node " << n << " dof:";
            for(int i=0; i<internal_nodes_dof[n].size(); i++){std::cout << " " << internal_nodes_dof[n][i];}
            std::cout << "\n";
        }
    }
        
    return;
}
    
void FEAModel::assemble_RHS_and_jacobian_matrix(){
    /*!==========================================
    |    assemble_RHS_and_jacobian_matrix    |
    ==========================================
       
    Assemble the global right hand side vector and 
    jacobian matrix if required for the solution 
    technique.
        
    */
       
    RHS = std::vector< double >(total_ndof,0.); //Zero the residual vector
    
    if(input.mms_fxn!=NULL){//Add the manufactured solution forcing function if required
        for(int i=0; i<RHS.size(); i++){
            RHS[i] = -F[i];   //Negative because the RHS is the residual for the nonlinear calculation
        }
    }
        
    std::cout << "=\n"<<
                "| Computing RHS and global stiffness matrix\n"<<
                "=\n";
        
    std::vector< double > element_coordinates(24,0.);      //!The coordinates of the nodes in a given element
    unsigned int internal_node_number;                     //!The number of the node as defined in the code
    std::vector< double > element_u(input.node_dof*8,0.);  //!The solution variable for the element
    std::vector< double > element_du(input.node_dof*8,0.); //!The change in solution variable for the element
        
    micro_element::Hex8 current_element;                   //!The current element
        
    for(int e=0; e<mapped_elements.size(); e++){//Iterate through the elements
        //Construct the reference coordinates of the element
        for(int n=0; n<8; n++){
                
            internal_node_number = mapped_elements[e].nodes[n];
                
            for(int i=0; i<input.nodes[internal_node_number].coordinates.size(); i++){
                element_coordinates[i+n*input.nodes[internal_node_number].coordinates.size()] = input.nodes[internal_node_number].coordinates[i];
            }
        }
            
        if(input.verbose){
            std::cout << "Element " << mapped_elements[e].number << " reference coords";
            for(int i=0; i<element_coordinates.size(); i++){std::cout << " " << element_coordinates[i];}
            std::cout << "\n";
        }
            
        //Construct the u vector and du for the element
        for(int n=0; n<8; n++){
            internal_node_number = mapped_elements[e].nodes[n];
            for(int i=0; i<input.node_dof; i++){
                element_u[i+n*input.node_dof]  =  u[internal_nodes_dof[internal_node_number][i]];
                element_du[i+n*input.node_dof] = du[internal_nodes_dof[internal_node_number][i]];
            }
        }
            
        //Construct the element
        current_element = micro_element::Hex8(element_coordinates, element_u, element_du,
                                              input.fprops, input.iprops);
            
        //Integrate the element
        current_element.integrate_element();
            
        //Update the RHS vector and jacobian matrix
        for(int n=0; n<8; n++){
            internal_node_number = mapped_elements[e].nodes[n];
            for(int i=0; i<input.node_dof; i++){
                RHS[internal_nodes_dof[internal_node_number][i]] += current_element.RHS[i+n*input.node_dof];
            }
        }
        
        //for(int n=0; n<8; n++){
        //    std::cout << "Element node " << n+1 << " RHS:";
        //    for(int i=0; i<input.node_dof; i++){
        //        std::cout << " " << current_element.RHS[i+n*input.node_dof];
        //    }
        //    std::cout << "\n";
        //}
    }
    //assert(1==0);
    return;
}
    
/*!=
|=> Manufactured solutions methods
=*/
    
void FEAModel::apply_manufactured_solution(){
    /*!=====================================
    |    apply_manufactured_solution    |
    =====================================
        
    Initialize the manufactured solution 
    parameters as required.
        
    */
        
    compute_mms_forcing_function();
    compute_mms_bc_values();
    return;
}
    
void FEAModel::compute_mms_forcing_function(){
    /*!======================================
    |    compute_mms_forcing_function    |
    ======================================
        
    Compute the forcing function for the method 
    of manufactured solutions.
        
    */
        
    std::cout << "\n|=> Computing the method of manufactured solutions forcing function.\n";
        
    F          = std::vector< double >(total_ndof,0.); //Initialize the forcing vector
        
    set_mms_dof_vector();                            //Set u to the manufactured solutions vector
    assemble_RHS_and_jacobian_matrix();              //Compute the residual value for the manufactured solution
    for(int i=0; i<RHS.size(); i++){F[i] = RHS[i];}  //Copy the residual vector to the forcing function vector
    
    for(int j=0; j<up.size(); j++){
        u[j]  = up[j];              //Reset u
        du[j] = 0;                  //Reset du
    }
        
    std::cout << "\n|=> Forcing function created\n";
        
    return;
}
    
void FEAModel::set_mms_dof_vector(){
    /*!============================
    |    set_mms_dof_vector    |
    ============================
        
    Set the manufactured solutions 
    degree of freedom vector.
        
    */
        
    std::vector<double> utmp(input.node_dof,0);  //!The temporary u value
    mms_u = std::vector< double >(total_ndof,0); //!Initialize the manufactured solution vector
        
    for(int n=0; n<input.nodes.size(); n++){
        utmp = input.mms_fxn(input.nodes[n].coordinates, input.t);
            
        if(input.verbose){
            std::cout << "Node " << n << " dof:";
            for(int i=0; i<utmp.size(); i++){std::cout << " " << utmp[i];}
            std::cout << "\n";
        }
            
        for(int j=0; j<input.node_dof; j++){
            mms_u[j+n*input.node_dof] = utmp[j];                         //!Update the manufactured solution vector
            du[j+n*input.node_dof]    = utmp[j] - u[j+n*input.node_dof]; //!Set the change in u for the manufactured solution
            u[j+n*input.node_dof]     = utmp[j];                         //!Set u to the manufactured solution vector
        }
    }
        
    return;
}
    
void FEAModel::compute_mms_bc_values(){
    /*!===============================
    |    compute_mms_bc_values    |
    ===============================
        
    Compute and apply the manufactured solution 
    on the boundary.
        
    */
        
    //!Initialize required values
    unsigned int node_number;                       //!The current node number (internal numbering)
    std::vector< double > utmp(input.node_dof,0.);  //!Degree of freedom vector at the node
    std::array< double, 3 > coordinates;            //!The coordinates of the node
    std::vector< unsigned int > internal_dof;       //!The internal numbering of the degree of freedom at the internal numbering of the nodes
        
    if(input.mms_dirichlet_set_number>0){
        for(int n=0; n<mapped_nodesets[input.mms_dirichlet_set_number].nodes.size(); n++){//Iterate through the boundary nodes (internal numbering)
            node_number = mapped_nodesets[input.mms_dirichlet_set_number].nodes[n];       //Set the internal node number
            
            coordinates      = input.nodes[node_number].coordinates;         //Set the coordinates of the node
            internal_dof     = internal_nodes_dof[node_number];              //Set the degrees of freedom of the node (internal numbering)
            
            utmp = input.mms_fxn(coordinates,input.t);                       //Compute the degree of freedom vector
            
            for(int i=0; i<internal_dof.size(); i++){     //Update the degree of freedom vector and the change in the 
                                                          //degree of freedom vector with the new boundary conditions
                du[internal_dof[i]] = utmp[i]-u[internal_dof[i]];
                u[internal_dof[i]]  = utmp[i];
            }
        } 
    }
    return;
}
    
void FEAModel::compare_manufactured_solution(){
    /*!=======================================
    |    compare_manufactured_solution    |
    =======================================
    
    Compare the manufactured solution to the output 
    from the code.
        
    */
        
    bool result=true;
    double max_error=0;
    double max_relative_error=0;
    double error;
    double relative_error;
        
    for(int i=0; i<mms_u.size(); i++){
        error          = fabs(mms_u[i]-u[i]);
        relative_error = error/std::max(fabs(mms_u[i]),fabs(u[i]));
        if(error>max_error){max_error = error;}
        if(relative_error>max_relative_error){max_relative_error = relative_error;}
    }
    
    if((max_error>mms_tol) && (max_relative_error>mms_tol)){result=false;}
    
    if(result){
        std::cout << "Manufactured solution passed\n";
    }
    else{
        std::cout << "Error: Manufactured solution did not pass\n";
        print_vector("u",u);
        print_vector("mms_u",mms_u);
    }
    
    std::cout << "\nMaximum error: " << max_error << "\n";
        
    
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
        FEAModel FM = FEAModel(IP);
        FM.solve();
    }
        
}
