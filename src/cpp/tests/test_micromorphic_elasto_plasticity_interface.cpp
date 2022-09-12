//Tests for constitutive_tools

#include<micromorphic_elasto_plasticity.h>
#include<micromorphic_material_library.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

#define BOOST_TEST_MODULE test_micromorphic_elasto_plasticity_interface
#include <boost/test/included/unit_test.hpp>

typedef micromorphicTools::constantType constantType;
typedef micromorphicTools::constantVector constantVector;
typedef micromorphicTools::constantMatrix constantMatrix;

typedef micromorphicTools::parameterType parameterType;
typedef micromorphicTools::parameterVector parameterVector;
typedef micromorphicTools::parameterMatrix parameterMatrix;

typedef micromorphicTools::variableType variableType;
typedef micromorphicTools::variableVector variableVector;
typedef micromorphicTools::variableMatrix variableMatrix;

typedef micromorphicTools::errorNode errorNode;
typedef micromorphicTools::errorOut errorOut;

BOOST_AUTO_TEST_CASE( testMaterialLibraryInterface ){
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    BOOST_CHECK( material );

    //Set up the inputs
    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                        {0.100, 0.001, 0.000 },
                                        {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                          { -0.18572739,  0.06847269,  0.22931628 },
                                          { -0.01829735, -0.48731265, -0.25277529 },
                                          {  0.26626212,  0.4844646 , -0.31965177 },
                                          {  0.49197846,  0.19051656, -0.0365349  },
                                          { -0.06607774, -0.33526875, -0.15803078 },
                                          {  0.09738707, -0.49482218, -0.39584868 },
                                          { -0.45599864,  0.08585038, -0.09432794 },
                                          {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    solverTools::homotopyMap DEBUG;
#endif

    solverTools::floatVector PK2_answer = { 172.484,   15.3785,   -0.917177,
                                             13.4848, 142.823,    -0.0214307,
                                             -1.7635,   1.77719, 141.069 };

    solverTools::floatVector SIGMA_answer = { 176.916,   15.8646,   -2.83731,
                                               15.8646, 144.538,     1.85836,
                                               -2.83731,  1.85836, 142.013 };

    solverTools::floatVector M_answer = { 0.598283, -0.512218,  0.620664,    3.22636,   1.16682,
                                          1.20593,   0.562825, -2.52317,     1.62616,  -2.61391,
                                         -0.60994,  -1.02147,   0.668187,    0.49348,  -0.23916,
                                         -2.77419,   0.760483,  1.71784,    -0.499389,  2.62828,
                                         -0.761044,  1.23369,  -0.00778206, -2.25643,  -0.729551,
                                          0.743204,  0.910521 };

    solverTools::floatVector SDVS_answer = { -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
                                             -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
                                              0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
                                              0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
                                              0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
                                             -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
                                             -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
                                             -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
                                              0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
                                             -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
                                              0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 };

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    //Evaluate the model
    int errorCode = material->evaluate_model( time, fparams,
                                              current_grad_u, current_phi, current_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
#ifdef DEBUG_MODE
                                              , DEBUG
#endif
                                            );

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVS_answer ) );

    //Check the Jacobian using the previously tested jacobian
    std::vector< std::vector< double > > DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
                                         DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
                                         DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer;

    std::vector< std::vector< std::vector< double > > > ADD_JACOBIANS;

    SDVS = SDVSDefault;

#ifdef DEBUG_MODE
    DEBUG.clear();
#endif

    errorCode = micromorphicElastoPlasticity::evaluate_model(
                                time, fparams,
                                current_grad_u, current_phi, current_grad_phi,
                                previous_grad_u, previous_phi, previous_grad_phi,
                                SDVS,
                                current_ADD_DOF, current_ADD_grad_DOF,
                                previous_ADD_DOF, previous_ADD_grad_DOF,
                                PK2_result, SIGMA_result, M_result,
                                DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
                                DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
                                DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer,
                                ADD_TERMS, ADD_JACOBIANS,
                                output_message
#ifdef DEBUG_MODE
                                , DEBUG
#endif
                              );

    BOOST_CHECK( errorCode <= 0 );

    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();

    SDVS = SDVSDefault;

    std::vector< std::vector< double > > DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                         DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                         DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result;

    errorCode = material->evaluate_model( time, fparams,
                                          current_grad_u, current_phi, current_grad_phi,
                                          previous_grad_u, previous_phi, previous_grad_phi,
                                          SDVS,
                                          current_ADD_DOF, current_ADD_grad_DOF,
                                          previous_ADD_DOF, previous_ADD_grad_DOF,
                                          PK2_result, SIGMA_result, M_result,
                                          DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                          DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                          DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
                                          ADD_TERMS, ADD_JACOBIANS,
                                          output_message
#ifdef DEBUG_MODE
                                          , DEBUG
#endif
                                        );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS_answer, SDVS ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer ) );

#ifdef DEBUG_MODE
    DEBUG.clear();
#endif

    //Test the computed numeric Jacobian values
    SDVS = SDVSDefault;
    errorCode = material->evaluate_model_numeric_gradients( time, fparams,
                                                            current_grad_u, current_phi, current_grad_phi,
                                                            previous_grad_u, previous_phi, previous_grad_phi,
                                                            SDVS,
                                                            current_ADD_DOF, current_ADD_grad_DOF,
                                                            previous_ADD_DOF, previous_ADD_grad_DOF,
                                                            PK2_result, SIGMA_result, M_result,
                                                            DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                                            DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                                            DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
                                                            ADD_TERMS, ADD_JACOBIANS,
                                                            output_message,
#ifdef DEBUG_MODE
                                                            DEBUG,
#endif
                                                            1e-6 );

    BOOST_CHECK( errorCode <= 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVS, SDVS_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer, 1e-4, 1e-5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer, 1e-4 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer, 1e-4 ) );

}

BOOST_AUTO_TEST_CASE( testMaterialLibraryInterface2 ){
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     * NOTE: This function mostly exists to perform debugging
     *       on the implementation of the function into a 
     *       larger solver code.
     *
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    //Set up the inputs
    //Initialize the time
    std::vector< double > time = { 0.045, 0.01 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 170, 15, 2, 140, 20, 2, 2, 27, 2, 0.56, 0.2, 2, 0.15, 0.3, 2, 0.82, 0.1, 2, 0.42, 0.3, 2, 0.05, 0.2, 2, 0.52, 0.4, 2, 29480, 25480, 5, 1000, 400, -1500, -1400, -3000, 11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0, 2, 400, -3000, 0.5, 0.5, 0.5, 1e-09, 1e-09 };

    //Initialize the gradient of the macro displacement
    double current_grad_u[ 3 ][ 3 ] =
    {
        { -0.00124343, -6.55319e-14, 3.99657e-13},
        { 0, 0.0045, 0},
        { -1.75135e-13, -1.35481e-13, -0.00124343 },
    };

    double previous_grad_u[ 3 ][ 3 ] =
    {
        { -0.00123858, -1.22379e-17, 5.04154e-18},
        { 0, 0.004, 0},
        { -1.47723e-18, 4.44523e-18, -0.00123858 },
    };

    //Initialize the micro displacement
    double current_phi[ 9 ] = { -0.00153489, -3.04626e-13, 5.16537e-13, 1.58771e-13, 0.00303407, 4.29828e-14, -4.38368e-13, -1.80694e-13, -0.00153489 };

    double previous_phi[ 9 ] = { -0.00164749, -2.63663e-17, 1.35603e-17, 8.65138e-19, 0.00325613, -2.13082e-20, -1.17433e-17, 2.24626e-18, -0.00164749 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    solverTools::homotopyMap DEBUG;
#endif

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    //Evaluate the model
    int errorCode = material->evaluate_model( time, fparams,
                                              current_grad_u, current_phi, current_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
#ifdef DEBUG_MODE
                                              , DEBUG
#endif
                                            );

    BOOST_CHECK( errorCode <= 0 );

//    std::cout << "SDVS:\n"; vectorTools::print( SDVS );
//    std::cout << "PK2:\n"; vectorTools::print( PK2_result );
//    std::cout << "SIGMA:\n"; vectorTools::print( SIGMA_result );
//    std::cout << "M:\n"; vectorTools::print( M_result );

#ifdef DEBUG_MODE
    for ( auto step = DEBUG.begin(); step != DEBUG.end(); step++ ){
        std::cout << step->first << "\n";
        for ( auto iteration = step->second.begin(); iteration != step->second.end(); iteration++ ){
            std::cout << "    " << iteration->first << "\n";
            for ( auto value = iteration->second.begin(); value != iteration->second.end(); value++ ){
                if ( value->second.size() <= 27 ) {
                    std::cout << "        " << value->first << "\n";
                    std::cout << "            "; vectorTools::print( value->second );
                }
            }
        }
    }
#endif

}

BOOST_AUTO_TEST_CASE( testEvaluate_model_history ){
    /*!
     * Test the material model undergoing a time history.
     *
     */

#ifdef DEBUG_MODE
    std::ofstream output_file;
    output_file.open( "output_file.txt" );
#endif

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    std::vector< std::vector< double > > grad_u_0 = { { 0, 0, 0 },
                                                      { 0, 0, 0 },
                                                      { 0, 0, 0 } };

    std::vector< std::vector< double > > grad_u_f = { { 0.5, 0, 0 },
                                                      { 0.0, 0, 0 },
                                                      { 0.0, 0, 0 } };

    std::vector< double > phi_0 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector< double > phi_f = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    std::vector< std::vector< double > > grad_phi_0 = { { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0} };

    std::vector< std::vector< double > > grad_phi_f = { { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0} };

    double dt = 0.05;
    double t0 = 0.;
    double tf = 0.25;

    double t = t0;

    //Set up the model parameters
    std::vector< double > fparams = { 2, 1e3, 1e2,
                                      2, 7e2, 1e4,
                                      2, 1e3, 1e4,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 29480, 25480,
                                      5, 1000, 400, -1500, -1400, -3000,
                                      11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0,
                                      2, 400, -3000,
                                      0.5, 0.5, 0.5, 1e-09, 1e-09 };
//                                      0.0, 0.0, 0.0, 1e-09, 1e-09 };

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    solverTools::homotopyMap DEBUG;
#endif

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    std::vector< double > PK2Answer = {  5.16027494e+03, -4.30718970e-18,  4.23714575e-18,
                                        -4.30349349e-18,  4.05236467e+03, -1.78599775e-18,
                                         4.25600330e-18, -1.80293236e-18,  4.05236467e+03 };

    std::vector< double > SIGMAAnswer = {  5051.82,     -4.38296e-18,  4.03864e-18,
                                          -4.38296e-18,  4116.16,     -1.6831e-18,
                                           4.03864e-18, -1.6831e-18,   4116.16 };

    std::vector< double > MAnswer = {  1.0019e-46,  -9.52574e-49,  1.00995e-48, -8.8467e-25,   8.11869e-46,
                                      -5.63453e-34, -1.10865e-32,  2.44433e-34,  1.6875e-49,  -9.38296e-38,
                                       4.14595e-47, -6.23442e-59, -1.42808e-16,  1.31056e-37, -9.48874e-38,
                                      -8.77566e-16,  8.05349e-37, -5.83091e-37, -1.10003e-37, -1.80712e-47,
                                       3.77743e-56, -1.06975e-37,  9.81714e-59, -1.63467e-56, -8.77566e-16,
                                       8.05349e-37, -1.05373e-36 };

    std::vector< double > SDVSAnswer = { 0.0536436,    0.0316833,    0,            0,            0,
                                         0.441627,     2.79534e-10,  0,            0,            0,
                                         0.0406588,   -4.90538e-23,  6.74381e-23,  2.31867e-22, -0.0197043,
                                        -2.04791e-23, -2.35533e-22,  8.57369e-23, -0.0197043,    0.0150881,
                                        -5.35887e-23, -5.0789e-24,   2.923e-22,   -0.00745964,  -6.81252e-24,
                                        -8.98761e-24,  1.5285e-24,  -0.00745964,   1.47911e-31, -4.70981e-22,
                                        -1.60422e-39, -4.82898e-51,  1.54956e-54,  5.15605e-62,  2.70794e-38,
                                         1.60531e-36,  1.84864e-40,  1.93673e-43, -1.471e-31,    2.93537e-22,
                                         1.85391e-40,  4.83968e-64,  1.1391e-40,   5.43161e-38, -1.48104e-62,
                                        -1.08569e-49, -2.16264e-24, -7.51219e-32,  1.97215e-31, -5.93004e-39,
                                        -9.83265e-57,  2.94509e-61,  1.67967e-41, -1.14207e-61, -7.96338e-64 };

    std::vector< std::vector< double > > grad_u_prev   = grad_u_0;
    std::vector< double > phi_prev                     = phi_0;
    std::vector< std::vector< double > > grad_phi_prev = grad_phi_0;

    std::vector< std::vector< double > > grad_u_curr;
    std::vector< double > phi_curr;
    std::vector< std::vector< double > > grad_phi_curr;

    std::vector< double > time;

    double current_grad_u[ 3 ][ 3 ], current_phi[ 9 ], current_grad_phi[ 9 ][ 3 ];
    double previous_grad_u[ 3 ][ 3 ], previous_phi[ 9 ], previous_grad_phi[ 9 ][ 3 ];

    //Initial state
    //Update the arrays
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        previous_phi[ i ] = phi_prev[ i ];

        for ( unsigned int j = 0; j < 3; j++ ){
            previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
        }
    }
    //Evaluate the model
    time = { 0., 0. };
    int errorCode = material->evaluate_model( time, fparams,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
#ifdef DEBUG_MODE
                                              , DEBUG
#endif
                                            );
   
    BOOST_CHECK( errorCode <= 0);
#ifdef DEBUG_MODE
    if( errorCode > 0 ){ 
        output_file.close();
    }
#endif

#ifdef DEBUG_MODE
    output_file << "NEW_INCREMENT\n";

    //Output the time
    output_file << time[ 0 ] << ", " << time[ 1 ] << "\n";

    //Output the current gradient of u
    output_file << current_grad_u[ 0 ][ 0 ] << ", " << current_grad_u[ 0 ][ 1 ] << ", " << current_grad_u[ 0 ][ 2 ] << ", "
                << current_grad_u[ 1 ][ 0 ] << ", " << current_grad_u[ 1 ][ 1 ] << ", " << current_grad_u[ 1 ][ 2 ] << ", "
                << current_grad_u[ 2 ][ 0 ] << ", " << current_grad_u[ 2 ][ 1 ] << ", " << current_grad_u[ 2 ][ 2 ] << "\n";

    //Output the current micro displacement
    output_file << current_phi[ 0 ] << ", " << current_phi[ 1 ] << ", " << current_phi[ 2 ] << ", "
                << current_phi[ 3 ] << ", " << current_phi[ 4 ] << ", " << current_phi[ 5 ] << ", "
                << current_phi[ 6 ] << ", " << current_phi[ 7 ] << ", " << current_phi[ 8 ] << "\n";

    //Output the current gradient of the micro displacement
    output_file << current_grad_phi[ 0 ][ 0 ] << ", " << current_grad_phi[ 0 ][ 1 ] << ", " << current_grad_phi[ 0 ][ 2 ] << ", "
                << current_grad_phi[ 1 ][ 0 ] << ", " << current_grad_phi[ 1 ][ 1 ] << ", " << current_grad_phi[ 1 ][ 2 ] << ", "
                << current_grad_phi[ 2 ][ 0 ] << ", " << current_grad_phi[ 2 ][ 1 ] << ", " << current_grad_phi[ 2 ][ 2 ] << ", "
                << current_grad_phi[ 3 ][ 0 ] << ", " << current_grad_phi[ 3 ][ 1 ] << ", " << current_grad_phi[ 3 ][ 2 ] << ", "
                << current_grad_phi[ 4 ][ 0 ] << ", " << current_grad_phi[ 4 ][ 1 ] << ", " << current_grad_phi[ 4 ][ 2 ] << ", "
                << current_grad_phi[ 5 ][ 0 ] << ", " << current_grad_phi[ 5 ][ 1 ] << ", " << current_grad_phi[ 5 ][ 2 ] << ", "
                << current_grad_phi[ 6 ][ 0 ] << ", " << current_grad_phi[ 6 ][ 1 ] << ", " << current_grad_phi[ 6 ][ 2 ] << ", "
                << current_grad_phi[ 7 ][ 0 ] << ", " << current_grad_phi[ 7 ][ 1 ] << ", " << current_grad_phi[ 7 ][ 2 ] << ", "
                << current_grad_phi[ 8 ][ 0 ] << ", " << current_grad_phi[ 8 ][ 1 ] << ", " << current_grad_phi[ 8 ][ 2 ] << "\n";

    //Output the PK2 stress
    output_file << PK2_result[ 0 ] << ", " << PK2_result[ 1 ] << ", " << PK2_result[ 2 ] << ", "
                << PK2_result[ 3 ] << ", " << PK2_result[ 4 ] << ", " << PK2_result[ 5 ] << ", "
                << PK2_result[ 6 ] << ", " << PK2_result[ 7 ] << ", " << PK2_result[ 8 ] << "\n";

    //Output the SIGMA stress
    output_file << SIGMA_result[ 0 ] << ", " << SIGMA_result[ 1 ] << ", " << SIGMA_result[ 2 ] << ", "
                << SIGMA_result[ 3 ] << ", " << SIGMA_result[ 4 ] << ", " << SIGMA_result[ 5 ] << ", "
                << SIGMA_result[ 6 ] << ", " << SIGMA_result[ 7 ] << ", " << SIGMA_result[ 8 ] << "\n";

    //Output the M stress
    output_file << M_result[  0 ] << ", " << M_result[  1 ] << ", " << M_result[  2 ] << ", "
                << M_result[  3 ] << ", " << M_result[  4 ] << ", " << M_result[  5 ] << ", "
                << M_result[  6 ] << ", " << M_result[  7 ] << ", " << M_result[  8 ] << ", "
                << M_result[  9 ] << ", " << M_result[ 10 ] << ", " << M_result[ 11 ] << ", "
                << M_result[ 12 ] << ", " << M_result[ 13 ] << ", " << M_result[ 14 ] << ", "
                << M_result[ 15 ] << ", " << M_result[ 16 ] << ", " << M_result[ 17 ] << ", "
                << M_result[ 18 ] << ", " << M_result[ 19 ] << ", " << M_result[ 20 ] << ", "
                << M_result[ 21 ] << ", " << M_result[ 22 ] << ", " << M_result[ 23 ] << ", "
                << M_result[ 24 ] << ", " << M_result[ 25 ] << ", " << M_result[ 26 ] << "\n";

    std::vector< double > PK2_intermediate, SIGMA_intermediate, M_intermediate;

    //Determine if there were non-linear iterations
    auto inc = DEBUG.find( "converged_values" );
    if ( inc != DEBUG.end() ){
        PK2_intermediate   = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ];
        SIGMA_intermediate = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ];
        M_intermediate     = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ];
    }
    else{
        PK2_intermediate   = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediatePK2Stress" ];
        SIGMA_intermediate = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceMicroStress" ];
        M_intermediate     = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceHigherOrderStress" ];
    }

    //Output the intermediate PK2 stress
    output_file << PK2_intermediate[ 0 ] << ", " << PK2_intermediate[ 1 ] << ", " << PK2_intermediate[ 2 ] << ", "
                << PK2_intermediate[ 3 ] << ", " << PK2_intermediate[ 4 ] << ", " << PK2_intermediate[ 5 ] << ", "
                << PK2_intermediate[ 6 ] << ", " << PK2_intermediate[ 7 ] << ", " << PK2_intermediate[ 8 ] << "\n";

    //Output the intermediate SIGMA stress
    output_file << SIGMA_intermediate[ 0 ] << ", " << SIGMA_intermediate[ 1 ] << ", " << SIGMA_intermediate[ 2 ] << ", "
                << SIGMA_intermediate[ 3 ] << ", " << SIGMA_intermediate[ 4 ] << ", " << SIGMA_intermediate[ 5 ] << ", "
                << SIGMA_intermediate[ 6 ] << ", " << SIGMA_intermediate[ 7 ] << ", " << SIGMA_intermediate[ 8 ] << "\n";

    //Output the intermediate M stress
    output_file << M_intermediate[  0 ] << ", " << M_intermediate[  1 ] << ", " << M_intermediate[  2 ] << ", "
                << M_intermediate[  3 ] << ", " << M_intermediate[  4 ] << ", " << M_intermediate[  5 ] << ", "
                << M_intermediate[  6 ] << ", " << M_intermediate[  7 ] << ", " << M_intermediate[  8 ] << ", "
                << M_intermediate[  9 ] << ", " << M_intermediate[ 10 ] << ", " << M_intermediate[ 11 ] << ", "
                << M_intermediate[ 12 ] << ", " << M_intermediate[ 13 ] << ", " << M_intermediate[ 14 ] << ", "
                << M_intermediate[ 15 ] << ", " << M_intermediate[ 16 ] << ", " << M_intermediate[ 17 ] << ", "
                << M_intermediate[ 18 ] << ", " << M_intermediate[ 19 ] << ", " << M_intermediate[ 20 ] << ", "
                << M_intermediate[ 21 ] << ", " << M_intermediate[ 22 ] << ", " << M_intermediate[ 23 ] << ", "
                << M_intermediate[ 24 ] << ", " << M_intermediate[ 25 ] << ", " << M_intermediate[ 26 ] << "\n";

    //Output the state variables
    for ( unsigned int i = 0; i < SDVS.size()-1; i++ ){
        output_file << SDVS[ i ] << ", ";
    }
    output_file << SDVS[ SDVS.size() - 1 ] << "\n";

#endif


    //Begin iteration
    while ( t + dt < tf ){

        time = { t + dt, dt };

        //Increment the displacements
        grad_u_curr   = grad_u_prev   + dt * ( grad_u_f - grad_u_0 );
        phi_curr      = phi_prev      + dt * ( phi_f - phi_0 );
        grad_phi_curr = grad_phi_prev + dt * ( grad_phi_f - grad_phi_0 );

        //Update the arrays
        for ( unsigned int i = 0; i < 3; i++ ){
            for ( unsigned int j = 0; j < 3; j++ ){
                current_grad_u[ i ][ j ]  = grad_u_curr[ i ][ j ];
                previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
            }
        }

        for ( unsigned int i = 0; i < 9; i++ ){
            current_phi[ i ] = phi_curr[ i ];
            previous_phi[ i ] = phi_prev[ i ];

            for ( unsigned int j = 0; j < 3; j++ ){
                current_grad_phi[ i ][ j ] = grad_phi_curr[ i ][ j ];
                previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
            }
        }

        //Evaluate the model
#ifdef DEBUG_MODE
        DEBUG.clear();
#endif
        int errorCode = material->evaluate_model( time, fparams,
                                                  current_grad_u, current_phi, current_grad_phi,
                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                  SDVS,
                                                  current_ADD_DOF, current_ADD_grad_DOF,
                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                  PK2_result, SIGMA_result, M_result,
                                                  ADD_TERMS,
                                                  output_message
#ifdef DEBUG_MODE
                                                  , DEBUG
#endif
                                                );

//        std::cout << "SDVS:\n"; vectorTools::print( SDVS );

#ifdef DEBUG_MODE

//        if ( ( fabs( SDVS[ 0 ] ) > 1e-8 ) || ( errorCode != 0 ) ){
//            for ( auto inc = DEBUG.begin(); inc != DEBUG.end(); inc++ ){
//                std::cout << inc->first << "\n";
//                for ( auto itr = inc->second.begin(); itr != inc->second.end(); itr++ ){
//                    if ( itr->first.compare( "converged_values" ) != 0 ){
//                        std::cout << "    " << itr->first << "\n";
//                        std::cout << "        currentMacroGamma:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroGamma" ] );
//                        std::cout << "        currentMicroGamma:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGamma" ] );
//                        std::cout << "        currentMicroGradientGamma:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientGamma" ] );
//                        std::cout << "        currentDeformationGradient:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentDeformationGradient" ] );
//                        std::cout << "        currentElasticDeformationGradient:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentElasticDeformationGradient" ] );
//                        std::cout << "        currentPlasticDeformationGradient:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentPlasticDeformationGradient" ] );
//                        std::cout << "        currentPK2Stress:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentPK2Stress" ] );
//                        std::cout << "        currentMacroCohesion:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroCohesion" ] );
//                        std::cout << "        dMacroCohesiondMacroStrainISV:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "dMacroCohesiondMacroStrainISV" ] );
//                        std::cout << "        previousMacroFlowDirection:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "previousMacroFlowDirection" ] );
//                        std::cout << "        previousMicroFlowDirection:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "previousMicroFlowDirection" ] );
//                        std::cout << "        previousMicroGradientFlowDirection:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "previousMicroGradientFlowDirection" ] );
//                        std::cout << "        currentMacroFlowDirection:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroFlowDirection" ] );
//                        std::cout << "        currentMicroFlowDirection:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroFlowDirection" ] );
//                        std::cout << "        currentMicroGradientFlowDirection:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientFlowDirection" ] );
//                        std::cout << "        currentMacroStrainISV\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMacroStrainISV" ] );
//                        std::cout << "        currentMicroStrainISV\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroStrainISV" ] );
//                        std::cout << "        currentMicroGradientStrainISV\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "currentMicroGradientStrainISV" ] );
//                        std::cout << "        macroYieldFunction:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "macroYieldFunction" ] );
//                        std::cout << "        microYieldFunction:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "microYieldFunction" ] );
//                        std::cout << "        microGradientYieldFunction:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "microGradientYieldFunction" ] );
//                    }
//                    else{
//                        std::cout << "        convergedPlasticDeformationGradient:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "convergedPlasticDeformationGradient" ] );
//                        std::cout << "        convergedPlasticMicroDeformation:\n";
//                        std::cout << "        "; vectorTools::print( itr->second[ "convergedPlasticMicroDeformation" ] );
//                    }
//                }
//            }
//        }
//
        output_file << "NEW_INCREMENT\n";

        //Output the time
        output_file << time[ 0 ] << ", " << time[ 1 ] << "\n";

        //Output the current gradient of u
        output_file << current_grad_u[ 0 ][ 0 ] << ", " << current_grad_u[ 0 ][ 1 ] << ", " << current_grad_u[ 0 ][ 2 ] << ", "
                    << current_grad_u[ 1 ][ 0 ] << ", " << current_grad_u[ 1 ][ 1 ] << ", " << current_grad_u[ 1 ][ 2 ] << ", "
                    << current_grad_u[ 2 ][ 0 ] << ", " << current_grad_u[ 2 ][ 1 ] << ", " << current_grad_u[ 2 ][ 2 ] << "\n";

        //Output the current micro displacement
        output_file << current_phi[ 0 ] << ", " << current_phi[ 1 ] << ", " << current_phi[ 2 ] << ", "
                    << current_phi[ 3 ] << ", " << current_phi[ 4 ] << ", " << current_phi[ 5 ] << ", "
                    << current_phi[ 6 ] << ", " << current_phi[ 7 ] << ", " << current_phi[ 8 ] << "\n";

        //Output the current gradient of the micro displacement
        output_file << current_grad_phi[ 0 ][ 0 ] << ", " << current_grad_phi[ 0 ][ 1 ] << ", " << current_grad_phi[ 0 ][ 2 ] << ", "
                    << current_grad_phi[ 1 ][ 0 ] << ", " << current_grad_phi[ 1 ][ 1 ] << ", " << current_grad_phi[ 1 ][ 2 ] << ", "
                    << current_grad_phi[ 2 ][ 0 ] << ", " << current_grad_phi[ 2 ][ 1 ] << ", " << current_grad_phi[ 2 ][ 2 ] << ", "
                    << current_grad_phi[ 3 ][ 0 ] << ", " << current_grad_phi[ 3 ][ 1 ] << ", " << current_grad_phi[ 3 ][ 2 ] << ", "
                    << current_grad_phi[ 4 ][ 0 ] << ", " << current_grad_phi[ 4 ][ 1 ] << ", " << current_grad_phi[ 4 ][ 2 ] << ", "
                    << current_grad_phi[ 5 ][ 0 ] << ", " << current_grad_phi[ 5 ][ 1 ] << ", " << current_grad_phi[ 5 ][ 2 ] << ", "
                    << current_grad_phi[ 6 ][ 0 ] << ", " << current_grad_phi[ 6 ][ 1 ] << ", " << current_grad_phi[ 6 ][ 2 ] << ", "
                    << current_grad_phi[ 7 ][ 0 ] << ", " << current_grad_phi[ 7 ][ 1 ] << ", " << current_grad_phi[ 7 ][ 2 ] << ", "
                    << current_grad_phi[ 8 ][ 0 ] << ", " << current_grad_phi[ 8 ][ 1 ] << ", " << current_grad_phi[ 8 ][ 2 ] << "\n";

        //Output the PK2 stress
        output_file << PK2_result[ 0 ] << ", " << PK2_result[ 1 ] << ", " << PK2_result[ 2 ] << ", "
                    << PK2_result[ 3 ] << ", " << PK2_result[ 4 ] << ", " << PK2_result[ 5 ] << ", "
                    << PK2_result[ 6 ] << ", " << PK2_result[ 7 ] << ", " << PK2_result[ 8 ] << "\n";

        //Output the SIGMA stress
        output_file << SIGMA_result[ 0 ] << ", " << SIGMA_result[ 1 ] << ", " << SIGMA_result[ 2 ] << ", "
                    << SIGMA_result[ 3 ] << ", " << SIGMA_result[ 4 ] << ", " << SIGMA_result[ 5 ] << ", "
                    << SIGMA_result[ 6 ] << ", " << SIGMA_result[ 7 ] << ", " << SIGMA_result[ 8 ] << "\n";

        //Output the M stress
        output_file << M_result[  0 ] << ", " << M_result[  1 ] << ", " << M_result[  2 ] << ", "
                    << M_result[  3 ] << ", " << M_result[  4 ] << ", " << M_result[  5 ] << ", "
                    << M_result[  6 ] << ", " << M_result[  7 ] << ", " << M_result[  8 ] << ", "
                    << M_result[  9 ] << ", " << M_result[ 10 ] << ", " << M_result[ 11 ] << ", "
                    << M_result[ 12 ] << ", " << M_result[ 13 ] << ", " << M_result[ 14 ] << ", "
                    << M_result[ 15 ] << ", " << M_result[ 16 ] << ", " << M_result[ 17 ] << ", "
                    << M_result[ 18 ] << ", " << M_result[ 19 ] << ", " << M_result[ 20 ] << ", "
                    << M_result[ 21 ] << ", " << M_result[ 22 ] << ", " << M_result[ 23 ] << ", "
                    << M_result[ 24 ] << ", " << M_result[ 25 ] << ", " << M_result[ 26 ] << "\n";

        //Determine if there were non-linear iterations
        auto inc = DEBUG.find( "converged_values" );
        if ( inc != DEBUG.end() ){
            PK2_intermediate   = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediatePK2Stress" ];
            SIGMA_intermediate = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceMicroStress" ];
            M_intermediate     = DEBUG[ "converged_values" ][ "converged_values" ][ "intermediateReferenceHigherOrderStress" ];
        }
        else{
            PK2_intermediate   = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediatePK2Stress" ];
            SIGMA_intermediate = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceMicroStress" ];
            M_intermediate     = DEBUG[ "pre_iteration_values" ][ "pre_iteration_values" ][ "intermediateReferenceHigherOrderStress" ];
        }

        //Output the intermediate PK2 stress
        output_file << PK2_intermediate[ 0 ] << ", " << PK2_intermediate[ 1 ] << ", " << PK2_intermediate[ 2 ] << ", "
                    << PK2_intermediate[ 3 ] << ", " << PK2_intermediate[ 4 ] << ", " << PK2_intermediate[ 5 ] << ", "
                    << PK2_intermediate[ 6 ] << ", " << PK2_intermediate[ 7 ] << ", " << PK2_intermediate[ 8 ] << "\n";

        //Output the intermediate SIGMA stress
        output_file << SIGMA_intermediate[ 0 ] << ", " << SIGMA_intermediate[ 1 ] << ", " << SIGMA_intermediate[ 2 ] << ", "
                    << SIGMA_intermediate[ 3 ] << ", " << SIGMA_intermediate[ 4 ] << ", " << SIGMA_intermediate[ 5 ] << ", "
                    << SIGMA_intermediate[ 6 ] << ", " << SIGMA_intermediate[ 7 ] << ", " << SIGMA_intermediate[ 8 ] << "\n";

        //Output the intermediate M stress
        output_file << M_intermediate[  0 ] << ", " << M_intermediate[  1 ] << ", " << M_intermediate[  2 ] << ", "
                    << M_intermediate[  3 ] << ", " << M_intermediate[  4 ] << ", " << M_intermediate[  5 ] << ", "
                    << M_intermediate[  6 ] << ", " << M_intermediate[  7 ] << ", " << M_intermediate[  8 ] << ", "
                    << M_intermediate[  9 ] << ", " << M_intermediate[ 10 ] << ", " << M_intermediate[ 11 ] << ", "
                    << M_intermediate[ 12 ] << ", " << M_intermediate[ 13 ] << ", " << M_intermediate[ 14 ] << ", "
                    << M_intermediate[ 15 ] << ", " << M_intermediate[ 16 ] << ", " << M_intermediate[ 17 ] << ", "
                    << M_intermediate[ 18 ] << ", " << M_intermediate[ 19 ] << ", " << M_intermediate[ 20 ] << ", "
                    << M_intermediate[ 21 ] << ", " << M_intermediate[ 22 ] << ", " << M_intermediate[ 23 ] << ", "
                    << M_intermediate[ 24 ] << ", " << M_intermediate[ 25 ] << ", " << M_intermediate[ 26 ] << "\n";
        
        //Output the state variables
        for ( unsigned int i = 0; i < SDVS.size()-1; i++ ){
            output_file << SDVS[ i ] << ", ";
        }
        output_file << SDVS[ SDVS.size() - 1 ] << "\n";

#endif

        BOOST_CHECK( errorCode <= 0 );
#ifdef DEBUG_MODE
        if ( errorCode > 0 ){
            output_file.close();
        }
#endif

        t += dt;

        grad_u_prev   = grad_u_curr;
        phi_prev      = phi_curr;
        grad_phi_prev = grad_phi_curr;

#ifdef DEBUG_MODE    
        if ( t > 0.25 ){ output_file.close(); return 1; }
#endif
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( SDVSAnswer, SDVS ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer, PK2_result ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMAAnswer, SIGMA_result ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( MAnswer, M_result ) );
    
#ifdef DEBUG_MODE
    output_file.close();
#endif
}
