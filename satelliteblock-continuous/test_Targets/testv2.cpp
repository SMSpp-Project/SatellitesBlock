/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing 
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_LEVEL 0
// 0 = only pass/fail
// 1 = result of each test
// 2 = + solver log

#if( LOG_LEVEL >= 1 )
 #define LOG1( x ) cout << x
 #define CLOG1( y , x ) if( y ) cout << x

 #if( LOG_LEVEL >= 2 )
  #define LOG_ON_COUT 1
  // if nonzero, the NDO Solver log is sent on cout rather than on a file
 #endif
#else
 #define LOG1( x )
 #define CLOG1( y , x )
#endif

/*--------------------------------------------------------------------------*/

#define USECOLORS 1
#if( USECOLORS )
 #define RED( x ) "\x1B[31m" #x "\033[0m"
 #define GREEN( x ) "\x1B[32m" #x "\033[0m"
#else
 #define RED( x ) #x
 #define GREEN( x ) #x
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <iomanip>

#include <random>
#include <filesystem>  

#include <chrono>

#include "BlockSolverConfig.h"
#include "SingleTargetBlock.h"
#include "MultiTargetBlock.h"
#include "CPXMILPSolver.h"
#include "GRBMILPSolver.h"
#include "MILPSolver.h"
#include "UpdateSolver.h"
#include "LagrangianDualSolver.h"
#include "BundleSolver.h"
#include "PrimalProximalHeur.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;
using namespace SMSpp_di_unipi_it;

MultiTargetBlock * oCONST = nullptr;    // original MultiTargetBlock
MultiTargetBlock * oCONST1 = nullptr;    // original MultiTargetBlock
SingleTargetBlock * oSAT = nullptr;    // original SingleTargetBlock

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
  namespace stdfs = std::filesystem;
  std::string file = argv[ 1 ];
  std::string solver = argv[ 2 ];

  std::ifstream ifile( file );

  if( ! ifile.is_open() )
   throw( std::invalid_argument( "can't open input file" ) );

  oCONST = dynamic_cast< MultiTargetBlock * >( Block::new_Block( "MultiTargetBlock" ) );
  assert( oCONST );

  oCONST->load( file );
  std::cout << "Instance charged!\n";

  auto bsc = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( solver ) );
    if( ! bsc ) {
    cerr << "Error: configuration file not a BlockSolverConfig" << endl;
    exit( 1 );    
    }

  bsc->apply( oCONST );
  bsc->clear();

  if( oCONST->get_registered_solvers().empty() ) {
      cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
      exit( 1 );    
    }

  Solver * slvr = oCONST->get_registered_solvers().front();
  auto c_start_chrono = std::chrono::high_resolution_clock::now();//std::clock();
  int rtrn = slvr->compute();
  auto c_end_chrono = std::chrono::high_resolution_clock::now();//std::clock();
  double time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end_chrono - c_start_chrono).count()/1e+9;
  std::cout << "MAXIMUM REVISIT TIME INTERVAL LB: " << slvr->get_lb() << "\n";
  std::cout << "MAXIMUM REVISIT TIME INTERVAL UB: " << slvr->get_ub() << "\n";
  std::cout << "TOTAL CPU TIME: " << time_elapsed << "\n";

  for(int j = 0; j < oCONST->get_numTarget()-1; j++){
    for(int i = 0; i < oCONST->get_numSat(); i++){
      for(int k = 0; k < oCONST->get_numOrbits(); k++){
        if(std::abs(oCONST->get_activations( j , i , k ) - oCONST->get_activations( j + 1 , i , k )) > 1e-4){
          std::cout << "SOLUTION INFEASIBLE!" << std::endl;
          break;
        }
      }
    }
  }

  oCONST1 = dynamic_cast< MultiTargetBlock * >( Block::new_Block( "MultiTargetBlock" ) );
    oCONST1->load( file );
    std::cout << "Instance charged!\n";

    auto bsc1 = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( "MILPPar3.txt" ) );
    if( ! bsc1 ) {
      cerr << "Error: configuration file not a BlockSolverConfig" << endl;
      exit( 1 );    
    }

    bsc1->apply( oCONST1 );
    //bsc->clear();
    //((oCONST1->get_registered_solvers()).front())->set_log( &cout );

    slvr->get_var_solution();
    for(int s = 0; s < oCONST1->get_numSat(); s++){
      double orbits_new = 0.0;
      //std::cout << oCONST->get_numOrbits( s ) << std::endl;
      for(int jj = 0; jj < oCONST1->get_numOrbits(); jj++){
        bool orbit_flag = true;
        for(int m = 0; m < oCONST1->get_numTarget(); m++){
          if(oCONST->get_activations(m,s,jj) > 1e-6){
            orbit_flag = false;
            break; 
          }
        }
        if(orbit_flag){
          for(int m = 0; m < oCONST1->get_numTarget(); m++){
            //std::cout << oCONST->get_activations(s,jj) << std::endl;
            oCONST1->set_activations(m,s,jj,0);
            //std::cout << oCONST1->get_activations(s,jj) << std::endl;
          }
          //if(oCONST->get_zs(s) == 0)
          //  oCONST1->set_zs(s,0);
        }
        orbit_flag = true;
        for(int m = 0; m < oCONST1->get_numTarget(); m++){
          if(oCONST->get_activations(m,s,jj) < 1-1e-6){
            orbit_flag = false;
            break; 
          }
        }
        if(orbit_flag){
          for(int m = 0; m < oCONST1->get_numTarget(); m++){
            //std::cout << oCONST->get_activations(s,jj) << std::endl;
            oCONST1->set_activations(m,s,jj,1);
            //std::cout << oCONST1->get_activations(s,jj) << std::endl;
          }
          //if(oCONST->get_zs(s) == 0)
          //  oCONST1->set_zs(s,0);
        }
        std::cout << "orbits_new: " << orbits_new << std::endl;
      }
    }

    Solver * slvr1 = oCONST1->get_registered_solvers().front();
    c_start_chrono = std::chrono::high_resolution_clock::now();//std::clock();


    rtrn = slvr1->compute();
    c_end_chrono = std::chrono::high_resolution_clock::now();//std::clock();
    double time_elapsed_exact = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end_chrono - c_start_chrono).count()/1e+9;
    std::cout << "REVISIT TIME LB: " << slvr1->get_lb() << "\n";
    std::cout << "REVISIT TIME UB: " << slvr1->get_ub() << "\n";
    std::cout << "CPU TIME: " << time_elapsed_exact << "\n";

    std::cout << "TOTAL CPU TIME: " << time_elapsed + time_elapsed_exact << "\n\n";


  return 0;
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
