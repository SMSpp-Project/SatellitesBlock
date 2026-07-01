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
#include "MultiTargetBlockv2.h"
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

MultiTargetBlockv2 * oCONSTv2 = nullptr;    // original MultiTargetBlock

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
  namespace stdfs = std::filesystem;
  std::string file = argv[ 1 ];
  std::string solver = argv[ 2 ];

  std::ifstream ifile( file );

  if( ! ifile.is_open() )
   throw( std::invalid_argument( "can't open input file" ) );

  oCONSTv2 = dynamic_cast< MultiTargetBlockv2 * >( Block::new_Block( "MultiTargetBlockv2" ) );
  assert( oCONSTv2 );

  oCONSTv2->load( file );
  std::cout << "Instance charged!\n";

  auto bsc = dynamic_cast< BlockSolverConfig * >(Configuration::deserialize( solver ) );
    if( ! bsc ) {
    cerr << "Error: configuration file not a BlockSolverConfig" << endl;
    exit( 1 );    
    }

  bsc->apply( oCONSTv2 );
  bsc->clear();

  if( oCONSTv2->get_registered_solvers().empty() ) {
      cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
      exit( 1 );    
    }

  Solver * slvr = oCONSTv2->get_registered_solvers().front();
  auto c_start_chrono = std::chrono::high_resolution_clock::now();//std::clock();
  int rtrn = slvr->compute();
  auto c_end_chrono = std::chrono::high_resolution_clock::now();//std::clock();
  double time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end_chrono - c_start_chrono).count()/1e+9;
  std::cout << "MAXIMUM REVISIT TIME INTERVAL LB: " << slvr->get_lb() << "\n";
  std::cout << "MAXIMUM REVISIT TIME INTERVAL UB: " << slvr->get_ub() << "\n";
  std::cout << "TOTAL CPU TIME: " << time_elapsed << "\n";

  return 0;
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
