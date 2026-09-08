/*--------------------------------------------------------------------------*/
/*---------------------------- File test_ld.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Lagrangian decomposition test for SatellitesBlock: loads a small SCDP
 * instance into a MultiTargetBlock, configures a LagrangianDualSolver
 * (with a BundleSolver on the Lagrangian dual and a *MILPSolver on the
 * subproblems) on it out of a BlockSolverConfig txt file and solves it,
 * checking that a finite optimum is found. The instance and the
 * configuration can be overridden on the command line.
 *
 * NOTE: this file is not wired into test/CMakeLists.txt nor into the plain
 * makefile in this directory (unlike test.cpp and test_milp.cpp), so it is
 * not built as part of the standard test suite; it must be compiled and
 * linked manually against MultiTargetBlock, LagrangianDualSolver and
 * BundleSolver.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>

#include <iostream>

#include "BlockSolverConfig.h"

#include "MultiTargetBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 std::string instance = argc > 1 ? argv[ 1 ] : "inputheur-target";
 std::string config = argc > 2 ? argv[ 2 ] : "BSMILPPar.txt";

 // construct a MultiTargetBlock via the factory and load the instance
 auto block = dynamic_cast< MultiTargetBlock * >(
                                    Block::new_Block( "MultiTargetBlock" ) );
 if( ! block ) {
  std::cerr << "MultiTargetBlock not present in Block factory" << std::endl;
  return( 1 );
  }

 block->load( instance );

 // configure a Solver on the Block out of the BlockSolverConfig
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                                     Configuration::deserialize( config ) );
 if( ! bsc ) {
  std::cerr << "'" << config << "' is not a BlockSolverConfig" << std::endl;
  delete block;
  return( 1 );
  }

 bsc->apply( block ); 

 if( block->get_registered_solvers().empty() ) {
  std::cerr << "the BlockSolverConfig did not register any Solver"
            << std::endl;
  delete bsc;
  delete block;
  return( 1 );
  }

 // solve
 auto slvr = block->get_registered_solvers().front();
 slvr->compute();

 auto lb = slvr->get_lb();
 auto ub = slvr->get_ub();
 std::cout << "lb = " << lb << ", ub = " << ub << std::endl;

 bool ok = std::isfinite( lb ) && std::isfinite( ub ) &&
           ( ub - lb <= 1e-6 * std::max( 1.0 , std::abs( ub ) ) );

 // cleanup: unregister the Solver and delete everything
 bsc->clear();
 bsc->apply( block );
 delete bsc;
 delete block;

 if( ok ) {
  std::cout << "SatellitesBlock LD: all tests passed" << std::endl;
  return( 0 );
  }

 std::cout << "SatellitesBlock LD: test failed" << std::endl;
 return( 1 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File test_ld.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
