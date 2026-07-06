/*--------------------------------------------------------------------------*/
/*--------------------------- File test_milp.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * MILP test for SatellitesBlock: loads a small SCDP instance into a
 * ConstellationBlock, configures a *MILPSolver on it out of a
 * BlockSolverConfig txt file and solves it, checking that a finite optimum
 * is found. The instance and the configuration can be overridden on the
 * command line.
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

#include "ConstellationBlock_discrete.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 std::string instance = argc > 1 ? argv[ 1 ] : "inputheur";
 std::string config = argc > 2 ? argv[ 2 ] : "MILPPar.txt";

 // construct a ConstellationBlock via the factory and load the instance
 auto block = dynamic_cast< ConstellationBlock_discrete * >(
                                    Block::new_Block( "ConstellationBlock_discrete" ) );
 if( ! block ) {
  std::cerr << "ConstellationBlock_discrete not present in Block factory" << std::endl;
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
  std::cout << "SatellitesBlock MILP: all tests passed" << std::endl;
  return( 0 );
  }

 std::cout << "SatellitesBlock MILP: test failed" << std::endl;
 return( 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test_milp.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
