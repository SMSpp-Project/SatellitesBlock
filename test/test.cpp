/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of SatellitesBlock, using nothing but the core SMS++ library.
 *
 * Each of the Block that are the root of a model of the module, i.e.,
 * ConstellationBlock, DiscreteConstellationBlock and MultiTargetBlock, is
 * constructed by the Block factory and loaded with the instance of the
 * module it is meant for [see data/txt], and then:
 *
 * - its abstract representation, and that of all its sub-Block, is
 *   generated with a FakeSolver registered to it, which has to receive no
 *   Modification at all: the abstract representation is part of the
 *   construction of the Block, and a Modification issued while generating
 *   it reaches whoever listens, e.g., the LagBFunction of a Lagrangian
 *   decomposition, which does not know the Block yet and throws;
 *
 * - the instance is loaded again into the same Block, and the abstract
 *   representation generated again has to have the same number of static
 *   Variable and Constraint, in the Block and in each sub-Block.
 *
 * The solution of the models, by a :MILPSolver and by the Lagrangian
 * decomposition with the Solver of this module on the sub-Block, is
 * cross-checked in the SatellitesBlock suite of the tests of the umbrella.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <string>
#include <vector>

#include "Block.h"
#include "FakeSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// generates the abstract representation of block and of all its sub-Block

static void generate_all( Block * block )
{
 block->generate_abstract_variables();
 for( auto sb : block->get_nested_Blocks() )
  generate_all( sb );
 block->generate_abstract_constraints();
 block->generate_objective();
 }

/*--------------------------------------------------------------------------*/
/// the number of static Variable and Constraint of block and its sub-Block

static std::vector< Block::Index > shape( Block * block )
{
 std::vector< Block::Index > s = { block->get_number_static_variables() ,
                                   block->get_number_static_constraints() };
 for( auto sb : block->get_nested_Blocks() ) {
  auto ss = shape( sb );
  s.insert( s.end() , ss.begin() , ss.end() );
  }
 return( s );
 }

/*--------------------------------------------------------------------------*/
/// the checks on the Block named classname, loaded with instance

static bool check( const std::string & classname ,
                   const std::string & instance )
{
 auto block = Block::new_Block( classname );
 if( ! block ) {
  std::cout << classname << ": not present in the Block factory"
            << std::endl;
  return( false );
  }

 bool ok = true;
 try {
  block->load( instance );
  if( block->get_nested_Blocks().empty() ) {
   std::cout << classname << ": no sub-Block after load()" << std::endl;
   ok = false;
   }

  auto fake = new FakeSolver();
  block->register_Solver( fake );
  generate_all( block );
  if( ! fake->get_Modification_list().empty() ) {
   std::cout << classname << ": " << fake->get_Modification_list().size()
             << " Modification issued while generating the abstract "
             << "representation" << std::endl;
   ok = false;
   }
  block->unregister_Solvers( true );

  const auto first = shape( block );
  block->load( instance );
  generate_all( block );
  if( shape( block ) != first ) {
   std::cout << classname << ": a second load() gives a different abstract "
             << "representation" << std::endl;
   ok = false;
   }
  }
 catch( std::exception & e ) {
  std::cout << classname << ": " << e.what() << std::endl;
  ok = false;
  }

 delete block;

 if( ok )
  std::cout << classname << ": OK" << std::endl;
 return( ok );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // the directory of the instances, data/txt of the module by default
 const std::string data = argc > 1 ? argv[ 1 ] : "../data/txt";

 bool ok = check( "ConstellationBlock" , data + "/inputheur-const" );
 ok = check( "DiscreteConstellationBlock" , data + "/inputheur-const" ) && ok;
 ok = check( "MultiTargetBlock" , data + "/inputheur-target" ) && ok;

 if( ok )
  std::cout << "SatellitesBlock: all tests passed" << std::endl;
 else
  std::cout << "SatellitesBlock: some test failed" << std::endl;

 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
