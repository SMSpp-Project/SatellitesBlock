/*--------------------------------------------------------------------------*/
/*-------------------------- File test_block.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Smoke test for SatellitesBlock: constructs a ConstellationBlock via the
 * Block factory, checking that the module links correctly and the class is
 * registered. Replace it with real tests exercising the module.
 *
 * NOTE: this file is currently NOT wired into test/CMakeLists.txt nor into
 * the plain makefile in this directory (both instead build test.cpp as the
 * "smoke"/default test); it is kept here but not built as part of the
 * standard test suite.
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

#include <iostream>

#include "ConstellationBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // construct a ConstellationBlock via the Block factory: this checks that
 // the class is registered and the library is linked in (whole-archive)
 auto block = Block::new_Block( "ConstellationBlock" );

 if( ! block ) {
  std::cerr << "ConstellationBlock not present in Block factory" << std::endl;
  return( 1 );
  }

 if( ! dynamic_cast< ConstellationBlock * >( block ) ) {
  std::cerr << "factory did not return a ConstellationBlock" << std::endl;
  delete block;
  return( 1 );
  }

 delete block;

 std::cout << "SatellitesBlock: all tests passed" << std::endl;

 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File test_block.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
