/*--------------------------------------------------------------------------*/
/*------------------------------ File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Smoke test for SatellitesBlock: constructs a SatellitesBlock via the Block
 * factory, checking that the module links correctly and the class is
 * registered. Replace it with real tests exercising the module.
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
 // construct a SatellitesBlock via the Block factory: this checks that the
 // class is registered and the library is linked in (whole-archive)
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
/*---------------------------- End File test.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
