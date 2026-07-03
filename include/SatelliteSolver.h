/*--------------------------------------------------------------------------*/
/*------------------------- File SatelliteSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SatelliteSolver class, implementing a
 * Solver for Delay-Constrained Routing problems (DCR) relative to a Single
 * Flow, as set by SatelliteBlock, via a "Benders with nested Lagrange"
 * approach.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SatelliteSolver
#define __SatelliteSolver
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "SatelliteBlock.h"

#include "BlockSolverConfig.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

//using namespace MCFClass_di_unipi_it;
//using Index = Block::Index;

//class SatelliteSolverState;  // forward declaration of SatelliteSolverState

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteSolver_CLASSES Classes in SatelliteSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SatelliteSolver --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for SatelliteBlock

class SatelliteSolver : public Solver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Public Types
 *  @{ */

 /** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING SatelliteSolver --------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructing and destructing SatelliteSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from BenBound. */

 SatelliteSolver( void ) : Solver() {}

/*--------------------------------------------------------------------------*/
 /// destructor: it has to release all the Modifications

 virtual ~SatelliteSolver() {}

 /** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Other initializations
 *
 * Parameter-wise, SatelliteSolver maps the parameters of [CDA]Solver
**/

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override
 {
  if( block == f_Block ) // actually doing nothing
   return; // cowardly and silently return

  Solver::set_Block( block ); // attach to the new Block

  if( block ) { // this is not just resetting everything
   auto SATB = dynamic_cast< SatelliteBlock * >( f_Block );
   if( !SATB )
    throw( std::invalid_argument(
     "SatelliteSolver:set_Block: block must be a SatelliteBlock" ) );

   bool owned = SATB->is_owned_by( f_id );
   if( ( !owned ) && ( !SATB->read_lock() ) )
    throw( std::logic_error( "cannot acquire read_lock on SatelliteBlock" ) );
   // load the new SatelliteBlock into the :BenBound object

   // once done, read_unlock the SatelliteBlock (if it was read-lock()-ed)
   if( !owned )
    SATB->read_unlock();

   // TODO: maybe log it
  }
 } // end( set_Block )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Solving the Satellite
 *  @{ */

 /// (try to) solve the Satellite

 int compute( bool changedvars = true ) override
 {
  lock(); // first of all, acquire self-lock

  if( !f_Block ) // there is no [SatelliteBlock] to solve
   return ( kBlockLocked ); // return error

  bool owned = f_Block->is_owned_by( f_id ); // check if already locked
  if( ( !owned ) && ( !f_Block->read_lock() ) ) // if not try to read_lock
   return ( kBlockLocked ); // return error on failure

  if( !owned ) // if the [Satellite]Block was actually read_locked
   f_Block->read_unlock(); // read_unlock it

  process_outstanding_Modification();

  /****** COMPUTE PROCEDURE ******/

  auto SATB = dynamic_cast< SatelliteBlock * >( f_Block );
  auto t = SATB->get_t();
  auto n = SATB->get_n();
  auto theta = SATB->get_theta();

  xi_new.resize( boost::extents[ n ][ t ] );

  zeta = 0.0;

  auto tot_periods = 0;
  sum_lambda = 0.0;
  int max_period = 0;

  for( int k = 0 ; k < n ; k++ ) {
   tot_periods += SATB->get_period( k );
   if( max_period < SATB->get_period( k ) )
    max_period = SATB->get_period( k );
  }

  lambda1.resize( boost::extents[ n ][ max_period ] );
  lambda2.resize( boost::extents[ n ][ t ] );

  for( int j = 0 ; j < n ; j++ ) {
   for( int i = 0 ; i < max_period ; i++ ) {
    lambda1[ j ][ i ] = 0.0;
   }
  }

  lambdaz =
   static_cast< LinearFunction * >(
    static_cast< FRealObjective * >( SATB->get_objective() )->get_function() )
    ->get_coefficient( 0 );
  lambda3 =
   static_cast< LinearFunction * >(
    static_cast< FRealObjective * >( SATB->get_objective() )->get_function() )
    ->get_coefficient( 1 );

  int index = 2;
  for( int j = 0 ; j < n ; j++ ) {
   for( int i = 0 ; i < SATB->get_period( j ) ; i++ ) {
    lambda1[ j ][ i ] =
     static_cast< LinearFunction * >(
      static_cast< FRealObjective * >( SATB->get_objective() )
       ->get_function() )
      ->get_coefficient( index );
    sum_lambda -= lambda1[ j ][ i ];
    index += 1;
   }
  }

  lambda11.resize( boost::extents[ n ][ t ] );

  index = 1;
  index += 1;
  for( int j = 0 ; j < n ; j++ ) {
   for( int tt = 0 ; tt < t ; tt++ ) {
    lambda11[ j ][ tt ] =
     static_cast< LinearFunction * >(
      static_cast< FRealObjective * >( SATB->get_objective() )
       ->get_function() )
      ->get_coefficient( index );
    index += 1;
   }
  }

  objective_opt = 0.0;

  int flag = 0;

  for( int j = 0 ; j < n ; ++j )
   for( int k = 0 ; k < t ; ++k )
    xi_new[ j ][ k ] = 0;

  boost::multi_array< double, 2 > xi_new1;
  xi_new1.resize( boost::extents[ n ][ t ] );

  boost::multi_array< double, 2 > xi_sum1p;
  xi_sum1p.resize( boost::extents[ n ][ max_period ] );

  boost::multi_array< double, 2 > xi_sump1;
  xi_sump1.resize( boost::extents[ n ][ max_period ] );

  double pp;

  double objective_opt1 = 0.0;
  zeta = 1.0;
  xi_sum = 0.0;
  auto xi_sum1 = 0.0;
  auto sum_p = 0.0;
  orbit_opt = -1;

  //std::cout << "NUM_ORBIT: " << SATB->get_numOrbit() << std::endl;

  for( int i = 0 ; i < SATB->get_numOrbit() ; i++ ) {
   xi_sum1 = 0.0;
   sum_p = 0.0;

   for( int j = 0 ; j < n ; j++ ) {
    pp = SATB->get_horizon() / SATB->get_timeStep() / SATB->get_period( j );
    for( int tt = 0 ; tt < t ; tt++ ) {
     if( SATB->get_Delta_lat( i , j , tt ) <= theta and
         SATB->get_Delta_long( i , j , tt ) <= theta ) {
      xi_sum1 += lambda11[ j ][ tt ];
      xi_new1[ j ][ tt ] = 1.0;
      zeta1 = 1.0;
     } else {
      xi_new1[ j ][ tt ] = 0.0;
     }
    }
   }

   if( xi_sum1 < objective_opt1 ) {
    objective_opt1 = xi_sum1;
    orbit_opt = i;
    zeta = zeta1;
    xi_new = xi_new1;
    xi_sum = xi_sum1;
    flag = 1;
   }
  }

  objective_opt = std::min( lambdaz + objective_opt1 + lambda3 * theta , 0.0 );

  unlock(); // unlock the mutex

  return ( kOK );
 }

 /** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Accessing the found solutions (if any)
 *  @{ */

/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override
 {
  if( orbit_opt == -1 )
   objective_opt = 0.0;

  return ( objective_opt );
 }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override
 {
  if( orbit_opt == -1 )
   objective_opt = 0.0;

  return ( objective_opt );
 }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_var_value( void ) override { return ( get_ub() ); }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override { return ( true ); }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration * solc = nullptr ) override
 {
  if( !f_Block ) // no [SingleFlowDCR]Block to write to
   return; // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  auto SATB = static_cast< SatelliteBlock * >( f_Block );
  auto t = SATB->get_t();
  auto n = SATB->get_n();

  if( objective_opt >= 0.0 ) {
   orbit_opt = -1;
  }

  SATB->set_zeta1( 1.0 );

  if( orbit_opt > -1 ) {
   SATB->set_zeta1( 1.0 );
  } else {
   SATB->set_zeta1( 0.0 );
  }

  for( int i = 0 ; i < SATB->get_numOrbit() ; i++ )
   SATB->set_activation1( i , 0 );

  if( orbit_opt > -1 ) {
   SATB->set_activation1( orbit_opt , 1 );
  }

  auto sumx = 0.0;
  for( int j = 0 ; j < n ; j++ ) {
   for( int k = 0 ; k < t ; k++ ) {
    if( orbit_opt > -1 ) {
     SATB->set_xi1( j , k , xi_new[ j ][ k ] );
     sumx += xi_new[ j ][ k ];
    } else {
     SATB->set_xi1( j , k , 0.0 );
    }
   }
  }
 }

 /** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
 /** @name Changing the data of the model
 *  @{ */

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void )
 {
  bool reload = false;

  // note: since processing the Modification is fast, we don't bother with
  // being nice to other processes and do it all with v_mod under lock
  // try to acquire lock, spin on failure
  while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
   ;

  // process all the Modifications
  for( auto mod : v_mod )
   if( auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod.get() ) ) {
    reload = true; // a reset must be done
    break; // ignore all the remaining Modifications
   }

  v_mod.clear(); // all Modifications tackled, clear the list

  f_mod_lock.clear( std::memory_order_release ); // release lock

  if( reload ) {
   auto SATB = dynamic_cast< SatelliteBlock * >( f_Block );
   auto n = SATB->get_n();

   int tot_periods = 0;

   for( int k = 0 ; k < n ; k++ ) {
    tot_periods += SATB->get_period( k );
    //std::cout << "TOT: " << SATB->get_period( k ) << std::endl;
   }
  }
 }

 /** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 double zeta;
 double zeta1;
 double objective_opt;
 double sum_lambda;
 double xi_sum;

 int orbit_opt;

 boost::multi_array< double, 2 > xi_new;
 boost::multi_array< double, 2 > lambda1;
 boost::multi_array< double, 2 > lambda2;
 boost::multi_array< double, 2 > lambda11;
 double lambda3;
 double lambdaz;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:
/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

}; // end( class( SatelliteSolver ) )
} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SatelliteSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File SatelliteSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
