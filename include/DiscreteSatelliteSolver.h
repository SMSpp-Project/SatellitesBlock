/*--------------------------------------------------------------------------*/
/*------------------ File DiscreteSatelliteSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the DiscreteSatelliteSolver class.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DiscreteSatelliteSolver
 #define __DiscreteSatelliteSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "DiscreteSatelliteBlock.h"

#include "BlockSolverConfig.h"

#include "DQuadFunction.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup DiscreteSatelliteSolver_CLASSES Classes in DiscreteSatelliteSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS DiscreteSatelliteSolver --------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for DiscreteSatelliteBlock
/** DiscreteSatelliteSolver plays, for a DiscreteSatelliteBlock, the same
 * role that SatelliteSolver [see SatelliteSolver.h] plays for a
 * SatelliteBlock: it is meant to be the pricing Solver attached to each
 * DiscreteSatelliteBlock while the corresponding DiscreteConstellationBlock
 * is tackled by a Lagrangian decomposition approach, and it solves the
 * subproblem "by inspection" rather than via a general-purpose solver.
 *
 * Since the only Variable of a DiscreteSatelliteBlock are the binary
 * y[ i ][ j ] (i indexing the candidate orbit, j the discretized
 * observability level, see DiscreteSatelliteBlock.h) subject to the single
 * constraint sum_{i,j} y[ i ][ j ] <= 1, and its (dualized) Objective is a
 * plain LinearFunction with one coefficient eta per y[ i ][ j ] (laid out
 * in row-major (i,j) order), the subproblem trivially reduces to picking
 * the single most negative coefficient eta_min (if any is negative; doing
 * nothing, i.e., all y == 0, otherwise), which is exactly what compute()
 * does by linear scan. */

class DiscreteSatelliteSolver : public Solver {

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
/*----------------- CONSTRUCTING AND DESTRUCTING DiscreteSatelliteSolver ---*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DiscreteSatelliteSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from BenBound. */

 DiscreteSatelliteSolver( void ) : Solver() { }

/*--------------------------------------------------------------------------*/
 /// destructor: it has to release all the Modifications

 virtual ~DiscreteSatelliteSolver() { }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, DiscreteSatelliteSolver maps the parameters of [CDA]Solver
**/

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override
 {

  if( block == f_Block )  // actually doing nothing
   return;                // cowardly and silently return

  Solver::set_Block( block );  // attach to the new Block

  if( block ) {  // this is not just resetting everything
   auto SATB = dynamic_cast< DiscreteSatelliteBlock * >( f_Block );
   if( ! SATB )
    throw( std::invalid_argument(
		         "DiscreteSatelliteSolver:set_Block: block must be a DiscreteSatelliteBlock" ) );

   bool owned = SATB->is_owned_by( f_id );
   if( ( ! owned ) && ( !SATB->read_lock() ) )
    throw( std::logic_error( "cannot acquire read_lock on DiscreteSatelliteBlock" ) );
   // load the new DiscreteSatelliteBlock into the :BenBound object

   // once done, read_unlock the DiscreteSatelliteBlock (if it was read-lock()-ed)
   if( ! owned )
    SATB->read_unlock();

   // TODO: maybe log it
   }
  }  // end( set_Block )

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

  lock();  // first of all, acquire self-lock

  if( ! f_Block )           // there is no [DiscreteSatelliteBlock] to solve
   return( kBlockLocked );  // return error 

  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure
  
  if( ! owned )             // if the [Satellite]Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it

  process_outstanding_Modification();

  auto SATB = dynamic_cast< DiscreteSatelliteBlock * >( f_Block );
  int index = 0;
  eta_min = 0.0;
  double eta;

  auto o = SATB->get_orbits();
  auto ell = SATB->get_ell();

  // linear scan of the o * ell coefficients of the (dense) Objective, one
  // per y[ j ][ tt ], in the same row-major order in which
  // DiscreteSatelliteBlock::generate_objective() laid them out; eta_min
  // (initialized to 0, i.e., the "do nothing" cost) and (o_min,ell_min)
  // track the best (y , orbit) pair found so far

  for(int j = 0; j < o ; j++){
    for(int tt = 0; tt < ell ; tt++){
      eta = static_cast< LinearFunction *>(static_cast< FRealObjective *>
                ( SATB->get_objective())->get_function())->get_coefficient(index);
      if(eta < eta_min){
        o_min = j;
        ell_min = tt;
        eta_min = eta;
      }
      index++;
    }
  }

  unlock();                  // unlock the mutex    
  
  return( kOK );
 }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */
 
/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override {  

  return( eta_min );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { 
 
  return( eta_min );  
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  OFValue get_var_value( void ) override { return( get_ub() ); }

/*--------------------------------------------------------------------------*/

bool has_var_solution( void ) override {
  return( true );
}

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration * solc = nullptr ) override
 {

  if( ! f_Block )  // no [SingleFlowDCR]Block to write to
   return;         // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  auto SATB = static_cast< DiscreteSatelliteBlock * >( f_Block );
  auto o = SATB->get_orbits();
  auto ell = SATB->get_ell();

  // resets all y[][] to 0, then activates only the single (orbit, level)
  // pair found by compute(), and only if it actually improved on "do
  // nothing" (eta_min < 0); otherwise the all-zero solution is left in place

  for(int j = 0; j < o ; j++)
    for(int tt = 0; tt < ell ; tt++)
      SATB->set_y(j , tt , 0.0);

  if(eta_min < 0.0)
    SATB->set_y(o_min , ell_min , 1.0);

}

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

/*--------------------------------------------------------------------------*/
  
  void process_outstanding_Modification( void ) {

    bool reload = false;

      // note: since processing the Modification is fast, we don't bother with
      // being nice to other processes and do it all with v_mod under lock
       // try to acquire lock, spin on failure
      while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
       ;
     
      // process all the Modifications
      for( auto mod : v_mod )
        if( auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod.get() ) ) {
          reload = true;  // a reset must be done
          break;          // ignore all the remaining Modifications
        }
     
      v_mod.clear();  // all Modifications tackled, clear the list
     
      f_mod_lock.clear( std::memory_order_release );  // release lock
     
      if( reload ){

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

  double eta_min;
  int o_min;
  int ell_min;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( DiscreteSatelliteSolver ) )
}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* DiscreteSatelliteSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File DiscreteSatelliteSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/





