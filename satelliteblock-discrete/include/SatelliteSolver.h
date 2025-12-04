/*--------------------------------------------------------------------------*/
/*------------------ File SatelliteSolver.h ---------------------*/
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

//#include "MCFClass.h"

//#include "MILPSolver.h"

#include "BlockSolverConfig.h"

#include "DQuadFunction.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
  
 //using namespace MCFClass_di_unipi_it;
 //using Index = Block::Index;
 
 //class SatelliteSolverState;  // forward declaration of SatelliteSolverState

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteSolver_CLASSES Classes in SatelliteSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SatelliteSolver --------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for SatelliteBlock

class SatelliteSolver : public Solver {

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
/*----------------- CONSTRUCTING AND DESTRUCTING SatelliteSolver */
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SatelliteSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special, except verifying that the
  * template argument derives from BenBound. */

 SatelliteSolver( void ) : Solver() { }

/*--------------------------------------------------------------------------*/
 /// destructor: it has to release all the Modifications

 virtual ~SatelliteSolver() { }

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

  if( block == f_Block )  // actually doing nothing
   return;                // cowardly and silently return

  Solver::set_Block( block );  // attach to the new Block

  if( block ) {  // this is not just resetting everything
   auto SATB = dynamic_cast< SatelliteBlock * >( f_Block );
   if( ! SATB )
    throw( std::invalid_argument(
		         "SatelliteSolver:set_Block: block must be a SatelliteBlock" ) );

   bool owned = SATB->is_owned_by( f_id );
   if( ( ! owned ) && ( !SATB->read_lock() ) )
    throw( std::logic_error( "cannot acquire read_lock on SatelliteBlock" ) );
   // load the new SatelliteBlock into the :BenBound object

   // once done, read_unlock the SatelliteBlock (if it was read-lock()-ed)
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

  if( ! f_Block )           // there is no [SatelliteBlock] to solve
   return( kBlockLocked );  // return error 

  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure
  
  if( ! owned )             // if the [Satellite]Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it

  process_outstanding_Modification();

  /****** COMPUTE PROCEDURE ******/

  auto SATB = dynamic_cast< SatelliteBlock * >( f_Block );
  int index = 0;
  eta_min = 0.0;
  double eta;

  auto o = SATB->get_orbits();
  auto ell = SATB->get_ell();

  for(int j = 0; j < o ; j++){
    for(int tt = 0; tt < ell ; tt++){
      eta = static_cast< LinearFunction *>(static_cast< FRealObjective *>( SATB->get_objective())->get_function())->get_coefficient(index);
      if(eta < eta_min){
        o_min = j;
        ell_min = tt;
        eta_min = eta;
      }
      index++;
    }
  }

  //std::cout << eta_min << " " << ell_min << " " << o_min << std::endl;

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

  auto SATB = static_cast< SatelliteBlock * >( f_Block );
  auto o = SATB->get_orbits();
  auto ell = SATB->get_ell();

  
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

};  // end( class( SatelliteSolver ) )
}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SatelliteSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File SatelliteSolver.h -------------------*/
/*--------------------------------------------------------------------------*/





