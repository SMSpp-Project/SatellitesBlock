/*--------------------------------------------------------------------------*/
/*--------------------- File DiscreteSatelliteSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the DiscreteSatelliteSolver class, a Solver that solves
 * a DiscreteSatelliteBlock "by inspection".
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DiscreteSatelliteSolver
 #define __DiscreteSatelliteSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "DiscreteSatelliteBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE & USING ----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup DiscreteSatelliteSolver_CLASSES Discrete satellite Solver
 *  @{ */

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteSatelliteSolver ----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for DiscreteSatelliteBlock
/** DiscreteSatelliteSolver plays, for a DiscreteSatelliteBlock, the role
 * that SatelliteSolver [see SatelliteSolver.h] plays for a SatelliteBlock:
 * it is the Solver attached to each DiscreteSatelliteBlock when the
 * DiscreteConstellationBlock is solved by Lagrangian decomposition, and it
 * solves the subproblem "by inspection" rather than via a general-purpose
 * MILP solver.
 *
 * The only Variable of a DiscreteSatelliteBlock are the binary y[ o ][ l ]
 * (o indexing the candidate orbit, l the level of the observability
 * threshold, see DiscreteSatelliteBlock.h), subject to the single
 * constraint sum_{ o , l } y[ o ][ l ] <= 1, and its Objective is a
 * LinearFunction of the y[][]. The subproblem therefore amounts to setting
 * to 1 the y with the smallest coefficient, if it is negative, and all the
 * y to 0 otherwise, which is what compute() does with a linear scan. Since
 * compute() always reads the whole Objective, the Modification received by
 * the Solver are simply discarded. */

class DiscreteSatelliteSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

 using Index = Block::Index; ///< the type of the indices

/** @} ---------------------------------------------------------------------*/
/*---------- CONSTRUCTING AND DESTRUCTING DiscreteSatelliteSolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DiscreteSatelliteSolver
 *  @{ */

 /// constructor: initializes the (empty) solution

 DiscreteSatelliteSolver( void )
  : Solver() , f_value( 0 ) , eta_min( 0 ) , o_min( 0 ) , ell_min( 0 ) ,
    f_has_sol( false ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: there is nothing to release

 virtual ~DiscreteSatelliteSolver() {}

/** @} ---------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** Sets the Block that the Solver has to solve, which must be a
  * DiscreteSatelliteBlock (exception is thrown otherwise); any solution
  * found for the previous Block is discarded. */

 void set_Block( Block * block ) override {
  if( block == f_Block ) // actually doing nothing
   return;               // cowardly and silently return

  if( block && ( ! dynamic_cast< DiscreteSatelliteBlock * >( block ) ) )
   throw( std::invalid_argument( "DiscreteSatelliteSolver::set_Block: "
                                 "block must be a DiscreteSatelliteBlock" ) );

  Solver::set_Block( block ); // attach to the new Block
  f_has_sol = false;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the DiscreteSatelliteBlock
 *  @{ */

 /// solve the DiscreteSatelliteBlock
 /** Solves the DiscreteSatelliteBlock by a linear scan of the coefficients
  * of its Objective. Returns kOK on success, kBlockLocked if the Block
  * cannot be read_lock()-ed, and kError if there is no Block or if its
  * Objective is not a FRealObjective with a LinearFunction on the
  * y[][] Variable. */

 int compute( bool changedvars = true ) override {
  lock(); // first of all, acquire self-lock

  if( ! f_Block ) { // there is no DiscreteSatelliteBlock to solve
   unlock();
   return( kError );
   }

  bool owned = f_Block->is_owned_by( f_id );      // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) ) { // if not try to lock
   unlock();
   return( kBlockLocked ); // return error on failure
   }

  process_outstanding_Modification();
  f_has_sol = false;

  int status = kError;
  auto SATB = dynamic_cast< DiscreteSatelliteBlock * >( f_Block );
  auto obj = SATB ? dynamic_cast< FRealObjective * >( SATB->get_objective() )
                  : nullptr;
  auto lf = obj ? dynamic_cast< LinearFunction * >( obj->get_function() )
                : nullptr;

  if( lf ) {
   // linear scan of the coefficients: the position of y[ o ][ l ] in the
   // (row-major) storage of the y[][] is o * ell + l; eta_min starts from
   // 0, the cost of "doing nothing"
   const Index ell = SATB->get_ell();
   const Index ny = SATB->get_orbits() * ell;
   const ColVariable * y0 = ny ? SATB->i2p_y( 0 , 0 ) : nullptr;
   eta_min = 0;
   o_min = ell_min = 0;
   status = kOK;
   for( const auto & [ var , eta ] : lf->get_v_var() ) {
    const auto pos = var - y0;
    if( ( ! y0 ) || ( pos < 0 ) || ( Index( pos ) >= ny ) ) {
     status = kError; // a Variable that is not one of the y[][]
     break;
     }
    if( eta < eta_min ) {
     o_min = Index( pos ) / ell;
     ell_min = Index( pos ) % ell;
     eta_min = eta;
     }
    }

   if( status == kOK ) {
    f_value = lf->get_constant_term() + eta_min;
    f_has_sol = true;
    }
   }

  if( ! owned )             // if the Block was actually read_locked
   f_Block->read_unlock(); // read_unlock it

  unlock(); // unlock the mutex

  return( status );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 /// returns the optimal value found by the last call to compute()
 /** Returns the optimal value found by the last call to compute(); since
  * the DiscreteSatelliteBlock is solved exactly, lower and upper bound
  * coincide. Before any successful compute() this is -Inf. */

 OFValue get_lb( void ) override {
  return( f_has_sol ? f_value : -Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns the optimal value found by the last call to compute()
 /** Same as get_lb(), but returning +Inf before any successful compute(). */

 OFValue get_ub( void ) override {
  return( f_has_sol ? f_value : Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns the value of the solution found by the last call to compute()

 OFValue get_var_value( void ) override { return( get_ub() ); }

/*--------------------------------------------------------------------------*/
 /// returns true if the last call to compute() has found a solution

 bool has_var_solution( void ) override { return( f_has_sol ); }

/*--------------------------------------------------------------------------*/
 /// writes the solution found by compute() in the DiscreteSatelliteBlock
 /** Writes the solution found by the last call to compute() in the y[][]
  * Variable of the DiscreteSatelliteBlock: all the y are set to 0, except
  * y[ o_min ][ ell_min ] which is set to 1 if its coefficient is negative.
  * Nothing is done if there is no solution, or if solc is a
  * SimpleConfiguration< int > with value 2. */

 void get_var_solution( Configuration * solc = nullptr ) override {
  if( ( ! f_Block ) || ( ! f_has_sol ) ) // no solution to write
   return;                             // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  auto SATB = static_cast< DiscreteSatelliteBlock * >( f_Block );

  for( Index o = 0 ; o < SATB->get_orbits() ; ++o )
   for( Index l = 0 ; l < SATB->get_ell() ; ++l )
    SATB->set_y( o , l , 0 );

  if( eta_min < 0 )
   SATB->set_y( o_min , ell_min , 1 );
  }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /// discards all the outstanding Modification
 /** Since compute() reads the whole Objective every time, the
  * Modification need not be processed: they are just discarded. */

 void process_outstanding_Modification( void ) {
  // try to acquire lock, spin on failure
  while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
   ;

  v_mod.clear(); // all Modification tackled, clear the list

  f_mod_lock.clear( std::memory_order_release ); // release lock
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 double f_value; ///< the optimal value found by compute()
 double eta_min; ///< the smallest coefficient, or 0 if none is negative
 Index o_min;    ///< the orbit o of the y[ o ][ l ] with smallest cost
 Index ell_min;  ///< the level l of the y[ o ][ l ] with smallest cost
 bool f_has_sol; ///< true if compute() has found a solution

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*----------------------------- PRIVATE FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class( DiscreteSatelliteSolver ) )

/** @} end( group( DiscreteSatelliteSolver_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* DiscreteSatelliteSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File DiscreteSatelliteSolver.h -------------------*/
/*--------------------------------------------------------------------------*/
