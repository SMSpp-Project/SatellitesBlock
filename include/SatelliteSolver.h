/*--------------------------------------------------------------------------*/
/*------------------------- File SatelliteSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SatelliteSolver class, a Solver that solves the
 * Lagrangian subproblem of a SatelliteBlock by enumeration of its orbits.
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

#ifndef __SatelliteSolver
 #define __SatelliteSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Solver.h"

#include "SatelliteBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE & USING ----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SatelliteSolver_CLASSES Classes in SatelliteSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SatelliteSolver --------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for SatelliteBlock
/** SatelliteSolver is a Solver that exploits the structure of a single
 * SatelliteBlock to solve it "by inspection", without any general-purpose
 * MILP solver. It is the Solver attached to each SatelliteBlock when the
 * ConstellationBlock is solved by Lagrangian decomposition (e.g., by
 * LagrangianDualSolver and BundleSolver): the constraints of
 * ConstellationBlock, the only ones linking different SatelliteBlock, are
 * dualized, so that the Lagrangian function decomposes into one subproblem
 * per SatelliteBlock, whose Objective is a LinearFunction with coefficient
 * lambdaz for zeta, lambda3 for thetaVar and lambda11[ m ][ t ] for every
 * xi[ m ][ t ] (any of them may be missing, i.e., zero).
 *
 * The subproblem is solved with the threshold fixed to its maximum value
 * \f$ \theta^{\max} \f$ (see SatelliteBlock::get_theta()). For every
 * candidate orbit c, the xi[][] are then fixed: xi[ m ][ t ] = 1 iff both
 * distances Delta_lat and Delta_long of orbit c from target m at time t are
 * within \f$ \theta^{\max} \f$, with cost sum_{ m , t } lambda11[ m ][ t ]
 * xi[ m ][ t ]. The orbit with the smallest cost is chosen if this cost is
 * negative, the satellite being inactive (all the Variable zero, cost 0)
 * otherwise; the value of the solution is the constant term of the
 * Objective plus min( lambdaz + cost + lambda3 \f$ \theta^{\max} \f$ , 0 ).
 * Since the threshold is fixed, the result is optimal only among the
 * solutions with \f$ \theta = \theta^{\max} \f$ for the active satellite.
 * Since compute() always reads the whole Objective, the Modification
 * received by the Solver are simply discarded. */

class SatelliteSolver : public Solver {

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
/*-------------- CONSTRUCTING AND DESTRUCTING SatelliteSolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SatelliteSolver
 *  @{ */

 /// constructor: initializes the (empty) solution

 SatelliteSolver( void )
  : Solver() , objective_opt( 0 ) , orbit_opt( -1 ) , lambda3( 0 ) ,
    lambdaz( 0 ) ,
    f_has_sol( false ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: there is nothing to release

 virtual ~SatelliteSolver() {}

/** @} ---------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** Sets the Block that the Solver has to solve, which must be a
  * SatelliteBlock (exception is thrown otherwise); any solution found for
  * the previous Block is discarded. */

 void set_Block( Block * block ) override {
  if( block == f_Block ) // actually doing nothing
   return;               // cowardly and silently return

  if( block && ( ! dynamic_cast< SatelliteBlock * >( block ) ) )
   throw( std::invalid_argument(
    "SatelliteSolver::set_Block: block must be a SatelliteBlock" ) );

  Solver::set_Block( block ); // attach to the new Block
  f_has_sol = false;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the SatelliteBlock
 *  @{ */

 /// solve the SatelliteBlock
 /** Solves the SatelliteBlock by enumeration of its candidate orbits [see
  * the class comments]. Returns kOK on success, kBlockLocked if the Block
  * cannot be read_lock()-ed, and kError if there is no Block or if its
  * Objective is not a FRealObjective with a LinearFunction of the Variable
  * zeta, thetaVar and xi[][]. */

 int compute( bool changedvars = true ) override {
  lock(); // first of all, acquire self-lock

  if( ! f_Block ) { // there is no SatelliteBlock to solve
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
  auto SATB = dynamic_cast< SatelliteBlock * >( f_Block );
  auto obj = SATB ? dynamic_cast< FRealObjective * >( SATB->get_objective() )
                  : nullptr;
  auto lf = obj ? dynamic_cast< LinearFunction * >( obj->get_function() )
                : nullptr;

  if( lf && read_multipliers( SATB , lf ) ) {
   const auto t = SATB->get_t();
   const auto n = SATB->get_n();
   const auto theta = SATB->get_theta();

   // enumeration of the candidate orbits: for orbit i, xi_new1[ j ][ tt ]
   // is 1 iff target j is within the threshold at time stamp tt, and
   // xi_sum1 is the corresponding cost; the best orbit is kept only if it
   // improves on the cost 0 of "doing nothing"

   boost::multi_array< double , 2 > xi_new1( boost::extents[ n ][ t ] );
   xi_new.resize( boost::extents[ n ][ t ] );
   std::fill_n( xi_new.data() , xi_new.num_elements() , 0.0 );

   double objective_opt1 = 0.0;
   orbit_opt = -1;

   for( Index i = 0 ; i < SATB->get_numOrbit() ; ++i ) {
    double xi_sum1 = 0.0;
    for( Index j = 0 ; j < n ; ++j )
     for( Index tt = 0 ; tt < t ; ++tt )
      if( ( SATB->get_Delta_lat( i , j , tt ) <= theta ) &&
          ( SATB->get_Delta_long( i , j , tt ) <= theta ) ) {
       xi_sum1 += lambda11[ j ][ tt ];
       xi_new1[ j ][ tt ] = 1.0;
       }
      else
       xi_new1[ j ][ tt ] = 0.0;

    if( xi_sum1 < objective_opt1 ) {
     objective_opt1 = xi_sum1;
     orbit_opt = i;
     xi_new = xi_new1;
     }
    }

   // the value also includes the costs of zeta and of the threshold;
   // "doing nothing" (cost 0) is always possible

   objective_opt =
    std::min( lambdaz + objective_opt1 + lambda3 * theta , 0.0 );
   if( ( orbit_opt == -1 ) || ( objective_opt >= 0.0 ) ) {
    objective_opt = 0.0;
    orbit_opt = -1;
    }
   objective_opt += lf->get_constant_term();

   f_has_sol = true;
   status = kOK;
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

 /// returns the value found by the last call to compute()
 /** Returns the value found by the last call to compute() [see the class
  * comments], which is used as both lower and upper bound; before any
  * successful compute() this is -Inf. Note that, the threshold being fixed
  * to \f$ \theta^{\max} \f$, the value is the optimum of the SatelliteBlock
  * restricted to the solutions with that threshold, hence an upper bound
  * on its optimum but not, in general, a lower bound: a Lagrangian dual of
  * ConstellationBlock computed with this Solver on the SatelliteBlock is
  * therefore not, in general, a lower bound on the optimum of the
  * ConstellationBlock. */

 OFValue get_lb( void ) override {
  return( f_has_sol ? objective_opt : -Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns the value found by the last call to compute()
 /** Same as get_lb(), but returning +Inf before any successful compute(). */

 OFValue get_ub( void ) override {
  return( f_has_sol ? objective_opt : Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --*/
 /// returns the value of the solution found by the last call to compute()

 OFValue get_var_value( void ) override { return( get_ub() ); }

/*--------------------------------------------------------------------------*/
 /// returns true if the last call to compute() has found a solution

 bool has_var_solution( void ) override { return( f_has_sol ); }

/*--------------------------------------------------------------------------*/
 /// writes the solution found by compute() in the SatelliteBlock
 /** Writes the solution found by the last call to compute() in the
  * Variable of the SatelliteBlock: if an orbit has been chosen, zeta is 1,
  * the activation of that orbit is 1 (all the others 0), the threshold is
  * \f$ \theta^{\max} \f$ and the xi[][] are those of the orbit; otherwise,
  * all the Variable are 0. Nothing is done if there is no solution, or if
  * solc is a SimpleConfiguration< int > with value 2. */

 void get_var_solution( Configuration * solc = nullptr ) override {
  if( ( ! f_Block ) || ( ! f_has_sol ) ) // no solution to write
   return;                             // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  auto SATB = static_cast< SatelliteBlock * >( f_Block );
  const bool active = orbit_opt > -1;

  SATB->set_zeta1( active ? 1 : 0 );
  SATB->set_thetaVar1( active ? SATB->get_theta() : 0 );

  for( Index i = 0 ; i < SATB->get_numOrbit() ; ++i )
   SATB->set_activation1( i , 0 );
  if( active )
   SATB->set_activation1( orbit_opt , 1 );

  for( Index j = 0 ; j < SATB->get_n() ; ++j )
   for( Index k = 0 ; k < SATB->get_t() ; ++k )
    SATB->set_xi1( j , k , active ? xi_new[ j ][ k ] : 0.0 );
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
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// reads lambdaz, lambda3 and lambda11[][] out of the Objective
 /** Reads the coefficients of zeta (lambdaz), of thetaVar (lambda3) and of
  * the xi[][] (lambda11[][]) out of the LinearFunction lf of the Objective
  * of SATB, those of the missing Variable being 0. Returns false if lf has
  * some other Variable. */

 bool read_multipliers( SatelliteBlock * SATB , LinearFunction * lf ) {
  const auto t = SATB->get_t();
  const auto n = SATB->get_n();

  lambdaz = lambda3 = 0;
  lambda11.resize( boost::extents[ n ][ t ] );
  std::fill_n( lambda11.data() , lambda11.num_elements() , 0.0 );

  const ColVariable * xi0 = ( n && t ) ? SATB->i2p_r( 0 , 0 ) : nullptr;
  for( const auto & [ var , coeff ] : lf->get_v_var() ) {
   if( var == SATB->i2p_z() )
    lambdaz = coeff;
   else if( var == SATB->i2p_theta() )
    lambda3 = coeff;
   else {
    const auto pos = xi0 ? var - xi0 : -1;
    if( ( pos < 0 ) || ( Index( pos ) >= n * t ) )
     return( false );
    lambda11[ pos / t ][ pos % t ] = coeff;
    }
   }

  return( true );
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 double objective_opt; ///< the value found by compute()
 int orbit_opt;        ///< the chosen orbit, -1 if the satellite is inactive

 boost::multi_array< double , 2 > xi_new;   ///< the xi[][] of orbit_opt
 boost::multi_array< double , 2 > lambda11; ///< the costs of the xi[][]
 double lambda3;                           ///< the cost of thetaVar
 double lambdaz;                           ///< the cost of zeta

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

 }; // end( class( SatelliteSolver ) )

/** @} end( group( SatelliteSolver_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SatelliteSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File SatelliteSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
