/*--------------------------------------------------------------------------*/
/*------------------------- File SatelliteBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SatelliteBlock class.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SatelliteBlock.h"
#include <iomanip>
#include <cmath>

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
#define CHECK_DS 0
/* Perform long and costly checks on the data structures representing the
  * abstract and the physical representations agree. */
#else
#define CHECK_DS 0
// never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using FNumber = SatelliteBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const double RAYON = 6378136.3;      // mean Earth radius, in meters
static const double M_limit = 10 * 3.14159265;  // unused generic big-M bound
static const auto PI = 3.14159265;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

static void print_UB( std::ostream & os , FNumber ub )
{
 if( ub == Inf< FNumber >() )
  os << "+Inf";
 else
  os << ub;
}

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ

template < typename T > static Index countdiff( T beg , T end , T cmp )
{
 Index ndiff = 0;
 for( ; beg != end ;)
  if( *( beg++ ) != *( cmp++ ) )
   ndiff++;

 return ( ndiff );
}

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template < typename T >
static bool is_equal( std::vector< T > & vec , c_Subset & nms ,
                      typename std::vector< T >::const_iterator cmp,
                      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *( cmp++ ) )
   return ( false );
 }

 return ( true );
}

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices

template < typename T >
static Index countdiff( std::vector< T > & vec , c_Subset & nms ,
                        typename std::vector< T >::const_iterator cmp,
                        Index n_max )
{
 Index ndiff = 0;
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *( cmp++ ) )
   ndiff++;
 }

 return ( ndiff );
}

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another

template < typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
                     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *( cpy++ );
}

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SatelliteBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SatelliteBlock );

// register SatelliteSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( SatelliteSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF SatelliteBlock --------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( FNumber num_targets , FNumber time_step ,
                           FNumber horizon, FNumber altValues,
                           FNumber thetaValues, FNumber indOrbit,
                           FNumber aHalf,
                           boost::multi_array< double, 3 > CoverageLat,
                           boost::multi_array< double, 3 > CoverageLong,
                           std::vector< double > periods_tgt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 n = num_targets;
 OrbitSet = indOrbit;  // number of candidate orbital configurations

 thetaVal = thetaValues;  // the theta^{\max} threshold, cf. generate_abstract_constraints()

 dt = time_step;
 T = horizon;
 t = T / dt;  // number of discrete time stamps in the horizon
 alphaHalf = aHalf;

 periods = periods_tgt;

 // CoverageSatLat/Long[ i ][ j ][ jj ] hold, for target i, time stamp j and
 // orbital configuration jj, the (precomputed) angular distance in
 // latitude/longitude between the satellite ground track and the target;
 // these are the \Delta lat[] and \Delta long[] terms of constraint (1) in
 // the class comments, and are copied over from the caller-provided arrays
 // (typically built by ConstellationBlock::load(), see there)

 CoverageSatLat.resize( boost::extents[ num_targets ][ t ][ OrbitSet ] );
 CoverageSatLong.resize( boost::extents[ num_targets ][ t ][ OrbitSet ] );

 for( Index jj = 0 ; jj < OrbitSet ; ++jj ) {
  for( Index j = 0 ; j < t ; ++j ) {
   for( Index i = 0 ; i < num_targets ; ++i ) {
    CoverageSatLat[ i ][ j ][ jj ] = CoverageLat[ i ][ j ][ jj ];
    CoverageSatLong[ i ][ j ][ jj ] = CoverageLong[ i ][ j ][ jj ];
   }
  }
 }

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 //if( anyone_there() )
 // add_Modification( std::make_shared< NBModification >( this ) );

} // end( SatelliteBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

} // end( SatelliteBlock::load( istream ) )

/*--------------------------------------------------------------------------*/
// creates the four groups of static Variable of the SatelliteBlock: the
// continuous thetaVar (the theta^{\max} threshold actually used, which
// depends on the selected orbit), the binary zeta (satellite activation),
// the binary activation[] (one per candidate orbital configuration, this
// is \pi[] in the class comments) and the binary xi[][] (one per target
// per time stamp, this is \xi[][] in the class comments)

void SatelliteBlock::generate_abstract_variables( Configuration * stvv )
{
 if( AR3 & HasVar ) // the variables are there already
  return; // nothing to do

 thetaVar.resize( 1 );
 for( auto & var : thetaVar )
  var.set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( thetaVar );

 zeta.resize( 1 );
 for( auto & var : zeta )
  var.set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( zeta );

 activation.resize( OrbitSet );
 for( auto & var : activation )
  var.set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( activation );

 xi.resize( boost::extents[ n ][ t ] );
 for( Index i = 0 ; i < n ; ++i )
  for( Index j = 0 ; j < t ; ++j )
   xi[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( xi );

 AR3 |= HasVar;

} // end( SatelliteBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( AR2 & HasCnst ) // the constraints are there already
  return; // nothing to do

 // thetaUB: thetaVar <= thetaVal * zeta, i.e., if the satellite is not
 // active (zeta == 0) then the threshold theta^{\max} is forced to 0, so
 // that no target can possibly be observed (cf. the obs2_cnst/obs4_cnst
 // constraints below)

 thetaUB.resize( 1 );
 LinearFunction::v_coeff_pair var_theta1;
 var_theta1.push_back( std::make_pair( &thetaVar[ 0 ] , 1.0 ) );
 var_theta1.push_back( std::make_pair( &zeta[ 0 ] , -thetaVal ) );
 LinearFunction * FunctTheta1 = new LinearFunction( std::move( var_theta1 ) );
 thetaUB[ 0 ].set_rhs( 0.0 , eNoBlck );
 thetaUB[ 0 ].set_lhs( -Inf< double >() , eNoBlck );
 thetaUB[ 0 ].set_function( FunctTheta1 , eNoBlck );

 add_static_constraint( thetaUB );

 // thetaLB: thetaVar >= ( thetaVal / 20 ) * zeta, i.e., if the satellite
 // is active then theta^{\max} is bounded away from 0 by a fixed fraction
 // of thetaVal (this keeps the linearized observability constraints below
 // from becoming numerically degenerate)

 thetaLB.resize( 1 );
 LinearFunction::v_coeff_pair var_theta2;
 var_theta2.push_back( std::make_pair( &thetaVar[ 0 ] , 1.0 ) );
 var_theta2.push_back( std::make_pair( &zeta[ 0 ] , -thetaVal / 20 ) );
 LinearFunction * FunctTheta2 = new LinearFunction( std::move( var_theta2 ) );
 thetaLB[ 0 ].set_rhs( Inf< double >() , eNoBlck );
 thetaLB[ 0 ].set_lhs( 0.0 , eNoBlck );
 thetaLB[ 0 ].set_function( FunctTheta2 , eNoBlck );

 add_static_constraint( thetaLB );

 // generate the orbitSelection constraint, which select exctly one orbit
 // configuration for the satellite: sum_{j \in OrbitSet} activation[ j ] == 1
 // remember: activation[ j ] = 1 iff. the j-th configuration is selected

 orbitSelection.resize( 1 );
 LinearFunction::v_coeff_pair orbit_var;
 for( Index j = 0 ; j < OrbitSet ; ++j ) {
  orbit_var.push_back( std::make_pair( &activation[ j ] , 1.0 ) );
 }
 orbit_var.push_back( std::make_pair( &zeta[ 0 ] , -1.0 ) );
 LinearFunction * FunctAnm = new LinearFunction( std::move( orbit_var ) );
 orbitSelection[ 0 ].set_rhs( 0.0 , eNoBlck );
 orbitSelection[ 0 ].set_lhs( 0.0 , eNoBlck );
 orbitSelection[ 0 ].set_function( FunctAnm , eNoBlck );

 add_static_constraint( orbitSelection );

 // generate the activationSat_cnst, which active the current satellite if
 // there exists a target m is observed by this satellite at a given time step:
 // zeta \leq xi[ i ][ j ], for all targets i's and time steps j's
 // remember: zeta = 1 iff. the current satellite is active in the constellation
 // and xi[ i ][ j ] = 1 iff. satellite observe target i at time step j

 activationSat_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ] );
 for( Index i = 0 ; i < n ; ++i ) {
  for( Index j = 0 ; j < t ; ++j ) {
   LinearFunction::v_coeff_pair v_var;
   v_var.push_back( std::make_pair( &xi[ i ][ j ] , -1.0 ) );
   v_var.push_back( std::make_pair( &zeta[ 0 ] , 1.0 ) );
   LinearFunction * FunctSat = new LinearFunction( std::move( v_var ) );
   activationSat_cnst[ i ][ j ].set_rhs( Inf< double >() , eNoBlck );
   activationSat_cnst[ i ][ j ].set_lhs( 0.0 , eNoBlck );
   activationSat_cnst[ i ][ j ].set_function( FunctSat , eNoBlck );
  }
 }
 add_static_constraint( activationSat_cnst );

 // generate the converse aggregate constraint, forcing the satellite to be
 // inactive if it observes no target at all:
 // zeta \leq sum_{i,j} xi[ i ][ j ]
 // together with activationSat_cnst above, this makes zeta == 1 iff. the
 // satellite observes at least one target at some time step

 activationSat_cnst_1.resize( 1 );

 LinearFunction::v_coeff_pair v_var_z;
 v_var_z.push_back( std::make_pair( &zeta[ 0 ] , 1.0 ) );
 for( Index i = 0 ; i < n ; ++i ) {
  for( Index j = 0 ; j < t ; ++j ) {
   v_var_z.push_back( std::make_pair( &xi[ i ][ j ] , -1.0 ) );
  }
 }
 LinearFunction * FunctSat_1 = new LinearFunction( std::move( v_var_z ) );
 activationSat_cnst_1[ 0 ].set_rhs( 0.0 , eNoBlck );
 activationSat_cnst_1[ 0 ].set_lhs( -Inf< double >() , eNoBlck );
 activationSat_cnst_1[ 0 ].set_function( FunctSat_1 , eNoBlck );

 add_static_constraint( activationSat_cnst_1 );

 // generate observability (linearize) constraints via big-M approac
 // these constraints traslate the fact that a target is observed by
 // the current satellite, i.e., xi[ i ][ j ] = 1, if the (scaled)
 // distance between the projection of the satellite position onto
 // the Earth surface and the position of the target is smaller than
 // a threshold theta^{\max} with respect to Latitude and Longitude.

 // theta^{\max} + (1 - xi[ i ][ j ]) * MLAT \geq
 // sum_{jj \in OrbitSet} activation[ jj ] * CoverageSatLat[ i ][ j ][ jj ]
 // for all targets m's and time steps j's

 // theta^{\max} + (1 - xi[ i ][ j ]) * MLONG \geq
 // sum_{jj \in OrbitSet} activation[ jj ] * CoverageSatLong[ i ][ j ][ jj ]
 // for all targets m's and time steps j's

 // MLAT and MLONG are two big-M parameters automatically computed
 // such that their numerical values are the smallest to guarantee
 // that constraint are valid (redundant when xi[ i ][ j ] = 0)

 obs2_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ] );
 obs4_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ] );

 double MLAT = PI;
 double MLONG = 2 * PI;

 for( Index i = 0 ; i < n ; ++i ) {
  for( Index j = 0 ; j < t ; ++j ) {
   // MLAT / MLONG are recomputed, for each (target, time stamp) pair, as
   // the largest coverage distance over all orbital configurations: since
   // exactly one activation[ jj ] is 1 (cf. orbitSelection), this is the
   // smallest big-M value for which the constraint is guaranteed to be
   // redundant whenever xi[ i ][ j ] == 0, whatever the selected orbit is
   MLAT = 0.0;
   MLONG = 0.0;
   LinearFunction::v_coeff_pair v_obs1, v_obs2;
   for( Index jj = 0 ; jj < OrbitSet ; ++jj ) {
    v_obs1.push_back(
     std::make_pair( &activation[ jj ] , -CoverageSatLat[ i ][ j ][ jj ] ) );
    v_obs2.push_back(
     std::make_pair( &activation[ jj ] , -CoverageSatLong[ i ][ j ][ jj ] ) );
    MLAT = std::max( MLAT , ( CoverageSatLat[ i ][ j ][ jj ] ) );
    MLONG = std::max( MLONG , ( CoverageSatLong[ i ][ j ][ jj ] ) );
   }

   v_obs1.push_back( std::make_pair( &xi[ i ][ j ] , -MLAT ) );
   v_obs2.push_back( std::make_pair( &xi[ i ][ j ] , -MLONG ) );

   v_obs1.push_back( std::make_pair( &thetaVar[ 0 ] , 1.0 ) );
   v_obs2.push_back( std::make_pair( &thetaVar[ 0 ] , 1.0 ) );

   // theta^{\max} - sum_{jj} activation[ jj ] * Coverage[ i ][ j ][ jj ]
   //   - MLAT * xi[ i ][ j ] >= -MLAT
   // i.e., theta^{\max} - Coverage(selected orbit) >= -MLAT * ( 1 - xi ),
   // which is vacuous when xi == 0 and forces the coverage distance to be
   // within the threshold when xi == 1 (same reasoning for MLONG)

   LinearFunction * Funct1 = new LinearFunction( std::move( v_obs1 ) );
   LinearFunction * Funct2 = new LinearFunction( std::move( v_obs2 ) );

   obs2_cnst[ i ][ j ].set_rhs( Inf< double >() , eNoBlck );
   obs2_cnst[ i ][ j ].set_lhs( -MLAT , eNoBlck );
   obs2_cnst[ i ][ j ].set_function( Funct1 , eNoBlck );

   obs4_cnst[ i ][ j ].set_rhs( Inf< double >() , eNoBlck );
   obs4_cnst[ i ][ j ].set_lhs( -MLONG , eNoBlck );
   obs4_cnst[ i ][ j ].set_function( Funct2 , eNoBlck );
  }
 }

 add_static_constraint( obs2_cnst );
 add_static_constraint( obs4_cnst );

 AR2 |= HasCnst;
} // end( SatelliteBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
// the objective of a single SatelliteBlock is just "minimize zeta"; since
// Block objectives sum up over a Block tree and ConstellationBlock does not
// define an Objective of its own (see ConstellationBlock.h), the overall
// problem obtained by putting several SatelliteBlock together under one
// ConstellationBlock amounts to minimizing the number of active satellites
// subject to every target being observed (the observability constraints
// defined by ConstellationBlock)

void SatelliteBlock::generate_objective( Configuration * objc )
{
 if( AR1 & HasObj ) // the objective is there already
  return; // cowardly (and silently) return

 LinearFunction::v_coeff_pair p( 1 );

 p[ 0 ].first = &zeta[ 0 ];
 p[ 0 ].second = 1.0;

 c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
 set_objective( &c , eNoMod );

 AR1 |= HasObj;

} // end( SatelliteBlock::generate_objective )

/*--------------------------------------------------------------------------*/
// TODO: the tolerance eps is extracted from fsbc / f_BlockConfig exactly as
// documented in the header, but the actual feasibility check against the
// abstract/physical representation (cf. MCFBlock::is_feasible() for the
// pattern) is not implemented yet: the method unconditionally reports the
// Block as infeasible

bool SatelliteBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 FNumber eps = 0;
 auto tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >( fsbc );

 if( ( !tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
   f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return ( 0 );

} // end( SatelliteBlock::is_feasible )

/*--------------------------------------------------------------------------*/
// TODO: the tolerances ceps and feps are extracted exactly as documented in
// the header, but, like is_feasible(), the actual optimality check is not
// implemented yet: the method unconditionally reports the Block as
// non-optimal

bool SatelliteBlock::is_optimal( bool useabstract , Configuration * optc )
{
 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  if( auto toptc =
       dynamic_cast< SimpleConfiguration< std::pair< CNumber, FNumber > > * >(
        optc ) ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
  } else {
   auto ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >( optc );

   if( ( !ttoptc ) && f_BlockConfig &&
       f_BlockConfig->f_is_optimal_Configuration )
    ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >(
     f_BlockConfig->f_is_optimal_Configuration );
   if( ttoptc )
    ceps = ttoptc->f_value;

   if( f_BlockConfig && f_BlockConfig->f_is_feasible_Configuration ) {
    auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
     f_BlockConfig->f_is_feasible_Configuration );
    if( fsbc )
     feps = fsbc->f_value;
   }
  }
 } else if( f_BlockConfig ) {
  if( f_BlockConfig->f_is_optimal_Configuration )
   if( auto csbc = dynamic_cast< SimpleConfiguration< CNumber > * >(
        f_BlockConfig->f_is_optimal_Configuration ) )
    ceps = csbc->f_value;

  if( f_BlockConfig->f_is_feasible_Configuration )
   if( auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
        f_BlockConfig->f_is_feasible_Configuration ) )
    feps = fsbc->f_value;
 }

 return ( 0 );

} //  end( SatelliteBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
// creates a new SatelliteSolution, read from the current status of the
// SatelliteBlock unless an empty one (emptys == true) is asked for; wsol is
// extracted from solc / f_BlockConfig but is currently unused since
// SatelliteSolution only stores zeta (see SatelliteSolution::read())

Solution * SatelliteBlock::get_Solution( Configuration * solc , bool emptys )
{
 int wsol = 0;
 if( ( !solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;

 auto * sol = new SatelliteSolution();

 if( !emptys )
  sol->read( this );

 return ( sol );

} // end( SatelliteBlock::get_Solution )

/*--------------------------------------------------------------------------*/
// sets the value of the zeta ColVariable (there is only one, zeta being a
// scalar) from the range [ rng.first , rng.second ) of fstrt; note that the
// loop variable is named xi purely as a local iterator, unrelated to the
// xi[][] member Variable

void SatelliteBlock::set_zeta( c_Vec_FNumber_it fstrt , Range rng )
{
 if( !( AR3 & HasVar ) ) // nowhere to put the value in
  return; // cowardly (and silently) return

 Index i = rng.first;

 for( auto xi = zeta.begin() + i ; i < 1 ; ++i )
  ( xi++ )->set_value( *( fstrt++ ) );

} // end( SatelliteBlock::set_zeta )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
// intercepts Modification concerning this Block (as opposed to some nested
// sub-Block) before forwarding them to the base Block::add_Modification();
// "concerns_Block( false )" is set so that Block::add_Modification() does
// not process it again as if it were still to be dispatched

void SatelliteBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
 }

 Block::add_Modification( mod , chnl );

} // end( SatelliteBlock::add_Modification )

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE SatelliteBlock ---*/
/*--------------------------------------------------------------------------*/
// TODO: printing the SatelliteBlock instance is not implemented yet

void SatelliteBlock::print( std::ostream & output , char vlvl ) const {

} // end( SatelliteBlock::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints
 Constraint::clear( orbitSelection ); // static
 Constraint::clear( activationSat_cnst ); // static
 Constraint::clear( obs2_cnst );
 Constraint::clear( obs4_cnst );

 c.clear(); // clear the Objective

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 //reset_dynamic_constraints();
 //reset_dynamic_variables();
 reset_objective();

} // end( SatelliteBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/
// dispatches "abstract" Modification originated by changes to the Variable
// / Constraint / Objective of this SatelliteBlock; currently no Modification
// is supported (all the Variable/Constraint of a SatelliteBlock are meant to
// be immutable once generated), so this is a no-op stub: TODO throw an
// exception here once the set of legal Modification (if any) is decided,
// mirroring what MCFBlock::guts_of_add_Modification() does

void SatelliteBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -

} // end( SatelliteBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF SatelliteSolution ---------------------*/
/*--------------------------------------------------------------------------*/
// a SatelliteSolution only stores the value of the zeta (satellite
// activation) Variable of a SatelliteBlock; deserialize()/serialize() are
// currently no-ops, i.e., loading/saving a SatelliteSolution to/from a
// netCDF::NcGroup is not implemented yet

void SatelliteSolution::deserialize( const netCDF::NcGroup & group ) {}

/*--------------------------------------------------------------------------*/
// reads the current value of zeta out of the given SatelliteBlock into
// v_zeta; note that this only happens if v_zeta is already non-empty (i.e.,
// this SatelliteSolution has previously been sized to hold zeta), and that
// the value returned by SATB->get_zeta() is presently discarded: TODO store
// it into v_zeta[ 0 ]

void SatelliteSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const SatelliteBlock * >( block );
 if( !SATB )
  throw( std::invalid_argument( "block is not a SatelliteBlock" ) );

 if( !v_zeta.empty() ) {
  v_zeta.resize( 1 );

  SATB->get_zeta();
 }
}

/*--------------------------------------------------------------------------*/
// writes the value of zeta stored in v_zeta into the given SatelliteBlock

void SatelliteSolution::write( Block * block )
{
 auto SATB = dynamic_cast< SatelliteBlock * >( block );
 if( !SATB )
  throw( std::invalid_argument( "block is not a SatelliteBlock" ) );

 if( !v_zeta.empty() ) {
  SATB->set_zeta( v_zeta.begin() );
 }
}

/*--------------------------------------------------------------------------*/

void SatelliteSolution::serialize( netCDF::NcGroup & group ) const {}

/*--------------------------------------------------------------------------*/
// zeta being a binary (0/1) activation flag, "scaling" it by factor makes
// little sense: an empty clone is returned instead, i.e., scale() amounts
// to discarding the stored value

SatelliteSolution * SatelliteSolution::scale( double factor ) const
{
 auto * sol = SatelliteSolution::clone( true );
 return ( sol );
}

/*--------------------------------------------------------------------------*/
// TODO: summing SatelliteSolution (e.g. for Lagrangian aggregation) is not
// implemented yet

void SatelliteSolution::sum( const Solution * solution , double multiplier ) {}

/*--------------------------------------------------------------------------*/
// if empty is true, only allocates v_zeta (if this SatelliteSolution has
// one) without copying its value; otherwise deep-copies v_zeta as well

SatelliteSolution * SatelliteSolution::clone( bool empty ) const
{
 auto * sol = new SatelliteSolution();

 if( empty ) {
  if( !v_zeta.empty() )
   sol->v_zeta.resize( 1 );
 } else {
  sol->v_zeta = v_zeta;
 }

 return ( sol );
}

/*--------------------------------------------------------------------------*/
/*------------------- End File SatelliteBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
