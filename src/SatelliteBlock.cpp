/*--------------------------------------------------------------------------*/
/*------------------------ File SatelliteBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SatelliteBlock and SatelliteSolution classes.
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Luca Mencarelli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SatelliteBlock.h"

#include <cmath>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

using Range = Block::Range;

using FNumber = SatelliteBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SatelliteBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SatelliteBlock );

// register SatelliteSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( SatelliteSolution );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SatelliteBlock ------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( Index num_targets , FNumber time_step ,
                           FNumber horizon , FNumber altValues ,
                           FNumber thetaValues , Index indOrbit ,
                           FNumber aHalf ,
                           const boost::multi_array< double , 3 > &
                            CoverageLat ,
                           const boost::multi_array< double , 3 > &
                            CoverageLong ,
                           const std::vector< double > & periods_tgt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 n = num_targets;
 OrbitSet = indOrbit;
 thetaVal = thetaValues;
 altitudeVal = altValues;
 dt = time_step;
 T = horizon;
 t = T / dt;
 alphaHalf = aHalf;
 periods = periods_tgt;

 if( ( CoverageLat.shape()[ 0 ] < n ) || ( CoverageLat.shape()[ 1 ] < t ) ||
     ( CoverageLat.shape()[ 2 ] < OrbitSet ) ||
     ( CoverageLong.shape()[ 0 ] < n ) || ( CoverageLong.shape()[ 1 ] < t ) ||
     ( CoverageLong.shape()[ 2 ] < OrbitSet ) )
  throw( std::invalid_argument(
   "SatelliteBlock::load: coverage arrays of the wrong size" ) );

 // CoverageSatLat[ i ][ j ][ jj ] is the latitude distance between the
 // ground track of orbit jj and target i at time stamp j, and the same for
 // CoverageSatLong: these are the Delta lat and Delta long of constraint
 // (1) in the class comments

 CoverageSatLat.resize( boost::extents[ n ][ t ][ OrbitSet ] );
 CoverageSatLong.resize( boost::extents[ n ][ t ][ OrbitSet ] );

 for( Index i = 0 ; i < n ; ++i )
  for( Index j = 0 ; j < t ; ++j )
   for( Index jj = 0 ; jj < OrbitSet ; ++jj ) {
    CoverageSatLat[ i ][ j ][ jj ] = CoverageLat[ i ][ j ][ jj ];
    CoverageSatLong[ i ][ j ][ jj ] = CoverageLong[ i ][ j ][ jj ];
    }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( SatelliteBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( std::istream & input , char frmt )
{
 throw( std::logic_error( "SatelliteBlock::load: loading out of an istream "
                          "is not supported" ) );

 } // end( SatelliteBlock::load( istream ) )

/*--------------------------------------------------------------------------*/
// the four groups of static Variable: thetaVar (the threshold), zeta (the
// activation), activation[] (pi[] in the class comments) and xi[][]

void SatelliteBlock::generate_abstract_variables( Configuration * stvv )
{
 if( AR3 & HasVar ) // the variables are there already
  return;           // nothing to do

 thetaVar.resize( 1 );
 for( auto & var : thetaVar )
  var.set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( thetaVar , "thetaVar" );

 zeta.resize( 1 );
 for( auto & var : zeta )
  var.set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( zeta , "zeta" );

 activation.resize( OrbitSet );
 for( auto & var : activation )
  var.set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( activation , "activation" );

 xi.resize( boost::extents[ n ][ t ] );
 for( Index i = 0 ; i < n ; ++i )
  for( Index j = 0 ; j < t ; ++j )
   xi[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( xi , "xi" );

 AR3 |= HasVar;

 } // end( SatelliteBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( AR2 & HasCnst ) // the constraints are there already
  return;            // nothing to do

 // thetaUB: thetaVar <= thetaVal * zeta, i.e., the threshold is at most
 // the maximum one, and it is 0 if the satellite is not active

 thetaUB.resize( 1 );
 LinearFunction::v_coeff_pair var_theta1;
 var_theta1.push_back( std::make_pair( &thetaVar[ 0 ] , 1.0 ) );
 var_theta1.push_back( std::make_pair( &zeta[ 0 ] , -thetaVal ) );
 LinearFunction * FunctTheta1 = new LinearFunction( std::move( var_theta1 ) );
 thetaUB[ 0 ].set_rhs( 0.0 , eNoBlck );
 thetaUB[ 0 ].set_lhs( -Inf< double >() , eNoBlck );
 thetaUB[ 0 ].set_function( FunctTheta1 , eNoBlck );

 add_static_constraint( thetaUB , "thetaUB" );

 // thetaLB: thetaVar >= ( thetaVal / 20 ) * zeta, i.e., the threshold of
 // an active satellite is at least a fixed fraction of the maximum one

 thetaLB.resize( 1 );
 LinearFunction::v_coeff_pair var_theta2;
 var_theta2.push_back( std::make_pair( &thetaVar[ 0 ] , 1.0 ) );
 var_theta2.push_back( std::make_pair( &zeta[ 0 ] , -thetaVal / 20 ) );
 LinearFunction * FunctTheta2 = new LinearFunction( std::move( var_theta2 ) );
 thetaLB[ 0 ].set_rhs( Inf< double >() , eNoBlck );
 thetaLB[ 0 ].set_lhs( 0.0 , eNoBlck );
 thetaLB[ 0 ].set_function( FunctTheta2 , eNoBlck );

 add_static_constraint( thetaLB , "thetaLB" );

 // orbitSelection, constraint (2) of the class comments: an active
 // satellite has exactly one orbit, sum_j activation[ j ] == zeta

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

 add_static_constraint( orbitSelection , "orbitSelection" );

 // activationSat_cnst, constraints (3) of the class comments: a satellite
 // observing some target is active, xi[ i ][ j ] <= zeta for all targets i
 // and time stamps j

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
 add_static_constraint( activationSat_cnst , "activationSat_cnst" );

 // activationSat_cnst_1, constraint (3) of the class comments: a
 // satellite observing no target is inactive, zeta <= sum_{ i , j }
 // xi[ i ][ j ]; with activationSat_cnst, zeta == 1 iff the satellite
 // observes some target at some time stamp

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

 add_static_constraint( activationSat_cnst_1 , "activationSat_cnst_1" );

 // obs2_cnst and obs4_cnst, the big-M linearization of constraint (1) of
 // the class comments: target i can be observed at time stamp j
 // (xi[ i ][ j ] = 1) only if the latitude and longitude distances of the
 // chosen orbit are within the threshold, i.e.,
 //
 //   theta + ( 1 - xi[ i ][ j ] ) MLAT >=
 //     sum_{ jj } activation[ jj ] CoverageSatLat[ i ][ j ][ jj ]
 //
 // and the same with MLONG and CoverageSatLong; MLAT and MLONG are the
 // smallest values making the constraints redundant when xi[ i ][ j ] = 0

 obs2_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ] );
 obs4_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ] );

 for( Index i = 0 ; i < n ; ++i ) {
  for( Index j = 0 ; j < t ; ++j ) {
   // the largest distance over all the orbits: since at most one
   // activation[ jj ] is 1, this is the smallest big-M making the
   // constraint redundant when xi[ i ][ j ] == 0, whatever the orbit
   double MLAT = 0.0;
   double MLONG = 0.0;
   LinearFunction::v_coeff_pair v_obs1 , v_obs2;
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

 add_static_constraint( obs2_cnst , "obs2_cnst" );
 add_static_constraint( obs4_cnst , "obs4_cnst" );

 AR2 |= HasCnst;
 } // end( SatelliteBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
// the objective is zeta: since ConstellationBlock has no Objective of its
// own, the objective of the problem is the number of active satellites

void SatelliteBlock::generate_objective( Configuration * objc )
{
 if( AR1 & HasObj ) // the objective is there already
  return;           // cowardly (and silently) return

 LinearFunction::v_coeff_pair p( 1 );

 p[ 0 ].first = &zeta[ 0 ];
 p[ 0 ].second = 1.0;

 c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
 set_objective( &c , eNoMod );

 AR1 |= HasObj;

 } // end( SatelliteBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR CHECKING THE Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool SatelliteBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 if( ! ( ( AR3 & HasVar ) && ( AR2 & HasCnst ) ) ) // nothing to check
  return( false );

 FNumber eps = 0;
 auto tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
   f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return( ColVariable::is_feasible( thetaVar , eps ) &&
          ColVariable::is_feasible( zeta , eps ) &&
          ColVariable::is_feasible( activation , eps ) &&
          ColVariable::is_feasible( xi , eps ) &&
          RowConstraint::is_feasible( thetaUB , eps ) &&
          RowConstraint::is_feasible( thetaLB , eps ) &&
          RowConstraint::is_feasible( orbitSelection , eps ) &&
          RowConstraint::is_feasible( activationSat_cnst , eps ) &&
          RowConstraint::is_feasible( activationSat_cnst_1 , eps ) &&
          RowConstraint::is_feasible( obs2_cnst , eps ) &&
          RowConstraint::is_feasible( obs4_cnst , eps ) );

 } // end( SatelliteBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool SatelliteBlock::is_optimal( bool useabstract , Configuration * optc )
{
 return( false );

 } // end( SatelliteBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/

Solution * SatelliteBlock::get_Solution( Configuration * solc , bool emptys )
{
 auto * sol = new SatelliteSolution();

 if( AR3 & HasVar ) {
  sol->v_zeta.assign( 1 , 0 );
  if( ! emptys )
   sol->read( this );
  }

 return( sol );

 } // end( SatelliteBlock::get_Solution )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::set_zeta( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR3 & HasVar ) ) // nowhere to put the value in
  return;                // cowardly (and silently) return

 if( ( rng.first == 0 ) && ( rng.second > 0 ) )
  zeta[ 0 ].set_value( *fstrt );

 } // end( SatelliteBlock::set_zeta )

/*--------------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 } // end( SatelliteBlock::add_Modification )

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR PRINTING & SAVING THE SatelliteBlock ------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::print( std::ostream & output , char vlvl ) const
{
 output << "SatelliteBlock: " << n << " targets, " << t << " time stamps, "
        << OrbitSet << " orbits, maximum threshold " << thetaVal << std::endl;

 } // end( SatelliteBlock::print )

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::guts_of_destructor( void )
{
 // clear() the Constraint and the Objective, so that they do not bother to
 // un-register themselves from Variable that are going to be deleted anyway
 Constraint::clear( orbitSelection );
 Constraint::clear( thetaUB );
 Constraint::clear( thetaLB );
 Constraint::clear( activationSat_cnst );
 Constraint::clear( activationSat_cnst_1 );
 Constraint::clear( obs2_cnst );
 Constraint::clear( obs4_cnst );
 c.clear();

 // explicitly reset all Constraint and Variable, so that a new abstract
 // representation is not added to the (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 thetaVar.clear();
 zeta.clear();
 activation.clear();
 xi.resize( boost::extents[ 0 ][ 0 ] );

 AR1 = AR2 = AR3 = 0;

 } // end( SatelliteBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/
// the SatelliteBlock has no physical representation that an abstract
// Modification may change, hence there is nothing to do here

void SatelliteBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl ) {
 } // end( SatelliteBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS OF SatelliteSolution ----------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteSolution::deserialize( const netCDF::NcGroup & group )
{
 throw(
  std::logic_error( "SatelliteSolution::deserialize: not implemented" ) );
 }

/*--------------------------------------------------------------------------*/

void SatelliteSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const SatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument(
   "SatelliteSolution::read: block is not a SatelliteBlock" ) );

 if( ! v_zeta.empty() )
  v_zeta[ 0 ] = SATB->get_zeta();
 }

/*--------------------------------------------------------------------------*/

void SatelliteSolution::write( Block * block )
{
 auto SATB = dynamic_cast< SatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument(
   "SatelliteSolution::write: block is not a SatelliteBlock" ) );

 if( ! v_zeta.empty() )
  SATB->set_zeta( v_zeta.begin() );
 }

/*--------------------------------------------------------------------------*/

void SatelliteSolution::serialize( netCDF::NcGroup & group ) const
{
 throw( std::logic_error( "SatelliteSolution::serialize: not implemented" ) );
 }

/*--------------------------------------------------------------------------*/

SatelliteSolution * SatelliteSolution::scale( double factor ) const
{
 auto * sol = clone();
 for( auto & val : sol->v_zeta )
  val *= factor;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void SatelliteSolution::sum( const Solution * solution , double multiplier )
{
 auto tsol = dynamic_cast< const SatelliteSolution * >( solution );
 if( ! tsol )
  throw( std::invalid_argument(
   "SatelliteSolution::sum: solution is not a SatelliteSolution" ) );

 if( tsol->v_zeta.size() != v_zeta.size() )
  throw( std::invalid_argument( "SatelliteSolution::sum: wrong size" ) );

 for( Index i = 0 ; i < v_zeta.size() ; ++i )
  v_zeta[ i ] += multiplier * tsol->v_zeta[ i ];
 }

/*--------------------------------------------------------------------------*/

SatelliteSolution * SatelliteSolution::clone( bool empty ) const
{
 auto * sol = new SatelliteSolution();

 if( empty )
  sol->v_zeta.assign( v_zeta.size() , 0 );
 else
  sol->v_zeta = v_zeta;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void SatelliteSolution::print( std::ostream & output ) const
{
 output << "SatelliteSolution:";
 for( auto val : v_zeta )
  output << " zeta = " << val;
 output << std::endl;
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File SatelliteBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
