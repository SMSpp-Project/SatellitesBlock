/*--------------------------------------------------------------------------*/
/*--------------------- File Satellite.cpp ------------------------*/
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
static const double RAYON = 6378136.3;
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

template< typename T >
static Index countdiff( T beg , T end , T cmp )
{
 Index ndiff = 0;
 for( ; beg != end ; )
  if( *(beg++) != *(cmp++) )
   ndiff++;

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector< T > & vec , c_Subset & nms ,
		      typename std::vector< T >::const_iterator cmp ,
		      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices

template< typename T >
static Index countdiff( std::vector< T > & vec , c_Subset & nms ,
			typename std::vector< T >::const_iterator cmp ,
			Index n_max )
{
 Index ndiff = 0;
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   ndiff++;
  }

 return( ndiff );
 }

/*--------------------------------------------------------------------------*/
// copys one vector to a given subset of another

template< typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SatelliteBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SatelliteBlock );

// register SatelliteSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( SatelliteSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF SatelliteBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( FNumber num_targets , FNumber time_step, FNumber horizon, FNumber altValues , FNumber thetaValues ,
                            FNumber indOrbit , FNumber aHalf ,
                            boost::multi_array< double , 3 > CoverageLat , boost::multi_array< double , 3 > CoverageLong )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

  n = num_targets;
  OrbitSet = indOrbit;

  thetaVal = thetaValues;
  
  dt = time_step;
  T = horizon;
  t = T / dt;
  alphaHalf = aHalf;

  CoverageSatLat.resize(boost::extents[num_targets][t][OrbitSet]);
  CoverageSatLong.resize(boost::extents[num_targets][t][OrbitSet]);

  int index1 = -1;
    for( Index jj = 0 ; jj < OrbitSet ; ++jj ) {
            for( Index j = 0 ; j < t ; ++j ) {
              for( Index i = 0 ; i < num_targets ; ++i ){
                CoverageSatLat[i][j][jj] = CoverageLat[i][j][jj];
                CoverageSatLong[i][j][jj] = CoverageLong[i][j][jj];
              }
            } 
         }

 }  // end( SatelliteBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( SatelliteBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_abstract_variables( Configuration *stvv )
{

  if( AR3 & HasVar )  // the variables are there already
    return;           // nothing to do

  zeta.resize( 1 );
  for( auto & var : zeta )
   var.set_type( ColVariable::kBinary );

  add_static_variable( zeta );

  activation.resize( OrbitSet );
  for( auto & var : activation )
   var.set_type( ColVariable::kBinary );

  add_static_variable( activation );

  xi.resize( boost::extents[ n ][ t ] );
  for( Index i = 0 ; i < n ; ++i )
    for( Index j = 0 ; j < t ; ++j )
      xi[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( xi );

  AR3 |= HasVar;

 }  // end( SatelliteBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_abstract_constraints( Configuration *stcc )
{

  if( AR2 & HasCnst )  // the constraints are there already
    return;           // nothing to do

  // generate the orbitSelection constraint, which select exctly one orbit
  // configuration for the satellite: sum_{j \in OrbitSet} activation[ j ] == 1
  // remember: activation[ j ] = 1 iff. the j-th configuration is selected 

  orbitSelection.resize( 1 );
  LinearFunction::v_coeff_pair orbit_var;
  for( Index j = 0 ; j < OrbitSet ; ++j ) {
    orbit_var.push_back( std::make_pair( &activation[ j ], 1.0 ));
    }
  LinearFunction* FunctAnm = new LinearFunction( std::move( orbit_var ));
  orbitSelection[ 0 ].set_rhs( 1.0 );
  orbitSelection[ 0 ].set_lhs( 1.0 ); 
  orbitSelection[ 0 ].set_function( FunctAnm );

  add_static_constraint( orbitSelection );

  // generate the activationSat_cnst, which active the current satellite if
  // there exists a target m is observed by this satellite at a given time step:
  // zeta \leq xi[ i ][ j ], for all targets i's and time steps j's
  // remember: zeta = 1 iff. the current satellite is active in the constellation
  // and xi[ i ][ j ] = 1 iff. satellite observe target i at time step j

  activationSat_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ] );
  for( Index i = 0 ; i < n ; ++i ){
    for( Index j = 0 ; j < t ; ++j ){
      LinearFunction::v_coeff_pair v_var;
      v_var.push_back( std::make_pair( &xi[ i ][ j ], -1.0 ));
      v_var.push_back( std::make_pair( &zeta[ 0 ] ,  1.0 ));
      LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
      activationSat_cnst[ i ][ j ].set_rhs( Inf< double >() );
      activationSat_cnst[ i ][ j ].set_lhs( 0.0 ); 
      activationSat_cnst[ i ][ j ].set_function( FunctSat );
    }
  }
  add_static_constraint( activationSat_cnst );
  
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
  // sum_{jj \in OrbitSet}  (activation[ jj ] * CoverageSatLong[ i ][ j ][ jj ] 
  // for all targets m's and time steps j's
  
  // MLAT and MLONG are two big-M parameters automatically computed 
  // such that their numerical values are the smallest to guarantee 
  // that constraint are valid (redundant when xi[ i ][ j ] = 0)

  obs2_cnst.resize(
    boost::multi_array_types::
    extent_gen()[ n ][ t ] );
  obs4_cnst.resize(
    boost::multi_array_types::
    extent_gen()[ n ][ t ] );

  double MLAT = PI;
  double MLONG = 2*PI;
  
  for( Index i = 0 ; i < n ; ++i )
  {
    for( Index j = 0 ; j < t ; ++j )
    {
      MLAT = 0.0;
      MLONG = 0.0;
      LinearFunction::v_coeff_pair v_obs1, v_obs2;
        for( Index jj = 0 ; jj < OrbitSet ; ++jj )
        {
          v_obs1.push_back( std::make_pair( &activation[ jj ], 
                                -CoverageSatLat[ i ][ j ][ jj ])); 
          v_obs2.push_back( std::make_pair( &activation[ jj ], 
                                -CoverageSatLong[ i ][ j ][ jj ])); 
          MLAT = std::max(MLAT, (CoverageSatLat[ i ][ j ][ jj ]-thetaVal));
          MLONG = std::max(MLONG, (CoverageSatLong[ i ][ j ][ jj ]-thetaVal));
        }

      v_obs1.push_back( std::make_pair( &xi[i][j], -MLAT )); 
      v_obs2.push_back( std::make_pair( &xi[i][j], -MLONG )); 
      LinearFunction* Funct1 = new LinearFunction( std::move( v_obs1 ));
      LinearFunction* Funct2 = new LinearFunction( std::move( v_obs2 ));

      obs2_cnst[ i ][ j ].set_rhs( Inf< double >() );
      obs2_cnst[ i ][ j ].set_lhs( -thetaVal-MLAT ); 
      obs2_cnst[ i ][ j ].set_function( Funct1 );

      obs4_cnst[ i ][ j ].set_rhs( Inf< double >() );
      obs4_cnst[ i ][ j ].set_lhs( -thetaVal-MLONG ); 
      obs4_cnst[ i ][ j ].set_function( Funct2 );
    }
  }

  add_static_constraint( obs2_cnst );
  add_static_constraint( obs4_cnst );

  AR2 |= HasCnst;
 }  // end( SatelliteBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::generate_objective( Configuration *objc )
{

 if( AR1 & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

  LinearFunction::v_coeff_pair p( 1 );

  p[ 0 ].first = &zeta[ 0 ];
  p[ 0 ].second = 1.0;
  
  c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
  set_objective( & c , eNoMod );

  AR1 |= HasObj;

 }  // end( SatelliteBlock::generate_objective )

/*--------------------------------------------------------------------------*/

 bool SatelliteBlock::is_feasible( bool useabstract , Configuration *fsbc )
{
 FNumber eps = 0;
 auto tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                        f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return( 0 );

 }  // end( SatelliteBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool SatelliteBlock::is_optimal( bool useabstract , Configuration *optc )
{
 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  if( auto toptc =
      dynamic_cast< SimpleConfiguration< std::pair< CNumber , FNumber > > * >(
								    optc ) ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
   }
  else {
   auto ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >( optc );

   if( ( ! ttoptc ) && f_BlockConfig &&
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
  }
 else
  if( f_BlockConfig ) {
   if( f_BlockConfig->f_is_optimal_Configuration )
    if( auto csbc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                              f_BlockConfig->f_is_optimal_Configuration ) )
     ceps = csbc->f_value;

   if( f_BlockConfig->f_is_feasible_Configuration )
    if( auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration ) )
     feps = fsbc->f_value;
   }

 return( 0 );

 }  //  end( SatelliteBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/


 Solution * SatelliteBlock::get_Solution( Configuration * solc , bool emptys )
{
   
 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;
 
 auto *sol = new SatelliteSolution();

 if( ! emptys )
  sol->read( this );
 
 return( sol );

 }  // end( SatelliteBlock::get_Solution )

 void SatelliteBlock::set_zeta( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR3 & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 Index i = rng.first;

 for( auto xi = zeta.begin() + i ;
       i < 1 ; ++i )
   (xi++)->set_value( *(fstrt++) );
}


/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE SatelliteBlock ---*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::print( std::ostream  & output , char vlvl ) const
{
 
 }  // end( SatelliteBlock::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteBlock::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints
 Constraint::clear( orbitSelection );   // static
 Constraint::clear( activationSat_cnst );   // static
 Constraint::clear( obs2_cnst );
 Constraint::clear( obs4_cnst );

 c.clear();  // clear the Objective

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 //reset_dynamic_constraints();
 //reset_dynamic_variables();
 reset_objective();

 }  // end( SatelliteBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void SatelliteBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since SatelliteBlock is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN SatelliteBlock WHEN THIS METHOD IS
  *   EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS ISSUED:
  *   NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some of the logic here.*/

  // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {

  auto lfo = static_cast< LinearFunction * const >( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Range in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble

  return;
  }
*/

 //throw( std::invalid_argument( "unsupported Modification to SatelliteBlock" ) );

 }  // end( SatelliteBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF SatelliteSolution ------------------*/
/*--------------------------------------------------------------------------*/

void SatelliteSolution::deserialize( const netCDF::NcGroup & group )
{}

void SatelliteSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const SatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a SatelliteBlock" ) );

 if( ! v_zeta.empty() ) {
  v_zeta.resize( 1 );

  SATB->get_zeta();
  }
}

void SatelliteSolution::write( Block * block ) 
{

 auto SATB = dynamic_cast<SatelliteBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a SatelliteBlock" ) );

 if( ! v_zeta.empty() ) {
  SATB->set_zeta( v_zeta.begin() );
  }

}

void SatelliteSolution::serialize( netCDF::NcGroup & group ) const
{}

SatelliteSolution * SatelliteSolution::scale( double factor ) const
{
  auto * sol = SatelliteSolution::clone( true );
  return( sol );
}

void SatelliteSolution::sum( const Solution * solution , double multiplier )
{}

SatelliteSolution * SatelliteSolution::clone( bool empty ) const
{
  auto * sol = new SatelliteSolution();

 if( empty ) {
  if( ! v_zeta.empty() )
   sol->v_zeta.resize( 1 );
  }
 else {
  sol->v_zeta = v_zeta;
  }
 
 return( sol );
}

/*--------------------------------------------------------------------------*/
/*------------------- End File SatelliteBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
