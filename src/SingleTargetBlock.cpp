/*--------------------------------------------------------------------------*/
/*--------------------- File SingleTarget.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SingleTargetBlock class.
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

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SingleTargetBlock.h"
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

using FNumber = SingleTargetBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const double RAYON = 6378136.3; //[m]
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

// register SingleTargetBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SingleTargetBlock );

// register SingleTargetSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( SingleTargetSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF SingleTargetBlock -----------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SingleTargetBlock::load( FNumber num_satellites , FNumber time_step , 
                            FNumber horizon , FNumber altValues , FNumber thetaValues ,
                            FNumber indOrbit , FNumber aHalf ,
                            boost::multi_array< double , 2 > CoverageLat , 
                            boost::multi_array< double , 2 > CoverageLong )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

  n = num_satellites;
  OrbitSet = indOrbit;

  thetaVal = thetaValues;
  std::cout << "thetaVal=" << thetaValues << "\n";
  
  dt = time_step;
  T = horizon;
  t = T / dt;
  alphaHalf = aHalf;

  CoverageSatLat.resize(boost::extents[ t ][ OrbitSet ]);
  CoverageSatLong.resize(boost::extents[ t ][ OrbitSet ]);

  int index1 = -1;
  //for( Index ii = 0 ; ii < altSet ; ++ii )
    for( Index j = 0 ; j < t ; ++j ) {
      for( Index jj = 0 ; jj < OrbitSet ; ++jj ) {
        CoverageSatLat[ j ][ jj ] = CoverageLat[ j ][ jj ];
        CoverageSatLong[ j ][ jj ] = CoverageLong[j][ jj ];
      }
    } 

 }  // end( SingleTargetBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void SingleTargetBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( SingleTargetBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void SingleTargetBlock::generate_abstract_variables( Configuration *stvv )
{

  if( AR3 & HasVar )  // the variables are there already
    return;           // nothing to do

  Deltat.resize( 1 );
  for( auto & var : Deltat )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat );
  
  Deltat_k1A.resize( t * (t-1)/2 );
  for( auto & var : Deltat_k1A )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k1A );

  Deltat_k1.resize( t-1 );
  for( auto & var : Deltat_k1 )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k1 );

  b1.resize( t-1 );
  for( auto & var : b1 )
   var.set_type( ColVariable::kBinary );

  add_static_variable( b1 );

  Deltat_k2A.resize( t*(t-1)/2 );
  for( auto & var : Deltat_k2A )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k2A );

  Deltat_k2.resize( t-1 );
  for( auto & var : Deltat_k2 )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k2 );

  b2.resize( t-1 );
  for( auto & var : b2 )
   var.set_type( ColVariable::kBinary );

  add_static_variable( b2 );

  zeta.resize( t );
  for( auto & var : zeta )
   var.set_type( ColVariable::kBinary );

  add_static_variable( zeta );

  activation.resize( boost::extents[ n ][ OrbitSet ] );
  for( Index i = 0 ; i < n ; ++i )
    for( Index j = 0 ; j < OrbitSet ; ++j )
      activation[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( activation );

  h.resize( t * (t-1)/2 );
  for( auto & var : h )
   var.set_type( ColVariable::kBinary );

  add_static_variable( h );

  xi.resize( boost::extents[ n ][ t ] );
  for( Index i = 0 ; i < n ; ++i )
    for( Index j = 0 ; j < t ; ++j )
      xi[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( xi );

  d1.resize( t * (t-1)/2 );
  for( auto & var : d1 )
   var.set_type( ColVariable::kBinary );

  add_static_variable( d1 );

  d2.resize( t * (t-1)/2 );
  for( auto & var : d2 )
   var.set_type( ColVariable::kBinary );

  add_static_variable( d2 );
  
  AR3 |= HasVar;

 }  // end( SingleTargetBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SingleTargetBlock::generate_abstract_constraints( Configuration *stcc )
{

  if( AR2 & HasCnst )  // the constraints are there already
    return;           // nothing to do
   
  // generate the orbitSelection constraint, which select exctly one orbit
  // configuration for the satellite: sum_{j \in OrbitSet} activation[ i ][ j ] == 1
  // \forall i in [n] where [n] := \{1,2,...,n\} and n is the number of satellites  
  // remember: activation[ i ][ j ] = 1 iff. the j-th configuration is selected 
  // for the i-th satellite active in the constellation

  orbitSelection.resize( n );
  for( Index i = 0 ; i < n ; ++i ) {
    LinearFunction::v_coeff_pair orbit_var;
    for( Index j = 0 ; j < OrbitSet ; ++j ) 
      orbit_var.push_back( std::make_pair( &activation[ i ][ j ], 1.0 ));
    LinearFunction* FunctAnm = new LinearFunction( std::move( orbit_var ));
    orbitSelection[ i ].set_rhs( 1.0 );
    orbitSelection[ i ].set_lhs( 1.0 ); 
    orbitSelection[ i ].set_function( FunctAnm );
  }

  add_static_constraint( orbitSelection , "orbitSelection" );

  activationSat_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ] );
  for( Index i = 0 ; i < n ; ++i ){
    for( Index j = 0 ; j < t ; ++j ){
      LinearFunction::v_coeff_pair v_var;
      v_var.push_back( std::make_pair( &xi[ i ][ j ], -1.0 ));
      v_var.push_back( std::make_pair( &zeta[ j ] ,  1.0 ));
      LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
      activationSat_cnst[ i ][ j ].set_rhs( Inf< double >() );
      activationSat_cnst[ i ][ j ].set_lhs( 0.0 ); 
      activationSat_cnst[ i ][ j ].set_function( FunctSat );
    }
  }
  add_static_constraint( activationSat_cnst , "activationSat_cnst" );

  activationSat1_cnst.resize( t );
  for( Index i = 0 ; i < t ; ++i ){
    LinearFunction::v_coeff_pair v_var;
    for( Index j = 0 ; j < n ; ++j ){
      v_var.push_back( std::make_pair( &xi[ j ][ i ], 1.0 ));
    }
    v_var.push_back( std::make_pair( &zeta[ i ] ,  -1.0 ));
    LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
    activationSat1_cnst[ i ].set_rhs( Inf< double >() );
    activationSat1_cnst[ i ].set_lhs( 0.0 ); 
    activationSat1_cnst[ i ].set_function( FunctSat );
  }
  add_static_constraint( activationSat1_cnst , "activationSat1_cnst" );
  

  h_cnst_1.resize( t*(t-1)/2 );
  Index ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ){
    for( Index j = i+1 ; j < t ; ++j ){
      //if ( j > i ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &h[ ii ], -1.0 ));
        v_var.push_back( std::make_pair( &zeta[ i ] ,  1.0 ));
        LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
        h_cnst_1[ ii ].set_rhs( Inf< double >() );
        h_cnst_1[ ii ].set_lhs( 0.0 ); 
        h_cnst_1[ ii ].set_function( FunctSat );
        ii += 1;
    }
  }
  add_static_constraint( h_cnst_1 , "h_cnst_1" );


  h_cnst_2.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ){
    for( Index j = i+1 ; j < t ; ++j ){
      //if ( j > i ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &h[ ii ], -1.0 ));
        v_var.push_back( std::make_pair( &zeta[ j ] ,  1.0 ));
        LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
        h_cnst_2[ ii ].set_rhs( Inf< double >() );
        h_cnst_2[ ii ].set_lhs( 0.0 ); 
        h_cnst_2[ ii ].set_function( FunctSat );
        ii += 1;
    }
  }
  add_static_constraint( h_cnst_2 , "h_cnst_2" );

  h_cnst_3.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ){
    for( Index j = i+1 ; j < t ; ++j ){
      //if ( j > i ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &h[ ii ], 1.0 ));
        v_var.push_back( std::make_pair( &zeta[ i ] ,  -1.0 ));
        v_var.push_back( std::make_pair( &zeta[ j ] ,  -1.0 ));
        LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
        h_cnst_3[ ii ].set_rhs( Inf< double >() );
        h_cnst_3[ ii ].set_lhs( -1.0 ); 
        h_cnst_3[ ii ].set_function( FunctSat );
        ii += 1;
    }
  }
  add_static_constraint( h_cnst_3 , "h_cnst_3" );

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
          v_obs1.push_back( std::make_pair( &activation[ i ] [ jj ], 
                                -CoverageSatLat[ j ][ jj ])); 
          v_obs2.push_back( std::make_pair( &activation[ i ][ jj ], 
                                -CoverageSatLong[ j ][ jj ])); 
          MLAT = std::max(MLAT, (CoverageSatLat[ j ][ jj ]-thetaVal));
          MLONG = std::max(MLONG, (CoverageSatLong[ j ][ jj ]-thetaVal));
        }

      v_obs1.push_back( std::make_pair( &xi[ i ][ j ], -MLAT )); 
      v_obs2.push_back( std::make_pair( &xi[ i ][ j ], -MLONG )); 
      LinearFunction* Funct1 = new LinearFunction( std::move( v_obs1 ));
      LinearFunction* Funct2 = new LinearFunction( std::move( v_obs2 ));

      obs2_cnst[ i ][ j ].set_rhs( Inf< double >() );
      obs2_cnst[ i ][ j ].set_lhs(( -thetaVal-MLAT )); 
      obs2_cnst[ i ][ j ].set_function( Funct1 );

      obs4_cnst[ i ][ j ].set_rhs( Inf< double >() );
      obs4_cnst[ i ][ j ].set_lhs(( -thetaVal-MLONG )); 
      obs4_cnst[ i ][ j ].set_function( Funct2 );
    }
  }
  
  add_static_constraint( obs2_cnst );
  add_static_constraint( obs4_cnst );
     
  Deltat_max_dt.resize( 1 );
  LinearFunction::v_coeff_pair v_vart;
  v_vart.push_back( std::make_pair( &Deltat[ 0 ], 1.0 ));
  LinearFunction* Functt = new LinearFunction( std::move( v_vart ));
  Deltat_max_dt[ 0 ].set_rhs( T/2.0 - 1.0 );
  Deltat_max_dt[ 0 ].set_lhs( -Inf< double >() );
  Deltat_max_dt[ 0 ].set_function( Functt );

  //add_static_constraint( Deltat_max_dt , "Deltat_max_dt" );

  Deltat_max_dt1.resize( 1 );
  LinearFunction::v_coeff_pair v_vart1;
  v_vart1.push_back( std::make_pair( &Deltat[ 0 ], 1.0 ));
  LinearFunction* Functt1 = new LinearFunction( std::move( v_vart1 ));
  Deltat_max_dt1[ 0 ].set_rhs( Inf< double >() );
  Deltat_max_dt1[ 0 ].set_lhs( dt );
  Deltat_max_dt1[ 0 ].set_function( Functt1 );

  add_static_constraint( Deltat_max_dt1 , "Deltat_max_dt1" );

  Deltat_max1.resize( t-1 );
  for( Index i = 0 ; i < t-1 ; ++i ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat[ 0 ], 1.0 ));
    v_var.push_back( std::make_pair( &Deltat_k1[ i ], -1.0 ));
    v_var.push_back( std::make_pair( &b1[ i ], 0.5*T ));
    LinearFunction* Funct = new LinearFunction( std::move( v_var ));
    Deltat_max1[ i ].set_rhs( Inf< double >() );
    Deltat_max1[ i ].set_lhs( 0.0 ); 
    Deltat_max1[ i ].set_function( Funct );
  }

  add_static_constraint( Deltat_max1 , "Deltat_max1" );

  Deltat_max11.resize( t-1 );
  for( Index i = 0 ; i < t-1 ; ++i ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat_k1[ i ], 1.0 ));
    v_var.push_back( std::make_pair( &b1[ i ], -0.5*T ));
    LinearFunction* Funct = new LinearFunction( std::move( v_var ));
    Deltat_max11[ i ].set_rhs( Inf< double >() );
    Deltat_max11[ i ].set_lhs( 0.0 ); 
    Deltat_max11[ i ].set_function( Funct );
  }

  add_static_constraint( Deltat_max11 , "Deltat_max11" );

  Deltat_max2.resize( t-1 );
  for( Index i = 0 ; i < t-1 ; ++i ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat[ 0 ], 1.0 ));
    v_var.push_back( std::make_pair( &Deltat_k2[ i ], -1.0 ));
    v_var.push_back( std::make_pair( &b2[ i ], 0.5*T ));
    LinearFunction* Funct = new LinearFunction( std::move( v_var ));
    Deltat_max2[ i ].set_rhs( Inf< double >() );
    Deltat_max2[ i ].set_lhs( 0.0 ); 
    Deltat_max2[ i ].set_function( Funct );
  }

  add_static_constraint( Deltat_max2 , "Deltat_max2" );

  Deltat_max22.resize( t-1 ); 
  for( Index i = 0 ; i < t-1 ; ++i ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat_k2[ i ], 1.0 ));
    v_var.push_back( std::make_pair( &b2[ i ], -0.5*T ));
    LinearFunction* Funct = new LinearFunction( std::move( v_var ));
    Deltat_max22[ i ].set_rhs( Inf< double >() );
    Deltat_max22[ i ].set_lhs( 0.0 ); 
    Deltat_max22[ i ].set_function( Funct );
  }

  add_static_constraint( Deltat_max22 , "Deltat_max22" );


  Deltat_min_k1_11.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &Deltat_k1A[ ii ], -1.0 ));
        v_var.push_back( std::make_pair( &h[ ii ], (( j - i)*dt - T )));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        Deltat_min_k1_11[ ii ].set_rhs( -T );
        Deltat_min_k1_11[ ii ].set_lhs( -T ); 
        Deltat_min_k1_11[ ii ].set_function( Funct );
        ii += 1;
    }
  }

  //add_static_constraint( Deltat_min_k1_11 , "Deltat_min_k1_1" );

  Deltat_min_k2_11.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &Deltat_k2A[ ii ], 1.0 ));
        v_var.push_back( std::make_pair( &h[ ii ], (( j - i)*dt )));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        Deltat_min_k2_11[ ii ].set_rhs( T );
        Deltat_min_k2_11[ ii ].set_lhs( T ); 
        Deltat_min_k2_11[ ii ].set_function( Funct );
        ii += 1;
    }
  }

  //add_static_constraint( Deltat_min_k2_11 , "Deltat_min_k2_1" );

  Deltat_min_k1_1.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &Deltat_k1[ i ], -1.0 ));
        v_var.push_back( std::make_pair( &h[ ii ], (( j - i)*dt - 0.5*T )));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        Deltat_min_k1_1[ ii ].set_rhs( Inf< double >() );
        Deltat_min_k1_1[ ii ].set_lhs( -0.5*T ); 
        Deltat_min_k1_1[ ii ].set_function( Funct );
        ii += 1;
    }
  }

  add_static_constraint( Deltat_min_k1_1 , "Deltat_min_k1_1" );

  Deltat_min_k1_2.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &Deltat_k1[ i ], -1.0 ));
        v_var.push_back( std::make_pair( &h[ ii ], (( j - i)*dt - 0.5*T )));
        v_var.push_back( std::make_pair( &d1[ ii ], 0.5*T));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        Deltat_min_k1_2[ ii ].set_rhs( 0.0 );
        Deltat_min_k1_2[ ii ].set_lhs( -Inf< double >() ); 
        Deltat_min_k1_2[ ii ].set_function( Funct );
        ii += 1;
    }
  }

  add_static_constraint( Deltat_min_k1_2 , "Deltat_min_k1_2" );

  d1_cnst.resize( t-1 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    LinearFunction::v_coeff_pair v_var;
    for( Index j = i+1 ; j < t ; ++j ) {
	v_var.push_back( std::make_pair( &d1[ ii ], 1.0));
        ii += 1;
    }
    LinearFunction* Funct = new LinearFunction( std::move( v_var ));
    d1_cnst[ i ].set_rhs( 1.0 );
    d1_cnst[ i ].set_lhs( 1.0 ); 
    d1_cnst[ i ].set_function( Funct );
  }

  add_static_constraint( d1_cnst , "d1_cnst" );

  Deltat_min_k2_1.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &Deltat_k2[ i ], 1.0 ));
        v_var.push_back( std::make_pair( &h[ ii ], (( j - i)*dt - 0.5*T)));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        Deltat_min_k2_1[ ii ].set_rhs( 0.5*T );
        Deltat_min_k2_1[ ii ].set_lhs( -Inf< double >() ); 
        Deltat_min_k2_1[ ii ].set_function( Funct );
        ii += 1;
    }
  }

  add_static_constraint( Deltat_min_k2_1 , "Deltat_min_k2_1" );

  Deltat_min_k2_2.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &Deltat_k2[ i ], 1.0 ));
        v_var.push_back( std::make_pair( &h[ ii ], (( j - i)*dt - 0.5*T)));
        v_var.push_back( std::make_pair( &d2[ ii ], -0.5*T));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        Deltat_min_k2_2[ ii ].set_rhs( Inf< double >() );
        Deltat_min_k2_2[ ii ].set_lhs( 0.0 ); 
        Deltat_min_k2_2[ ii ].set_function( Funct );
        ii += 1;
    }
  }

  add_static_constraint( Deltat_min_k2_2 , "Deltat_min_k2_2" );

  d2_cnst.resize( t-1 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    LinearFunction::v_coeff_pair v_var;
    for( Index j = i+1 ; j < t ; ++j ) {
        v_var.push_back( std::make_pair( &d2[ ii ], 1.0));
        ii += 1;
    }
    LinearFunction* Funct = new LinearFunction( std::move( v_var ));
    d2_cnst[ i ].set_rhs( 1.0 );
    d2_cnst[ i ].set_lhs( 1.0 ); 
    d2_cnst[ i ].set_function( Funct );
  }

  add_static_constraint( d2_cnst , "d2_cnst" );

  d2a_cnst.resize( 1 );
  LinearFunction::v_coeff_pair v_varz;
  v_varz.push_back( std::make_pair( &zeta[ t-1 ], 1.0));
  LinearFunction* Functz = new LinearFunction( std::move( v_varz ));
  d2a_cnst[ 0 ].set_rhs( 0.0 );
  d2a_cnst[ 0 ].set_lhs( 0.0 ); 
  d2a_cnst[ 0 ].set_function( Functz );

  //add_static_constraint( d2a_cnst , "d2a_cnst" );
  
  obs_cnst_h.resize( 1 );
  ii = 0;
  LinearFunction::v_coeff_pair v_var1;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
      v_var1.push_back( std::make_pair( &h[ ii ], 1.0));
      ii += 1;
    }
  }
  LinearFunction* FunctA = new LinearFunction( std::move( v_var1 ));
  obs_cnst_h[ 0 ].set_rhs( Inf< double >() );
  obs_cnst_h[ 0 ].set_lhs( 1.0 );
  obs_cnst_h[ 0 ].set_function( FunctA );

  add_static_constraint( obs_cnst_h , "obs_cnst_h" );
  

  obs_cnst_xi.resize( n );
  for( Index i = 0 ; i < n ; ++i ) {
    LinearFunction::v_coeff_pair v_varxi;
    for( Index j = 0 ; j < t ; ++j ) {
      v_varxi.push_back( std::make_pair( &xi[ i ][ j ], 1.0));
    }
    LinearFunction* FunctB = new LinearFunction( std::move( v_varxi ));
    obs_cnst_xi[ i ].set_rhs( Inf< double >() );
    obs_cnst_xi[ i ].set_lhs( 3.0 );
    obs_cnst_xi[ i ].set_function( FunctB );
  }

  add_static_constraint( obs_cnst_xi , "obs_cnst_xi" );


  d1A_cnst.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &d1[ ii ], 1.0));
        v_var.push_back( std::make_pair( &h[ ii ], -1.0));
        LinearFunction* Funct = new LinearFunction( std::move( v_var ));
        d1A_cnst[ ii ].set_rhs( 0.0 );
        d1A_cnst[ ii ].set_lhs( -Inf< double >() );
        d1A_cnst[ ii ].set_function( Funct );
    }
  }

  //add_static_constraint( d1A_cnst , "d1A_cnst" );
  
  d2A_cnst.resize( t*(t-1)/2 );
  ii = 0;
  for( Index i = 0 ; i < t-1 ; ++i ) {
    for( Index j = i+1 ; j < t ; ++j ) {
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &d2[ ii ], 1.0));
        v_var.push_back( std::make_pair( &h[ ii ], -1.0));
	LinearFunction* Funct = new LinearFunction( std::move( v_var ));
	d2A_cnst[ ii ].set_rhs( 0.0 );
	d2A_cnst[ ii ].set_lhs( -Inf< double >() );
	d2A_cnst[ ii ].set_function( Funct );
    }
  }

  //add_static_constraint( d2A_cnst , "d2A_cnst" );

  AR2 |= HasCnst;
 }  // end( SingleTargetBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void SingleTargetBlock::generate_objective( Configuration *objc )
{

 if( AR1 & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return
 
  LinearFunction::v_coeff_pair p( 1 );

  p[ 0 ].first = &Deltat[ 0 ];
  p[ 0 ].second = 1.0;

  LinearFunction* Functobj = new LinearFunction( std::move( p ));
  
  c.set_function( Functobj );
  set_objective( &c );
  
  AR1 |= HasObj;

 }  // end( SingleTargetBlock::generate_objective )

/*--------------------------------------------------------------------------*/

 bool SingleTargetBlock::is_feasible( bool useabstract , Configuration *fsbc )
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

 }  // end( SingleTargetBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool SingleTargetBlock::is_optimal( bool useabstract , Configuration *optc )
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

 }  //  end( SingleTargetBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/


 Solution * SingleTargetBlock::get_Solution( Configuration * solc , bool emptys )
{
   
 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;
 
 auto *sol = new SingleTargetSolution();

 if( ! emptys )
  sol->read( this );
 
 return( sol );

 }  // end( SingleTargetBlock::get_Solution )

 void SingleTargetBlock::set_zeta( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR3 & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 Index i = rng.first;

 for( auto xi = Deltat.begin() + i ;
       i < 1 ; ++i )
   (xi++)->set_value( *(fstrt++) );
}


/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void SingleTargetBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE SingleTargetBlock ---*/
/*--------------------------------------------------------------------------*/

void SingleTargetBlock::print( std::ostream  & output , char vlvl ) const
{
  // TO DO: implement print() method for printing instance data to output file 
 }  // end( SingleTargetBlock::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SingleTargetBlock::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints
  
 Constraint::clear( orbitSelection ); 
 Constraint::clear( Deltat_max1 );
 Constraint::clear( Deltat_max2 );
 Constraint::clear( Deltat_max3 );
 Constraint::clear( Deltat_min_k1_1 );
 Constraint::clear( Deltat_min_k1_2 );
 Constraint::clear( Deltat_min_k2_1 );
 Constraint::clear( Deltat_min_k2_2 );
 Constraint::clear( d1_cnst );
 Constraint::clear( d2_cnst );
 //Constraint::clear( dA_cnst );
 //Constraint::clear( dB_cnst );
 Constraint::clear( h_cnst_1 );
 Constraint::clear( h_cnst_2 );
 Constraint::clear( h_cnst_3 );
 Constraint::clear( activationSat_cnst );
 //Constraint::clear( obs1_cnst );
 Constraint::clear( obs2_cnst );
 //Constraint::clear( obs3_cnst );
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

 }  // end( SingleTargetBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void SingleTargetBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since SingleTargetBlock is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN SingleTargetBlock WHEN THIS METHOD IS
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

 //throw( std::invalid_argument( "unsupported Modification to SingleTargetBlock" ) );

 }  // end( SingleTargetBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF DCRSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void SingleTargetSolution::deserialize( const netCDF::NcGroup & group )
{
  // TO DO: implement deserialize() method for dealign with netCDF input files
}

void SingleTargetSolution::serialize( netCDF::NcGroup & group ) const
{
  // TO DO: implement serialize() method for dealing with netCDF input files
}

void SingleTargetSolution::sum( const Solution * solution , double multiplier )
{
  // TO DO: implement sum() method
}

void SingleTargetSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const SingleTargetBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a SingleTargetBlock" ) );

 if( ! v_zeta.empty() ) {
  v_zeta.resize( 1 );

  SATB->get_zeta();
  }
}

void SingleTargetSolution::write( Block * block ) 
{
 auto SATB = dynamic_cast<SingleTargetBlock * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a SingleTargetBlock" ) );

 if( ! v_zeta.empty() ) {
  SATB->set_zeta( v_zeta.begin() );
  }
}

SingleTargetSolution * SingleTargetSolution::scale( double factor ) const
{
  auto * sol = SingleTargetSolution::clone( true );
  return( sol );
}

SingleTargetSolution * SingleTargetSolution::clone( bool empty ) const
{
  auto * sol = new SingleTargetSolution();

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
/*------------------- End File SingleTargetBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
