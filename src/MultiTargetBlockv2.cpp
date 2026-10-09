/*--------------------------------------------------------------------------*/
/*---------------------- File MultiTargetBlockv2.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MultiTargetBlockv2 and MultiTargetSolution classes.
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

#include "MultiTargetBlockv2.h"

#include <cmath>

#include <fstream>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

using Range = Block::Range;

using FNumber = MultiTargetBlockv2::FNumber;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

static const auto RAYON = 6378136.3; // mean Earth radius [m]
static const auto PI = 3.14159265;
static const auto MU = 3.986004418e14; // Earth gravity constant [m^3/s^2]
static const auto WE = 7.2921e-5;      // Earth rotation speed [rad/s]
static const auto FACTOR = 1.2;        // safety factor on Theta_min
static const auto angle0 = -1.3882860164509252; // Greenwich angle at t = 0

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register MultiTargetBlockv2 to the Block factory

SMSpp_insert_in_factory_cpp_1( MultiTargetBlockv2 );

// register MultiTargetSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( MultiTargetSolution );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF MultiTargetBlockv2 ----------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::load( const std::string & input , char frmt )
{
 std::ifstream iFile( input );
 if( ! iFile.is_open() )
  throw( std::invalid_argument(
   "MultiTargetBlockv2::load: cannot open file " + input ) );

 load( iFile , frmt );

 } // end( MultiTargetBlockv2::load( std::string ) )

/*--------------------------------------------------------------------------*/
// the candidate orbits are built as in MultiTargetBlock::load() [see there]

void MultiTargetBlockv2::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // read the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double horizon;
 if( ! ( input >> eatcomments >> horizon ) || ( horizon <= 0 ) )
  throw( std::invalid_argument(
   "MultiTargetBlockv2::load: error reading the horizon" ) );
 horizon *= 3600.0; // from hours to seconds

 double time_step;
 if( ! ( input >> eatcomments >> time_step ) || ( time_step <= 0 ) )
  throw( std::invalid_argument(
   "MultiTargetBlockv2::load: error reading the time step" ) );

 if( ! ( input >> eatcomments >> targets ) || ( ! targets ) )
  throw( std::invalid_argument(
   "MultiTargetBlockv2::load: error reading the number of targets" ) );

 std::vector< double > Latitude( targets );
 std::vector< double > Longitude( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
  if( ! ( input >> eatcomments >> Latitude[ i ] >> Longitude[ i ] ) )
   throw( std::invalid_argument(
    "MultiTargetBlockv2::load: error reading the targets" ) );
  if( std::abs( Latitude[ i ] ) > 90 )
   throw( std::invalid_argument(
    "MultiTargetBlockv2::load: latitude out of [ -90 , 90 ]" ) );
  Latitude[ i ] *= PI / 180;
  Longitude[ i ] *= PI / 180;
  }

 if( ! ( input >> eatcomments >> n ) || ( ! n ) )
  throw( std::invalid_argument(
   "MultiTargetBlockv2::load: error reading the number of satellites" ) );

 // the altitudes - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the circular orbits making j revolutions within the horizon, for every
 // integer j, have period horizon / j and altitude (by Kepler's third law)
 // cbrt( MU ( horizon / j )^2 / ( 4 pi^2 ) ) - RAYON: only those in the
 // range [ 400 , 1400 ] km are kept

 std::vector< double > altitude;
 std::vector< double > periodSat;
 for( Index i = 0 ; i < horizon / 3600.0 ; ++i ) {
  const double j1 = i + 1;
  const double altitudeSetVal =
   cbrt( ( MU * pow( horizon / j1 , 2.0 ) ) / ( 4.0 * pow( PI , 2.0 ) ) ) -
   RAYON;
  if( ( altitudeSetVal >= 400000.0 ) && ( altitudeSetVal <= 1400000.0 ) ) {
   altitude.push_back( altitudeSetVal );
   periodSat.push_back( horizon / j1 );
   }
  }

 if( altitude.empty() )
  throw(
   std::invalid_argument( "MultiTargetBlockv2::load: the horizon "
                          "allows no orbit altitude in [ 400 , 1400 ] km" ) );

 const Index altSet = altitude.size();

 // Theta_min is the angle swept by the satellite on the lowest orbit in one
 // time step, times FACTOR; aHalf is the corresponding half-aperture of the
 // sensor, seen from the satellite

 const double Theta_min =
  ( ( 2 * PI * time_step ) / ( 2 * periodSat[ altSet - 1 ] ) ) * FACTOR;
 double aHalf =
  atan( sin( Theta_min ) /
        ( ( RAYON + altitude[ altSet - 1 ] ) / RAYON - cos( Theta_min ) ) );

 // if the sensor cone with half-aperture aHalf does not intersect the
 // Earth at some altitude, aHalf is reduced to reach the Earth's limb

 for( Index ii = 0 ; ii < altSet ; ++ii )
  if( ( ( RAYON + altitude[ ii ] ) / RAYON ) * sin( aHalf ) > 1 ) {
   aHalf = asin( ( RAYON / ( RAYON + altitude[ ii ] ) ) );
   break;
   }

 // mean motion (t_p), orbital speed (t_u) and sqrt( a / MU ) (t_GM) of
 // the circular orbit at the (single) altitude used

 const Index ii = altSet - 1;
 const double t_p = sqrt( MU / ( pow( RAYON + altitude[ ii ] , 3 ) ) );
 const double t_u = sqrt( MU / ( RAYON + altitude[ ii ] ) );
 const double t_GM = sqrt( ( RAYON + altitude[ ii ] ) / MU );

 // inclination is sampled uniformly in [ 0 , PI ], RAAN and mean anomaly
 // in [ 0 , 2 PI ], each with numbOfDiscretize points

 const Index numbOfDiscretize = ceil( PI / Theta_min );
 if( numbOfDiscretize < 2 )
  throw(
   std::invalid_argument( "MultiTargetBlockv2::load: time step too long "
                          "w.r.t. the orbital period" ) );

 const Index incSet = numbOfDiscretize;
 const Index ascSet = numbOfDiscretize;
 const Index anmSet = numbOfDiscretize;
 const Index nOrbits = incSet * ascSet * anmSet;

 std::vector< double > inclination( numbOfDiscretize );
 std::vector< double > nodeAscendant( numbOfDiscretize );
 for( Index i = 0 ; i < numbOfDiscretize ; ++i ) {
  inclination[ i ] = i * ( PI / ( numbOfDiscretize - 1 ) );
  nodeAscendant[ i ] = i * ( 2 * PI / ( numbOfDiscretize - 1 ) );
  }
 const std::vector< double > & meanAnomaly = nodeAscendant;

 thetaValFinal = Theta_min;
 alphaHalf = aHalf;
 dt = time_step;
 T = horizon;
 t = T / dt;

 CoverageSatLat.resize( boost::extents[ targets ][ t ][ nOrbits ] );
 CoverageSatLong.resize( boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLat1(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLong1(
  boost::extents[ targets ][ t ][ nOrbits ] );

 // enumerate the candidate orbits (inclination jj, RAAN k, mean anomaly l),
 // with index index1, and propagate the ground track of each one

 Index index1 = 0;
 OrbitSet = 0;
 for( Index jj = 0 ; jj < incSet ; ++jj )
  for( Index k = 0 ; k < ascSet ; ++k )
   for( Index l = 0 ; l < anmSet ; ++l , ++index1 ) {
    bool observes = false;
    for( Index j = 0 ; j < t ; ++j ) {
     // latitude of the projection of the satellite on the Earth surface
     const double lat_Sat =
      asin( ( ( sin( inclination[ jj ] ) * ( altitude[ ii ] + RAYON ) *
                sin( meanAnomaly[ l ] ) ) *
              cos( t_p * ( j * time_step ) ) / ( altitude[ ii ] + RAYON ) ) +
            ( ( sin( inclination[ jj ] ) * t_u * cos( meanAnomaly[ l ] ) ) *
              sin( t_p * ( j * time_step ) ) * t_GM ) );

     // longitude of the projection of the satellite on the Earth surface
     double long_Sat = fmod(
      -( angle0 + ( WE * ( j * time_step ) ) ) +
       atan2(
        ( ( ( sin( nodeAscendant[ k ] ) * ( altitude[ ii ] + RAYON ) *
              cos( meanAnomaly[ l ] ) ) +
            ( cos( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
              ( altitude[ ii ] + RAYON ) * sin( meanAnomaly[ l ] ) ) ) *
          cos( t_p * ( j * time_step ) ) / ( altitude[ ii ] + RAYON ) ) +
         ( ( -( sin( nodeAscendant[ k ] ) * t_u * sin( meanAnomaly[ l ] ) ) +
             ( cos( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) * t_u *
               cos( meanAnomaly[ l ] ) ) ) *
           sin( t_p * ( j * time_step ) ) * t_GM ) ,
        ( ( ( ( cos( nodeAscendant[ k ] ) * ( altitude[ ii ] + RAYON ) *
                cos( meanAnomaly[ l ] ) ) -
              ( sin( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
                ( altitude[ ii ] + RAYON ) * sin( meanAnomaly[ l ] ) ) ) *
            cos( t_p * ( j * time_step ) ) / ( altitude[ ii ] + RAYON ) ) +
          ( ( -( cos( nodeAscendant[ k ] ) * t_u * sin( meanAnomaly[ l ] ) ) -
              ( sin( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) * t_u *
                cos( meanAnomaly[ l ] ) ) ) *
            sin( ( t_p * ( j * time_step ) ) ) * t_GM ) ) ) ,
      ( 2 * PI ) );
     if( long_Sat <= 0 )
      long_Sat += 2 * PI;

     // the latitude and longitude distances to every target
     for( Index i = 0 ; i < targets ; ++i ) {
      CoverageSatLat1[ i ][ j ][ index1 ] =
       std::abs( Latitude[ i ] - lat_Sat );
      CoverageSatLong1[ i ][ j ][ index1 ] =
       std::abs( Longitude[ i ] - long_Sat ) * cos( Latitude[ i ] );
      if( ( CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValFinal ) &&
          ( CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValFinal ) )
       observes = true;
      }
     }

    // the orbits that observe no target at any time stamp are discarded,
    // the others are compacted into CoverageSatLat / CoverageSatLong

    if( ! observes )
     continue;

    for( Index j = 0 ; j < t ; ++j )
     for( Index i = 0 ; i < targets ; ++i ) {
      CoverageSatLat[ i ][ j ][ OrbitSet ] =
       CoverageSatLat1[ i ][ j ][ index1 ];
      CoverageSatLong[ i ][ j ][ OrbitSet ] =
       CoverageSatLong1[ i ][ j ][ index1 ];
      }

    ++OrbitSet;
    }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( MultiTargetBlockv2::load( std::istream ) )

/*--------------------------------------------------------------------------*/
// the Variable of SingleTargetBlock.h, each with a further index, the last
// one, for the target, save activation[][] and theta[], of which there is
// a single copy for all the targets

void MultiTargetBlockv2::generate_abstract_variables( Configuration * stvv )
{
 if( AR3 & HasVar ) // the variables are there already
  return;           // nothing to do

 theta.resize( n );
 for( auto & var : theta )
  var.set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( theta );

 Deltat.resize( targets );
 for( auto & var : Deltat )
  var.set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( Deltat );

 Deltat_k1A.resize( boost::extents[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index i = 0 ; i < t * ( t - 1 ) / 2 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   Deltat_k1A[ i ][ j ].set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( Deltat_k1A );

 Deltat_k1.resize( boost::extents[ t - 1 ][ targets ] );
 for( Index i = 0 ; i < t - 1 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   Deltat_k1[ i ][ j ].set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( Deltat_k1 );

 b1.resize( boost::extents[ t - 1 ][ targets ] );
 for( Index i = 0 ; i < t - 1 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   b1[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( b1 );

 Deltat_k2A.resize( boost::extents[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index i = 0 ; i < t * ( t - 1 ) / 2 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   Deltat_k2A[ i ][ j ].set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( Deltat_k2A );

 Deltat_k2.resize( boost::extents[ t - 1 ][ targets ] );
 for( Index i = 0 ; i < t - 1 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   Deltat_k2[ i ][ j ].set_type( ColVariable::kNonNegative , eNoBlck );

 add_static_variable( Deltat_k2 );

 b2.resize( boost::extents[ t - 1 ][ targets ] );
 for( Index i = 0 ; i < t - 1 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   b2[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( b2 );

 zeta.resize( boost::extents[ t ][ targets ] );
 for( Index i = 0 ; i < t ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   zeta[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( zeta );

 activation.resize( boost::extents[ n ][ OrbitSet ] );
 for( Index i = 0 ; i < n ; ++i )
  for( Index j = 0 ; j < OrbitSet ; ++j )
   activation[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( activation );

 h.resize( boost::extents[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index i = 0 ; i < t * ( t - 1 ) / 2 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   h[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( h );

 xi.resize( boost::extents[ n ][ t ][ targets ] );
 for( Index i = 0 ; i < n ; ++i )
  for( Index j = 0 ; j < t ; ++j )
   for( Index k = 0 ; k < targets ; ++k )
    xi[ i ][ j ][ k ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( xi );

 d1.resize( boost::extents[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index i = 0 ; i < t * ( t - 1 ) / 2 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   d1[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( d1 );

 d2.resize( boost::extents[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index i = 0 ; i < t * ( t - 1 ) / 2 ; ++i )
  for( Index j = 0 ; j < targets ; ++j )
   d2[ i ][ j ].set_type( ColVariable::kBinary , eNoBlck );

 add_static_variable( d2 );

 AR3 |= HasVar;

 } // end( MultiTargetBlockv2::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

// the constraints of SingleTargetBlock::generate_abstract_constraints()
// [see there], one copy per target k, save those involving only
// activation[][] and theta[], of which there is a single copy

void MultiTargetBlockv2::generate_abstract_constraints( Configuration * stcc )
{
 if( AR2 & HasCnst ) // the constraints are there already
  return;            // nothing to do
 // orbitSelection[ i ]: sum_j activation[ i ][ j ] == 1 for satellite i,
 // a single copy for all the targets

 orbitSelection.resize( n );
 for( Index i = 0 ; i < n ; ++i ) {
  LinearFunction::v_coeff_pair orbit_var;
  for( Index j = 0 ; j < OrbitSet ; ++j )
   orbit_var.push_back( std::make_pair( &activation[ i ][ j ] , 1.0 ) );
  LinearFunction * FunctAnm = new LinearFunction( std::move( orbit_var ) );
  orbitSelection[ i ].set_rhs( 1.0 , eNoBlck );
  orbitSelection[ i ].set_lhs( 1.0 , eNoBlck );
  orbitSelection[ i ].set_function( FunctAnm , eNoBlck );
  }

 add_static_constraint( orbitSelection , "orbitSelection" );

 activationSat_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < n ; ++i ) {
   for( Index j = 0 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &xi[ i ][ j ][ k ] , -1.0 ) );
    v_var.push_back( std::make_pair( &zeta[ j ][ k ] , 1.0 ) );
    LinearFunction * FunctSat = new LinearFunction( std::move( v_var ) );
    activationSat_cnst[ i ][ j ][ k ].set_rhs( Inf< double >() , eNoBlck );
    activationSat_cnst[ i ][ j ][ k ].set_lhs( 0.0 , eNoBlck );
    activationSat_cnst[ i ][ j ][ k ].set_function( FunctSat , eNoBlck );
    }
   }
  }
 add_static_constraint( activationSat_cnst , "activationSat_cnst" );

 activationSat1_cnst.resize(
  boost::multi_array_types::extent_gen()[ t ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < t ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index j = 0 ; j < n ; ++j ) {
    v_var.push_back( std::make_pair( &xi[ j ][ i ][ k ] , 1.0 ) );
    }
   v_var.push_back( std::make_pair( &zeta[ i ][ k ] , -1.0 ) );
   LinearFunction * FunctSat = new LinearFunction( std::move( v_var ) );
   activationSat1_cnst[ i ][ k ].set_rhs( Inf< double >() , eNoBlck );
   activationSat1_cnst[ i ][ k ].set_lhs( 0.0 , eNoBlck );
   activationSat1_cnst[ i ][ k ].set_function( FunctSat , eNoBlck );
   }
  }
 add_static_constraint( activationSat1_cnst , "activationSat1_cnst" );

 h_cnst_1.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &h[ ii ][ k ] , -1.0 ) );
    v_var.push_back( std::make_pair( &zeta[ i ][ k ] , 1.0 ) );
    LinearFunction * FunctSat = new LinearFunction( std::move( v_var ) );
    h_cnst_1[ ii ][ k ].set_rhs( Inf< double >() , eNoBlck );
    h_cnst_1[ ii ][ k ].set_lhs( 0.0 , eNoBlck );
    h_cnst_1[ ii ][ k ].set_function( FunctSat , eNoBlck );
    ii += 1;
    }
   }
  }
 add_static_constraint( h_cnst_1 , "h_cnst_1" );

 h_cnst_2.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &h[ ii ][ k ] , -1.0 ) );
    v_var.push_back( std::make_pair( &zeta[ j ][ k ] , 1.0 ) );
    LinearFunction * FunctSat = new LinearFunction( std::move( v_var ) );
    h_cnst_2[ ii ][ k ].set_rhs( Inf< double >() , eNoBlck );
    h_cnst_2[ ii ][ k ].set_lhs( 0.0 , eNoBlck );
    h_cnst_2[ ii ][ k ].set_function( FunctSat , eNoBlck );
    ii += 1;
    }
   }
  }
 add_static_constraint( h_cnst_2 , "h_cnst_2" );

 h_cnst_3.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &h[ ii ][ k ] , 1.0 ) );
    v_var.push_back( std::make_pair( &zeta[ i ][ k ] , -1.0 ) );
    v_var.push_back( std::make_pair( &zeta[ j ][ k ] , -1.0 ) );
    LinearFunction * FunctSat = new LinearFunction( std::move( v_var ) );
    h_cnst_3[ ii ][ k ].set_rhs( Inf< double >() , eNoBlck );
    h_cnst_3[ ii ][ k ].set_lhs( -1.0 , eNoBlck );
    h_cnst_3[ ii ][ k ].set_function( FunctSat , eNoBlck );
    ii += 1;
    }
   }
  }
 add_static_constraint( h_cnst_3 , "h_cnst_3" );

 // obs2_cnst and obs4_cnst, the big-M linearization of constraint (5) of
 // SingleTargetBlock.h for every target: target k can be observed by
 // satellite i at time stamp j only if the latitude and longitude
 // distances of the chosen orbit are within the threshold theta[ i ]

 obs2_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );
 obs4_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );

 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < n ; ++i ) {
   for( Index j = 0 ; j < t ; ++j ) {
    // the largest distance over all the orbits, i.e., the smallest big-M
    // making the constraint redundant when xi == 0, whatever the orbit
    double MLAT = 0.0;
    double MLONG = 0.0;
    LinearFunction::v_coeff_pair v_obs1 , v_obs2;
    for( Index jj = 0 ; jj < OrbitSet ; ++jj ) {
     v_obs1.push_back( std::make_pair( &activation[ i ][ jj ] ,
                                       -CoverageSatLat[ k ][ j ][ jj ] ) );
     v_obs2.push_back( std::make_pair( &activation[ i ][ jj ] ,
                                       -CoverageSatLong[ k ][ j ][ jj ] ) );
     MLAT = std::max( MLAT , ( CoverageSatLat[ k ][ j ][ jj ] ) );
     MLONG = std::max( MLONG , ( CoverageSatLong[ k ][ j ][ jj ] ) );
     }

    v_obs1.push_back( std::make_pair( &xi[ i ][ j ][ k ] , -MLAT ) );
    v_obs2.push_back( std::make_pair( &xi[ i ][ j ][ k ] , -MLONG ) );

    v_obs1.push_back( std::make_pair( &theta[ i ] , 1.0 ) );
    v_obs2.push_back( std::make_pair( &theta[ i ] , 1.0 ) );

    LinearFunction * Funct1 = new LinearFunction( std::move( v_obs1 ) );
    LinearFunction * Funct2 = new LinearFunction( std::move( v_obs2 ) );

    obs2_cnst[ i ][ j ][ k ].set_rhs( Inf< double >() , eNoBlck );
    obs2_cnst[ i ][ j ][ k ].set_lhs( ( -MLAT ) , eNoBlck );
    obs2_cnst[ i ][ j ][ k ].set_function( Funct1 , eNoBlck );

    obs4_cnst[ i ][ j ][ k ].set_rhs( Inf< double >() , eNoBlck );
    obs4_cnst[ i ][ j ][ k ].set_lhs( ( -MLONG ) , eNoBlck );
    obs4_cnst[ i ][ j ][ k ].set_function( Funct2 , eNoBlck );
    }
   }
  }

 add_static_constraint( obs2_cnst , "obs2_cnst" );
 add_static_constraint( obs4_cnst , "obs4_cnst" );

 theta_UB.resize( n );
 for( Index i = 0 ; i < n ; ++i ) {
  LinearFunction::v_coeff_pair v_var_theta;
  v_var_theta.push_back( std::make_pair( &theta[ i ] , 1.0 ) );
  LinearFunction * Funct_theta =
   new LinearFunction( std::move( v_var_theta ) );
  theta_UB[ i ].set_rhs( thetaValFinal , eNoBlck );
  theta_UB[ i ].set_lhs( -Inf< double >() , eNoBlck );
  theta_UB[ i ].set_function( Funct_theta , eNoBlck );
  }

 add_static_constraint( theta_UB , "theta_UB" );

 Deltat_max_dt1.resize( targets );
 for( Index k = 0 ; k < targets ; ++k ) {
  LinearFunction::v_coeff_pair v_vart1;
  v_vart1.push_back( std::make_pair( &Deltat[ k ] , 1.0 ) );
  LinearFunction * Functt1 = new LinearFunction( std::move( v_vart1 ) );
  Deltat_max_dt1[ k ].set_rhs( Inf< double >() , eNoBlck );
  Deltat_max_dt1[ k ].set_lhs( dt , eNoBlck );
  Deltat_max_dt1[ k ].set_function( Functt1 , eNoBlck );
  }

 add_static_constraint( Deltat_max_dt1 , "Deltat_max_dt1" );

 Deltat_max1.resize(
  boost::multi_array_types::extent_gen()[ t - 1 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   v_var.push_back( std::make_pair( &Deltat[ k ] , 1.0 ) );
   v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ] , -1.0 ) );
   v_var.push_back( std::make_pair( &b1[ i ][ k ] , 0.5 * T ) );
   LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
   Deltat_max1[ i ][ k ].set_rhs( Inf< double >() , eNoBlck );
   Deltat_max1[ i ][ k ].set_lhs( 0.0 , eNoBlck );
   Deltat_max1[ i ][ k ].set_function( Funct , eNoBlck );
   }
  }

 add_static_constraint( Deltat_max1 , "Deltat_max1" );

 Deltat_max11.resize(
  boost::multi_array_types::extent_gen()[ t - 1 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ] , 1.0 ) );
   v_var.push_back( std::make_pair( &b1[ i ][ k ] , -0.5 * T ) );
   LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
   Deltat_max11[ i ][ k ].set_rhs( Inf< double >() , eNoBlck );
   Deltat_max11[ i ][ k ].set_lhs( 0.0 , eNoBlck );
   Deltat_max11[ i ][ k ].set_function( Funct , eNoBlck );
   }
  }

 add_static_constraint( Deltat_max11 , "Deltat_max11" );

 Deltat_max2.resize(
  boost::multi_array_types::extent_gen()[ t - 1 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   v_var.push_back( std::make_pair( &Deltat[ k ] , 1.0 ) );
   v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ] , -1.0 ) );
   v_var.push_back( std::make_pair( &b2[ i ][ k ] , 0.5 * T ) );
   LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
   Deltat_max2[ i ][ k ].set_rhs( Inf< double >() , eNoBlck );
   Deltat_max2[ i ][ k ].set_lhs( 0.0 , eNoBlck );
   Deltat_max2[ i ][ k ].set_function( Funct , eNoBlck );
   }
  }

 add_static_constraint( Deltat_max2 , "Deltat_max2" );

 Deltat_max22.resize(
  boost::multi_array_types::extent_gen()[ t - 1 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ] , 1.0 ) );
   v_var.push_back( std::make_pair( &b2[ i ][ k ] , -0.5 * T ) );
   LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
   Deltat_max22[ i ][ k ].set_rhs( Inf< double >() , eNoBlck );
   Deltat_max22[ i ][ k ].set_lhs( 0.0 , eNoBlck );
   Deltat_max22[ i ][ k ].set_function( Funct , eNoBlck );
   }
  }

 add_static_constraint( Deltat_max22 , "Deltat_max22" );

 Deltat_min_k1_1.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ] , -1.0 ) );
    v_var.push_back(
     std::make_pair( &h[ ii ][ k ] , ( ( j - i ) * dt - 0.5 * T ) ) );
    LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
    Deltat_min_k1_1[ ii ][ k ].set_rhs( Inf< double >() , eNoBlck );
    Deltat_min_k1_1[ ii ][ k ].set_lhs( -0.5 * T , eNoBlck );
    Deltat_min_k1_1[ ii ][ k ].set_function( Funct , eNoBlck );
    ii += 1;
    }
   }
  }

 add_static_constraint( Deltat_min_k1_1 , "Deltat_min_k1_1" );

 Deltat_min_k1_2.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ] , -1.0 ) );
    v_var.push_back(
     std::make_pair( &h[ ii ][ k ] , ( ( j - i ) * dt - 0.5 * T ) ) );
    v_var.push_back( std::make_pair( &d1[ ii ][ k ] , 0.5 * T ) );
    LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
    Deltat_min_k1_2[ ii ][ k ].set_rhs( 0.0 , eNoBlck );
    Deltat_min_k1_2[ ii ][ k ].set_lhs( -Inf< double >() , eNoBlck );
    Deltat_min_k1_2[ ii ][ k ].set_function( Funct , eNoBlck );
    ii += 1;
    }
   }
  }

 add_static_constraint( Deltat_min_k1_2 , "Deltat_min_k1_2" );

 d1_cnst.resize( boost::multi_array_types::extent_gen()[ t - 1 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index j = i + 1 ; j < t ; ++j ) {
    v_var.push_back( std::make_pair( &d1[ ii ][ k ] , 1.0 ) );
    ii += 1;
    }
   LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
   d1_cnst[ i ][ k ].set_rhs( 1.0 , eNoBlck );
   d1_cnst[ i ][ k ].set_lhs( 1.0 , eNoBlck );
   d1_cnst[ i ][ k ].set_function( Funct , eNoBlck );
   }
  }

 add_static_constraint( d1_cnst , "d1_cnst" );

 Deltat_min_k2_1.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ] , 1.0 ) );
    v_var.push_back(
     std::make_pair( &h[ ii ][ k ] , ( ( j - i ) * dt - 0.5 * T ) ) );
    LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
    Deltat_min_k2_1[ ii ][ k ].set_rhs( 0.5 * T , eNoBlck );
    Deltat_min_k2_1[ ii ][ k ].set_lhs( -Inf< double >() , eNoBlck );
    Deltat_min_k2_1[ ii ][ k ].set_function( Funct , eNoBlck );
    ii += 1;
    }
   }
  }

 add_static_constraint( Deltat_min_k2_1 , "Deltat_min_k2_1" );

 Deltat_min_k2_2.resize(
  boost::multi_array_types::extent_gen()[ t * ( t - 1 ) / 2 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var;
    v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ] , 1.0 ) );
    v_var.push_back(
     std::make_pair( &h[ ii ][ k ] , ( ( j - i ) * dt - 0.5 * T ) ) );
    v_var.push_back( std::make_pair( &d2[ ii ][ k ] , -0.5 * T ) );
    LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
    Deltat_min_k2_2[ ii ][ k ].set_rhs( Inf< double >() , eNoBlck );
    Deltat_min_k2_2[ ii ][ k ].set_lhs( 0.0 , eNoBlck );
    Deltat_min_k2_2[ ii ][ k ].set_function( Funct , eNoBlck );
    ii += 1;
    }
   }
  }

 add_static_constraint( Deltat_min_k2_2 , "Deltat_min_k2_2" );

 d2_cnst.resize( boost::multi_array_types::extent_gen()[ t - 1 ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index j = i + 1 ; j < t ; ++j ) {
    v_var.push_back( std::make_pair( &d2[ ii ][ k ] , 1.0 ) );
    ii += 1;
    }
   LinearFunction * Funct = new LinearFunction( std::move( v_var ) );
   d2_cnst[ i ][ k ].set_rhs( 1.0 , eNoBlck );
   d2_cnst[ i ][ k ].set_lhs( 1.0 , eNoBlck );
   d2_cnst[ i ][ k ].set_function( Funct , eNoBlck );
   }
  }

 add_static_constraint( d2_cnst , "d2_cnst" );

 obs_cnst_h.resize( targets );
 for( Index k = 0 ; k < targets ; ++k ) {
  Index ii = 0;
  LinearFunction::v_coeff_pair v_var1;
  for( Index i = 0 ; i < t - 1 ; ++i ) {
   for( Index j = i + 1 ; j < t ; ++j ) {
    v_var1.push_back( std::make_pair( &h[ ii ][ k ] , 1.0 ) );
    ii += 1;
    }
   }
  LinearFunction * FunctA = new LinearFunction( std::move( v_var1 ) );
  obs_cnst_h[ k ].set_rhs( Inf< double >() , eNoBlck );
  obs_cnst_h[ k ].set_lhs( 1.0 , eNoBlck );
  obs_cnst_h[ k ].set_function( FunctA , eNoBlck );
  }

 add_static_constraint( obs_cnst_h , "obs_cnst_h" );

 obs_cnst_xi.resize( boost::multi_array_types::extent_gen()[ n ][ targets ] );
 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < n ; ++i ) {
   LinearFunction::v_coeff_pair v_varxi;
   for( Index j = 0 ; j < t ; ++j ) {
    v_varxi.push_back( std::make_pair( &xi[ i ][ j ][ k ] , 1.0 ) );
    }
   LinearFunction * FunctB = new LinearFunction( std::move( v_varxi ) );
   obs_cnst_xi[ i ][ k ].set_rhs( Inf< double >() , eNoBlck );
   obs_cnst_xi[ i ][ k ].set_lhs( 3.0 , eNoBlck );
   obs_cnst_xi[ i ][ k ].set_function( FunctB , eNoBlck );
   }
  }

 add_static_constraint( obs_cnst_xi , "obs_cnst_xi" );

 observation1.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()[ targets ][ t ] );

 for( Index i = 0 ; i < targets ; ++i ) {
  for( Index j = 0 ; j < t ; ++j ) {
   LinearFunction::v_coeff_pair v_var1;
   for( Index k = 0 ; k < n ; ++k )
    v_var1.push_back( std::make_pair( &xi[ k ][ j ][ i ] , 1.0 ) );

   observation1[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var1 ) ) , eNoBlck );
   observation1[ i ][ j ].set_rhs( 1.0 , eNoBlck );
   observation1[ i ][ j ].set_lhs( -Inf< double >() , eNoBlck );
   }
  }

 add_static_constraint( observation1 , "observation1" );

 AR2 |= HasCnst;
 } // end( MultiTargetBlockv2::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

// the objective is the average of the maximum revisit times Deltat[ k ] of
// the targets

void MultiTargetBlockv2::generate_objective( Configuration * objc )
{
 if( AR1 & HasObj ) // the objective is there already
  return;           // cowardly (and silently) return

 LinearFunction::v_coeff_pair v_obj;

 for( Index k = 0 ; k < targets ; ++k ) {
  v_obj.push_back( std::make_pair( &Deltat[ k ] , 1.0 / targets ) );
  }

 LinearFunction * Functobj = new LinearFunction( std::move( v_obj ) );

 c.set_function( Functobj , eNoBlck );
 set_objective( &c , eNoMod );

 AR1 |= HasObj;

 } // end( MultiTargetBlockv2::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR CHECKING THE Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool MultiTargetBlockv2::is_feasible( bool useabstract ,
                                      Configuration * fsbc )
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

 return( ColVariable::is_feasible( theta , eps ) &&
          ColVariable::is_feasible( Deltat , eps ) &&
          ColVariable::is_feasible( Deltat_k1 , eps ) &&
          ColVariable::is_feasible( Deltat_k2 , eps ) &&
          ColVariable::is_feasible( Deltat_k1A , eps ) &&
          ColVariable::is_feasible( Deltat_k2A , eps ) &&
          ColVariable::is_feasible( zeta , eps ) &&
          ColVariable::is_feasible( b1 , eps ) &&
          ColVariable::is_feasible( b2 , eps ) &&
          ColVariable::is_feasible( d1 , eps ) &&
          ColVariable::is_feasible( d2 , eps ) &&
          ColVariable::is_feasible( h , eps ) &&
          ColVariable::is_feasible( xi , eps ) &&
          ColVariable::is_feasible( activation , eps ) &&
          RowConstraint::is_feasible( orbitSelection , eps ) &&
          RowConstraint::is_feasible( theta_UB , eps ) &&
          RowConstraint::is_feasible( Deltat_max_dt1 , eps ) &&
          RowConstraint::is_feasible( Deltat_max1 , eps ) &&
          RowConstraint::is_feasible( Deltat_max11 , eps ) &&
          RowConstraint::is_feasible( Deltat_max2 , eps ) &&
          RowConstraint::is_feasible( Deltat_max22 , eps ) &&
          RowConstraint::is_feasible( Deltat_min_k1_1 , eps ) &&
          RowConstraint::is_feasible( Deltat_min_k1_2 , eps ) &&
          RowConstraint::is_feasible( d1_cnst , eps ) &&
          RowConstraint::is_feasible( Deltat_min_k2_1 , eps ) &&
          RowConstraint::is_feasible( Deltat_min_k2_2 , eps ) &&
          RowConstraint::is_feasible( d2_cnst , eps ) &&
          RowConstraint::is_feasible( activationSat_cnst , eps ) &&
          RowConstraint::is_feasible( activationSat1_cnst , eps ) &&
          RowConstraint::is_feasible( observation1 , eps ) &&
          RowConstraint::is_feasible( h_cnst_1 , eps ) &&
          RowConstraint::is_feasible( h_cnst_2 , eps ) &&
          RowConstraint::is_feasible( h_cnst_3 , eps ) &&
          RowConstraint::is_feasible( obs2_cnst , eps ) &&
          RowConstraint::is_feasible( obs4_cnst , eps ) &&
          RowConstraint::is_feasible( obs_cnst_h , eps ) &&
          RowConstraint::is_feasible( obs_cnst_xi , eps ) );

 } // end( MultiTargetBlockv2::is_feasible )

/*--------------------------------------------------------------------------*/

bool MultiTargetBlockv2::is_optimal( bool useabstract , Configuration * optc )
{
 return( false );

 } // end( MultiTargetBlockv2::is_optimal )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/

Solution * MultiTargetBlockv2::get_Solution( Configuration * solc ,
                                             bool emptys )
{
 auto * sol = new MultiTargetSolution();

 if( AR3 & HasVar ) {
  sol->v_Deltat.assign( Deltat.size() , 0 );
  if( ! emptys )
   sol->read( this );
  }

 return( sol );

 } // end( MultiTargetBlockv2::get_Solution )

/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::set_Deltat( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR3 & HasVar ) ) // nowhere to put the value in
  return;                // cowardly (and silently) return

 const Index stop = std::min( Index( rng.second ) , Index( Deltat.size() ) );
 for( Index k = rng.first ; k < stop ; ++k )
  Deltat[ k ].set_value( *( fstrt++ ) );

 } // end( MultiTargetBlockv2::set_Deltat )

/*--------------------------------------------------------------------------*/
/*------------------- Methods for handling Modification --------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 } // end( MultiTargetBlockv2::add_Modification )

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE MultiTargetBlockv2 ----------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::print( std::ostream & output , char vlvl ) const
{
 output << "MultiTargetBlockv2: " << targets << " targets, " << n
        << " satellites, " << t << " time stamps, " << OrbitSet << " orbits"
        << std::endl;

 } // end( MultiTargetBlockv2::print )

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::guts_of_destructor( void )
{
 // clear() the Constraint and the Objective, so that they do not bother to
 // un-register themselves from Variable that are going to be deleted anyway
 Constraint::clear( orbitSelection );
 Constraint::clear( theta_UB );
 Constraint::clear( Deltat_max_dt1 );
 Constraint::clear( Deltat_max1 );
 Constraint::clear( Deltat_max11 );
 Constraint::clear( Deltat_max2 );
 Constraint::clear( Deltat_max22 );
 Constraint::clear( Deltat_min_k1_1 );
 Constraint::clear( Deltat_min_k1_2 );
 Constraint::clear( d1_cnst );
 Constraint::clear( Deltat_min_k2_1 );
 Constraint::clear( Deltat_min_k2_2 );
 Constraint::clear( d2_cnst );
 Constraint::clear( activationSat_cnst );
 Constraint::clear( activationSat1_cnst );
 Constraint::clear( observation1 );
 Constraint::clear( h_cnst_1 );
 Constraint::clear( h_cnst_2 );
 Constraint::clear( h_cnst_3 );
 Constraint::clear( obs2_cnst );
 Constraint::clear( obs4_cnst );
 Constraint::clear( obs_cnst_h );
 Constraint::clear( obs_cnst_xi );
 c.clear();

 // explicitly reset all Constraint and Variable, so that a new abstract
 // representation is not added to the (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 reset_objective();

 theta.clear();
 Deltat.clear();
 Deltat_k1.resize( boost::extents[ 0 ][ 0 ] );
 Deltat_k2.resize( boost::extents[ 0 ][ 0 ] );
 Deltat_k1A.resize( boost::extents[ 0 ][ 0 ] );
 Deltat_k2A.resize( boost::extents[ 0 ][ 0 ] );
 zeta.resize( boost::extents[ 0 ][ 0 ] );
 b1.resize( boost::extents[ 0 ][ 0 ] );
 b2.resize( boost::extents[ 0 ][ 0 ] );
 d1.resize( boost::extents[ 0 ][ 0 ] );
 d2.resize( boost::extents[ 0 ][ 0 ] );
 h.resize( boost::extents[ 0 ][ 0 ] );
 xi.resize( boost::extents[ 0 ][ 0 ][ 0 ] );
 activation.resize( boost::extents[ 0 ][ 0 ] );

 AR1 = AR2 = AR3 = 0;

 } // end( MultiTargetBlockv2::guts_of_destructor )

/*--------------------------------------------------------------------------*/
// the MultiTargetBlockv2 has no physical representation that an abstract
// Modification may change, hence there is nothing to do here

void MultiTargetBlockv2::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{} // end( MultiTargetBlockv2::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF MultiTargetSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetSolution::deserialize( const netCDF::NcGroup & group )
{
 throw(
  std::logic_error( "MultiTargetSolution::deserialize: not implemented" ) );
 }

/*--------------------------------------------------------------------------*/

void MultiTargetSolution::read( const Block * block )
{
 auto MTB = dynamic_cast< const MultiTargetBlockv2 * >( block );
 if( ! MTB )
  throw( std::invalid_argument(
   "MultiTargetSolution::read: block is not a MultiTargetBlockv2" ) );

 if( v_Deltat.empty() )
  return;

 if( v_Deltat.size() != MTB->Deltat.size() )
  throw(
   std::invalid_argument( "MultiTargetSolution::read: wrong Block size" ) );

 for( Index k = 0 ; k < v_Deltat.size() ; ++k )
  v_Deltat[ k ] = MTB->get_Deltat( k );
 }

/*--------------------------------------------------------------------------*/

void MultiTargetSolution::write( Block * block )
{
 auto MTB = dynamic_cast< MultiTargetBlockv2 * >( block );
 if( ! MTB )
  throw( std::invalid_argument(
   "MultiTargetSolution::write: block is not a MultiTargetBlockv2" ) );

 if( ! v_Deltat.empty() )
  MTB->set_Deltat( v_Deltat.begin() , Range( 0 , v_Deltat.size() ) );
 }

/*--------------------------------------------------------------------------*/

void MultiTargetSolution::serialize( netCDF::NcGroup & group ) const
{
 throw(
  std::logic_error( "MultiTargetSolution::serialize: not implemented" ) );
 }

/*--------------------------------------------------------------------------*/

MultiTargetSolution * MultiTargetSolution::scale( double factor ) const
{
 auto * sol = clone();
 for( auto & val : sol->v_Deltat )
  val *= factor;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void MultiTargetSolution::sum( const Solution * solution , double multiplier )
{
 auto tsol = dynamic_cast< const MultiTargetSolution * >( solution );
 if( ! tsol )
  throw( std::invalid_argument(
   "MultiTargetSolution::sum: solution is not a MultiTargetSolution" ) );

 if( tsol->v_Deltat.size() != v_Deltat.size() )
  throw( std::invalid_argument( "MultiTargetSolution::sum: wrong size" ) );

 for( Index i = 0 ; i < v_Deltat.size() ; ++i )
  v_Deltat[ i ] += multiplier * tsol->v_Deltat[ i ];
 }

/*--------------------------------------------------------------------------*/

MultiTargetSolution * MultiTargetSolution::clone( bool empty ) const
{
 auto * sol = new MultiTargetSolution();

 if( empty )
  sol->v_Deltat.assign( v_Deltat.size() , 0 );
 else
  sol->v_Deltat = v_Deltat;

 return( sol );
 }

/*--------------------------------------------------------------------------*/

void MultiTargetSolution::print( std::ostream & output ) const
{
 output << "MultiTargetSolution: Deltat =";
 for( auto val : v_Deltat )
  output << " " << val;
 output << std::endl;
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File MultiTargetBlockv2.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
