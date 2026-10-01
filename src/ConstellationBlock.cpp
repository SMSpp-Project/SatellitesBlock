/*--------------------------------------------------------------------------*/
/*---------------------- File ConstellationBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ConstellationBlock class.
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

#include "ConstellationBlock.h"

#include <algorithm>

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

// register ConstellationBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ConstellationBlock );

// the symbol that whoever links the module asks for, so that the linker keeps
// the module, and with it the registration of all its classes in the factory

SMSpp_define_force_load( SatellitesBlock )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF ConstellationBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

 void ConstellationBlock::load( const std::string & input , char frmt )
{
 std::ifstream iFile( input );
 if( ! iFile.is_open() )
  throw( std::invalid_argument(
   "ConstellationBlock::load: cannot open file " + input ) );

 load( iFile , frmt );

 } // end( ConstellationBlock::load( std::string ) )

/*--------------------------------------------------------------------------*/
// after reading the instance, for each satellite a discretized set of
// candidate circular orbits (varying inclination, right ascension of the
// ascending node and mean anomaly) is built, together with the distances
// between their ground tracks and the targets, which are used to load()
// the SatelliteBlock of the satellite. The steps are:
//
// 1) among the circular orbits making an integer number of revolutions
//    within the time horizon, those whose altitude (by Kepler's third law)
//    is in the range [ 400 , 1400 ] km are selected;
//
// 2) the angular resolution Theta_min, the angle swept in one time step on
//    the lowest of these orbits (times FACTOR), gives the half-aperture
//    aHalf of the sensor and the threshold thetaValF = Theta_min / 3;
//
// 3) inclination in [ 0 , PI ], and RAAN and mean anomaly in [ 0 , 2 PI ],
//    are each discretized into numbOfDiscretize points, so that consecutive
//    points are at most Theta_min apart;
//
// 4) for each satellite (the satellites are split into three groups, each
//    using a different altitude) and each candidate orbit, the ground track
//    of the circular orbit on the rotating Earth is computed at all the
//    time stamps, together with its geodesic latitude and longitude
//    distances to every target; the orbits that never bring any target
//    within the threshold are discarded, so as to keep the number of
//    candidate orbits [C] of the SatelliteBlock small.

void ConstellationBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // read the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( input >> eatcomments >> horizon ) || ( horizon <= 0 ) )
  throw( std::invalid_argument(
   "ConstellationBlock::load: error reading the horizon" ) );
 horizon *= 3600.0; // from hours to seconds

 if( ! ( input >> eatcomments >> time_step ) || ( time_step <= 0 ) )
  throw( std::invalid_argument(
   "ConstellationBlock::load: error reading the time step" ) );

 if( ! ( input >> eatcomments >> targets ) || ( ! targets ) )
  throw( std::invalid_argument(
   "ConstellationBlock::load: error reading the number of targets" ) );

 std::vector< double > Latitude( targets );
 std::vector< double > Longitude( targets );
 periods.resize( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
  if( ! ( input >> eatcomments >> periods[ i ] >> Latitude[ i ] >>
         Longitude[ i ] ) ||
      ( periods[ i ] < 1 ) )
   throw( std::invalid_argument(
    "ConstellationBlock::load: error reading the targets" ) );
  if( std::abs( Latitude[ i ] ) > 90 )
   throw( std::invalid_argument(
    "ConstellationBlock::load: latitude out of [ -90 , 90 ]" ) );
  Latitude[ i ] *= PI / 180;
  Longitude[ i ] *= PI / 180;
  }

 if( ! ( input >> eatcomments >> satellites ) || ( ! satellites ) )
  throw( std::invalid_argument(
   "ConstellationBlock::load: error reading the number of satellites" ) );

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

 // three altitudes are needed, one for each group of satellites
 if( altitude.size() < 3 )
  throw( std::invalid_argument( "ConstellationBlock::load: the horizon "
                                "allows less than three orbit altitudes" ) );

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
 // the circular orbit at every altitude, used in the propagation below

 std::vector< double > t_p( altSet );
 std::vector< double > t_u( altSet );
 std::vector< double > t_GM( altSet );
 for( Index ii = 0 ; ii < altSet ; ++ii ) {
  t_p[ ii ] = sqrt( MU / ( pow( RAYON + altitude[ ii ] , 3 ) ) );
  t_u[ ii ] = sqrt( MU / ( RAYON + altitude[ ii ] ) );
  t_GM[ ii ] = sqrt( ( RAYON + altitude[ ii ] ) / MU );
  }

 // inclination is sampled uniformly in [ 0 , PI ], RAAN and mean anomaly
 // in [ 0 , 2 PI ], each with numbOfDiscretize points

 const Index numbOfDiscretize = ceil( PI / Theta_min );
 if( numbOfDiscretize < 2 )
  throw( std::invalid_argument( "ConstellationBlock::load: time step too "
                                "long w.r.t. the orbital period" ) );

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

 thetaValF = Theta_min / 3.0;
 const double altitudeFinal = altitude[ altSet - 1 ];

 const Index t = horizon / time_step;

 boost::multi_array< double , 3 > CoverageSatLat(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLong(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLat1(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLong1(
  boost::extents[ targets ][ t ][ nOrbits ] );

 v_Block.resize( satellites );

 // the satellites are split into three groups of (almost) the same size,
 // using the altitude with index 2, 1 and 0, respectively

 const Index ialt = satellites / 3;

 for( Index isat = 0 ; isat < satellites ; ++isat ) {
  Index ii = 1;
  if( isat < ialt )
   ii = 2;
  if( isat > 2 * ialt )
   ii = 0;

  // enumerate the candidate orbits (inclination jj, RAAN k, mean anomaly
  // l), with index index1, and propagate the ground track of each one

  Index index1 = 0;
  Index indexOrbit = 0;
  for( Index jj = 0 ; jj < incSet ; ++jj )
   for( Index k = 0 ; k < ascSet ; ++k )
    for( Index l = 0 ; l < anmSet ; ++l , ++index1 ) {
     bool observes = false;
     for( Index j = 0 ; j < t ; ++j ) {
      // latitude of the projection of the satellite on the Earth surface
      const double lat_Sat = asin(
       ( ( sin( inclination[ jj ] ) * ( altitude[ ii ] + RAYON ) *
           sin( meanAnomaly[ l ] ) ) *
         cos( t_p[ ii ] * ( j * time_step ) ) / ( altitude[ ii ] + RAYON ) ) +
       ( ( sin( inclination[ jj ] ) * t_u[ ii ] * cos( meanAnomaly[ l ] ) ) *
         sin( t_p[ ii ] * ( j * time_step ) ) * t_GM[ ii ] ) );

      // longitude of the projection of the satellite on the Earth surface
      double long_Sat = fmod(
       -( angle0 + ( WE * ( j * time_step ) ) ) +
        atan2(
         ( ( ( sin( nodeAscendant[ k ] ) * ( altitude[ ii ] + RAYON ) *
               cos( meanAnomaly[ l ] ) ) +
             ( cos( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
               ( altitude[ ii ] + RAYON ) * sin( meanAnomaly[ l ] ) ) ) *
           cos( t_p[ ii ] * ( j * time_step ) ) /
           ( altitude[ ii ] + RAYON ) ) +
          ( ( -( sin( nodeAscendant[ k ] ) * t_u[ ii ] *
                 sin( meanAnomaly[ l ] ) ) +
              ( cos( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
                t_u[ ii ] * cos( meanAnomaly[ l ] ) ) ) *
            sin( t_p[ ii ] * ( j * time_step ) ) * t_GM[ ii ] ) ,
         ( ( ( ( cos( nodeAscendant[ k ] ) * ( altitude[ ii ] + RAYON ) *
                 cos( meanAnomaly[ l ] ) ) -
               ( sin( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
                 ( altitude[ ii ] + RAYON ) * sin( meanAnomaly[ l ] ) ) ) *
             cos( t_p[ ii ] * ( j * time_step ) ) /
             ( altitude[ ii ] + RAYON ) ) +
           ( ( -( cos( nodeAscendant[ k ] ) * t_u[ ii ] *
                  sin( meanAnomaly[ l ] ) ) -
               ( sin( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
                 t_u[ ii ] * cos( meanAnomaly[ l ] ) ) ) *
             sin( ( t_p[ ii ] * ( j * time_step ) ) ) * t_GM[ ii ] ) ) ) ,
       ( 2 * PI ) );
      if( long_Sat <= 0 )
       long_Sat += 2 * PI;

      // the (geodesic) latitude and longitude distances to every target
      for( Index i = 0 ; i < targets ; ++i ) {
       CoverageSatLat1[ i ][ j ][ index1 ] = std::abs(
        2 * asin( 0.5 * sqrt( 1 - cos( Latitude[ i ] - lat_Sat ) ) ) );
       CoverageSatLong1[ i ][ j ][ index1 ] =
        std::abs(
         2 *
         asin( 0.5 * sqrt( cos( Latitude[ i ] ) * cos( Latitude[ i ] ) *
                           ( 1 - cos( Longitude[ i ] - long_Sat ) ) ) ) ) *
        cos( Latitude[ i ] );
       if( ( CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValF ) &&
           ( CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValF ) )
        observes = true;
       }
      }

     // the orbits that observe no target at any time stamp are discarded,
     // the others are compacted into CoverageSatLat / CoverageSatLong

     if( ! observes )
      continue;

     for( Index j = 0 ; j < t ; ++j )
      for( Index i = 0 ; i < targets ; ++i ) {
       CoverageSatLat[ i ][ j ][ indexOrbit ] =
        CoverageSatLat1[ i ][ j ][ index1 ];
       CoverageSatLong[ i ][ j ][ indexOrbit ] =
        CoverageSatLong1[ i ][ j ][ index1 ];
       }

     ++indexOrbit;
     }

  // the indexOrbit kept orbits are the candidate orbits [C] of the
  // SatelliteBlock of the satellite

  auto SB = new SatelliteBlock( this );
  SB->load( targets , time_step , horizon , altitudeFinal , thetaValF ,
            indexOrbit , aHalf , CoverageSatLat , CoverageSatLong ,
            periods );
  v_Block[ isat ] = SB;
  }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( ConstellationBlock::load( std::istream ) )

/*--------------------------------------------------------------------------*/

void ConstellationBlock::generate_abstract_variables( Configuration * stvv )
{
 for( auto blck : v_Block )
  blck->generate_abstract_variables();

 } // end( ConstellationBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ConstellationBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( AR ) // the constraints are there already
  return; // nothing to do

 // the constraints below use the Variable of the sub-Block
 for( auto blck : v_Block ) {
  blck->generate_abstract_variables();
  blck->generate_abstract_constraints();
  }

 auto SB = [ this ]( Index k ) {
  return( static_cast< SatelliteBlock * >( v_Block[ k ] ) );
  };

 // thetaM: sum_k thetaVar[ k ] <= 0.9 thetaValF sum_k zeta[ k ], i.e.,
 // the average threshold of the active satellites is at most 0.9 times
 // the reference threshold thetaValF

 thetaM.resize( 1 );

 LinearFunction::v_coeff_pair v_var12;
 for( Index k = 0 ; k < satellites ; ++k ) {
  v_var12.push_back( std::make_pair( SB( k )->i2p_theta() , 1.0 ) );
  v_var12.push_back( std::make_pair( SB( k )->i2p_z() , -thetaValF * 0.9 ) );
  }
 thetaM[ 0 ].set_function( new LinearFunction( std::move( v_var12 ) ) ,
                           eNoBlck );
 thetaM[ 0 ].set_rhs( 0.0 , eNoBlck );
 thetaM[ 0 ].set_lhs( -Inf< double >() , eNoBlck );

 add_static_constraint( thetaM , "thetaM" );

 // observation, constraints (1) of the class comments: target i is
 // observed at least once in each of its revisit windows j, made of the pp
 // time stamps [ j pp , ( j + 1 ) pp ). Rows are allocated up to the
 // largest number of periods, those with j >= periods[ i ] being no
 // revisit window of target i: they get the void constraint
 // sum_k 0 * zeta[ k ] == 0

 const Index maxPeriods =
  *std::max_element( periods.begin() , periods.end() );
 observation.resize( boost::extents[ targets ][ maxPeriods ] );
 const Index t = horizon / time_step; // the number of time stamps

 for( Index i = 0 ; i < targets ; ++i ) {
  const double pp = horizon / time_step / periods[ i ];
  for( Index j = 0 ; j < periods[ i ] ; ++j ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index k = 0 ; k < satellites ; ++k )
    for( Index tt = j * pp ; ( tt < ( j + 1 ) * pp ) && ( tt < t ) ; ++tt )
     v_var.push_back( std::make_pair( SB( k )->i2p_r( i , tt ) , 1.0 ) );

   observation[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var ) ) , eNoBlck );
   observation[ i ][ j ].set_rhs( Inf< double >() , eNoBlck );
   observation[ i ][ j ].set_lhs( 1.0 , eNoBlck );
   }

  for( Index j = periods[ i ] ; j < maxPeriods ; ++j ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index k = 0 ; k < satellites ; ++k )
    v_var.push_back( std::make_pair( SB( k )->i2p_z() , 0.0 ) );

   observation[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var ) ) , eNoBlck );
   observation[ i ][ j ].set_rhs( 0.0 , eNoBlck );
   observation[ i ][ j ].set_lhs( 0.0 , eNoBlck );
   }
  }

 add_static_constraint( observation , "observation" );

 // observation1: no two satellites observe target i at time stamp j, i.e.,
 // sum_k xi_k[ i ][ j ] <= 1

 observation1.resize( boost::extents[ targets ][ t ] );

 for( Index i = 0 ; i < targets ; ++i )
  for( Index j = 0 ; j < t ; ++j ) {
   LinearFunction::v_coeff_pair v_var1;
   for( Index k = 0 ; k < satellites ; ++k )
    v_var1.push_back( std::make_pair( SB( k )->i2p_r( i , j ) , 1.0 ) );

   observation1[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var1 ) ) , eNoBlck );
   observation1[ i ][ j ].set_rhs( 1.0 , eNoBlck );
   observation1[ i ][ j ].set_lhs( -Inf< double >() , eNoBlck );
   }

 add_static_constraint( observation1 , "observation1" );

 AR = true;

 } // end( ConstellationBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR CHECKING THE ConstellationBlock ---------------*/
/*--------------------------------------------------------------------------*/

bool ConstellationBlock::is_feasible( bool useabstract ,
                                      Configuration * fsbc )
{
 if( ! AR ) // the constraints are not there
  return( false );

 // the tolerance and the type of violation
 double tol = 1e-1;
 bool rel_viol = true;

 // extract, from "c", the parameters that determine feasibility: if it
 // succeeds, it sets the values of the parameters and returns true,
 // otherwise it returns false
 auto extract_parameters = [ &tol , &rel_viol ]( Configuration * c ) -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
   }
  if( auto tc =
       dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >(
        c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
   }
  return( false );
  };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return( RowConstraint::is_feasible( thetaM , tol , rel_viol ) &&
          RowConstraint::is_feasible( observation , tol , rel_viol ) &&
          RowConstraint::is_feasible( observation1 , tol , rel_viol ) );

 } // end( ConstellationBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR PRINTING & SAVING THE ConstellationBlock ----------*/
/*--------------------------------------------------------------------------*/

void ConstellationBlock::print( std::ostream & output , char vlvl ) const
{
 output << "ConstellationBlock: horizon " << horizon << " s, time step "
        << time_step << " s, " << targets << " targets, " << satellites
        << " satellites" << std::endl;
 for( Index k = 0 ; k < v_Block.size() ; ++k )
  output << "satellite " << k << ": " << get_numOrbits( k ) << " orbits"
         << std::endl;

 } // end( ConstellationBlock::print )

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void ConstellationBlock::guts_of_destructor( void )
{
 // clear() all the Constraint, so that they do not bother to un-register
 // themselves from the Variable of the sub-Block, deleted below
 Constraint::clear( thetaM );
 Constraint::clear( observation );
 Constraint::clear( observation1 );

 // explicitly reset all Constraint and Variable, so that a new abstract
 // representation is not added to the (no longer current) one
 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 for( auto blck : v_Block )
  delete blck;
 v_Block.clear();

 AR = false;

 } // end( ConstellationBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*-------------------- End File ConstellationBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
