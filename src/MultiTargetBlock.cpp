/*--------------------------------------------------------------------------*/
/*----------------------- File MultiTargetBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MultiTargetBlock class.
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

#include "MultiTargetBlock.h"

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

// register MultiTargetBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( MultiTargetBlock );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS OF MultiTargetBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlock::load( const std::string & input , char frmt )
{
 std::ifstream iFile( input );
 if( ! iFile.is_open() )
  throw( std::invalid_argument( "MultiTargetBlock::load: cannot open file " +
                                input ) );

 load( iFile , frmt );

 } // end( MultiTargetBlock::load( std::string ) )

/*--------------------------------------------------------------------------*/
// the candidate orbits are built with the same ground-track propagation as
// ConstellationBlock::load() [see there], with two differences: a single
// altitude, the lowest one, is used for all the satellites, since the set
// [C] of candidate orbits must be the same for all of them, and the
// latitude and longitude distances are the plain (not geodesic) ones

void MultiTargetBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // read the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( input >> eatcomments >> horizon ) || ( horizon <= 0 ) )
  throw( std::invalid_argument(
   "MultiTargetBlock::load: error reading the horizon" ) );
 horizon *= 3600.0; // from hours to seconds

 if( ! ( input >> eatcomments >> time_step ) || ( time_step <= 0 ) )
  throw( std::invalid_argument(
   "MultiTargetBlock::load: error reading the time step" ) );

 if( ! ( input >> eatcomments >> targets ) || ( ! targets ) )
  throw( std::invalid_argument(
   "MultiTargetBlock::load: error reading the number of targets" ) );

 std::vector< double > Latitude( targets );
 std::vector< double > Longitude( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
  if( ! ( input >> eatcomments >> Latitude[ i ] >> Longitude[ i ] ) )
   throw( std::invalid_argument(
    "MultiTargetBlock::load: error reading the targets" ) );
  if( std::abs( Latitude[ i ] ) > 90 )
   throw( std::invalid_argument(
    "MultiTargetBlock::load: latitude out of [ -90 , 90 ]" ) );
  Latitude[ i ] *= PI / 180;
  Longitude[ i ] *= PI / 180;
  }

 if( ! ( input >> eatcomments >> satellites ) || ( ! satellites ) )
  throw( std::invalid_argument(
   "MultiTargetBlock::load: error reading the number of satellites" ) );

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
   std::invalid_argument( "MultiTargetBlock::load: the horizon "
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
  throw( std::invalid_argument( "MultiTargetBlock::load: time step too long "
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

 const double thetaValFinal = Theta_min;

 const Index t = horizon / time_step;

 boost::multi_array< double , 3 > CoverageSatLat(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLong(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLat1(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLong1(
  boost::extents[ targets ][ t ][ nOrbits ] );

 // enumerate the candidate orbits (inclination jj, RAAN k, mean anomaly l),
 // with index index1, and propagate the ground track of each one

 Index index1 = 0;
 indexOrbit = 0;
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
      CoverageSatLat[ i ][ j ][ indexOrbit ] =
       CoverageSatLat1[ i ][ j ][ index1 ];
      CoverageSatLong[ i ][ j ][ indexOrbit ] =
       CoverageSatLong1[ i ][ j ][ index1 ];
      }

    ++indexOrbit;
    }

 // one SingleTargetBlock per target, with the (time stamp , orbit) slice
 // of the distances of that target; the set [C] of candidate orbits and
 // the number of satellites are the same for all, as required by the
 // "duplicate" constraints of generate_abstract_constraints()

 v_Block.resize( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
  boost::multi_array< double , 2 > CoverageSatLatThis(
   boost::extents[ t ][ indexOrbit ] );
  boost::multi_array< double , 2 > CoverageSatLongThis(
   boost::extents[ t ][ indexOrbit ] );

  for( Index j = 0 ; j < t ; ++j )
   for( Index k = 0 ; k < indexOrbit ; ++k ) {
    CoverageSatLatThis[ j ][ k ] = CoverageSatLat[ i ][ j ][ k ];
    CoverageSatLongThis[ j ][ k ] = CoverageSatLong[ i ][ j ][ k ];
    }

  auto SB = new SingleTargetBlock( this );
  SB->load( targets , satellites , time_step , horizon , thetaValFinal ,
            indexOrbit , aHalf , CoverageSatLatThis , CoverageSatLongThis );
  v_Block[ i ] = SB;
  }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( MultiTargetBlock::load( std::istream ) )

/*--------------------------------------------------------------------------*/

void MultiTargetBlock::generate_abstract_variables( Configuration * stvv )
{
 for( auto blck : v_Block )
  blck->generate_abstract_variables();

 } // end( MultiTargetBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void MultiTargetBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( AR ) // the constraints are there already
  return; // nothing to do

 // the constraints below use the Variable of the sub-Block
 for( auto blck : v_Block ) {
  blck->generate_abstract_variables();
  blck->generate_abstract_constraints();
  }

 auto SB = [ this ]( Index i ) {
  return( static_cast< SingleTargetBlock * >( v_Block[ i ] ) );
  };

 // duplicate_pi[ i ][ j ][ k ]: activation_i[ j ][ k ] ==
 // activation_{i+1}[ j ][ k ] for every satellite j and candidate orbit k,
 // linking target i to target i + 1

 duplicate_pi.resize(
  boost::extents[ targets - 1 ][ satellites ][ indexOrbit ] );
 for( Index i = 0 ; i + 1 < targets ; ++i )
  for( Index j = 0 ; j < satellites ; ++j )
   for( Index k = 0 ; k < indexOrbit ; ++k ) {
    LinearFunction::v_coeff_pair v_vars;
    v_vars.push_back( std::make_pair( SB( i )->i2p_pi( j , k ) , 1.0 ) );
    v_vars.push_back( std::make_pair( SB( i + 1 )->i2p_pi( j , k ) , -1.0 ) );
    duplicate_pi[ i ][ j ][ k ].set_function(
     new LinearFunction( std::move( v_vars ) ) , eNoBlck );
    duplicate_pi[ i ][ j ][ k ].set_rhs( 0.0 , eNoBlck );
    duplicate_pi[ i ][ j ][ k ].set_lhs( 0.0 , eNoBlck );
    }

 add_static_constraint( duplicate_pi , "duplicate_pi" );

 // duplicate_theta[ i ][ j ]: theta_i[ j ] == theta_{i+1}[ j ] for every
 // satellite j, linking target i to target i + 1

 duplicate_theta.resize( boost::extents[ targets - 1 ][ satellites ] );
 for( Index i = 0 ; i + 1 < targets ; ++i )
  for( Index j = 0 ; j < satellites ; ++j ) {
   LinearFunction::v_coeff_pair v_vars;
   v_vars.push_back( std::make_pair( SB( i )->i2p_theta( j ) , 1.0 ) );
   v_vars.push_back( std::make_pair( SB( i + 1 )->i2p_theta( j ) , -1.0 ) );
   duplicate_theta[ i ][ j ].set_function(
    new LinearFunction( std::move( v_vars ) ) , eNoBlck );
   duplicate_theta[ i ][ j ].set_rhs( 0.0 , eNoBlck );
   duplicate_theta[ i ][ j ].set_lhs( 0.0 , eNoBlck );
   }

 add_static_constraint( duplicate_theta , "duplicate_theta" );

 AR = true;

 } // end( MultiTargetBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
/*----------- METHODS FOR PRINTING & SAVING THE MultiTargetBlock -----------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlock::print( std::ostream & output , char vlvl ) const
{
 output << "MultiTargetBlock: horizon " << horizon << " s, time step "
        << time_step << " s, " << targets << " targets, " << satellites
        << " satellites, " << indexOrbit << " orbits" << std::endl;

 } // end( MultiTargetBlock::print )

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlock::guts_of_destructor( void )
{
 // clear() all the Constraint, so that they do not bother to un-register
 // themselves from the Variable of the sub-Block, deleted below
 Constraint::clear( duplicate_pi );
 Constraint::clear( duplicate_theta );

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

 } // end( MultiTargetBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*--------------------- End File MultiTargetBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
