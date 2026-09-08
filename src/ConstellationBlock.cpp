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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ConstellationBlock.h"

#include <math.h>

#include <ctype.h>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC MEMBERS -------------------------------*/
/*--------------------------------------------------------------------------*/

// register ConstellationBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( ConstellationBlock );

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const auto RAYON = 6378136.3; //[m]
static const auto PI = 3.14159265;
static const auto MU = 3.986004418e14; //[m^3/s^2]
static const auto WE = 7.2921e-5;
static const auto FACTOR = 1.2;
static const auto angle0 = -1.3882860164509252;

/*--------------------------------------------------------------------------*/
/*------------------------ OTHER INITIALIZATIONS ---------------------------*/
/*--------------------------------------------------------------------------*/

void ConstellationBlock::load( std::istream & input , char frmt )
{
 // TO DO: implement load() method for loading instance data
 // from input file (for the file format, see next load() method)
}

/*--------------------------------------------------------------------------*/
// this load() reads the instance description (time horizon, targets with
// their revisit periods and coordinates, number of satellites) from file
// and then, for each satellite, builds a discretized set of candidate
// circular orbits (varying inclination, right ascension of the ascending
// node and mean anomaly) together with the corresponding target coverage
// distances, and uses them to construct and load() the satellite's
// SatelliteBlock (see SatelliteBlock::load()). The overall algorithm is:
//
// 1) determine, among all orbital periods compatible with an integer number
//    of orbits within the time horizon, those whose corresponding altitude
//    (via Kepler's third law) falls in the "usable" LEO band [400,1400] Km;
//
// 2) compute aHalf (half the Earth central angle subtended by the sensor
//    swath) and the corresponding observability threshold theta^{\max} for
//    each candidate altitude, out of the (worst-case) angular resolution
//    Theta_min imposed by the revisit requirement time_step;
//
// 3) discretize inclination in [0,PI], and RAAN / mean anomaly in
//    [0,2*PI], each into numbOfDiscretize samples, so that the Cartesian
//    product of the three angles (together with the altitude "band" ii,
//    see below) defines the candidate orbital configurations [C];
//
// 4) for each satellite (partitioned into three altitude bands so as to
//    diversify the constellation) and for each candidate configuration,
//    propagate the (simplified, J2-perturbed) ground track over all time
//    stamps and compute the geodesic latitude/longitude distance to every
//    target (CoverageSatLat1/CoverageSatLong1); configurations that never
//    bring any target within the threshold are discarded so as to keep the
//    number of columns [C] of the resulting SatelliteBlock manageable, and
//    the surviving ones are compacted into CoverageSatLat/CoverageSatLong.

void ConstellationBlock::load( const std::string & input , char frmt )
{
 /* The structure of the input file to load the instance data should be the
 * following (see, as example, input file in testConstellation directory):
 *
 * < time horizon (in hours) >
 * < time step for time discretization (in seconds) >
 * < number m of targets >
 * < number of revisit period of the target 1 > < latitude of the target 1 > < longitude of the target 1 >
 * < ... >
 * < number of revisit period of the target m > < latitude of the target m > < longitude of the target m >
 * < maximum number of satellites active in the constellation >
 *
 */

 // ensure starting from clean slate
 guts_of_destructor();

 std::ifstream iFile( input );

 iFile >> horizon;
 horizon *= 3600.0;  // from hours to seconds
 iFile >> time_step;
 iFile >> targets;

 std::cout << "TARGETS: " << targets << "\n";
 std::cout << "TIME STEP: " << time_step << "\n";

 Vec_CNumber Latitude;
 Vec_CNumber Longitude;

 double altitudeSetVal;
 Vec_CNumber altitude;
 Vec_CNumber periodSat;
 int indexLen = 0;
 double j = 0.0;

 // scan all the (integer) numbers j of orbits that fit within the horizon
 // and, via Kepler's third law, back out the corresponding altitude
 // altitudeSetVal = cbrt( MU * ( horizon / j )^2 / ( 4 pi^2 ) ) - RAYON;
 // only altitudes in the usable LEO range [400,1400] Km are kept

 for( Index i = 0 ; i < horizon / 3600.0 ; ++i ) {
  j++;
  altitudeSetVal =
   cbrt( ( MU * pow( ( horizon ) / j , 2.0 ) ) / ( 4.0 * pow( PI , 2.0 ) ) ) -
   RAYON;
  if( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 ) {
   indexLen++;
  }
 }

 // same loop as above, this time actually storing the altitude / orbital
 // period pairs that were only counted before

 FNumber altSet = indexLen;
 altitude.resize( indexLen );
 periodSat.resize( indexLen );
 indexLen = 0;

 j = 0.0;
 for( Index i = 0 ; i < horizon / 3600.0 ; ++i ) {
  j++;
  altitudeSetVal =
   cbrt( ( MU * pow( ( horizon ) / j , 2.0 ) ) / ( 4.0 * pow( PI , 2.0 ) ) ) -
   RAYON;
  if( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 ) {
   altitude[ indexLen ] = altitudeSetVal;
   periodSat[ indexLen ] = horizon / j;
   indexLen++;
  }
 }

 Latitude.resize( targets );
 Longitude.resize( targets );
 periods.resize( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
  iFile >> periods[ i ];
  iFile >> Latitude[ i ];
  iFile >> Longitude[ i ];
  Latitude[ i ] *= PI / 180;
  Longitude[ i ] *= PI / 180;
 }

 Vec_CNumber t_p;
 Vec_CNumber t_u;
 Vec_CNumber t_GM;
 Vec_CNumber thetaVal;

 t_p.resize( altSet );
 t_u.resize( altSet );
 t_GM.resize( altSet );
 thetaVal.resize( altSet );

 // Theta_min is the angular resolution (scaled by FACTOR as a safety
 // margin) that the constellation must achieve to guarantee the revisit
 // time_step at the largest computed orbital period; aHalf is the
 // corresponding half-cone angle of the sensor as seen from the satellite,
 // obtained by inverting the (planar) visibility geometry

 double alphalim;
 double Theta_min =
  ( ( 2 * PI * time_step ) / ( 2 * periodSat[ indexLen - 1 ] ) ) * FACTOR;
 std::cout << "theta_min: " << Theta_min << "\n";
 double aHalf =
  atan( sin( Theta_min ) /
        ( ( RAYON + altitude[ indexLen - 1 ] ) / RAYON - cos( Theta_min ) ) );

 // guard against aHalf being geometrically infeasible (asin argument > 1)
 // at some candidate altitude: if so, aHalf is capped to the largest value
 // for which the sensor cone still reaches the Earth's limb at that altitude

 for( Index ii = 0 ; ii < altSet ; ++ii ) {
  if( ( ( RAYON + altitude[ ii ] ) / RAYON ) * sin( aHalf ) > 1 ) {
   aHalf = asin( ( RAYON / ( RAYON + altitude[ ii ] ) ) );
   std::cout << "WARNING: computed alpha_lim\n";
   break;
  }
 }

 // for every candidate altitude, precompute the mean motion (t_p), orbital
 // velocity (t_u) and its inverse scaled by sqrt(a) (t_GM) used below in
 // the ground-track propagation, as well as the Earth central angle
 // theta^{\max} (thetaVal) subtended by the sensor swath at that altitude

 double alt = altitude[ altSet - 1 ];
 for( Index ii = 0 ; ii < altSet ; ++ii ) {
  t_p[ ii ] = sqrt( MU / ( pow( RAYON + altitude[ ii ] , 3 ) ) );
  t_u[ ii ] = sqrt( MU / ( RAYON + altitude[ ii ] ) );
  t_GM[ ii ] = sqrt( ( RAYON + altitude[ ii ] ) / MU );
  thetaVal[ ii ] =
   -aHalf + asin( ( ( RAYON + altitude[ ii ] ) / RAYON ) * sin( aHalf ) );
 }

 double Theta_max = thetaVal[ altSet - 1 ];
 std::cout << "theta_min: " << Theta_min << "\n";
 Theta_max = Theta_min;

 // the number of discretization points for each orbital angle (inclination,
 // RAAN, mean anomaly) is chosen so that a full turn (2*PI) is sampled at a
 // resolution no coarser than Theta_min, i.e., consecutive samples are at
 // most Theta_min apart

 double numbOfDiscretize = ceil( PI / Theta_min );

 std::cout << "numbOfDiscretize: " << numbOfDiscretize << "\n";

 FNumber incSet = numbOfDiscretize;
 FNumber ascSet = numbOfDiscretize;
 FNumber anmSet = numbOfDiscretize;

 // inclination is sampled uniformly in [0,PI] ...

 double start_in = 0;
 double end_in = PI;
 size_t num_in = numbOfDiscretize;

 double dx = ( end_in - start_in ) / ( num_in - 1 );
 std::vector< double > x( num_in );
 int iter = 0;
 std::generate( x.begin() , x.end() ,
                [ & ] { return start_in + ( iter++ ) * dx; } );

 Vec_CNumber inclination = x;

 // ... while RAAN and mean anomaly are (independently) sampled uniformly
 // in [0,2*PI], reusing the same number of discretization points

 end_in = 2 * PI;
 dx = ( end_in - start_in ) / ( num_in - 1 );
 iter = 0;
 std::generate( x.begin() , x.end() ,
                [ & ] { return start_in + ( iter++ ) * dx; } );

 Vec_CNumber nodeAscendant = x;
 Vec_CNumber meanAnomaly = x;

 double thetaValFinal = Theta_min;
 double altitudeFinal = altitude[ indexLen - 1 ];

 std::cout << "alpha_half: " << aHalf * 180 / PI << "\n";
 std::cout << "theta: " << thetaValFinal << "\n";
 std::cout << "altitude: " << altitudeFinal << "\n";

 FNumber t = horizon / time_step;

 boost::multi_array< double, 3 > CoverageSatLat(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );
 boost::multi_array< double, 3 > CoverageSatLong(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );
 boost::multi_array< double, 3 > CoverageSatLat1(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );
 boost::multi_array< double, 3 > CoverageSatLong1(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );

 double lat_Sat;
 double long_Sat;

 int indexOrbit;
 int indexOrbit1;

 iFile >> satellites;

 v_Block.resize( satellites );
 std::cout << "SATELLITES: " << satellites << "\n";

 int satellites1 = satellites;

 // satellites are split into three equally-sized groups, each assigned a
 // different altitude index ii (2, 1 or 0, i.e., roughly high/mid/low
 // altitude among the candidates found above), so that the constellation
 // is not confined to a single altitude band

 for( Index isat = 0 ; isat < satellites ; ++isat ) {
  int ialt = std::floor( satellites / 3 );
  int ii = 0;
  if( isat < ialt )
   ii = 2;
  if( isat > 2 * ialt )
   ii = 0;
  if( isat >= ialt and isat <= 2 * ialt )
   ii = 1;

  int index1 = -1;
  indexOrbit = 0;
  int addOrbit = 0;
  indexOrbit1 = 0;
  thetaValF = thetaValFinal / 3.0;
  std::cout << thetaValF << std::endl;

  // enumerate the full Cartesian product of the discretized inclination
  // (jj), RAAN (k) and mean anomaly (l): each triple, combined with the
  // altitude band ii fixed above, is one candidate orbital configuration
  // for the current satellite; index1 numbers them all (kept or discarded)

  for( Index jj = 0 ; jj < incSet ; ++jj ) {
   for( Index k = 0 ; k < ascSet ; ++k ) {
    for( Index l = 0 ; l < anmSet ; ++l ) {
     index1 += 1;
     indexOrbit1 = 0;
     for( Index j = 0 ; j < t ; ++j ) {
      // formula to compute the latitude of the projection of the position
      // of satellite onto the Earth surface corresponding to a given configuration
      lat_Sat = asin(
       ( ( sin( inclination[ jj ] ) * ( altitude[ ii ] + RAYON ) *
           sin( meanAnomaly[ l ] ) ) *
         cos( t_p[ ii ] * ( (j)*time_step ) ) / ( altitude[ ii ] + RAYON ) ) +
       ( ( sin( inclination[ jj ] ) * t_u[ ii ] * cos( meanAnomaly[ l ] ) ) *
         sin( t_p[ ii ] * ( (j)*time_step ) ) * t_GM[ ii ] ) );

      // formula to compute the longitude of the projection of the position
      // of satellite onto the Earth surface corresponding to a given configuration
      long_Sat = fmod(
       -( angle0 + ( WE * ( (j)*time_step ) ) ) +
        atan2(
         ( ( ( sin( nodeAscendant[ k ] ) * ( altitude[ ii ] + RAYON ) *
               cos( meanAnomaly[ l ] ) ) +
             ( cos( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
               ( altitude[ ii ] + RAYON ) * sin( meanAnomaly[ l ] ) ) ) *
           cos( t_p[ ii ] * ( (j)*time_step ) ) /
           ( altitude[ ii ] + RAYON ) ) +
          ( ( -( sin( nodeAscendant[ k ] ) * t_u[ ii ] *
                 sin( meanAnomaly[ l ] ) ) +
              ( cos( nodeAscendant[ k ] ) * cos( inclination[ jj ] ) *
                t_u[ ii ] * cos( meanAnomaly[ l ] ) ) ) *
            sin( t_p[ ii ] * ( (j)*time_step ) ) * t_GM[ ii ] ),
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
             sin( ( t_p[ ii ] * ( j * time_step ) ) ) * t_GM[ ii ] ) ) ),
       ( 2 * PI ) );
      if( long_Sat <= 0 )
       long_Sat += 2 * PI;

      for( Index i = 0 ; i < targets ; ++i ) {
       // compute the (geodedical) difference between the latitude of the target and the lat_Sat
       CoverageSatLat1[ i ][ j ][ index1 ] = std::abs(
        2 * asin( 0.5 * sqrt( 1 - cos( Latitude[ i ] - lat_Sat ) ) ) );
       CoverageSatLong1[ i ][ j ][ index1 ] =
        std::abs(
         2 *
         asin( 0.5 * sqrt( cos( Latitude[ i ] ) * cos( Latitude[ i ] ) *
                           ( 1 - cos( Longitude[ i ] - long_Sat ) ) ) ) ) *
        cos( Latitude[ i ] );
       // compute the scaled (geodedical) difference between the latitude of the target and the long_Sat
       if( cos( Latitude[ i ] ) < 0 )
        std::cout << "ERROR!" << "\n";
       if( CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValF and
           CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValF )
        indexOrbit1 += 1;
      }
     }
     // the orbital configuration that do not observe any satellite in any time-step
     // are discarded so that the solution space is maintened reasonably "small"
     if( indexOrbit1 >= 1 ) {
      for( Index j = 0 ; j < t ; ++j ) {
       for( Index i = 0 ; i < targets ; ++i ) {
        CoverageSatLat[ i ][ j ][ indexOrbit ] =
         CoverageSatLat1[ i ][ j ][ index1 ];
        CoverageSatLong[ i ][ j ][ indexOrbit ] =
         CoverageSatLong1[ i ][ j ][ index1 ];
       }
      }
      indexOrbit += 1;
     }
    }
   }
  }
  std::cout << "number Orbits: " << indexOrbit << "\n";

  // indexOrbit surviving configurations (out of incSet*ascSet*anmSet) are
  // passed to the isat-th SatelliteBlock as its set [C] of candidate
  // orbits, together with their precomputed coverage distances

  auto SB = new SatelliteBlock( this );
  SB->load( targets , time_step , horizon , altitudeFinal , thetaValF , indexOrbit ,
            aHalf, CoverageSatLat, CoverageSatLong, periods );
  v_Block[ isat ] = SB;
 }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

} // end( ConstellationBlock::load( const std::string & input )

/*--------------------------------------------------------------------------*/
// simply delegates to each nested SatelliteBlock

void ConstellationBlock::generate_abstract_variables( Configuration * stvv )
{
 for( auto blck : v_Block )
  blck->generate_abstract_variables();
}

/*--------------------------------------------------------------------------*/
// first delegates to each nested SatelliteBlock (its own per-satellite
// constraints), then adds the two families of ConstellationBlock-level
// constraints that link them together: a global bound on theta (thetaM)
// and the target observability constraints (observation / observation1)

void ConstellationBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( !AR ) {
  for( auto blck : v_Block )
   blck->generate_abstract_constraints();

  // thetaM: for every satellite k, thetaVar[ k ] <= 0.9 * thetaValF * zeta[ k ]
  // i.e., caps the per-satellite observability threshold, when the
  // satellite is active, to (a safety-scaled fraction of) the reference
  // value thetaValF used while generating the candidate orbits

  thetaM.resize( 1 );

  LinearFunction::v_coeff_pair v_var12;
  for( Index k = 0 ; k < satellites ; ++k ) {
   v_var12.push_back( std::make_pair(
    static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_theta(), 1.0 ) );
   v_var12.push_back(
    std::make_pair( static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_z() ,
                    -thetaValF * 0.9 ) );
  }
  thetaM[ 0 ].set_function( new LinearFunction( std::move( v_var12 ) ) ,
                            eNoBlck );
  thetaM[ 0 ].set_rhs( 0.0 , eNoBlck );
  thetaM[ 0 ].set_lhs( -Inf< double >() , eNoBlck );

  add_static_constraint( thetaM , "thetaM" );


  // generate the observability constraints  - - - - - - - - - - - - - - -
  /* sum_{t \in T(k, \Delta t[ m ], dt), i \in [s]} \xi[ i ][ t ][ m ] \geq 1,
  * \forall k \in [\lfloor T/\Delta t_m \rfloor], \forall m \in \mathcal{X}
  * where: T is the time horizon (here: horizon), Delta t[ m ] is the revisit
  * period associated with the target m (here: period[ m ]) and dt is the time
  * step for time discretization (here: time_step) T(k, \Delta t[ m ], dt) is
  * the set of the time step corresponding to the interval, in which target m
  * should be observed by the constellation.
  */

  // observation[ i ][ j ] enforces the revisit-time requirement (1) of the
  // class comments for target i in its j-th revisit window: since each
  // target may have a different period[ i ], the number of time stamps
  // pp = t / periods[ i ] making up one revisit window also differs across
  // targets, hence the per-target loop bound "j < periods[ i ]" below;
  // rows are allocated up to maxPeriods (the largest periods[] over all
  // targets) so that boost::multi_array can hold them all, and rows with
  // j >= periods[ i ] (which do not correspond to an actual revisit window
  // for target i) are filled in with the harmless "fake" constraint below

  double maxPeriods = *max_element( periods.begin() , periods.end() );
  observation.resize(
   boost::multi_array< FRowConstraint, 2 >::extent_gen()[ targets ]
                                                        [ maxPeriods ] );
  double pp;

  for( Index i = 0 ; i < targets ; ++i ) {
   pp = horizon / time_step / periods[ i ];
   for( Index j = 0 ; j < periods[ i ] ; ++j ) {
    LinearFunction::v_coeff_pair v_var;

    for( Index k = 0 ; k < satellites ; ++k ) {
     for( Index tt = j * pp ; tt < ( j + 1 ) * pp ; ++tt ) {
      // retrieve observation variable \xi[ i ][ t ][ m ] for SatelliteBlock i
      v_var.push_back( std::make_pair(
       static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_r( i , tt ),
       1.0 ) );
     }
    }
    observation[ i ][ j ].set_function(
     new LinearFunction( std::move( v_var ) ), eNoBlck );
    observation[ i ][ j ].set_rhs( Inf< double >() , eNoBlck );
    observation[ i ][ j ].set_lhs( 1.0 , eNoBlck );
   }

   for( Index j = periods[ i ] ; j < maxPeriods ; ++j ) {
    LinearFunction::v_coeff_pair v_var;

    for( Index k = 0 ; k < satellites ; ++k ) {
     v_var.push_back( std::make_pair(
      static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_z(), 0.0 ) );
    }

    // fake constraints when j \geq periods[ i ]: we simply set 0 * z[ i ] == 0
    // for SatelliteBlock i \in [s], where s is the total number of the satellite
    observation[ i ][ j ].set_function(
     new LinearFunction( std::move( v_var ) ), eNoBlck );
    observation[ i ][ j ].set_rhs( 0.0 , eNoBlck );
    observation[ i ][ j ].set_lhs( 0.0 , eNoBlck );
   }
  }

  add_static_constraint( observation , "observation" );

  // observation1[ i ][ j ] is a *stronger*, per-time-stamp version of the
  // same idea: at most one satellite observes target i at time stamp j
  // (sum_k xi[ k ][ i ][ j ] <= 1), which rules out redundant simultaneous
  // observations of the same target by several satellites

  FNumber t = horizon / time_step;

  observation1.resize(
   boost::multi_array< FRowConstraint, 2 >::extent_gen()[ targets ][ t ] );

  for( Index i = 0 ; i < targets ; ++i ) {
   for( Index j = 0 ; j < t ; ++j ) {
    LinearFunction::v_coeff_pair v_var1;
    for( Index k = 0 ; k < satellites ; ++k )
     v_var1.push_back( std::make_pair(
      static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_r( i , j ), 1.0 ) );

    observation1[ i ][ j ].set_function(
     new LinearFunction( std::move( v_var1 ) ), eNoBlck );
    observation1[ i ][ j ].set_rhs( 1.0 , eNoBlck );
    observation1[ i ][ j ].set_lhs( -Inf< double >() , eNoBlck );
   }
  }

  add_static_constraint( observation1 , "observation1" );

  std::cout << "Constraints charged!\n";
 }
 AR = true;

} // end( ConstellationBlock::generate_abstract_constraints() )

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR CHECKING THE ConstellationBlock ----------*/
/*--------------------------------------------------------------------------*/
// checks feasibility of the ConstellationBlock-level constraints only
// (thetaM, observation and observation1); feasibility of each nested
// SatelliteBlock is not checked here, as it is assumed each of them is
// individually verified via its own is_feasible()

bool ConstellationBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 1e-1;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ &tol, &rel_viol ]( Configuration * c ) -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return ( true );
  }
  if( auto tc =
       dynamic_cast< SimpleConfiguration< std::pair< double, int > > * >(
        c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return ( true );
  }
  return ( false );
 };

 if( ( !extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return (
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  RowConstraint::is_feasible( thetaM , tol , rel_viol ) &&
  RowConstraint::is_feasible( observation , tol , rel_viol ) &&
  RowConstraint::is_feasible( observation1 , tol , rel_viol ) );

} // end( ConstellationBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
// TODO: printing the ConstellationBlock instance is not implemented yet

void ConstellationBlock::print( std::ostream & output , char vlvl ) const {}

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void ConstellationBlock::guts_of_destructor( void )
{
 /* clear() all Constraint to ensure that they do not bother to un-register
    themselves from Variable that are going to be deleted anyway. Then
    deletes all the "abstract representation", if any. */

 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

} // end( guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*---------------------- End File ConstellationBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
