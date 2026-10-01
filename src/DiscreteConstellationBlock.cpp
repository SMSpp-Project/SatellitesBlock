/*--------------------------------------------------------------------------*/
/*------------------ File DiscreteConstellationBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DiscreteConstellationBlock class.
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

#include "DiscreteConstellationBlock.h"

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

static constexpr Index ell = 3; // the levels of the threshold

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DiscreteConstellationBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( DiscreteConstellationBlock );

/*--------------------------------------------------------------------------*/
/*----------------- METHODS OF DiscreteConstellationBlock ------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteConstellationBlock::load( const std::string & input , char frmt )
{
 std::ifstream iFile( input );
 if( ! iFile.is_open() )
  throw( std::invalid_argument(
   "DiscreteConstellationBlock::load: cannot open file " + input ) );

 load( iFile , frmt );

 } // end( DiscreteConstellationBlock::load( std::string ) )

/*--------------------------------------------------------------------------*/
// the orbital-mechanics computations are those of ConstellationBlock::load()
// [see there], up to the ground-track propagation; what changes is that
// the threshold is discretized into the ell levels thetaValF[] and that,
// for every satellite / time stamp / target / orbit / level, the 0/1
// observability outcome is stored in obs[][][][][]

void DiscreteConstellationBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // read the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( input >> eatcomments >> horizon ) || ( horizon <= 0 ) )
  throw( std::invalid_argument(
   "DiscreteConstellationBlock::load: error reading the horizon" ) );
 horizon *= 3600.0; // from hours to seconds

 if( ! ( input >> eatcomments >> time_step ) || ( time_step <= 0 ) )
  throw( std::invalid_argument(
   "DiscreteConstellationBlock::load: error reading the time step" ) );

 if( ! ( input >> eatcomments >> targets ) || ( ! targets ) )
  throw( std::invalid_argument( "DiscreteConstellationBlock::load: error "
                                "reading the number of targets" ) );

 std::vector< double > Latitude( targets );
 std::vector< double > Longitude( targets );
 periods.resize( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
  if( ! ( input >> eatcomments >> periods[ i ] >> Latitude[ i ] >>
         Longitude[ i ] ) ||
      ( ! periods[ i ] ) )
   throw( std::invalid_argument(
    "DiscreteConstellationBlock::load: error reading the targets" ) );
  if( std::abs( Latitude[ i ] ) > 90 )
   throw( std::invalid_argument(
    "DiscreteConstellationBlock::load: latitude out of [ -90 , 90 ]" ) );
  Latitude[ i ] *= PI / 180;
  Longitude[ i ] *= PI / 180;
  }

 if( ! ( input >> eatcomments >> satellites ) || ( ! satellites ) )
  throw( std::invalid_argument( "DiscreteConstellationBlock::load: error "
                                "reading the number of satellites" ) );

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
  throw( std::invalid_argument(
   "DiscreteConstellationBlock::load: the "
   "horizon allows less than three orbit altitudes" ) );

 const Index altSet = altitude.size();

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

 // Theta_min is the angle swept by the satellite on the lowest orbit in one
 // time step, times FACTOR; it is both the resolution of the discretization
 // of the orbits and the loosest threshold

 const double Theta_min =
  ( ( 2 * PI * time_step ) / ( 2 * periodSat[ altSet - 1 ] ) ) * FACTOR;

 // inclination is sampled uniformly in [ 0 , PI ], RAAN and mean anomaly
 // in [ 0 , 2 PI ], each with numbOfDiscretize points

 const Index numbOfDiscretize = ceil( PI / Theta_min );
 if( numbOfDiscretize < 2 )
  throw( std::invalid_argument(
   "DiscreteConstellationBlock::load: time step too long w.r.t. the "
   "orbital period" ) );

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

 // the ell levels of the threshold, from the loosest (level 0) to the
 // tightest (level ell - 1), the same for all the satellites

 std::vector< double > thetaValF( ell );
 for( Index l = 0 ; l < ell ; ++l )
  thetaValF[ l ] = double( ell - l ) / ell * Theta_min / 10.0;

 const Index t = horizon / time_step;

 boost::multi_array< double , 3 > CoverageSatLat1(
  boost::extents[ targets ][ t ][ nOrbits ] );
 boost::multi_array< double , 3 > CoverageSatLong1(
  boost::extents[ targets ][ t ][ nOrbits ] );

 obs.resize( boost::extents[ satellites ][ t ][ targets ][ nOrbits ][ ell ] );
 v_Block.resize( satellites );
 indexOrbitSat.resize( satellites );
 thetaValSat.resize( boost::extents[ satellites ][ ell ] );
 ellSat.resize( satellites );

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
       if( ( CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValF[ 0 ] ) &&
           ( CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValF[ 0 ] ) )
        observes = true;
       }
      }

     // an orbit is kept only if it observes some target at some time
     // stamp with the loosest threshold thetaValF[ 0 ]; for a kept orbit,
     // obs[ isat ][ j ][ i ][ indexOrbit ][ l1 ] is 1 iff both distances
     // are within the threshold thetaValF[ l1 ]

     if( ! observes )
      continue;

     for( Index j = 0 ; j < t ; ++j )
      for( Index i = 0 ; i < targets ; ++i )
       for( Index l1 = 0 ; l1 < ell ; ++l1 )
        obs[ isat ][ j ][ i ][ indexOrbit ][ l1 ] =
         ( ( CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValF[ l1 ] ) &&
           ( CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValF[ l1 ] ) )
          ? 1.0
          : 0.0;

     ++indexOrbit;
     }

  // record the number of kept orbits and the levels of the satellite, and
  // create its DiscreteSatelliteBlock with the (orbit , level) grid

  indexOrbitSat[ isat ] = indexOrbit;
  for( Index l = 0 ; l < ell ; ++l )
   thetaValSat[ isat ][ l ] = thetaValF[ l ];
  ellSat[ isat ] = ell;

  auto SB = new DiscreteSatelliteBlock( this );
  SB->load( indexOrbit , ell );
  v_Block[ isat ] = SB;
  }

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( DiscreteConstellationBlock::load( std::istream ) )

/*--------------------------------------------------------------------------*/

void DiscreteConstellationBlock::generate_abstract_variables(
 Configuration * stvv )
{
 for( auto blck : v_Block )
  blck->generate_abstract_variables();

 } // end( DiscreteConstellationBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void DiscreteConstellationBlock::generate_abstract_constraints(
 Configuration * stcc )
{
 if( AR ) // the constraints are there already
  return; // nothing to do

 // the constraints below use the Variable of the sub-Block
 for( auto blck : v_Block ) {
  blck->generate_abstract_variables();
  blck->generate_abstract_constraints();
  }

 auto y = [ this ]( Index k , Index o , Index l ) {
  return(
   static_cast< DiscreteSatelliteBlock * >( v_Block[ k ] )->i2p_y( o , l ) );
  };

 // thetaM, constraint (3) of the class comments: the average threshold of
 // the active satellites is at most 0.9 times the loosest level

 LinearFunction::v_coeff_pair v_var12;
 for( Index k = 0 ; k < satellites ; ++k )
  for( Index o = 0 ; o < indexOrbitSat[ k ] ; ++o )
   for( Index l = 0 ; l < ellSat[ k ] ; ++l )
    v_var12.push_back( std::make_pair(
     y( k , o , l ) , thetaValSat[ k ][ l ] - 0.9 * thetaValSat[ k ][ 0 ] ) );

 thetaM.set_function( new LinearFunction( std::move( v_var12 ) ) , eNoMod );
 thetaM.set_rhs( 0.0 , eNoMod );
 thetaM.set_lhs( -Inf< double >() , eNoMod );

 add_static_constraint( thetaM , "thetaM" );

 // observation, constraints (1) of the class comments: target i is
 // observed at least once in each of its revisit windows j, made of the pp
 // time stamps [ j pp , ( j + 1 ) pp ); the coefficient of y[ o ][ l ] of
 // satellite k is the number of time stamps of the window in which it
 // observes target i. Rows are allocated up to the largest number of
 // periods, those with j >= periods[ i ] being no revisit window of target
 // i: they get the void constraint 0 == 0, with every y appearing once
 // with coefficient 0

 const Index maxPeriods =
  *std::max_element( periods.begin() , periods.end() );
 observation.resize( boost::extents[ targets ][ maxPeriods ] );
 const Index nt = obs.shape()[ 1 ]; // the number of time stamps

 for( Index i = 0 ; i < targets ; ++i ) {
  const double pp = horizon / time_step / periods[ i ];
  for( Index j = 0 ; j < maxPeriods ; ++j ) {
   const bool window = j < periods[ i ];
   LinearFunction::v_coeff_pair v_var;
   for( Index k = 0 ; k < satellites ; ++k )
    for( Index o = 0 ; o < indexOrbitSat[ k ] ; ++o )
     for( Index l = 0 ; l < ellSat[ k ] ; ++l ) {
      double obb = 0.0;
      if( window )
       for( Index tt = j * pp ; ( tt < ( j + 1 ) * pp ) && ( tt < nt ) ;
            ++tt )
        obb += obs[ k ][ tt ][ i ][ o ][ l ];
      v_var.push_back( std::make_pair( y( k , o , l ) , obb ) );
      }

   observation[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var ) ) , eNoBlck );
   observation[ i ][ j ].set_rhs( window ? Inf< double >() : 0.0 , eNoBlck );
   observation[ i ][ j ].set_lhs( window ? 1.0 : 0.0 , eNoBlck );
   }
  }

 add_static_constraint( observation , "observation" );

 // observation1, constraints (2) of the class comments: no two satellites
 // observe target i at time stamp j

 observation1.resize( boost::extents[ targets ][ nt ] );

 for( Index i = 0 ; i < targets ; ++i )
  for( Index j = 0 ; j < nt ; ++j ) {
   LinearFunction::v_coeff_pair v_var;
   for( Index k = 0 ; k < satellites ; ++k )
    for( Index o = 0 ; o < indexOrbitSat[ k ] ; ++o )
     for( Index l = 0 ; l < ellSat[ k ] ; ++l )
      v_var.push_back(
       std::make_pair( y( k , o , l ) , obs[ k ][ j ][ i ][ o ][ l ] ) );

   observation1[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var ) ) , eNoBlck );
   observation1[ i ][ j ].set_rhs( 1.0 , eNoBlck );
   observation1[ i ][ j ].set_lhs( -Inf< double >() , eNoBlck );
   }

 add_static_constraint( observation1 , "observation1" );

 AR = true;

 } // end( DiscreteConstellationBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR CHECKING THE DiscreteConstellationBlock -----------*/
/*--------------------------------------------------------------------------*/

bool DiscreteConstellationBlock::is_feasible( bool useabstract ,
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

 } // end( DiscreteConstellationBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------ METHODS FOR PRINTING & SAVING THE DiscreteConstellationBlock ------*/
/*--------------------------------------------------------------------------*/

void DiscreteConstellationBlock::print( std::ostream & output ,
                                        char vlvl ) const
{
 output << "DiscreteConstellationBlock: horizon " << horizon
        << " s, time step " << time_step << " s, " << targets << " targets, "
        << satellites << " satellites" << std::endl;
 for( Index k = 0 ; k < satellites ; ++k )
  output << "satellite " << k << ": " << indexOrbitSat[ k ] << " orbits, "
         << ellSat[ k ] << " threshold levels" << std::endl;

 } // end( DiscreteConstellationBlock::print )

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteConstellationBlock::guts_of_destructor( void )
{
 // clear() all the Constraint, so that they do not bother to un-register
 // themselves from the Variable of the sub-Block, deleted below
 thetaM.clear();
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

 } // end( DiscreteConstellationBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*---------------- End File DiscreteConstellationBlock.cpp -----------------*/
/*--------------------------------------------------------------------------*/
