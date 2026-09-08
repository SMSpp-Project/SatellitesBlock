/*--------------------------------------------------------------------------*/
/*------------------------ File MultiTargetBlockv2.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MultiTargetBlockv2 class.
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

#include "MultiTargetBlockv2.h"
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

using FNumber = MultiTargetBlockv2::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();
static const double RAYON = 6378136.3;
static const double M_limit = 10 * 3.14159265;
static const auto PI = 3.14159265;
static const auto MU = 3.986004418e14; //[m^3/s^2]
static const auto WE = 7.2921e-5;
static const auto FACTOR = 1.2;
static const auto angle0 = -1.3882860164509252;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register MultiTargetBlockv2 to the Block factory

SMSpp_insert_in_factory_cpp_1( MultiTargetBlockv2 );

// register MultiTargetSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( MultiTargetSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF MultiTargetBlockv2 ----------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

// reads the instance and builds the candidate orbits exactly as
// MultiTargetBlock::load() does [see there for the details of the shared,
// single-altitude ground-track propagation]: the only difference is that,
// having no per-target sub-Block to create, this just stores the resulting
// CoverageSatLat/CoverageSatLong and the problem sizes (n, t, OrbitSet,
// targets, ...) directly as members, deferring Variable/Constraint
// creation to generate_abstract_variables()/generate_abstract_constraints()

void MultiTargetBlockv2::load( const std::string & input , char frmt )
{
 // ensure starting from clean slate
 guts_of_destructor();

 std::ifstream iFile( input );

 iFile >> horizon;
 horizon *= 3600.0;
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

 for( Index i = 0 ; i < horizon / 3600.0 ; ++i ) {
  j++;
  altitudeSetVal =
   cbrt( ( MU * pow( ( horizon ) / j , 2.0 ) ) / ( 4.0 * pow( PI , 2.0 ) ) ) -
   RAYON;
  if( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 ) {
   indexLen++;
  }
 }

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
  //std::cout << altitudeSetVal << "\n";
  if( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 ) {
   altitude[ indexLen ] = altitudeSetVal;
   periodSat[ indexLen ] = horizon / j;
   indexLen++;
  }
 }

 Latitude.resize( targets );
 Longitude.resize( targets );
 for( Index i = 0 ; i < targets ; ++i ) {
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

 double alphalim;
 double Theta_min =
  ( ( 2 * PI * time_step ) / ( 2 * periodSat[ indexLen - 1 ] ) ) * FACTOR;
 std::cout << "theta_min: " << Theta_min << "\n";
 double aHalf =
  atan( sin( Theta_min ) /
        ( ( RAYON + altitude[ indexLen - 1 ] ) / RAYON - cos( Theta_min ) ) );

 int alpha_lim_flag = 0;

 for( Index ii = 0 ; ii < altSet ; ++ii ) {
  if( ( ( RAYON + altitude[ ii ] ) / RAYON ) * sin( aHalf ) > 1 ) {
   aHalf = asin( ( RAYON / ( RAYON + altitude[ ii ] ) ) );
   std::cout << "WARNING: computed alpha_lim\n";
   alpha_lim_flag = 1;
   break;
  }
 }

 double alt = altitude[ altSet - 1 ];
 //double alt = altitude[0];
 for( Index ii = 0 ; ii < altSet ; ++ii ) {
  t_p[ ii ] = sqrt( MU / ( pow( RAYON + altitude[ ii ] , 3 ) ) );
  t_u[ ii ] = sqrt( MU / ( RAYON + altitude[ ii ] ) );
  t_GM[ ii ] = sqrt( ( RAYON + altitude[ ii ] ) / MU );
  thetaVal[ ii ] =
   -aHalf + asin( ( ( RAYON + altitude[ ii ] ) / RAYON ) * sin( aHalf ) );
 }

 //double Theta_max = thetaVal[altSet-1];
 std::cout << "alpha_flag=" << alpha_lim_flag << "\n";

 /***********************/
 if( alpha_lim_flag == 1 ) {
  /////Theta_min = thetaVal[altSet-1];//*1.2;
  std::cout << "theta_min: " << Theta_min << "\n";
 }
 /***********************/

 //Theta_max = Theta_min;
 double numbOfDiscretize = ceil( PI / Theta_min );

 std::cout << "numbOfDiscretize: " << numbOfDiscretize << "\n";

 FNumber incSet = numbOfDiscretize;
 FNumber ascSet = numbOfDiscretize;
 FNumber anmSet = numbOfDiscretize;

 double start_in = 0;
 double end_in = PI;
 size_t num_in = numbOfDiscretize;

 double dx = ( end_in - start_in ) / ( num_in - 1 );
 std::vector< double > x( num_in );
 int iter = 0;
 std::generate( x.begin() , x.end() ,
                [ & ] { return start_in + ( iter++ ) * dx; } );

 Vec_CNumber inclination = x;
 //std::cout << inclination;

 end_in = 2 * PI;
 dx = ( end_in - start_in ) / ( num_in - 1 );
 iter = 0;
 std::generate( x.begin() , x.end() ,
                [ & ] { return start_in + ( iter++ ) * dx; } );

 Vec_CNumber nodeAscendant = x;
 Vec_CNumber meanAnomaly = x;

 thetaValFinal = Theta_min;
 double altitudeFinal = altitude[ indexLen - 1 ];
 //double altitudeFinal = altitude[0];

 //std::cout << "theta_min: " << Theta_min << "\n";
 std::cout << "alpha_half: " << aHalf * 180 / PI << "\n";
 std::cout << "theta: " << thetaValFinal << "\n";
 std::cout << "altitude: " << altitudeFinal << "\n";

 t = horizon / time_step;

 CoverageSatLat.resize(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );
 CoverageSatLong.resize(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );
 boost::multi_array< double, 3 > CoverageSatLat1(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );
 boost::multi_array< double, 3 > CoverageSatLong1(
  boost::extents[ targets ][ t ][ incSet * ascSet * anmSet ] );

 double lat_Sat;
 double long_Sat;

 //for( Index ii = 0 ; ii < altSet ; ++ii )
 //{
 int ii = indexLen - 1;
 //int ii = 0;
 int index1 = -1;
 indexOrbit = 0;
 int addOrbit = 0;
 int indexOrbit1 = 0;
 for( Index jj = 0 ; jj < incSet ; ++jj ) {
  for( Index k = 0 ; k < ascSet ; ++k ) {
   for( Index l = 0 ; l < anmSet ; ++l ) {
    index1 += 1;
    indexOrbit1 = 0;
    for( Index j = 0 ; j < t ; ++j ) {
     lat_Sat = asin(
      ( ( sin( inclination[ jj ] ) * ( altitude[ ii ] + RAYON ) *
          sin( meanAnomaly[ l ] ) ) *
        cos( t_p[ ii ] * ( (j)*time_step ) ) / ( altitude[ ii ] + RAYON ) ) +
      ( ( sin( inclination[ jj ] ) * t_u[ ii ] * cos( meanAnomaly[ l ] ) ) *
        sin( t_p[ ii ] * ( (j)*time_step ) ) * t_GM[ ii ] ) );

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
     /*
             for( Index i = 0 ; i < targets ; ++i )
             {
                CoverageSatLat[ i ][ j ][ index1 ] = std::abs(Latitude[i] - lat_Sat);
                CoverageSatLong[ i ][ j ][ index1 ] = std::abs(Longitude[i] - long_Sat) * cos(Latitude[i]);
             }}
             indexOrbit += 1;
             */
     for( Index i = 0 ; i < targets ; ++i ) {
      CoverageSatLat1[ i ][ j ][ index1 ] =
       //        std::abs(2 * asin( 0.5 * sqrt(1 - cos(Latitude[i]-lat_Sat))));
       std::abs( Latitude[ i ] - lat_Sat );
      CoverageSatLong1[ i ][ j ][ index1 ] =
       //        std::abs(2 * asin( 0.5 * sqrt(cos(Latitude[i]) * cos(Latitude[i]) * (1 - cos(Longitude[i] - long_Sat))))) * cos(Latitude[i]);
       std::abs( Longitude[ i ] - long_Sat ) * cos( Latitude[ i ] );
      if( cos( Latitude[ i ] ) < 0 )
       std::cout << "ERROR!" << "\n";
      //std::cout << CoverageSatLat1[ i ][ j ][ index1 ] << " " << CoverageSatLong1[ i ][ j ][ index1 ] << std::endl;
      //std::cout << i << " " << j << " " << ii << " " << jj << " " << k << " " << l << " " << "\n";
      if( CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValFinal and
          CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValFinal ) {
       //std::cout << i << " " << j << " " << index1 << " " << "\n";
       indexOrbit1 += 1;
      }
     }
    }
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
 //}

 iFile >> satellites;
 std::cout << "number Orbits: " << indexOrbit << "\n";
 std::cout << "SATELLITES: " << satellites << "\n";

 n = satellites;
 dt = time_step;
 T = horizon;
 t = T / dt;
 alphaHalf = aHalf;

 OrbitSet = indexOrbit;

 //thetaVal = thetaValues;
 //std::cout << "thetaVal=" << thetaValues << "\n";

} // end( MultiTargetBlockv2::load( memory ) )

/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::load( std::istream & input , char frmt )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 guts_of_destructor();

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // allocate observability variables - - - - - - - - - - - - - - - - - - - - -

 //generate_abstract_variables();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

} // end( MultiTargetBlockv2::load( istream ) )

/*--------------------------------------------------------------------------*/

// creates the static Variable of the model: these are exactly the ones of
// SingleTargetBlock.h (theta, Deltat, Deltat_k1/k2 and their b1/b2/d1/d2
// linearization helpers, zeta, h, xi), each carrying an extra "targets"
// dimension (the last index) since a single MultiTargetBlockv2 covers all
// targets at once, with the sole exception of activation[ i ][ j ] (the
// per-satellite orbit selection), which is *not* replicated per target,
// being shared by construction (cf. the class comments in
// MultiTargetBlockv2.h on why no "duplicate" constraints are needed here)

void MultiTargetBlockv2::generate_abstract_variables( Configuration * stvv )
{
 if( AR3 & HasVar ) // the variables are there already
  return; // nothing to do

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
   //for( Index k = 0 ; k < targets ; ++k )
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

// generates the same static Constraint families as
// SingleTargetBlock::generate_abstract_constraints() [see there for the
// full explanation, referring to the numbered constraints (1)-(10) of
// SingleTargetBlock.h], each now built once per target (the extra "[ k ]"
// / "targets" loop) rather than once per SingleTargetBlock instance; since
// activation[][] is shared (not duplicated) across targets, no counterpart
// to MultiTargetBlock's duplicate_pi/duplicate_theta constraints is needed
// (the dead code below is what such a duplicate_pi would have looked like,
// kept only for reference, cf. MultiTargetBlock::generate_abstract_constraints())

void MultiTargetBlockv2::generate_abstract_constraints( Configuration * stcc )
{
 if( AR2 & HasCnst ) // the constraints are there already
  return; // nothing to do
 /*
  duplicate_pi.resize( boost::multi_array_types::extent_gen()[ targets - 1  ][ satellites ][ OrbitSet ] );
    for( Index i = 0 ; i < targets - 1  ; ++i ) {
       for( Index j = 0 ; j < satellites; ++j ) {
          for( Index k = 0 ; k < OrbitSet; ++k ) {
             LinearFunction::v_coeff_pair v_vars_pi;
             v_vars_pi.push_back( std::make_pair( &activation[j][k][i], 1.0 ));
             v_vars_pi.push_back( std::make_pair( &activation[j][k][i+1], -1.0 ));
             duplicate_pi[i][j][k].set_function( new LinearFunction( std::move( v_vars_pi )) , eNoBlck );
             duplicate_pi[i][j][k].set_rhs( 0.0 , eNoBlck );
             duplicate_pi[i][j][k].set_lhs( 0.0 , eNoBlck );
          }
       }
    }
*/
 // orbitSelection[ i ]: sum_{j} activation[ i ][ j ] == 1 for satellite i;
 // since activation[][] has no target dimension, this constraint does not
 // depend on k and is therefore built only once (not once per target)

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
    //if ( j > i ) {
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
    //if ( j > i ) {
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
    //if ( j > i ) {
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

 // generate observability constraints via big-M approach

 //Thetamax[i]+(1-Xi[i,p,j])*M_limit >= (sum(act4[i,a,k,l,s]*CoverageSatLat[j,p,s,l,k,a]
 //      for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

 //Thetamax[i]+(1-Xi[i,p,j])*M_limit >= (sum(act4[i,a,k,l,s]*CoverageSatLong[j,p,s,l,k,a]
 //      for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

 //Thetamax[i]-(Xi[i,p,j])*M_limit <= (sum(act4[i,a,k,l,s]*CoverageSatLat[j,p,s,l,k,a]
 //      for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

 //Thetamax[i]-(Xi[i,p,j])*M_limit <= (sum(act4[i,a,k,l,s]*CoverageSatLong[j,p,s,l,k,a]
 //      for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

 obs1_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );
 obs3_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );

 obs2_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );
 obs4_cnst.resize(
  boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );

 double MLAT = PI;
 double MLONG = 2 * PI;

 for( Index k = 0 ; k < targets ; ++k ) {
  for( Index i = 0 ; i < n ; ++i ) {
   for( Index j = 0 ; j < t ; ++j ) {
    MLAT = 0.0;
    MLONG = 0.0;
    LinearFunction::v_coeff_pair v_obs1, v_obs2;
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

    obs1_cnst[ i ][ j ][ k ].set_rhs( ( -thetaValFinal ) , eNoBlck );
    obs1_cnst[ i ][ j ][ k ].set_lhs( -Inf< double >() , eNoBlck );
    obs1_cnst[ i ][ j ][ k ].set_function( Funct1 , eNoBlck );

    obs3_cnst[ i ][ j ][ k ].set_rhs( ( -thetaValFinal ) , eNoBlck );
    obs3_cnst[ i ][ j ][ k ].set_lhs( -Inf< double >() , eNoBlck );
    obs3_cnst[ i ][ j ][ k ].set_function( Funct2 , eNoBlck );
   }
  }
 }

 //add_static_constraint( obs1_cnst );
 //add_static_constraint( obs3_cnst );

 add_static_constraint( obs2_cnst );
 add_static_constraint( obs4_cnst );

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
    //if ( j > i ){
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
    //if ( j > i ){
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
    //if ( j > i ){
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
    //if ( j > i ){
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
  boost::multi_array< FRowConstraint, 2 >::extent_gen()[ targets ][ t ] );

 for( Index i = 0 ; i < targets ; ++i ) {
  for( Index j = 0 ; j < t ; ++j ) {
   LinearFunction::v_coeff_pair v_var1;
   for( Index k = 0 ; k < n ; ++k )
    v_var1.push_back( std::make_pair( &xi[ k ][ j ][ i ] , 1.0 ) );

   observation1[ i ][ j ].set_function(
    new LinearFunction( std::move( v_var1 ) ), eNoBlck );
   observation1[ i ][ j ].set_rhs( 1.0 , eNoBlck ); //Inf< double >()
   observation1[ i ][ j ].set_lhs( -Inf< double >() , eNoBlck );
  }
 }

 add_static_constraint( observation1 , "observation1" );

 AR2 |= HasCnst;
} // end( MultiTargetBlockv2::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

// the objective is the average of the per-target maximum revisit times
// Deltat[ k ] (weight 1/targets each), directly summed here since there
// are no per-target sub-Block to sum objectives over, unlike MultiTargetBlock

void MultiTargetBlockv2::generate_objective( Configuration * objc )
{
 if( AR1 & HasObj ) // the objective is there already
  return; // cowardly (and silently) return

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

// TODO: not implemented yet, cf. SatelliteBlock::is_feasible()

bool MultiTargetBlockv2::is_feasible( bool useabstract , Configuration * fsbc )
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

} // end( MultiTargetBlockv2::is_feasible )

/*--------------------------------------------------------------------------*/

// TODO: not implemented yet, cf. SatelliteBlock::is_optimal()

bool MultiTargetBlockv2::is_optimal( bool useabstract , Configuration * optc )
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

} //  end( MultiTargetBlockv2::is_optimal )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
// creates a new MultiTargetSolution, read from the current status of the
// MultiTargetBlockv2 unless an empty one (emptys == true) is asked for;
// wsol is extracted but currently unused, cf. SatelliteBlock::get_Solution()

Solution * MultiTargetBlockv2::get_Solution( Configuration * solc ,
                                             bool emptys )
{
 int wsol = 0;
 if( ( !solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;

 auto * sol = new MultiTargetSolution();

 if( !emptys )
  sol->read( this );

 return ( sol );

} // end( MultiTargetBlockv2::get_Solution )

// sets the values of the Deltat ColVariable-s (one per target, unlike the
// single scalar Deltat/zeta of SingleTargetBlock/SatelliteBlock) from the
// range [ rng.first , rng.second ) of fstrt

void MultiTargetBlockv2::set_zeta( c_Vec_FNumber_it fstrt , Range rng )
{
 if( !( AR3 & HasVar ) ) // nowhere to put the value in
  return; // cowardly (and silently) return

 Index i = rng.first;
 Index stop = std::min( Index( rng.second ) , Index( Deltat.size() ) );

 for( auto xi = Deltat.begin() + i ; i < stop ; ++i )
  ( xi++ )->set_value( *( fstrt++ ) );

} // end( MultiTargetBlockv2::set_zeta )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
// intercepts Modification concerning this Block before forwarding them to
// Block::add_Modification(), cf. SatelliteBlock::add_Modification()

void MultiTargetBlockv2::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
 }

 Block::add_Modification( mod , chnl );

} // end( MultiTargetBlockv2::add_Modification )

/*--------------------------------------------------------------------------*/
/*----- METHODS FOR LOADING, PRINTING & SAVING THE MultiTargetBlockv2 ------*/
/*--------------------------------------------------------------------------*/
// TODO: printing the MultiTargetBlockv2 instance is not implemented yet

void MultiTargetBlockv2::print( std::ostream & output , char vlvl ) const {

} // end( MultiTargetBlockv2::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::guts_of_destructor( void )
{
 // clear() all Constraint to ensure that they do not bother to un-register
 // themselves from Variable that are going to be deleted anyway

 // clear the bound constraints
 // note: as in SingleTargetBlock::guts_of_destructor(), Deltat_max3 is
 // never populated by generate_abstract_constraints() and is therefore
 // intentionally not cleared here; obs1_cnst and obs3_cnst are likewise
 // declared and filled in but never actually added as static constraints
 // (see the commented-out add_static_constraint() calls there), so they
 // are not registered against any Variable and do not need clear() either

 Constraint::clear( orbitSelection );
 Constraint::clear( Deltat_max1 );
 Constraint::clear( Deltat_max11 );
 Constraint::clear( Deltat_max2 );
 Constraint::clear( Deltat_max22 );
 Constraint::clear( Deltat_max_dt1 );
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
 Constraint::clear( activationSat1_cnst );
 Constraint::clear( theta_UB );
 //Constraint::clear( obs1_cnst );
 Constraint::clear( obs2_cnst );
 //Constraint::clear( obs3_cnst );
 Constraint::clear( obs4_cnst );
 Constraint::clear( obs_cnst_h );
 Constraint::clear( obs_cnst_xi );
 Constraint::clear( observation1 );

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

} // end( MultiTargetBlockv2::guts_of_destructor )

/*--------------------------------------------------------------------------*/

// no Modification is currently supported on a MultiTargetBlockv2, hence
// this is a no-op stub, cf. SatelliteBlock::guts_of_add_Modification()

void MultiTargetBlockv2::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -

} // end( MultiTargetBlockv2::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF MultiTargetSolution --------------------*/
/*--------------------------------------------------------------------------*/
// a MultiTargetSolution stores in v_zeta the value of Deltat[ k ] (the
// revisit time) for each of the targets k; deserialize()/serialize() are
// currently no-ops, i.e., loading/saving a MultiTargetSolution to/from a
// netCDF::NcGroup is not implemented yet

void MultiTargetSolution::deserialize( const netCDF::NcGroup & group ) {}

/*--------------------------------------------------------------------------*/
// reads the current value of Deltat[] out of the given MultiTargetBlockv2
// into v_zeta; note that this only happens if v_zeta is already non-empty
// (i.e., this MultiTargetSolution has previously been sized to hold it)

void MultiTargetSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const MultiTargetBlockv2 * >( block );
 if( !SATB )
  throw( std::invalid_argument( "block is not a MultiTargetBlockv2" ) );

 if( !v_zeta.empty() ) {
  v_zeta.resize( SATB->targets );

  for( Index k = 0 ; k < v_zeta.size() ; ++k )
   v_zeta[ k ] = SATB->get_zeta( k );
 }
}

/*--------------------------------------------------------------------------*/
// writes the values of Deltat[] stored in v_zeta into the given
// MultiTargetBlockv2

void MultiTargetSolution::write( Block * block )
{
 auto SATB = dynamic_cast< MultiTargetBlockv2 * >( block );
 if( !SATB )
  throw( std::invalid_argument( "block is not a MultiTargetBlockv2" ) );

 if( !v_zeta.empty() )
  SATB->set_zeta( v_zeta.begin() );
}

/*--------------------------------------------------------------------------*/

void MultiTargetSolution::serialize( netCDF::NcGroup & group ) const {}

/*--------------------------------------------------------------------------*/

MultiTargetSolution * MultiTargetSolution::scale( double factor ) const
{
 auto * sol = MultiTargetSolution::clone( true );
 return ( sol );
}

/*--------------------------------------------------------------------------*/
// TODO: not implemented yet

void MultiTargetSolution::sum( const Solution * solution , double multiplier )
{
}

/*--------------------------------------------------------------------------*/

MultiTargetSolution * MultiTargetSolution::clone( bool empty ) const
{
 auto * sol = new MultiTargetSolution();

 if( empty ) {
  if( !v_zeta.empty() )
   sol->v_zeta.resize( v_zeta.size() );
 } else {
  sol->v_zeta = v_zeta;
 }

 return ( sol );
}

/*--------------------------------------------------------------------------*/
/*-------------------- End File MultiTargetBlockv2.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
