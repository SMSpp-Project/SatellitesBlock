/*--------------------------------------------------------------------------*/
/*--------------------- File MultiTarget.cpp ------------------------*/
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
static const double M_limit = 10*3.14159265;
static const auto PI = 3.14159265;
static const auto MU = 3.986004418e14; //[m^3/s^2]
static const auto WE = 7.2921e-5;
static const auto FACTOR = 1.2;
static const auto angle0 = -1.3882860164509252;

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

// register MultiTargetBlockv2 to the Block factory

SMSpp_insert_in_factory_cpp_1( MultiTargetBlockv2 );

// register MultiTargetSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( MultiTargetSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF MultiTargetBlockv2 -----------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

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

 for( Index i = 0 ; i < horizon/3600.0 ; ++i ){
   j++;
	altitudeSetVal = cbrt( ( MU* pow((horizon)/j,2.0) ) / (4.0*pow(PI,2.0))) - RAYON;
	if ( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 )
   {
		indexLen++;
   }
 }

 FNumber altSet = indexLen;
 altitude.resize( indexLen );
 periodSat.resize( indexLen );
 indexLen = 0;

 j = 0.0;
 for( Index i = 0 ; i < horizon/3600.0 ; ++i )
 {
   j++;
	altitudeSetVal = cbrt( ( MU * pow((horizon)/j, 2.0) ) / (4.0*pow(PI, 2.0) )) - RAYON;
   //std::cout << altitudeSetVal << "\n";
	if ( altitudeSetVal >= 400000.0 && altitudeSetVal <= 1400000.0 )
   {
		altitude[indexLen] = altitudeSetVal;
		periodSat[indexLen] = horizon/j;
      indexLen++;
   }
 }

 Latitude.resize( targets );
 Longitude.resize( targets );
 for( Index i = 0 ; i < targets ; ++i )
 {
   iFile >> Latitude[i];
   iFile >> Longitude[i];
   Latitude[i] *= PI/180;
   Longitude[i] *= PI/180;
 }

 Vec_CNumber t_p;
 Vec_CNumber t_u;
 Vec_CNumber t_GM;
 Vec_CNumber thetaVal;

 t_p.resize(altSet);
 t_u.resize(altSet);
 t_GM.resize(altSet);
 thetaVal.resize(altSet);

 double alphalim;
 double Theta_min = ((2*PI*time_step)/(2*periodSat[indexLen-1]))*FACTOR;
 std::cout << "theta_min: " << Theta_min << "\n";
 double aHalf = atan(sin(Theta_min)/((RAYON+altitude[indexLen-1])/RAYON - cos(Theta_min)));

 int alpha_lim_flag = 0;
 
 for( Index ii = 0 ; ii < altSet ; ++ii ){
   if(((RAYON+altitude[ ii ])/RAYON)*sin(aHalf) > 1){
      aHalf = asin((RAYON/(RAYON+altitude[ ii ])));
      std::cout << "WARNING: computed alpha_lim\n";
      alpha_lim_flag = 1;
      break;
   } 
 }

 double alt = altitude[altSet-1];
 //double alt = altitude[0];
 for( Index ii = 0 ; ii < altSet ; ++ii ){
   t_p[ ii ] = sqrt(MU / ( pow(RAYON+altitude[ii],3)));
   t_u[ ii ] = sqrt(MU / (RAYON+altitude[ ii ]));
   t_GM[ ii ] = sqrt((RAYON+altitude[ ii ]) / MU);
   thetaVal[ ii ] = -aHalf + asin(((RAYON+altitude[ii])/RAYON)*sin(aHalf));
 }

 //double Theta_max = thetaVal[altSet-1];
 std::cout << "alpha_flag=" << alpha_lim_flag << "\n";
 
 /***********************/
 if (alpha_lim_flag == 1) {
   /////Theta_min = thetaVal[altSet-1];//*1.2;
   std::cout << "theta_min: " << Theta_min << "\n";
 }
 /***********************/
 
 //Theta_max = Theta_min;
 double numbOfDiscretize = ceil(PI/Theta_min);
 
 std::cout << "numbOfDiscretize: " << numbOfDiscretize << "\n";

 FNumber incSet = numbOfDiscretize;
 FNumber ascSet = numbOfDiscretize;
 FNumber anmSet = numbOfDiscretize;

 double start_in = 0;
 double end_in = PI;
 size_t num_in = numbOfDiscretize;

 double dx = (end_in - start_in) / (num_in - 1);
 std::vector<double> x(num_in);
 int iter = 0;
 std::generate(x.begin(), x.end(), [&] { return start_in + (iter++) * dx; });

 Vec_CNumber inclination = x;
 //std::cout << inclination;

 end_in = 2*PI;
 dx = (end_in - start_in) / (num_in - 1);
 iter = 0;
 std::generate(x.begin(), x.end(), [&] { return start_in + (iter++) * dx; });

 Vec_CNumber nodeAscendant = x;
 Vec_CNumber meanAnomaly = x;

 thetaValFinal = Theta_min;
 double altitudeFinal = altitude[indexLen-1]; 
 //double altitudeFinal = altitude[0];
 
 //std::cout << "theta_min: " << Theta_min << "\n";
 std::cout << "alpha_half: " << aHalf*180/PI << "\n";
 std::cout << "theta: " << thetaValFinal << "\n";
 std::cout << "altitude: " << altitudeFinal << "\n";

 t = horizon / time_step;

 CoverageSatLat.resize(boost::extents[targets][t][incSet*ascSet*anmSet]);
 CoverageSatLong.resize(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLat1(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong1(boost::extents[targets][t][incSet*ascSet*anmSet]);

 double lat_Sat;
 double long_Sat;

 //for( Index ii = 0 ; ii < altSet ; ++ii )
 //{
 int ii = indexLen-1;
 //int ii = 0; 
 int index1 = -1;
 indexOrbit = 0;
 int addOrbit = 0;
 int indexOrbit1 = 0;
    for( Index jj = 0 ; jj < incSet ; ++jj ) 
    {
       for( Index k = 0 ; k < ascSet ; ++k ) 
       {
         for( Index l = 0 ; l < anmSet ; ++l )
          {
            index1 += 1;
            indexOrbit1 = 0;
            for( Index j = 0 ; j < t ; ++j ) 
            {
             lat_Sat = asin(((sin(inclination[jj])*(altitude[ii]+RAYON)*sin(meanAnomaly[l]))*cos(t_p[ii]*((j)*time_step))/(altitude[ii]+RAYON)) 
		                        + ((sin(inclination[jj])*t_u[ii]*cos(meanAnomaly[l]))*sin(t_p[ii]*((j)*time_step))*t_GM[ii]));
             
             long_Sat = fmod(-(angle0 + (WE*((j)*time_step))) + 
		                        atan2((((sin(nodeAscendant[k])*(altitude[ii]+RAYON)*cos(meanAnomaly[l])) + (cos(nodeAscendant[k])*cos(inclination[jj])*
		                        (altitude[ii]+RAYON)*sin(meanAnomaly[l])))*cos(t_p[ii]*((j)*time_step))/(altitude[ii]+RAYON))	+ ((-(sin(nodeAscendant[k])*t_u[ii]*sin(meanAnomaly[l]))
		                        + (cos(nodeAscendant[k])*cos(inclination[jj])*t_u[ii]*cos(meanAnomaly[l])))*sin(t_p[ii]*((j)*time_step))*t_GM[ii]),
		                        ((((cos(nodeAscendant[k])*(altitude[ii]+RAYON)*cos(meanAnomaly[l])) - (sin(nodeAscendant[k])*cos(inclination[jj])*
		                        (altitude[ii]+RAYON)*sin(meanAnomaly[l])))* cos(t_p[ii]*(j*time_step))/(altitude[ii]+RAYON)) + ((-(cos(nodeAscendant[k])*t_u[ii]*sin(meanAnomaly[l]))
		                        -(sin(nodeAscendant[k])*cos(inclination[jj])*t_u[ii]*cos(meanAnomaly[l])))*sin((t_p[ii]*(j*time_step)))*t_GM[ii]))), (2*PI));
             if( long_Sat <= 0 )
			      long_Sat += 2*PI;
        /*
	     for( Index i = 0 ; i < targets ; ++i )                                                                                                                                                        
             {                                                                                                                                                                                              
                CoverageSatLat[ i ][ j ][ index1 ] = std::abs(Latitude[i] - lat_Sat);                                                                                                                      
                CoverageSatLong[ i ][ j ][ index1 ] = std::abs(Longitude[i] - long_Sat) * cos(Latitude[i]);
	     }}
	     indexOrbit += 1;
	     */
	     for( Index i = 0 ; i < targets ; ++i ) 
             {
                CoverageSatLat1[ i ][ j ][ index1 ] = 
		  //        std::abs(2 * asin( 0.5 * sqrt(1 - cos(Latitude[i]-lat_Sat))));
                std::abs(Latitude[i] - lat_Sat); 
                CoverageSatLong1[ i ][ j ][ index1 ] = 
		  //        std::abs(2 * asin( 0.5 * sqrt(cos(Latitude[i]) * cos(Latitude[i]) * (1 - cos(Longitude[i] - long_Sat))))) * cos(Latitude[i]);
                std::abs(Longitude[i] - long_Sat) * cos(Latitude[i]); 
                if(cos(Latitude[i]) < 0)
                  std::cout << "ERROR!" << "\n";
                  //std::cout << CoverageSatLat1[ i ][ j ][ index1 ] << " " << CoverageSatLong1[ i ][ j ][ index1 ] << std::endl;
                //std::cout << i << " " << j << " " << ii << " " << jj << " " << k << " " << l << " " << "\n";
                if (CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValFinal and CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValFinal){
                  //std::cout << i << " " << j << " " << index1 << " " << "\n";
                  indexOrbit1 += 1;
                }
             }
            }
            if (indexOrbit1>=1){
               for( Index j = 0 ; j < t ; ++j ) 
               {
                  for( Index i = 0 ; i < targets ; ++i ) 
                  {
                  CoverageSatLat[ i ][ j ][ indexOrbit ] = CoverageSatLat1[ i ][ j ][ index1 ];
                  CoverageSatLong[ i ][ j ][ indexOrbit ] = CoverageSatLong1[ i ][ j ][ index1 ];
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

 }  // end( MultiTargetBlockv2::load( memory ) )

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

 }  // end( MultiTargetBlockv2::load( istream ) )

/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::generate_abstract_variables( Configuration *stvv )
{

  if( AR3 & HasVar )  // the variables are there already
    return;           // nothing to do

  theta.resize( n );
  for( auto & var : theta )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( theta );

  Deltat.resize( targets );
  for( auto & var : Deltat )
   var.set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat );

  Deltat_k1A.resize( boost::extents[ t * (t-1)/2 ][ targets ] );
  for( Index i = 0 ; i < t * (t-1)/2 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
    Deltat_k1A[ i ][ j ].set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k1A );

  Deltat_k1.resize( boost::extents[ t-1 ][ targets ] );
  for( Index i = 0 ; i < t-1 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
    Deltat_k1[ i ][ j ].set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k1 );

  b1.resize( boost::extents[ t-1 ][ targets ] );
  for( Index i = 0 ; i < t-1 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
     b1[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( b1 );

  Deltat_k2A.resize( boost::extents[ t * (t-1)/2 ][ targets ] );
  for( Index i = 0 ; i < t * (t-1)/2 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
   Deltat_k2A[ i ][ j ].set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k2A );

  Deltat_k2.resize( boost::extents[ t-1 ][ targets ] );
  for( Index i = 0 ; i < t-1 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
    Deltat_k2[ i ][ j ].set_type( ColVariable::kNonNegative );

  add_static_variable( Deltat_k2 );

  b2.resize( boost::extents[ t-1 ][ targets ] );
  for( Index i = 0 ; i < t-1 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
     b2[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( b2 );

  zeta.resize( boost::extents[ t ][ targets ] );
   for( Index i = 0 ; i < t ; ++i )
    for( Index j = 0 ; j < targets ; ++j )
      zeta[ i ][ j ].set_type( ColVariable::kBinary );
 
  add_static_variable( zeta );

  activation.resize( boost::extents[ n ][ OrbitSet ] );
  for( Index i = 0 ; i < n ; ++i )
    for( Index j = 0 ; j < OrbitSet ; ++j )
      //for( Index k = 0 ; k < targets ; ++k )
        activation[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( activation );

  h.resize( boost::extents[ t * (t-1)/2 ][ targets ] );
  for( Index i = 0 ; i < t * (t-1)/2 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
     h[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( h );

  xi.resize( boost::extents[ n ][ t ][ targets ] );
  for( Index i = 0 ; i < n ; ++i )
    for( Index j = 0 ; j < t ; ++j )
      for( Index k = 0 ; k < targets ; ++k )
        xi[ i ][ j ][ k ].set_type( ColVariable::kBinary );

  add_static_variable( xi );

  d1.resize( boost::extents[ t * (t-1)/2 ][ targets ] );
  for( Index i = 0 ; i < t * (t-1)/2 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
     d1[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( d1 );

  d2.resize( boost::extents[ t * (t-1)/2 ][ targets ] );
  for( Index i = 0 ; i < t * (t-1)/2 ; ++i )
   for( Index j = 0 ; j < targets ; ++j )
     d2[ i ][ j ].set_type( ColVariable::kBinary );

  add_static_variable( d2 );

  AR3 |= HasVar;

 }  // end( MultiTargetBlockv2::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::generate_abstract_constraints( Configuration *stcc )
{

  if( AR2 & HasCnst )  // the constraints are there already
    return;           // nothing to do
/*
  duplicate_pi.resize( boost::multi_array_types::extent_gen()[ targets - 1  ][ satellites ][ OrbitSet ] );
    for( Index i = 0 ; i < targets - 1  ; ++i ) {
       for( Index j = 0 ; j < satellites; ++j ) {
          for( Index k = 0 ; k < OrbitSet; ++k ) {
             LinearFunction::v_coeff_pair v_vars_pi;
             v_vars_pi.push_back( std::make_pair( &activation[j][k][i], 1.0 ));
             v_vars_pi.push_back( std::make_pair( &activation[j][k][i+1], -1.0 ));
             duplicate_pi[i][j][k].set_function( new LinearFunction( std::move( v_vars_pi )));
             duplicate_pi[i][j][k].set_rhs( 0.0 ); 
             duplicate_pi[i][j][k].set_lhs( 0.0 );
          }
       }
    }
*/
  orbitSelection.resize( n );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < n ; ++i ) {
      LinearFunction::v_coeff_pair orbit_var;
      for( Index j = 0 ; j < OrbitSet ; ++j ) 
        orbit_var.push_back( std::make_pair( &activation[ i ][ j ], 1.0 ));
      LinearFunction* FunctAnm = new LinearFunction( std::move( orbit_var ));
      orbitSelection[ i ].set_rhs( 1.0 );
      orbitSelection[ i ].set_lhs( 1.0 ); 
      orbitSelection[ i ].set_function( FunctAnm );
    }
  }

  add_static_constraint( orbitSelection , "orbitSelection" );

  activationSat_cnst.resize( boost::multi_array_types::extent_gen()[ n ][ t ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < n ; ++i ){
      for( Index j = 0 ; j < t ; ++j ){
        LinearFunction::v_coeff_pair v_var;
        v_var.push_back( std::make_pair( &xi[ i ][ j ][ k ], -1.0 ));
        v_var.push_back( std::make_pair( &zeta[ j ][ k ] ,  1.0 ));
        LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
        activationSat_cnst[ i ][ j ][ k ].set_rhs( Inf< double >() );
        activationSat_cnst[ i ][ j ][ k ].set_lhs( 0.0 ); 
        activationSat_cnst[ i ][ j ][ k ].set_function( FunctSat );
      }
    }
  }
  add_static_constraint( activationSat_cnst , "activationSat_cnst" );

  activationSat1_cnst.resize( boost::multi_array_types::extent_gen()[ t ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < t ; ++i ){
      LinearFunction::v_coeff_pair v_var;
      for( Index j = 0 ; j < n ; ++j ){
        v_var.push_back( std::make_pair( &xi[ j ][ i ][ k ], 1.0 ));
      }
      v_var.push_back( std::make_pair( &zeta[ i ][ k ] ,  -1.0 ));
      LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
      activationSat1_cnst[ i ][ k ].set_rhs( Inf< double >() );
      activationSat1_cnst[ i ][ k ].set_lhs( 0.0 ); 
      activationSat1_cnst[ i ][ k ].set_function( FunctSat );
    }
  }
  add_static_constraint( activationSat1_cnst , "activationSat1_cnst" );


  h_cnst_1.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ){
      for( Index j = i+1 ; j < t ; ++j ){
        //if ( j > i ) {
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &h[ ii ][ k ], -1.0 ));
          v_var.push_back( std::make_pair( &zeta[ i ][ k ] ,  1.0 ));
          LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
          h_cnst_1[ ii ][ k ].set_rhs( Inf< double >() );
          h_cnst_1[ ii ][ k ].set_lhs( 0.0 ); 
          h_cnst_1[ ii ][ k ].set_function( FunctSat );
          ii += 1;
      }
    }
  }
  add_static_constraint( h_cnst_1 , "h_cnst_1" );


  h_cnst_2.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ){
      for( Index j = i+1 ; j < t ; ++j ){
        //if ( j > i ) {
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &h[ ii ][ k ], -1.0 ));
          v_var.push_back( std::make_pair( &zeta[ j ][ k ] ,  1.0 ));
          LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
          h_cnst_2[ ii ][ k ].set_rhs( Inf< double >() );
          h_cnst_2[ ii ][ k ].set_lhs( 0.0 ); 
          h_cnst_2[ ii ][ k ].set_function( FunctSat );
          ii += 1;
      }
    }
  }
  add_static_constraint( h_cnst_2 , "h_cnst_2" );

  h_cnst_3.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ){
      for( Index j = i+1 ; j < t ; ++j ){
        //if ( j > i ) {
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &h[ ii ][ k ], 1.0 ));
          v_var.push_back( std::make_pair( &zeta[ i ][ k ] ,  -1.0 ));
          v_var.push_back( std::make_pair( &zeta[ j ][ k ] ,  -1.0 ));
          LinearFunction* FunctSat = new LinearFunction( std::move( v_var ));
          h_cnst_3[ ii ][ k ].set_rhs( Inf< double >() );
          h_cnst_3[ ii ][ k ].set_lhs( -1.0 ); 
          h_cnst_3[ ii ][ k ].set_function( FunctSat );
          ii += 1;
      }
    }
  }
  add_static_constraint( h_cnst_3 , "h_cnst_3" );

  // generate observability constraints via big-M approach

  //Thetamax[i]+(1-Xi[i,p,j])*M_limit >= (sum(act4[i,a,k,l,s]*CoverageSatLat[j,p,s,l,k,a]
	//	for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

  //Thetamax[i]+(1-Xi[i,p,j])*M_limit >= (sum(act4[i,a,k,l,s]*CoverageSatLong[j,p,s,l,k,a]
	//	for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

  //Thetamax[i]-(Xi[i,p,j])*M_limit <= (sum(act4[i,a,k,l,s]*CoverageSatLat[j,p,s,l,k,a]
	//	for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

  //Thetamax[i]-(Xi[i,p,j])*M_limit <= (sum(act4[i,a,k,l,s]*CoverageSatLong[j,p,s,l,k,a]
	//	for s=1:n_inclinaison for l=1:n_noeudAscendant for k in 1:n_meanAnomaly for a=1:length(altitudeSet))))

  obs1_cnst.resize(
    boost::multi_array_types::
    extent_gen()[ n ][ t ][ targets ] );
  obs3_cnst.resize(
    boost::multi_array_types::
    extent_gen()[ n ][ t ][ targets ] );

  obs2_cnst.resize(
    boost::multi_array_types::
    extent_gen()[ n ][ t ][ targets ] );
  obs4_cnst.resize(
    boost::multi_array_types::
    extent_gen()[ n ][ t ][ targets ] );

  double MLAT = PI;
  double MLONG = 2*PI;

  for( Index k = 0 ; k < targets ; ++k ){

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
                                  -CoverageSatLat[ k ][ j ][ jj ])); 
            v_obs2.push_back( std::make_pair( &activation[ i ][ jj ], 
                                  -CoverageSatLong[ k ][ j ][ jj ])); 
            MLAT = std::max(MLAT, (CoverageSatLat[ k ][ j ][ jj ]));
            MLONG = std::max(MLONG, (CoverageSatLong[ k ][ j ][ jj ]));
          }

        v_obs1.push_back( std::make_pair( &xi[i][j][ k ], -MLAT )); 
        v_obs2.push_back( std::make_pair( &xi[i][j][ k ], -MLONG )); 

        v_obs1.push_back( std::make_pair( &theta[i], 1.0 )); 
        v_obs2.push_back( std::make_pair( &theta[i], 1.0 )); 

        LinearFunction* Funct1 = new LinearFunction( std::move( v_obs1 ));
        LinearFunction* Funct2 = new LinearFunction( std::move( v_obs2 ));

        obs2_cnst[ i ][ j ][ k ].set_rhs( Inf< double >() );
        obs2_cnst[ i ][ j ][ k ].set_lhs( (-MLAT )); 
        obs2_cnst[ i ][ j ][ k ].set_function( Funct1 );

        obs4_cnst[ i ][ j ][ k ].set_rhs( Inf< double >() );
        obs4_cnst[ i ][ j ][ k ].set_lhs( (-MLONG )); 
        obs4_cnst[ i ][ j ][ k ].set_function( Funct2 );

        obs1_cnst[ i ][ j ][ k ].set_rhs( (-thetaValFinal ));
        obs1_cnst[ i ][ j ][ k ].set_lhs( -Inf< double >() ); 
        obs1_cnst[ i ][ j ][ k ].set_function( Funct1 );

        obs3_cnst[ i ][ j ][ k ].set_rhs( (-thetaValFinal ));
        obs3_cnst[ i ][ j ][ k ].set_lhs( -Inf< double >() ); 
        obs3_cnst[ i ][ j ][ k ].set_function( Funct2 );

      }
    }
  }

  //add_static_constraint( obs1_cnst );
  //add_static_constraint( obs3_cnst );
  
  add_static_constraint( obs2_cnst );
  add_static_constraint( obs4_cnst );

  theta_UB.resize( n );
  for( Index i = 0 ; i < n ; ++i ){
    LinearFunction::v_coeff_pair v_var_theta;
    v_var_theta.push_back( std::make_pair( &theta[ i ], 1.0 ));
    LinearFunction* Funct_theta = new LinearFunction( std::move( v_var_theta ));
    theta_UB[ i ].set_rhs( thetaValFinal );
    theta_UB[ i ].set_lhs( -Inf< double >() );
    theta_UB[ i ].set_function( Funct_theta );
  }

  add_static_constraint( theta_UB , "theta_UB" );

  Deltat_max_dt1.resize( targets );
  for( Index k = 0 ; k < targets ; ++k ){
    LinearFunction::v_coeff_pair v_vart1;
    v_vart1.push_back( std::make_pair( &Deltat[ k ], 1.0 ));
    LinearFunction* Functt1 = new LinearFunction( std::move( v_vart1 ));
    Deltat_max_dt1[ k ].set_rhs( Inf< double >() );
    Deltat_max_dt1[ k ].set_lhs( dt );
    Deltat_max_dt1[ k ].set_function( Functt1 );
  }

  add_static_constraint( Deltat_max_dt1 , "Deltat_max_dt1" );

  Deltat_max1.resize( boost::multi_array_types::extent_gen()[ t-1 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < t-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_var;
      v_var.push_back( std::make_pair( &Deltat[ k ], 1.0 ));
      v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ], -1.0 ));
      v_var.push_back( std::make_pair( &b1[ i ][ k ], 0.5*T ));
      LinearFunction* Funct = new LinearFunction( std::move( v_var ));
      Deltat_max1[ i ][ k ].set_rhs( Inf< double >() );
      Deltat_max1[ i ][ k ].set_lhs( 0.0 ); 
      Deltat_max1[ i ][ k ].set_function( Funct );
    }
  }

  add_static_constraint( Deltat_max1 , "Deltat_max1" );

  Deltat_max11.resize( boost::multi_array_types::extent_gen()[ t-1 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < t-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_var;
      v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ], 1.0 ));
      v_var.push_back( std::make_pair( &b1[ i ][ k ], -0.5*T ));
      LinearFunction* Funct = new LinearFunction( std::move( v_var ));
      Deltat_max11[ i ][ k ].set_rhs( Inf< double >() );
      Deltat_max11[ i ][ k ].set_lhs( 0.0 ); 
      Deltat_max11[ i ][ k ].set_function( Funct );
    }
  }

  add_static_constraint( Deltat_max11 , "Deltat_max11" );

  Deltat_max2.resize( boost::multi_array_types::extent_gen()[ t-1 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < t-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_var;
      v_var.push_back( std::make_pair( &Deltat[ k ], 1.0 ));
      v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ], -1.0 ));
      v_var.push_back( std::make_pair( &b2[ i ][ k ], 0.5*T ));
      LinearFunction* Funct = new LinearFunction( std::move( v_var ));
      Deltat_max2[ i ][ k ].set_rhs( Inf< double >() );
      Deltat_max2[ i ][ k ].set_lhs( 0.0 ); 
      Deltat_max2[ i ][ k ].set_function( Funct );
    }
  }

  add_static_constraint( Deltat_max2 , "Deltat_max2" );

  Deltat_max22.resize( boost::multi_array_types::extent_gen()[ t-1 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < t-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_var;
      v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ], 1.0 ));
      v_var.push_back( std::make_pair( &b2[ i ][ k ], -0.5*T ));
      LinearFunction* Funct = new LinearFunction( std::move( v_var ));
      Deltat_max22[ i ][ k ].set_rhs( Inf< double >() );
      Deltat_max22[ i ][ k ].set_lhs( 0.0 ); 
      Deltat_max22[ i ][ k ].set_function( Funct );
    }
  }

  add_static_constraint( Deltat_max22 , "Deltat_max22" );


  Deltat_min_k1_1.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      for( Index j = i+1 ; j < t ; ++j ) {
        //if ( j > i ){
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ], -1.0 ));
          v_var.push_back( std::make_pair( &h[ ii ][ k ], (( j - i)*dt - 0.5*T )));
          LinearFunction* Funct = new LinearFunction( std::move( v_var ));
          Deltat_min_k1_1[ ii ][ k ].set_rhs( Inf< double >() );
          Deltat_min_k1_1[ ii ][ k ].set_lhs( -0.5*T ); 
          Deltat_min_k1_1[ ii ][ k ].set_function( Funct );
          ii += 1;
      }
    }
  }

  add_static_constraint( Deltat_min_k1_1 , "Deltat_min_k1_1" );

  Deltat_min_k1_2.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      for( Index j = i+1 ; j < t ; ++j ) {
        //if ( j > i ){
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &Deltat_k1[ i ][ k ], -1.0 ));
          v_var.push_back( std::make_pair( &h[ ii ][ k ], (( j - i)*dt - 0.5*T )));
          v_var.push_back( std::make_pair( &d1[ ii ][ k ], 0.5*T));
          LinearFunction* Funct = new LinearFunction( std::move( v_var ));
          Deltat_min_k1_2[ ii ][ k ].set_rhs( 0.0 );
          Deltat_min_k1_2[ ii ][ k ].set_lhs( -Inf< double >() ); 
          Deltat_min_k1_2[ ii ][ k ].set_function( Funct );
          ii += 1;
      }
    }
  }

  add_static_constraint( Deltat_min_k1_2 , "Deltat_min_k1_2" );

  d1_cnst.resize( boost::multi_array_types::extent_gen()[ t-1 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_var;
      for( Index j = i+1 ; j < t ; ++j ) {
    v_var.push_back( std::make_pair( &d1[ ii ][ k ], 1.0));
          ii += 1;
      }
      LinearFunction* Funct = new LinearFunction( std::move( v_var ));
      d1_cnst[ i ][ k ].set_rhs( 1.0 );
      d1_cnst[ i ][ k ].set_lhs( 1.0 ); 
      d1_cnst[ i ][ k ].set_function( Funct );
    }
  }

  add_static_constraint( d1_cnst , "d1_cnst" );

  Deltat_min_k2_1.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      for( Index j = i+1 ; j < t ; ++j ) {
        //if ( j > i ){
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ], 1.0 ));
          v_var.push_back( std::make_pair( &h[ ii ][ k ], (( j - i)*dt - 0.5*T)));
          LinearFunction* Funct = new LinearFunction( std::move( v_var ));
          Deltat_min_k2_1[ ii ][ k ].set_rhs( 0.5*T );
          Deltat_min_k2_1[ ii ][ k ].set_lhs( -Inf< double >() ); 
          Deltat_min_k2_1[ ii ][ k ].set_function( Funct );
          ii += 1;
      }
    }
  }

  add_static_constraint( Deltat_min_k2_1 , "Deltat_min_k2_1" );

  Deltat_min_k2_2.resize( boost::multi_array_types::extent_gen()[ t*(t-1)/2 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      for( Index j = i+1 ; j < t ; ++j ) {
        //if ( j > i ){
          LinearFunction::v_coeff_pair v_var;
          v_var.push_back( std::make_pair( &Deltat_k2[ i ][ k ], 1.0 ));
          v_var.push_back( std::make_pair( &h[ ii ][ k ], (( j - i)*dt - 0.5*T)));
          v_var.push_back( std::make_pair( &d2[ ii ][ k ], -0.5*T));
          LinearFunction* Funct = new LinearFunction( std::move( v_var ));
          Deltat_min_k2_2[ ii ][ k ].set_rhs( Inf< double >() );
          Deltat_min_k2_2[ ii ][ k ].set_lhs( 0.0 ); 
          Deltat_min_k2_2[ ii ][ k ].set_function( Funct );
          ii += 1;
      }
    }
  }

  add_static_constraint( Deltat_min_k2_2 , "Deltat_min_k2_2" );

  d2_cnst.resize( boost::multi_array_types::extent_gen()[ t-1 ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      LinearFunction::v_coeff_pair v_var;
      for( Index j = i+1 ; j < t ; ++j ) {
          v_var.push_back( std::make_pair( &d2[ ii ][ k ], 1.0));
          ii += 1;
      }
      LinearFunction* Funct = new LinearFunction( std::move( v_var ));
      d2_cnst[ i ][ k ].set_rhs( 1.0 );
      d2_cnst[ i ][ k ].set_lhs( 1.0 ); 
      d2_cnst[ i ][ k ].set_function( Funct );
    }
  }

  add_static_constraint( d2_cnst , "d2_cnst" );
  
  obs_cnst_h.resize( targets );
  for( Index k = 0 ; k < targets ; ++k ){
    Index ii = 0;
    LinearFunction::v_coeff_pair v_var1;
    for( Index i = 0 ; i < t-1 ; ++i ) {
      for( Index j = i+1 ; j < t ; ++j ) {
        v_var1.push_back( std::make_pair( &h[ ii ][ k ], 1.0));
        ii += 1;
      }
    }
    LinearFunction* FunctA = new LinearFunction( std::move( v_var1 ));
    obs_cnst_h[ k ].set_rhs( Inf< double >() );
    obs_cnst_h[ k ].set_lhs( 1.0 );
    obs_cnst_h[ k ].set_function( FunctA );
  }

  add_static_constraint( obs_cnst_h , "obs_cnst_h" );
  

  obs_cnst_xi.resize( boost::multi_array_types::extent_gen()[ n ][ targets ] );
  for( Index k = 0 ; k < targets ; ++k ){
    for( Index i = 0 ; i < n ; ++i ) {
      LinearFunction::v_coeff_pair v_varxi;
      for( Index j = 0 ; j < t ; ++j ) {
        v_varxi.push_back( std::make_pair( &xi[ i ][ j ][ k ], 1.0));
      }
      LinearFunction* FunctB = new LinearFunction( std::move( v_varxi ));
      obs_cnst_xi[ i ][ k ].set_rhs( Inf< double >() );
      obs_cnst_xi[ i ][ k ].set_lhs( 3.0 );
      obs_cnst_xi[ i ][ k ].set_function( FunctB );
    }
  }

  add_static_constraint( obs_cnst_xi , "obs_cnst_xi" );

  observation1.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()[ targets ][ t ]);

  for( Index i = 0 ; i < targets ; ++i ){
      for( Index j = 0 ; j < t ; ++j ){         
         LinearFunction::v_coeff_pair v_var1;
         for( Index k = 0 ; k < n ; ++k )
            v_var1.push_back( std::make_pair( &xi[ k ][ j ][ i ], 1.0));

         observation1[i][j].set_function( new LinearFunction( std::move( v_var1 )));
         observation1[i][j].set_rhs( 1.0 );//Inf< double >() 
         observation1[i][j].set_lhs( -Inf< double >()  );
      }
   }

   add_static_constraint( observation1 , "observation1" );

  AR2 |= HasCnst;
 }  // end( MultiTargetBlockv2::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::generate_objective( Configuration *objc )
{

 if( AR1 & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return
 
  LinearFunction::v_coeff_pair v_obj;

  for( Index k = 0 ; k < targets ; ++k ) {
    v_obj.push_back( std::make_pair( &Deltat[ k ], 1.0/targets));
  }

  LinearFunction* Functobj = new LinearFunction( std::move( v_obj ));
  
  c.set_function( Functobj );
  set_objective( &c );
  
  AR1 |= HasObj;

 }  // end( MultiTargetBlockv2::generate_objective )

/*--------------------------------------------------------------------------*/

 bool MultiTargetBlockv2::is_feasible( bool useabstract , Configuration *fsbc )
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

 }  // end( MultiTargetBlockv2::is_feasible )

/*--------------------------------------------------------------------------*/

bool MultiTargetBlockv2::is_optimal( bool useabstract , Configuration *optc )
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

 }  //  end( MultiTargetBlockv2::is_optimal )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/


 Solution * MultiTargetBlockv2::get_Solution( Configuration * solc , bool emptys )
{
   
 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;
 
 auto *sol = new MultiTargetSolution();

 if( ! emptys )
  sol->read( this );
 
 return( sol );

 }  // end( MultiTargetBlockv2::get_Solution )

 /*
 void MultiTargetBlockv2::set_zeta( c_Vec_FNumber_it fstrt , Range rng )
{
 if( ! ( AR3 & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 Index i = rng.first;

 for( auto xi = Deltat.begin() + i ;
       i < 1 ; ++i )
   (xi++)->set_value( *(fstrt++) );
}
 */


/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );

 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE MultiTargetBlockv2 ---*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::print( std::ostream  & output , char vlvl ) const
{
 
 }  // end( MultiTargetBlockv2::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::guts_of_destructor( void )
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

 }  // end( MultiTargetBlockv2::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void MultiTargetBlockv2::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since MultiTargetBlockv2 is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN MultiTargetBlockv2 WHEN THIS METHOD IS
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

 //throw( std::invalid_argument( "unsupported Modification to MultiTargetBlockv2" ) );

 }  // end( MultiTargetBlockv2::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF DCRSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void MultiTargetSolution::deserialize( const netCDF::NcGroup & group )
{}

void MultiTargetSolution::read( const Block * block )
{
 auto SATB = dynamic_cast< const MultiTargetBlockv2 * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a MultiTargetBlockv2" ) );

 if( ! v_zeta.empty() ) {
  v_zeta.resize( 1 );

  //SATB->get_zeta();
  }
}

void MultiTargetSolution::write( Block * block ) 
{

 auto SATB = dynamic_cast<MultiTargetBlockv2 * >( block );
 if( ! SATB )
  throw( std::invalid_argument( "block is not a MultiTargetBlockv2" ) );

 if( ! v_zeta.empty() ) {
  //SATB->set_zeta( v_zeta.begin() );
  }

}

void MultiTargetSolution::serialize( netCDF::NcGroup & group ) const
{}

MultiTargetSolution * MultiTargetSolution::scale( double factor ) const
{
  auto * sol = MultiTargetSolution::clone( true );
  return( sol );
}

void MultiTargetSolution::sum( const Solution * solution , double multiplier )
{}

MultiTargetSolution * MultiTargetSolution::clone( bool empty ) const
{
  auto * sol = new MultiTargetSolution();

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
/*------------------- End File MultiTargetBlockv2.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
