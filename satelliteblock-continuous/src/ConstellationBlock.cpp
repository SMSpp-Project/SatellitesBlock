/*--------------------------------------------------------------------------*/
/*------------------------- File ConstellationBlock.cpp ---------------------*/
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

void  ConstellationBlock::load( std::istream & input , char frmt ){
   // TO DO: implement load() method for loading instance data 
   // from input file (for the file format, see next load() method)
}

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
 periods.resize( targets );
 for( Index i = 0 ; i < targets ; ++i )
 {
   iFile >> periods[i];
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

 for( Index ii = 0 ; ii < altSet ; ++ii ){
   if(((RAYON+altitude[ ii ])/RAYON)*sin(aHalf) > 1){
      aHalf = asin((RAYON/(RAYON+altitude[ ii ])));
      std::cout << "WARNING: computed alpha_lim\n";
      break;
   } 
 }

 double alt = altitude[altSet-1];
 for( Index ii = 0 ; ii < altSet ; ++ii ){
   t_p[ ii ] = sqrt(MU / ( pow(RAYON+altitude[ii],3)));
   t_u[ ii ] = sqrt(MU / (RAYON+altitude[ ii ]));
   t_GM[ ii ] = sqrt((RAYON+altitude[ ii ]) / MU);
   thetaVal[ ii ] = -aHalf + asin(((RAYON+altitude[ii])/RAYON)*sin(aHalf));
 }

 double Theta_max = thetaVal[altSet-1];
 std::cout << "theta_min: " << Theta_min << "\n";
 Theta_max = Theta_min;
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

 double thetaValFinal = Theta_min;
 double altitudeFinal = altitude[indexLen-1]; 
 ///double altitudeFinal = altitude[0];
 
 //std::cout << "theta_min: " << Theta_min << "\n";
 std::cout << "alpha_half: " << aHalf*180/PI << "\n";
 std::cout << "theta: " << thetaValFinal << "\n";
 std::cout << "altitude: " << altitudeFinal << "\n";

 FNumber t = horizon / time_step;

 boost::multi_array< double , 3 > CoverageSatLat(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLat1(boost::extents[targets][t][incSet*ascSet*anmSet]);
 boost::multi_array< double , 3 > CoverageSatLong1(boost::extents[targets][t][incSet*ascSet*anmSet]);

 double lat_Sat;
 double long_Sat;

 int indexOrbit;
 int indexOrbit1;

 iFile >> satellites;
 //std::cout << "number Orbits: " << indexOrbit << "\n";
 /*
 satellites = 0;
 for( Index i = 0 ; i < targets ; ++i )
 {
   satellites += periods[i];
 }
 */
 v_Block.resize( satellites );
 std::cout << "SATELLITES: " << satellites << "\n";

 int satellites1 = satellites;

 for(Index isat = 0; isat < satellites; ++isat){
 
 int ialt = std::floor(satellites/3);
 int ii = 0;
 if( isat < ialt )
   ii = 2;
 if( isat > 2*ialt)
   ii = 0;
 if(isat >= ialt and isat <= 2*ialt)
   ii = 1;

 int index1 = -1;
 indexOrbit = 0;
 int addOrbit = 0;
 indexOrbit1 = 0;
 thetaValF = thetaValFinal/3.0;
 std::cout << thetaValF << std::endl;
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
             // formula to compute the latitude of the projection of the position 
             // of satellite onto the Earth surface corresponding to a given configuration
             lat_Sat = asin(((sin(inclination[jj])*(altitude[ii]+RAYON)*sin(meanAnomaly[l]))*cos(t_p[ii]*((j)*time_step))/(altitude[ii]+RAYON)) 
		                        + ((sin(inclination[jj])*t_u[ii]*cos(meanAnomaly[l]))*sin(t_p[ii]*((j)*time_step))*t_GM[ii]));
             
             // formula to compute the longitude of the projection of the position 
             // of satellite onto the Earth surface corresponding to a given configuration                  
             long_Sat = fmod(-(angle0 + (WE*((j)*time_step))) + 
		                        atan2((((sin(nodeAscendant[k])*(altitude[ii]+RAYON)*cos(meanAnomaly[l])) + (cos(nodeAscendant[k])*cos(inclination[jj])*
		                        (altitude[ii]+RAYON)*sin(meanAnomaly[l])))*cos(t_p[ii]*((j)*time_step))/(altitude[ii]+RAYON))	+ ((-(sin(nodeAscendant[k])*t_u[ii]*sin(meanAnomaly[l]))
		                        + (cos(nodeAscendant[k])*cos(inclination[jj])*t_u[ii]*cos(meanAnomaly[l])))*sin(t_p[ii]*((j)*time_step))*t_GM[ii]),
		                        ((((cos(nodeAscendant[k])*(altitude[ii]+RAYON)*cos(meanAnomaly[l])) - (sin(nodeAscendant[k])*cos(inclination[jj])*
		                        (altitude[ii]+RAYON)*sin(meanAnomaly[l])))* cos(t_p[ii]*(j*time_step))/(altitude[ii]+RAYON)) + ((-(cos(nodeAscendant[k])*t_u[ii]*sin(meanAnomaly[l]))
		                        -(sin(nodeAscendant[k])*cos(inclination[jj])*t_u[ii]*cos(meanAnomaly[l])))*sin((t_p[ii]*(j*time_step)))*t_GM[ii]))), (2*PI));
             if( long_Sat <= 0 )
			      long_Sat += 2*PI;

	     for( Index i = 0 ; i < targets ; ++i ) 
             {
                // compute the (geodedical) difference between the latitude of the target and the lat_Sat
                CoverageSatLat1[ i ][ j ][ index1 ] = 
                        std::abs(2 * asin( 0.5 * sqrt(1 - cos(Latitude[i]-lat_Sat))));
                CoverageSatLong1[ i ][ j ][ index1 ] = std::abs(2 * asin( 0.5 * sqrt(cos(Latitude[i]) * cos(Latitude[i]) * (1 - cos(Longitude[i] - long_Sat))))) * cos(Latitude[i]);
                // compute the scaled (geodedical) difference between the latitude of the target and the long_Sat
                if(cos(Latitude[i]) < 0)
                  std::cout << "ERROR!" << "\n";
                if (CoverageSatLat1[ i ][ j ][ index1 ] <= thetaValF and CoverageSatLong1[ i ][ j ][ index1 ] <= thetaValF)
                  indexOrbit1 += 1;
             }
            }
            // the orbital configuration that do not observe any satellite in any time-step 
            // are discarded so that the solution space is maintened reasonably "small"
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
      std::cout << "number Orbits: " << indexOrbit << "\n";
      auto SB = new SatelliteBlock( this );
      SB->load(targets , time_step, horizon, altitudeFinal , thetaValF ,
                  indexOrbit , aHalf ,
                  CoverageSatLat , CoverageSatLong , periods );
      v_Block[ isat ] = SB;
 }

 //std::cout << satellites1 << std::endl;
 //satellites = satellites1;
    
 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( ConstellationBlock::load( const std::string & input )

/*--------------------------------------------------------------------------*/

void ConstellationBlock::generate_abstract_variables( Configuration * stvv )
{

  for( auto blck : v_Block )
    blck->generate_abstract_variables();
 }

/*--------------------------------------------------------------------------*/

void ConstellationBlock::generate_abstract_constraints( Configuration * stcc )
{

 if(!AR){
   for( auto blck : v_Block )
    blck->generate_abstract_constraints();

   thetaM.resize( 1 );

   LinearFunction::v_coeff_pair v_var12;
   for( Index k = 0 ; k < satellites ; ++k ){
      v_var12.push_back( std::make_pair( static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_theta(), 1.0 ));
      v_var12.push_back( std::make_pair( static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_z(), -thetaValF * 0.9 ));
   }
   thetaM[0].set_function( new LinearFunction( std::move( v_var12 )));
   thetaM[0].set_rhs( 0.0 );//Inf< double >()
   thetaM[0].set_lhs( -Inf< double >() );

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

   double maxPeriods = *max_element(periods.begin(), periods.end());
   observation.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()[ targets ][ maxPeriods ]);
   double pp;

   for( Index i = 0 ; i < targets ; ++i ) {
    pp = horizon/time_step/periods[i];
    for( Index j = 0 ; j < periods[i] ; ++j ) {
      LinearFunction::v_coeff_pair v_var;
      
      for( Index k = 0 ; k < satellites ; ++k ) {
         for( Index tt = j*pp ; tt < (j+1)*pp ; ++tt ) {
            // retrieve observation variable \xi[ i ][ t ][ m ] for SatelliteBlock i
            v_var.push_back( std::make_pair( static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_r( i , tt ), 1.0 ));
         }
      }
      observation[i][j].set_function( new LinearFunction( std::move( v_var )));
      observation[i][j].set_rhs( Inf< double >() ); 
      observation[i][j].set_lhs( 1.0 );
    }

    for( Index j = periods[i] ; j < maxPeriods ; ++j ) {
       LinearFunction::v_coeff_pair v_var;
      
       for( Index k = 0 ; k < satellites ; ++k ) {
          v_var.push_back( std::make_pair( static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_z(), 0.0 ));
       }

      // fake constraints when j \geq periods[ i ]: we simply set 0 * z[ i ] == 0 
      // for SatelliteBlock i \in [s], where s is the total number of the satellite
      observation[i][j].set_function( new LinearFunction( std::move( v_var )));
      observation[i][j].set_rhs( 0.0 ); 
      observation[i][j].set_lhs( 0.0 );
    }

   }

   add_static_constraint( observation , "observation" );

   FNumber t = horizon / time_step;

   observation1.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()[ targets ][ t ]);

   for( Index i = 0 ; i < targets ; ++i ){
      for( Index j = 0 ; j < t ; ++j ){         
         LinearFunction::v_coeff_pair v_var1;
         for( Index k = 0 ; k < satellites ; ++k )
            v_var1.push_back( std::make_pair( static_cast< SatelliteBlock * >( v_Block[ k ] )->i2p_r(i,j), 1.0 ));

         observation1[i][j].set_function( new LinearFunction( std::move( v_var1 )));
         observation1[i][j].set_rhs( 1.0 );//Inf< double >() 
         observation1[i][j].set_lhs( -Inf< double >()  );
      }
   }

   add_static_constraint( observation1 , "observation1" );

   std::cout << "Constraints charged!\n";
}
AR = true;

 }  // end( ConstellationBlock::generate_abstract_constraints() )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void ConstellationBlock::print( std::ostream & output , char vlvl ) const
{
 }

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

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) one
 reset_static_constraints();
 // not needed, there isn't any - reset_static_variables();
 // not needed, there isn't any - reset_dynamic_constraints();
 // not needed, there isn't any - reset_dynamic_variables();
 // not needed, there isn't any - reset_objective();

 }  // end( guts_of_destructor )

/*--------------------------------------------------------------------------*/
/*---------------------- End File ConstellationBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
